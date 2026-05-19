"""PhageFlow Module 03 — De novo assembly.

Tools and parameters (literature-based for purified phage sequencing):

SPAdes v3.15+ (Bankevich et al. 2012, Journal of Computational Biology):
    --meta
        Metagenomic mode: handles variable coverage better than --isolate
        for mixed preparations (phage + contamination residual).
        Disables repeat resolving that can collapse phage terminal repeats.
        Recommended by Prjibelski et al. 2020 (Current Protocols) for viral.

    --only-assembler
        Skip BayesHammer error correction — required at >200× coverage
        (BayesHammer collapses true variants into consensus, losing
        biological diversity in tail fiber / receptor binding proteins).
        Roux et al. 2019 (eLife): explicitly recommended for viral datasets.

    -k 21,33,55,77,99,127
        Multiple k-mers ensure coverage of both short repeats (k=21)
        and long unique regions (k=127). Standard for PE150 phage data
        (Bankevich et al. 2012; SPAdes manual recommendation for PE150).

    Post-assembly length filter: contigs < min_length (default 200 bp) are
        discarded. Phage CDSs average 600–800 bp (Casjens 2005, Curr Opin);
        contigs <200 bp cannot encode a full gene and are likely assembly
        artefacts at high coverage. 200 bp is more permissive than the
        500 bp default — appropriate for phage where compact genomes make
        every contig relevant. SPAdes does not have a built-in --min-contig-len
        flag, so filtering is applied post-hoc with a simple FASTA parser.

    Note: --cov-cutoff auto is INCOMPATIBLE with --meta (SPAdes ≥3.15)
    Note: contigs.fasta used (NOT scaffolds.fasta) — scaffolds introduce
          artificial N-gaps that break ORF prediction in Pharokka/Phold
          (geNomad, CheckV, Pharokka all use contigs as standard input)

MEGAHIT v1.2+ (Li et al. 2015, Bioinformatics):
    --k-list 21,29,39,59,79,99,119,141
        Progressive k-mer list — covers a wider range than SPAdes default,
        recovering contigs missed at intermediate k values.
        Li et al. 2015: recommended for PE150 metagenomic data.

    --min-count 2
        Minimum solid k-mer multiplicity — removes sequencing errors
        without losing low-coverage genuine phage k-mers.
        Li et al. 2016 (Methods): min-count=2 optimal for >50× coverage.

    --no-mercy
        Disable mercy k-mers (used to recover low-coverage tips).
        At phage coverage (>200×), mercy k-mers add noise rather than
        signal. Li et al. 2015: disable for high-coverage datasets.

    --min-contig-len 200
        Lowered from 500 to 200 bp — consistent with the SPAdes post-filter.
        Phage genomes are compact; a 200–499 bp contig can encode a partial
        or small CDS (e.g. holin, spanin subunit) with functional relevance.
        CheckV and geNomad handle short contigs gracefully.

cd-hit-est (Fu et al. 2012, Bioinformatics):
    -c 1.00 (100% identity): removes exact duplicates between assemblers
    while preserving biological variants (e.g. tail fiber diversity).
    At 95%, true SNP-level variants collapse — not appropriate for phage.
"""

from __future__ import annotations
from pathlib import Path
from typing import List

from phageflow.utils.config import Config
from phageflow.utils.logger import log_step, log_info, log_ok, log_warn
from phageflow.utils.tools import require_tools, run_silent, human_size, mkdirs, fasta_stats

STEP  = "03_assembly"
TOOLS = ["spades.py", "megahit", "cd-hit-est"]


def run(
    cfg:       Config,
    sample_id: str,
    r1:        Path,
    r2:        Path,
    force:     bool = False,
) -> Path:
    """
    Assemble a single sample with SPAdes + MEGAHIT, then dereplicate.

    Returns
    -------
    nr_fa : Path to non-redundant contigs FASTA (input for geNomad)
    """
    require_tools(*TOOLS)

    r1 = Path(r1); r2 = Path(r2)
    out_dir = cfg.results(STEP)
    rpt_dir = cfg.reports(STEP)
    mkdirs(out_dir / "spades", out_dir / "megahit", out_dir / "combined", rpt_dir)

    min_len = cfg.assembly.min_length   # 200 bp

    log_step(f"Module 03 — assembly [{sample_id}] · {cfg.threads} threads")
    log_info(f"  R1 : {human_size(r1)}")
    log_info("  SPAdes: --meta --only-assembler (Roux et al. 2019)")
    log_info("  MEGAHIT: --no-mercy --min-count 2 (Li et al. 2015)")
    log_info(f"  Min contig length: {min_len} bp (post-SPAdes filter + MEGAHIT flag)")
    log_info("  Dereplication: cd-hit-est 100% identity (Fu et al. 2012)")

    # ── SPAdes ────────────────────────────────────────────────────────────────
    sp_dir     = out_dir / "spades" / sample_id
    sp_raw     = sp_dir / "contigs.fasta"        # raw SPAdes output
    sp_fa      = sp_dir / "contigs_filtered.fasta"  # after length filter
    mkdirs(sp_dir)

    if sp_fa.exists() and sp_fa.stat().st_size > 0 and not force:
        log_info("  [SPAdes] already exists — skipping")
    else:
        log_info("  [SPAdes] assembling...")
        try:
            run_silent([
                "spades.py",
                "--pe1-1", str(r1), "--pe1-2", str(r2),
                "--meta",
                "--only-assembler",
                "-k", cfg.assembly.kmers,
                "--threads", str(cfg.threads),
                "--memory",  str(cfg.memory_gb),
                "-o", str(sp_dir),
            ], log_file=rpt_dir / f"{sample_id}_spades.log")
        except Exception as e:
            log_warn(f"  [SPAdes] warning: {e}")

        # Post-assembly length filter
        if sp_raw.exists() and sp_raw.stat().st_size > 0:
            n_raw, n_kept = _filter_fasta(sp_raw, sp_fa, min_len)
            n_disc = n_raw - n_kept
            log_info(
                f"  [SPAdes] length filter ≥{min_len}bp: "
                f"{n_raw} → {n_kept} contigs (discarded {n_disc})"
            )
        else:
            log_warn("  [SPAdes] no raw output to filter")

    if sp_fa.exists() and sp_fa.stat().st_size > 0:
        s = fasta_stats(sp_fa)
        log_ok(f"  [SPAdes] {s['n']} contigs | "
               f"largest={s['largest_bp']}bp | N50={s['n50']}bp | "
               f"total={s['total_bp']}bp")
    else:
        log_warn("  [SPAdes] no output after filtering")

    # ── MEGAHIT ───────────────────────────────────────────────────────────────
    mh_dir = out_dir / "megahit" / sample_id
    mh_fa  = mh_dir / "final.contigs.fa"

    if mh_fa.exists() and mh_fa.stat().st_size > 0 and not force:
        log_info("  [MEGAHIT] already exists — skipping")
    else:
        import shutil
        if mh_dir.exists(): shutil.rmtree(mh_dir)
        mh_dir.parent.mkdir(parents=True, exist_ok=True)
        log_info("  [MEGAHIT] assembling...")
        try:
            run_silent([
                "megahit",
                "-1", str(r1), "-2", str(r2),
                "-o", str(mh_dir),
                "--k-list",         "21,29,39,59,79,99,119,141",
                "--min-count",      "2",
                "--no-mercy",
                "--min-contig-len", str(min_len),
                "-t",               str(cfg.threads),
            ], log_file=rpt_dir / f"{sample_id}_megahit.log")
        except Exception as e:
            log_warn(f"  [MEGAHIT] warning: {e}")

    if mh_fa.exists() and mh_fa.stat().st_size > 0:
        s = fasta_stats(mh_fa)
        log_ok(f"  [MEGAHIT] {s['n']} contigs | "
               f"largest={s['largest_bp']}bp | N50={s['n50']}bp | "
               f"total={s['total_bp']}bp")
    else:
        log_warn("  [MEGAHIT] no output")

    # ── cd-hit-est 100% dereplication ─────────────────────────────────────────
    nr_fa  = out_dir / "combined" / f"{sample_id}_contigs_nr.fasta"
    tmp_fa = out_dir / "combined" / f"{sample_id}_all_tmp.fasta"

    if nr_fa.exists() and nr_fa.stat().st_size > 0 and not force:
        log_info("  [NR] already exists — skipping")
    else:
        n_sp = n_mh = 0
        with open(tmp_fa, "w") as out:
            # Use filtered SPAdes output
            if sp_fa.exists():
                with open(sp_fa) as f:
                    for line in f:
                        if line.startswith(">"):
                            out.write(f">{sample_id}_sp_{line[1:]}")
                            n_sp += 1
                        else:
                            out.write(line)
            if mh_fa.exists():
                with open(mh_fa) as f:
                    for line in f:
                        if line.startswith(">"):
                            out.write(f">{sample_id}_mh_{line[1:]}")
                            n_mh += 1
                        else:
                            out.write(line)

        n_input = n_sp + n_mh
        if n_input > 0:
            log_info(f"  [cd-hit-est] {n_input} contigs → dereplicating at 100% identity...")
            run_silent([
                "cd-hit-est",
                "-i", str(tmp_fa), "-o", str(nr_fa),
                "-c", "1.00", "-aS", "1.00",
                "-n", "10", "-d", "0",
                "-T", str(cfg.threads),
                "-M", str(cfg.memory_gb * 1000),
            ], log_file=rpt_dir / f"{sample_id}_cdhit.log")
            tmp_fa.unlink(missing_ok=True)
            Path(str(nr_fa)+".clstr").unlink(missing_ok=True)

            n_nr   = fasta_stats(nr_fa)["n"]
            redund = n_input - n_nr
            pct    = f"{redund/n_input*100:.1f}%"
            s      = fasta_stats(nr_fa)
            log_ok(f"  [NR] {n_input} → {n_nr} contigs "
                   f"| −{redund} redundant ({pct}) "
                   f"| largest={s['largest_bp']}bp")
        else:
            log_warn("  [NR] no contigs from either assembler")

    _save_tsv(sample_id, sp_fa, mh_fa, nr_fa,
              rpt_dir / "assembly_summary.tsv")

    log_step(f"Module 03 completed ✓  [{sample_id}]")
    log_info("Next: phageflow viral-id  (activate genomad env first)")
    return nr_fa


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _filter_fasta(src: Path, dst: Path, min_len: int) -> tuple[int, int]:
    """
    Write sequences from src to dst keeping only those >= min_len bp.

    Returns (n_total, n_kept).
    """
    n_total = n_kept = 0
    current_header = ""
    current_seq: List[str] = []

    def _flush(out):
        nonlocal n_kept
        seq = "".join(current_seq)
        if len(seq) >= min_len:
            out.write(current_header + "\n")
            out.write(seq + "\n")
            n_kept += 1

    with open(src) as f_in, open(dst, "w") as f_out:
        for line in f_in:
            line = line.rstrip()
            if line.startswith(">"):
                if current_header:
                    _flush(f_out)
                current_header = line
                current_seq = []
                n_total += 1
            else:
                current_seq.append(line)
        if current_header:
            _flush(f_out)

    return n_total, n_kept


def _save_tsv(sample_id, sp_fa, mh_fa, nr_fa, path):
    def _s(fa):
        if fa and fa.exists():
            s = fasta_stats(fa)
            return str(s["n"]), str(s["largest_bp"]), str(s["n50"]), str(s["total_bp"])
        return "0","0","0","0"

    headers = ["sample","assembler","contigs","largest_bp","n50_bp","total_bp"]
    rows = {}
    if path.exists():
        with open(path) as f:
            f.readline()
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if len(cols) >= 2: rows[f"{cols[0]}_{cols[1]}"] = cols

    for asm, fa in [("spades",sp_fa),("megahit",mh_fa),("NR(cdhit)",nr_fa)]:
        n,lg,n50,tot = _s(fa)
        key = f"{sample_id}_{asm}"
        rows[key] = [sample_id, asm, n, lg, n50, tot]

    with open(path,"w") as f:
        f.write("\t".join(headers)+"\n")
        for row in rows.values():
            f.write("\t".join(row)+"\n")

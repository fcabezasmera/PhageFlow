"""PhageFlow Module 02 — Host read removal.

Goal: remove host-derived reads from paired-end Illumina reads, retaining
reads of non-host origin for downstream assembly. No assumption is made
about the composition of the non-host fraction at this stage.

Two modes are available:

bwa-mem2 (RECOMMENDED)
  Maps reads against one or more host reference genome(s).
  Reads that do NOT map are retained as candidate non-host reads.
  Singletons (one mate maps, one does not) are retained separately
  for use as --s1 input to SPAdes, which is critical for recovering
  terminal repeats (DTR/ITR) of circular viral genomes.

  Alignment uses permissive parameters (-A 1 -B 2 -O 2,2) to maximise
  sensitivity for divergent host strains. Using default BWA parameters
  risks retaining divergent host reads as false non-host signal.
  Schmieder & Edwards 2011 (Bioinformatics 27:863) — DeconSeq;
  Li 2013 (arXiv:1303.3997) — BWA-MEM parameter rationale.

  Level A post-check: subsamples 50,000 reads and screens with Kraken2
  to detect residual bacterial contamination. Diagnostic only — no reads
  are removed. Warns if bacterial fraction exceeds contamination_warn_pct.

Kraken2 (alternative)
  Taxonomic classification against a k-mer database.
  Unclassified reads + viral-classified reads are retained.
  Appropriate when no host genome is available.

  Parameters: --confidence 0.5 --minimum-hit-groups 3
  Strict confidence prevents false-positive host classification of
  divergent viral reads with partial k-mer matches to bacterial sequences.
  Wood et al. 2019 (Genome Biol 20:257) — Kraken2.

  Optional Level B post-filter: bwa-mem2 alignment against bacterial
  species detected at >postfilter_min_pct by Kraken2. Recovers singletons
  lost by Kraken2 paired-end classification of terminal repeat reads.

Singleton recovery rationale
  In paired-end sequencing of circular genomes with DTR/ITR, read pairs
  that span the terminal repeat junction may have one mate map to host
  and one mate fail to map (or map to the virus). Retaining singletons
  ensures these reads contribute to assembly of genome termini.
  Antipov et al. 2016 (Bioinformatics 32:i60) — SPAdes plasmidSPAdes.
"""
from __future__ import annotations
import subprocess
import sys
import zipfile
import shutil
from pathlib import Path
from typing import List, Optional

from rich.panel import Panel
from rich.progress import (
    Progress, SpinnerColumn, TextColumn,
    BarColumn, MofNCompleteColumn, TimeElapsedColumn,
)

from phageflow.utils.config import Config
from phageflow.utils.logger import (
    log_step, log_info, log_ok, log_warn, log_error, console,
)
from phageflow.utils.tools import require_tools, run_silent, human_size, mkdirs

STEP          = "02_host_removal"
TOOLS_BWA     = ["bwa-mem2", "samtools"]
TOOLS_K2      = ["kraken2", "seqtk"]
FASTA_EXTS    = {".fasta", ".fa", ".fna", ".fas"}
TAXID_VIRUSES = 10239

# Non-host retention thresholds
# These are agnostic to sample composition — the non-host fraction
# may be viral, plasmid, or other non-host sequence.
_NONHOST_GOOD     = 50.0   # green: >50% retained
_NONHOST_WARN     = 20.0   # yellow: <20% may indicate reference mismatch
_NONHOST_CRITICAL =  5.0   # red: <5% almost certainly a reference error

# Kraken2: warn if >80% unclassified (DB may lack the host genome)
_K2_UNCLASS_WARN  = 0.80
_MIN_READS       = 50_000
_SAMTOOLS_MIN_VERSION = (1, 15)
_LEVEL_A_SUBSAMPLE    = 50_000

def run(
    cfg:             Config,
    sample_id:       str,
    r1:              Path,
    r2:              Path,
    force:           bool                = False,
    host_file:       Optional[Path]      = None,
    accessions:      Optional[List[str]] = None,
    accessions_file: Optional[Path]      = None,
    kraken_db:       Optional[Path]      = None,
) -> tuple:
    r1 = Path(r1); r2 = Path(r2)
    out_dir  = cfg.results(STEP)
    rpt_dir  = cfg.reports(STEP) / sample_id
    host_dir = rpt_dir / "host_genomes"
    tmp_dir  = out_dir / "tmp" / sample_id

    use_bwa = bool(host_file or accessions or accessions_file)
    if use_bwa:
        mkdirs(out_dir, rpt_dir, host_dir, tmp_dir)
    else:
        mkdirs(out_dir, rpt_dir, tmp_dir)

    r1_out        = out_dir / f"{sample_id}_R1.fastq.gz"
    r2_out        = out_dir / f"{sample_id}_R2.fastq.gz"
    singleton_out = out_dir / f"{sample_id}_singletons.fastq.gz"
    log_f         = rpt_dir / f"{sample_id}_host_removal.log"

    log_step(f"Module 02 — host removal  [{sample_id}]")
    log_info(f"  {sample_id}  ·  R1: {human_size(r1)}  ·  R2: {human_size(r2)}")

    if r1_out.exists() and r1_out.stat().st_size > 0 and not force:
        log_info("  Already processed — skipping  (--force to re-run)")
        return r1_out, r2_out, singleton_out

    active_warnings: list[str] = []
    always = list(cfg.host_removal.always_include_accessions or [])

    if use_bwa:
        require_tools(*TOOLS_BWA)
        _check_samtools_version()
        combined = host_dir / "combined_hosts.fasta"
        _host_label = ", ".join(accessions) if accessions else (
            host_file.name if host_file else "config defaults"
        )
        log_info(f"  mode: bwa-mem2  ·  host: {_host_label}")
        n_steps = 3
    else:
        require_tools(*TOOLS_K2)
        db = kraken_db or cfg.databases.kraken2
        if not db or not Path(db).exists():
            log_error(
                "  No host reference provided and no Kraken2 database found.\n"
                "  Provide at least one of:\n"
                "    --accessions GCF_XXXXXXXX.X\n"
                "    --host-file /path/to/host.fasta\n"
                "    --kraken-db /path/to/k2_db\n"
                "  Or set databases.kraken2 in config.yaml"
            )
            sys.exit(1)
        log_info(f"  mode: Kraken2  ·  db: {Path(str(db)).name}")
        if cfg.host_removal.kraken2_postfilter:
            log_info("  Kraken2 postfilter: bwa-mem2 active (recovers DTR/ITR singletons)")
        n_steps = 3 if cfg.host_removal.kraken2_postfilter else 2

    stats: dict = {}

    with Progress(
        SpinnerColumn(spinner_name="dots2"),
        TextColumn("[bold cyan]{task.description:<60}"),
        BarColumn(bar_width=20, complete_style="cyan", finished_style="green"),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        console=console,
        transient=True,
    ) as progress:
        task = progress.add_task("Initializing...", total=n_steps)

        if use_bwa:
            progress.update(task, description="[1/3] bwa-mem2 — reference + index")
            idx = Path(str(combined) + ".bwt.2bit.64")
            if idx.exists() and not force:
                log_info("  [bwa-mem2] index cached — skipping  (--force to refresh)")
            else:
                fastas = _resolve_fastas(host_file, accessions, accessions_file,
                                         host_dir, rpt_dir, always)
                if not fastas:
                    log_error("  No host FASTAs resolved.")
                    raise FileNotFoundError("No host reference FASTAs found.")
                _build_index(fastas, combined, rpt_dir, force)
                _cleanup_ncbi_download(host_dir)
            progress.advance(task)
            log_ok("  [1/3] bwa-mem2 — index ready")

            progress.update(task, description="[2/3] bwa-mem2 — streaming alignment")
            stats = _run_bwa_pipeline(
                sample_id, r1, r2, combined, r1_out, r2_out, singleton_out,
                tmp_dir, log_f, cfg.threads,
            )
            progress.advance(task)
            log_ok("  [2/3] bwa-mem2 — complete")

            progress.update(task, description="[3/3] Level A  — contamination check")
            _level_a_check(r1_out, r2_out, rpt_dir, sample_id, cfg, active_warnings)
            progress.advance(task)
            log_ok("  [3/3] Level A  — complete")

        else:
            progress.update(task, description="[1/2] Kraken2  — classifying reads")
            report, k2_out = _classify_kraken2(
                sample_id, r1, r2, Path(db), tmp_dir, rpt_dir, log_f, cfg.threads
            )
            progress.advance(task)
            log_ok("  [1/2] Kraken2  — complete")

            if cfg.host_removal.kraken2_postfilter:
                r1_k2 = tmp_dir / f"{sample_id}_k2_R1.fastq.gz"
                r2_k2 = tmp_dir / f"{sample_id}_k2_R2.fastq.gz"
                progress.update(task, description="[2/3] seqtk    — extracting Kraken2 phage reads")
                stats_k2 = _extract_phage_kraken2(
                    sample_id, r1, r2, r1_k2, r2_k2, report, k2_out, tmp_dir, log_f,
                )
                progress.advance(task)
                progress.update(task, description="[3/3] Level B  — bwa-mem2 post-filter + singletons")
                stats_bwa = _level_b_postfilter(
                    sample_id, r1_k2, r2_k2, r1_out, r2_out, singleton_out,
                    report, tmp_dir, host_dir, rpt_dir, log_f, cfg, force,
                )
                progress.advance(task)
                stats = {**stats_k2,
                         "reads_phage":     stats_bwa["reads_phage"],
                         "reads_singleton": stats_bwa["reads_singleton"],
                         "pct_phage":       stats_bwa["pct_phage"],
                         "pct_host":        stats_bwa["pct_host"],
                         "pct_singleton":   stats_bwa["pct_singleton"]}
            else:
                progress.update(task, description="[2/2] seqtk    — extracting phage reads")
                stats = _extract_phage_kraken2(
                    sample_id, r1, r2, r1_out, r2_out, report, k2_out, tmp_dir, log_f,
                )
                progress.advance(task)
                log_ok("  [2/2] seqtk    — complete")

    shutil.rmtree(tmp_dir, ignore_errors=True)
    try:
        tmp_dir.parent.rmdir()
    except OSError:
        pass

    mode = "bwa-mem2" if use_bwa else "kraken2"
    _validate_output(r1_out, r2_out, singleton_out, stats, active_warnings)
    _save_tsv({**stats, "sample": sample_id, "mode": mode},
              rpt_dir / "host_removal_summary.tsv")
    _check_warnings(stats, mode, active_warnings)
    _print_completion_panel(sample_id, r1_out, r2_out, singleton_out,
                            rpt_dir, stats, mode, active_warnings)
    log_step(f"Module 02 completed ✓  [{sample_id}]")
    return r1_out, r2_out, singleton_out

# ── samtools version check ────────────────────────────────────────────────────

def _check_samtools_version() -> None:
    try:
        r = subprocess.run(["samtools", "--version"],
                           capture_output=True, text=True, timeout=10)
        parts = r.stdout.strip().split("\n")[0].split()
        if len(parts) >= 2:
            ver_str = parts[1].split("-")[0]
            major, minor = (int(x) for x in ver_str.split(".")[:2])
            if (major, minor) < _SAMTOOLS_MIN_VERSION:
                log_warn(
                    f"  CRITICAL: samtools {ver_str} < 1.15. "
                    "Fix: mamba install -n phageflow 'samtools>=1.15'"
                )
    except Exception:
        log_warn("  Could not verify samtools version — ensure ≥1.15 is active.")

# ── Reference resolution ──────────────────────────────────────────────────────

def _resolve_fastas(
    host_file, accessions, accessions_file, host_dir, rpt_dir, always_include,
) -> List[Path]:
    fastas: List[Path] = []
    if host_file:
        p = Path(host_file)
        if p.is_dir():
            fastas = sorted(f for f in p.iterdir() if f.suffix.lower() in FASTA_EXTS)
            log_info(f"  Reference : folder  ({len(fastas)} FASTAs)")
        elif p.suffix.lower() in FASTA_EXTS:
            fastas = [p]
            log_info(f"  Reference : {p.name}  ({human_size(p)})")
        else:
            fastas = [Path(line.strip()) for line in open(p)
                      if line.strip() and not line.startswith("#")
                      and Path(line.strip()).suffix.lower() in FASTA_EXTS
                      and Path(line.strip()).exists()]
            log_info(f"  Reference : FASTA list  ({len(fastas)} files)")

    user_accs: List[str] = []
    if accessions:
        user_accs = accessions
    elif accessions_file and Path(accessions_file).exists():
        user_accs = [l.strip() for l in open(accessions_file)
                     if l.strip() and not l.startswith("#")]

    all_accs = list(dict.fromkeys(user_accs + always_include))
    if all_accs:
        log_info(f"  Accessions: {len(user_accs)} user + {len(always_include)} config "
                 f"= {len(all_accs)} total")
        fastas += _download_accessions(all_accs, host_dir, rpt_dir)
    return fastas

def _download_accessions(accs: List[str], host_dir: Path, rpt_dir: Path) -> List[Path]:
    require_tools("datasets")
    dl_dir   = host_dir / "ncbi_dataset"
    zip_out  = host_dir / "ncbi_download.zip"
    acc_file = host_dir / "accessions.txt"
    acc_file.write_text("\n".join(accs) + "\n")
    log_info(f"  Downloading : {', '.join(accs)}")
    try:
        run_silent(
            ["datasets", "download", "genome", "accession",
             "--inputfile", str(acc_file), "--include", "genome",
             "--filename", str(zip_out), "--no-progressbar"],
            log_file=rpt_dir / "ncbi_download.log",
        )
    except Exception as e:
        log_warn(f"  datasets warning: {e}")
    fastas: List[Path] = []
    if zip_out.exists():
        with zipfile.ZipFile(zip_out) as z:
            z.extractall(dl_dir)
        zip_out.unlink()
        fastas = sorted(dl_dir.rglob("*.fna"))
        log_ok(f"  datasets · downloaded {len(fastas)} genome(s)")
    return fastas

def _build_index(fastas: List[Path], combined: Path, rpt_dir: Path, force: bool) -> None:
    idx = Path(str(combined) + ".bwt.2bit.64")
    if idx.exists() and not force:
        log_info("  [bwa-mem2] index already exists — skipping")
        return
    log_info(f"  [bwa-mem2] building combined reference ({len(fastas)} genome(s))")
    with open(combined, "w") as out:
        for fna in fastas:
            n = 0
            with open(fna) as f:
                for line in f:
                    out.write(f">{fna.stem}_{line[1:]}" if line.startswith(">") else line)
                    n += line.startswith(">")
            log_info(f"    + {fna.name}  ({n} seq(s), {human_size(fna)})")
    run_silent(["bwa-mem2", "index", str(combined)], log_file=rpt_dir / "bwa_index.log")
    run_silent(["samtools", "faidx", str(combined)], log_file=rpt_dir / "bwa_index.log")
    log_ok("  [bwa-mem2] index ready")

def _cleanup_ncbi_download(host_dir: Path) -> None:
    """Remove raw NCBI download folder after combined FASTA is built."""
    ncbi_dir = host_dir / "ncbi_dataset"
    if ncbi_dir.exists():
        shutil.rmtree(ncbi_dir, ignore_errors=True)
        log_info("  [cleanup] ncbi_dataset/ removed — combined FASTA retained")

# ── bwa-mem2 streaming pipeline ───────────────────────────────────────────────

def _run_bwa_pipeline(
    sample_id, r1, r2, combined, r1_out, r2_out, singleton_out,
    tmp_dir, log_f, threads,
) -> dict:
    sort_prefix = tmp_dir / f"{sample_id}.sort"

    cmd = (
        # Permissive alignment parameters for host removal sensitivity:
    # -A 1 -B 2 -O 2,2: lower mismatch/gap penalties capture divergent host strains
    # that default parameters would fail to map, preventing false non-host retention.
    # Schmieder & Edwards 2011 (Bioinformatics 27:863)
    f"bwa-mem2 mem -t {threads} -A 1 -B 2 -O 2,2 -M {combined} {r1} {r2} "
        f"| samtools view -@ {min(threads,8)} -bF 2304 -u "
        f"| samtools sort -@ {min(threads,8)} -n -u -T {sort_prefix} - "
        f"| samtools fastq -@ {min(threads,4)} -N "
        f"  -f 4 -F 256 "
        f"  -1 {r1_out} -2 {r2_out} "
        f"  -s {singleton_out} "
        f"  -0 /dev/null"
    )
    run_silent(cmd, log_file=log_f, shell=True)

    n_in  = _count_reads_fastq(r1,           threads)
    n_out = _count_reads_fastq(r1_out,        threads)
    n_sin = _count_reads_fastq(singleton_out, threads)
    total = n_in * 2
    phage = n_out * 2 + n_sin
    host  = max(0, total - phage)

    pct_p = f"{phage / total * 100:.2f}%" if total else "N/A"
    pct_h = f"{host  / total * 100:.2f}%" if total else "N/A"
    pct_s = f"{n_sin / total * 100:.2f}%" if total else "N/A"

    if n_sin:
        log_ok(f"  [singletons] {n_sin:,} retained ({pct_s}) — pass to assembly --s1")

    return {
        "reads_in": str(total), "reads_phage": str(phage), "reads_singleton": str(n_sin),
        "pct_host": pct_h, "pct_phage": pct_p, "pct_singleton": pct_s,
    }

def _count_reads_fastq(path: Path, threads: int = 4) -> int:
    if not path.exists() or path.stat().st_size == 0:
        return 0
    try:
        decomp = f"pigz -dc -p {threads}" if shutil.which("pigz") else "gzip -dc"
        r = subprocess.run(["bash", "-c", f"{decomp} {path} | wc -l"],
                           capture_output=True, text=True, timeout=600)
        return int(r.stdout.strip()) // 4
    except Exception:
        return 0

# ── Level A: post-extraction contamination check ──────────────────────────────

def _level_a_check(
    r1_phage, r2_phage, rpt_dir, sample_id, cfg, active_warnings,
) -> None:
    """Kraken2 contamination screen on a subsample of extracted phage reads.

    Diagnostic only — no reads removed.
    Only runs when databases.kraken2 is configured.
    """
    k2_db = cfg.databases.kraken2
    if not k2_db or not k2_db.exists():
        log_info("  [Level A] Kraken2 DB not configured — contamination check skipped")
        return
    if not r1_phage.exists() or r1_phage.stat().st_size == 0:
        return

    warn_pct = cfg.host_removal.contamination_warn_pct
    sub_r1   = rpt_dir / f"{sample_id}_la_sub_R1.fq"
    sub_r2   = rpt_dir / f"{sample_id}_la_sub_R2.fq"
    k2_rep   = rpt_dir / f"{sample_id}_levelA_contamination.report"

    for src, dst in [(r1_phage, sub_r1), (r2_phage, sub_r2)]:
        try:
            subprocess.run(
                ["bash", "-c", f"seqtk sample -s42 {src} {_LEVEL_A_SUBSAMPLE} > {dst}"],
                capture_output=True, timeout=120,
            )
        except Exception:
            pass

    if not sub_r1.exists() or sub_r1.stat().st_size == 0:
        for p in (sub_r1, sub_r2): p.unlink(missing_ok=True)
        return

    try:
        subprocess.run(
            ["kraken2", "--db", str(k2_db), "--report", str(k2_rep),
             "--output", "/dev/null", "--paired",
             "--confidence", "0.3", "--minimum-hit-groups", "2",
             "--threads", str(min(cfg.threads, 8)),
             str(sub_r1), str(sub_r2)],
            capture_output=True, timeout=300,
        )
    except Exception as e:
        log_info(f"  [Level A] Kraken2 check failed: {e}")
    finally:
        for p in (sub_r1, sub_r2): p.unlink(missing_ok=True)

    if not k2_rep.exists():
        return

    bacteria_pct = 0.0
    top_contaminants: list[tuple[str, float]] = []
    in_bacteria = False
    bact_depth  = 0

    with open(k2_rep) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            try:
                pct      = float(parts[0])
                name_raw = parts[5]
                depth    = len(name_raw) - len(name_raw.lstrip())
                name     = name_raw.strip()
            except (ValueError, IndexError):
                continue
            if name == "Bacteria":
                bacteria_pct = pct; in_bacteria = True; bact_depth = depth
                continue
            if in_bacteria:
                if depth <= bact_depth:
                    in_bacteria = False
                elif pct >= warn_pct and depth == bact_depth + 2:
                    top_contaminants.append((name, pct))

    if bacteria_pct > warn_pct:
        top_str = (" · ".join(f"{n} ({p:.1f}%)" for n, p in
                   sorted(top_contaminants, key=lambda x: -x[1])[:3])
                   or "none at genus level")
        msg = (
            f"Level A: {bacteria_pct:.1f}% bacterial in phage reads "
            f"(threshold {warn_pct}%). Top: {top_str}. "
            "Consider re-running with --accessions including the contaminant."
        )
        log_warn(f"  {msg}")
        active_warnings.append(msg)
    else:
        log_ok(f"  [Level A] {bacteria_pct:.1f}% bacterial (threshold: {warn_pct}%) — OK")

# ── Level B: adaptive post-filter ────────────────────────────────────────────

def _level_b_postfilter(
    sample_id, r1_k2, r2_k2, r1_out, r2_out, singleton_out,
    k2_report, tmp_dir, host_dir, rpt_dir, log_f, cfg, force,
) -> dict:
    """bwa-mem2 post-filter after Kraken2 extraction (opt-in via kraken2_postfilter)."""
    min_pct = cfg.host_removal.postfilter_min_pct
    taxids  = _parse_contaminants_from_report(k2_report, min_pct)

    if not taxids:
        log_info(f"  [Level B] No bacterial species > {min_pct}% — using Kraken2 output directly")
        for src, dst in [(r1_k2, r1_out), (r2_k2, r2_out)]:
            if src.exists(): src.rename(dst)
        singleton_out.write_bytes(b"")
        n = _count_reads_fastq(r1_out, cfg.threads)
        return {"reads_phage": str(n * 2), "reads_singleton": "0",
                "pct_phage": "N/A", "pct_host": "0.00%", "pct_singleton": "0.00%"}

    pf_dir = host_dir / "postfilter"
    mkdirs(pf_dir)
    fastas = _download_by_taxids(taxids, pf_dir, rpt_dir)
    if not fastas:
        log_warn("  [Level B] No genomes downloaded — skipping bwa-mem2 post-filter")
        for src, dst in [(r1_k2, r1_out), (r2_k2, r2_out)]:
            if src.exists(): src.rename(dst)
        return {"reads_phage": str(_count_reads_fastq(r1_out, cfg.threads) * 2),
                "reads_singleton": "0", "pct_phage": "N/A",
                "pct_host": "N/A", "pct_singleton": "N/A"}

    combined_pf = pf_dir / "combined_postfilter.fasta"
    _build_index(fastas, combined_pf, rpt_dir, force)
    _cleanup_ncbi_download(pf_dir)
    log_info(f"  [Level B] bwa-mem2 post-filter vs {len(fastas)} genome(s) "
             f"({len(taxids)} bacterial taxon(s) > {min_pct}%)")

    stats = _run_bwa_pipeline(
        f"{sample_id}_pf", r1_k2, r2_k2, combined_pf,
        r1_out, r2_out, singleton_out, tmp_dir, log_f, cfg.threads,
    )
    for p in (r1_k2, r2_k2): p.unlink(missing_ok=True)

    n_out = _count_reads_fastq(r1_out, cfg.threads)
    n_sin = _count_reads_fastq(singleton_out, cfg.threads)
    if n_sin:
        log_ok(f"  [Level B] {n_sin:,} DTR/ITR singletons recovered from Kraken2 loss")
    return stats

def _parse_contaminants_from_report(report: Path, min_pct: float) -> list[int]:
    taxids: list[int] = []
    if not report.exists(): return taxids
    in_bacteria = False; bact_depth = 0
    with open(report) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 6: continue
            try:
                pct = float(parts[0]); taxid = int(parts[4].strip())
                name_raw = parts[5]; depth = len(name_raw) - len(name_raw.lstrip())
                name = name_raw.strip()
            except (ValueError, IndexError): continue
            if name == "Bacteria":
                in_bacteria = True; bact_depth = depth; continue
            if in_bacteria:
                if depth <= bact_depth: break
                if pct >= min_pct: taxids.append(taxid)
    return taxids

def _download_by_taxids(taxids: list[int], dest_dir: Path, rpt_dir: Path) -> List[Path]:
    require_tools("datasets")
    fastas: List[Path] = []
    for taxid in taxids:
        zip_out = dest_dir / f"taxid_{taxid}.zip"
        dl_dir  = dest_dir / f"taxid_{taxid}"
        try:
            run_silent(
                ["datasets", "download", "genome", "taxon", str(taxid),
                 "--reference", "--include", "genome",
                 "--filename", str(zip_out), "--no-progressbar"],
                log_file=rpt_dir / "ncbi_download.log", check=False,
            )
            if zip_out.exists():
                with zipfile.ZipFile(zip_out) as z: z.extractall(dl_dir)
                zip_out.unlink()
                found = sorted(dl_dir.rglob("*.fna"))
                if found:
                    fastas.extend(found)
                    log_info(f"  [Level B] taxid {taxid}: {found[0].name}")
                else:
                    log_warn(f"  [Level B] taxid {taxid}: no genome in download")
        except Exception as e:
            log_warn(f"  [Level B] taxid {taxid}: {e}")
    return fastas

# ── Kraken2 pipeline ──────────────────────────────────────────────────────────

def _classify_kraken2(sample_id, r1, r2, db, tmp_dir, rpt_dir, log_f, threads) -> tuple:
    report = rpt_dir / f"{sample_id}_k2.report"
    k2_out = tmp_dir  / f"{sample_id}_k2.tsv"
    run_silent([
        "kraken2", "--db", str(db), "--report", str(report),
        "--output", str(k2_out),
        "--unclassified-out", f"{tmp_dir}/{sample_id}_unclassified#.fastq",
        "--paired", "--confidence", "0.5", "--minimum-hit-groups", "3",
        "--threads", str(threads), str(r1), str(r2),
    ], log_file=log_f, check=False)
    return report, k2_out

def _extract_phage_kraken2(sample_id, r1, r2, r1_out, r2_out,
                            report, k2_out, tmp_dir, log_f) -> dict:
    un_r1 = tmp_dir / f"{sample_id}_unclassified_1.fastq"
    un_r2 = tmp_dir / f"{sample_id}_unclassified_2.fastq"
    viral_ids = tmp_dir / f"{sample_id}_viral_ids.txt"

    viral_taxids  = _parse_viral_taxids(report)
    n_viral_reads = _extract_viral_ids(k2_out, viral_taxids, viral_ids)
    log_info(f"  Kraken2 : {len(viral_taxids)} viral taxid(s) | {n_viral_reads:,} viral reads")
    _merge_reads(un_r1, un_r2, r1, r2, r1_out, r2_out, viral_ids, n_viral_reads, log_f)

    n_total = n_unclass = n_viral = 0
    if report.exists():
        with open(report) as f:
            for line in f:
                p = line.strip().split("\t")
                if len(p) < 6: continue
                try:
                    count = int(p[1]); name = p[5].strip()
                    if name == "unclassified": n_unclass = count
                    if name == "root":         n_total   = count + n_unclass
                    if name == "Viruses":      n_viral   = count
                except ValueError: pass

    for p in (k2_out, viral_ids, un_r1, un_r2):
        if p and Path(p).exists(): Path(p).unlink()

    n_kept = n_unclass + n_viral
    return {
        "reads_in": str(n_total), "reads_phage": str(n_kept), "reads_singleton": "0",
        "pct_host":  f"{(n_total - n_kept) / n_total * 100:.2f}%" if n_total else "N/A",
        "pct_phage": f"{n_kept / n_total * 100:.2f}%" if n_total else "N/A",
        "pct_singleton": "N/A", "k2_unclass": str(n_unclass), "k2_viral": str(n_viral),
    }

def _parse_viral_taxids(report: Path) -> set:
    taxids: set = set()
    if not report.exists(): return taxids
    in_viral = False; viral_depth: Optional[int] = None
    with open(report) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6: continue
            try:
                taxid = int(p[4].strip()); name_raw = p[5]
                depth = len(name_raw) - len(name_raw.lstrip()); name = name_raw.strip()
            except (ValueError, IndexError): continue
            if name == "Viruses" or taxid == TAXID_VIRUSES:
                in_viral = True; viral_depth = depth; taxids.add(taxid); continue
            if in_viral:
                if depth > viral_depth: taxids.add(taxid)
                else: in_viral = False; viral_depth = None
    return taxids

def _extract_viral_ids(k2_out: Path, viral_taxids: set, out_file: Path) -> int:
    n = 0
    if not k2_out.exists() or not viral_taxids: return n
    viral_taxids_str = {str(t) for t in viral_taxids}
    with (open(k2_out, "r", buffering=1 << 20) as f_in,
          open(out_file, "w", buffering=1 << 20) as f_out):
        for line in f_in:
            if not line.startswith("C\t"): continue
            parts = line.split("\t")
            if len(parts) >= 3 and parts[2].strip() in viral_taxids_str:
                f_out.write(parts[1] + "\n"); n += 1
    return n

def _merge_reads(un_r1, un_r2, src_r1, src_r2, r1_out, r2_out,
                 viral_ids, n_viral, log_f) -> None:
    has_pigz = bool(shutil.which("pigz"))
    zip_c    = "pigz -c" if has_pigz else "gzip -c"
    zip_pipe = "pigz"    if has_pigz else "gzip"
    has_viral = viral_ids.exists() and n_viral > 0
    cmds: list[str] = []
    for un, src, out in [(un_r1, src_r1, r1_out), (un_r2, src_r2, r2_out)]:
        if un.exists() and has_viral:
            cmds.append(f"{{ {zip_c} {un}; seqtk subseq {src} {viral_ids} | {zip_pipe}; }} > {out}")
        elif un.exists():
            cmds.append(f"{zip_c} {un} > {out}")
        elif has_viral:
            cmds.append(f"seqtk subseq {src} {viral_ids} | {zip_pipe} > {out}")
    if cmds:
        run_silent(" & ".join(cmds) + " & wait", log_file=log_f, shell=True, check=False)

# ── Validation & warnings ─────────────────────────────────────────────────────

def _validate_output(r1_out, r2_out, singleton_out, stats, active_warnings) -> None:
    for label, path in [("R1", r1_out), ("R2", r2_out)]:
        if not path.exists() or path.stat().st_size == 0:
            msg = f"{label} output missing or empty: {path}"
            log_warn(f"  [validate] {msg}"); active_warnings.append(msg)
    try:
        n = int(stats.get("reads_phage", 0))
        if n < _MIN_READS:
            msg = (f"Only {n:,} phage reads retained — < {_MIN_READS:,} "
                   "may be insufficient for assembly of large viral genomes (>100 kb at 50×).")
            log_warn(f"  [validate] {msg}"); active_warnings.append(msg)
        else:
            log_ok(f"  [validate] {n:,} phage reads retained — OK")
    except (ValueError, TypeError): pass

def _check_warnings(stats, mode, active_warnings) -> None:
    try:
        ph = float(stats.get("pct_phage", "0").rstrip("%"))
        if ph < _NONHOST_CRITICAL:
            msg = (
                f"CRITICAL: only {stats.get('pct_phage')} of reads passed host removal "
                f"(threshold {_NONHOST_CRITICAL}%). Possible reference genome mismatch "
                "or near-complete host contamination. Verify --accessions or --host-file."
            )
            log_warn(f"  {msg}"); active_warnings.append(msg)
        elif ph < _NONHOST_WARN:
            msg = (
                f"Low non-host fraction: {stats.get('pct_phage')} "
                f"(threshold {_NONHOST_WARN}%). "
                "Check that the host reference matches the propagation strain."
            )
            log_warn(f"  {msg}"); active_warnings.append(msg)
    except (ValueError, AttributeError): pass
    if mode == "kraken2":
        try:
            n_t = int(stats.get("reads_in", 0)); n_u = int(stats.get("k2_unclass", 0))
            if n_t and n_u / n_t > _K2_UNCLASS_WARN:
                msg = (
                    f"> {_K2_UNCLASS_WARN * 100:.0f}% reads unclassified by Kraken2 — "
                    "host genome may not be in the database. "
                    "Consider bwa-mem2 mode with --accessions for accurate host removal."
                )
                log_warn(f"  [Kraken2] {msg}"); active_warnings.append(msg)
        except (ValueError, TypeError): pass

# ── Completion panel ──────────────────────────────────────────────────────────

def _print_completion_panel(
    sample_id, r1_out, r2_out, singleton_out, rpt_dir, stats, mode, active_warnings,
) -> None:
    def _c(val, good, warn):
        try:
            v = float(str(val).rstrip("%"))
            if v >= good: return f"[bold green]{val}[/bold green]"
            if v >= warn: return f"[bold yellow]{val}[/bold yellow]"
            return f"[bold red]{val}[/bold red]"
        except: return str(val)

    def _f(val):
        s = str(val)
        if s in ("N/A", "", "None"):
            return "[dim]N/A[/dim]"
        try:
            return f"{int(s):,}"
        except ValueError:
            return s

    n_in  = _f(stats.get("reads_in",        "N/A"))
    n_out = _f(stats.get("reads_phage",      "N/A"))
    n_sin = _f(stats.get("reads_singleton",  "0"))
    ph    = stats.get("pct_phage", "N/A")
    ho    = stats.get("pct_host",  "N/A")
    unit  = "read pairs" if mode == "kraken2" else "reads"

    lines = [
        "[bold]Key metrics[/bold]",
        f"  [cyan]{unit.capitalize():12s}:[/cyan] {n_in} → {n_out}"
        f"  phage={_c(ph, _NONHOST_GOOD, _NONHOST_WARN)}  host={_f(ho)}",
        f"  [cyan]Singletons  :[/cyan] {n_sin} retained  [dim](pass to assembly --s1)[/dim]",
        f"  [cyan]Mode        :[/cyan] {mode}",
    ]
    if mode == "kraken2":
        lines.append(f"  [dim]Kraken2: unclassified={_f(stats.get('k2_unclass','N/A'))}"
                     f"  viral={_f(stats.get('k2_viral','N/A'))}[/dim]")
    lines += [
        "", "[bold]Output files[/bold]",
        f"  Phage R1   : {r1_out}",
        f"  Phage R2   : {r2_out}",
        f"  Singletons : {singleton_out}",
        f"  Summary    : {rpt_dir / 'host_removal_summary.tsv'}",
    ]
    if active_warnings:
        lines += ["", "[bold yellow]Warnings[/bold yellow]"]
        lines += [f"  [yellow]⚠ {w}[/yellow]" for w in active_warnings]

    console.print(Panel("\n".join(lines),
                        title=f"[bold cyan]Host removal complete — {sample_id}[/bold cyan]",
                        border_style="cyan", padding=(0, 2), width=120))

# ── TSV summary ───────────────────────────────────────────────────────────────

_TSV_HEADERS = [
    "sample", "mode", "reads_in", "reads_phage", "reads_singleton",
    "pct_host", "pct_phage", "pct_singleton", "k2_unclass", "k2_viral",
]

def _save_tsv(row: dict, path: Path) -> None:
    rows: dict[str, dict] = {}
    if path.exists():
        with open(path) as f:
            old_hdrs = f.readline().rstrip("\n").split("\t")
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if cols and cols[0]: rows[cols[0]] = dict(zip(old_hdrs, cols))
    rows[row["sample"]] = {h: str(row.get(h, "")) for h in _TSV_HEADERS}
    with open(path, "w") as f:
        f.write("\t".join(_TSV_HEADERS) + "\n")
        for r in rows.values():
            f.write("\t".join(r.get(h, "") for h in _TSV_HEADERS) + "\n")

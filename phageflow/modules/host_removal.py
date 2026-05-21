"""PhageFlow Module 02 — Host read removal.

Goal: maximise retention of phage reads — including reads at DTR/ITR boundaries —
while eliminating host contamination. Complete/HQ CheckV classification requires
uniform phage coverage across the entire genome, including terminal regions
(Nayfach et al. 2021, Nat Biotechnol 39:578).

Modes (priority order):

  1. --host-file  path
     Local FASTA (.fasta / .fna), folder of FASTAs, or text file of FASTA paths.
     Uses bwa-mem2 alignment — most precise when the host genome is available.

  2. --accessions GCF_1,GCF_2,...
     Comma-separated NCBI accessions (GCF / GCA). Downloads via datasets CLI,
     then aligns with bwa-mem2.

  3. --accessions-file hosts.txt
     Text file with one GCF/GCA accession per line.

  4. (default — no reference provided)
     Kraken2 filter: classify reads and retain unclassified + Viruses (taxid 10239).
     Parameters: --confidence 0.5 --minimum-hit-groups 3
     (Wood et al. 2019, Genome Biology 20:257)

bwa-mem2 pipeline notes
-----------------------
  * samtools sort -n REMOVED: samtools fastq with -N does not require
    name-sorted input when reads are properly paired in the BAM
    (Li et al. 2009, Bioinformatics 25:2078). Removing the sort eliminates
    a disk-bound O(n log n) bottleneck at no cost to output correctness.

  * Singletons are RETAINED: reads whose pair mapped to the host but which
    themselves are unmapped (-f 4 -F 8) are written to a separate FASTQ and
    passed to SPAdes as --s1. These represent phage fragments at DTR/ITR
    boundaries where the mate may span a repeat and align ambiguously to the
    host. Discarding them reduces terminal coverage and lowers the probability
    of CheckV 'Complete' classification (Nayfach et al. 2021).

samtools version requirement
----------------------------
  samtools >= 1.15 required for -N flag and --singleton in samtools fastq.
  samtools 0.1.x (legacy) does NOT implement these flags and will fail.
  Verify with: samtools --version
  Update with: mamba install -n phageflow "samtools>=1.15"

Thresholds
----------
  _PHAGE_GOOD     = 90% : expected in a well-prepared purified phage
                          preparation (CsCl gradient or 0.22 µm filtration)
                          with correct host reference (Roux et al. 2019, eLife).
  _PHAGE_WARN     = 50% : indicates a real problem — wrong reference, heavy
                          contamination, or incorrect sample. Below this value
                          the assembled phage genome is unlikely to be complete.
  _PHAGE_CRITICAL = 20% : assembly failure highly probable; downstream
                          Complete/HQ classification will be severely affected.

Kraken2 parameters
------------------
  --confidence 0.5       : higher than the common default of 0.2 to maximise
                           specificity (fewer phage reads misclassified as host).
                           In purified phage mode the priority is retaining phage
                           reads, not exhaustive host removal.
                           (Wood et al. 2019, Genome Biology 20:257)
  --minimum-hit-groups 3 : requires ≥3 distinct k-mer hit groups for a
                           classification call; reduces false positives compared
                           to the default of 2 (Wood et al. 2019).
"""

from __future__ import annotations
import subprocess
import zipfile
from pathlib import Path
from typing import List, Optional

from rich.panel import Panel
from rich.text import Text
from rich.progress import (
    Progress, SpinnerColumn, TextColumn,
    BarColumn, MofNCompleteColumn, TimeElapsedColumn,
)

from phageflow.utils.config import Config
from phageflow.utils.logger import log_step, log_info, log_ok, log_warn, log_error, console
from phageflow.utils.tools import require_tools, run_silent, human_size, mkdirs

STEP          = "02_host_removal"
TOOLS_BWA     = ["bwa-mem2", "samtools"]
TOOLS_K2      = ["kraken2", "seqtk"]
FASTA_EXTS    = {".fasta", ".fa", ".fna", ".fas"}
TAXID_VIRUSES = 10239

# ── Thresholds ────────────────────────────────────────────────────────────────
# Expected phage fraction in purified preparations (Roux et al. 2019, eLife)
_PHAGE_GOOD     = 90.0   # ≥90%: expected for good purified phage prep
_PHAGE_WARN     = 50.0   # <50%: real problem — wrong reference or contamination
_PHAGE_CRITICAL = 20.0   # <20%: assembly failure highly probable

_K2_UNCLASS_WARN = 0.30  # <30% unclassified → phage may be in Kraken2 database
_MIN_READS       = 50_000  # minimum retained reads for 50× on 150 kb genome
                            # (Nayfach et al. 2021; Bankevich et al. 2012)

# ── Minimum samtools version (legacy 0.1.x lacks -N and --singleton) ──────────
_SAMTOOLS_MIN_VERSION = (1, 15)


# ── Public entry point ────────────────────────────────────────────────────────

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
    """
    Remove host reads for a single sample.

    Singletons (reads whose pair mapped to host but are themselves unmapped)
    are written to {sample_id}_singletons.fastq.gz and should be passed to
    SPAdes as --s1 in Module 03 to preserve terminal phage coverage.

    Returns
    -------
    (r1_out, r2_out, singleton_out) : paths to host-filtered reads
    """
    r1 = Path(r1); r2 = Path(r2)
    out_dir  = cfg.results(STEP)
    rpt_dir  = cfg.reports(STEP)
    host_dir = cfg.workdir / "data" / "host_genomes"
    tmp_dir  = out_dir / "tmp"
    mkdirs(out_dir, rpt_dir, host_dir, tmp_dir)

    r1_out       = out_dir / f"{sample_id}_R1.fastq.gz"
    r2_out       = out_dir / f"{sample_id}_R2.fastq.gz"
    singleton_out = out_dir / f"{sample_id}_singletons.fastq.gz"
    log_f        = rpt_dir / f"{sample_id}_host_removal.log"

    log_step(f"Module 02 — host removal  [{sample_id}]")
    _print_input_table(sample_id, r1, r2)

    if r1_out.exists() and r1_out.stat().st_size > 0 and not force:
        log_info("  Already processed — skipping  (--force to re-run)")
        return r1_out, r2_out, singleton_out

    # ── Mode selection ────────────────────────────────────────────────────────
    use_bwa = bool(host_file or accessions or accessions_file)

    if use_bwa:
        require_tools(*TOOLS_BWA)
        _check_samtools_version()
        combined = host_dir / "combined_hosts.fasta"
        log_info("  Mode     : bwa-mem2 alignment  (bwa-mem2 + samtools ≥1.15)")
        log_info("  Singletons retained → {sample_id}_singletons.fastq.gz")
        log_info("  Vasimuddin et al. 2019 (bwa-mem2) · Li et al. 2009 (samtools)")
        n_steps = 3
    else:
        require_tools(*TOOLS_K2)
        db = kraken_db or cfg.databases.kraken2
        if not db or not Path(db).exists():
            log_error(
                "  Kraken2 database not found. "
                "Set databases.kraken2 in config.yaml or pass --host-file / --accessions."
            )
            return r1_out, r2_out, singleton_out
        log_info("  Mode     : kraken2 filter  (confidence=0.5 · minimum-hit-groups=3)")
        log_info(f"  Keeping  : unclassified + Viruses (taxid {TAXID_VIRUSES})")
        log_info("  Wood et al. 2019 (Genome Biology 20:257)")
        n_steps = 2

    # ── Pipeline steps ────────────────────────────────────────────────────────
    stats: dict = {}

    with Progress(
        SpinnerColumn(spinner_name="dots2"),
        TextColumn("[bold cyan]{task.description:<52}"),
        BarColumn(bar_width=20, complete_style="cyan", finished_style="green"),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        console=console,
        transient=True,
    ) as progress:
        task = progress.add_task("Initializing...", total=n_steps)

        if use_bwa:
            # 1/3 — resolve FASTAs + validate + build index
            progress.update(task, description="[1/3] bwa-mem2 — preparing reference")
            fastas = _resolve_fastas(host_file, accessions, accessions_file,
                                     host_dir, rpt_dir)
            if not fastas:
                log_error(
                    "  No host FASTAs resolved — check --host-file path, "
                    "accepted extensions (.fasta/.fa/.fna/.fas), or --accessions."
                )
                raise FileNotFoundError(
                    "No host reference FASTAs found. "
                    "Provide a valid --host-file or --accessions."
                )
            _build_index(fastas, combined, rpt_dir, force)
            progress.advance(task)

            # 2/3 — align reads (no intermediate name-sort)
            progress.update(task, description="[2/3] bwa-mem2 — aligning reads")
            bam = _align_bwa(sample_id, r1, r2, combined, tmp_dir, log_f,
                             cfg.threads)
            progress.advance(task)

            # 3/3 — extract phage reads + singletons
            progress.update(task, description="[3/3] samtools — extracting phage reads + singletons")
            stats = _extract_phage_bwa(bam, r1_out, r2_out, singleton_out, log_f)
            progress.advance(task)

        else:
            # 1/2 — classify reads
            progress.update(task, description="[1/2] kraken2 — classifying reads")
            report, k2_out = _classify_kraken2(
                sample_id, r1, r2, Path(db), tmp_dir, rpt_dir, log_f, cfg.threads
            )
            progress.advance(task)

            # 2/2 — extract + merge phage reads
            progress.update(task, description="[2/2] seqtk   — extracting phage reads")
            stats = _extract_phage_kraken2(
                sample_id, r1, r2, r1_out, r2_out,
                report, k2_out, tmp_dir, log_f,
            )
            progress.advance(task)

    import shutil as _sh
    _sh.rmtree(tmp_dir, ignore_errors=True)

    mode = "bwa-mem2" if use_bwa else "kraken2"

    # ── Metrics, display, save ────────────────────────────────────────────────
    _validate_output(sample_id, r1_out, r2_out, singleton_out, stats)
    _save_tsv({**stats, "sample": sample_id, "mode": mode},
              rpt_dir / "host_removal_summary.tsv")
    _print_summary_table(sample_id, stats, mode)
    _check_warnings(stats, mode)
    _print_completion_panel(sample_id, r1_out, r2_out, singleton_out, rpt_dir, stats)

    log_step(f"Module 02 completed ✓  [{sample_id}]")
    log_info(
        f"  Next: phageflow assembly --sample-id {sample_id} "
        f"--r1 {r1_out} --r2 {r2_out} --s1 {singleton_out}"
    )
    return r1_out, r2_out, singleton_out


# ── samtools version check ────────────────────────────────────────────────────

def _check_samtools_version() -> None:
    """
    Verify samtools >= 1.15 is available.

    samtools 0.1.x (legacy) does not implement the -N flag or --singleton
    option required by _extract_phage_bwa(). The conda environment may
    install 0.1.19 as a dependency of abricate/blast-legacy, silently
    shadowing the modern version.

    Reference: Li et al. (2009) Bioinformatics 25:2078 (samtools 1.x API).
    """
    try:
        result = subprocess.run(
            ["samtools", "--version"],
            capture_output=True, text=True, timeout=10,
        )
        first_line = result.stdout.strip().split("\n")[0]
        # Expected: "samtools 1.17" or "samtools 1.15-..."
        parts = first_line.split()
        if len(parts) >= 2:
            version_str = parts[1].split("-")[0]
            major, minor = (int(x) for x in version_str.split(".")[:2])
            if (major, minor) < _SAMTOOLS_MIN_VERSION:
                log_warn(
                    f"  samtools {version_str} detected — "
                    f"version ≥{'.'.join(str(x) for x in _SAMTOOLS_MIN_VERSION)} "
                    f"required for -N and --singleton flags. "
                    f"Update with: mamba install -n phageflow 'samtools>=1.15'"
                )
    except Exception:
        log_warn("  Could not determine samtools version — ensure ≥1.15 is active.")


# ── bwa-mem2 pipeline ─────────────────────────────────────────────────────────

def _resolve_fastas(
    host_file:       Optional[Path],
    accessions:      Optional[List[str]],
    accessions_file: Optional[Path],
    host_dir:        Path,
    rpt_dir:         Path,
) -> List[Path]:
    if host_file:
        p = Path(host_file)
        if p.is_dir():
            fastas = sorted(f for f in p.iterdir() if f.suffix.lower() in FASTA_EXTS)
            log_info(f"  Reference : folder  ({len(fastas)} FASTAs)")
            return fastas
        if p.suffix.lower() in FASTA_EXTS:
            log_info(f"  Reference : {p.name}  ({human_size(p)})")
            return [p]
        # treat as path-list file
        valid = [
            Path(line.strip()) for line in open(p)
            if line.strip() and not line.startswith("#")
            and Path(line.strip()).suffix.lower() in FASTA_EXTS
            and Path(line.strip()).exists()
        ]
        log_info(f"  Reference : FASTA list  ({len(valid)} files)")
        return valid

    if accessions:
        log_info(f"  Reference : NCBI accessions  ({len(accessions)} genome(s))")
        return _download_accessions(accessions, host_dir, rpt_dir)

    if accessions_file and Path(accessions_file).exists():
        accs = [
            line.strip() for line in open(accessions_file)
            if line.strip() and not line.startswith("#")
        ]
        log_info(f"  Reference : accessions file  ({len(accs)} genome(s))")
        return _download_accessions(accs, host_dir, rpt_dir)

    return []


def _download_accessions(
    accs: List[str], host_dir: Path, rpt_dir: Path
) -> List[Path]:
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
        log_warn(f"  datasets · download warning: {e}")
    fastas = []
    if zip_out.exists():
        with zipfile.ZipFile(zip_out) as z:
            z.extractall(dl_dir)
        zip_out.unlink()
        fastas = sorted(dl_dir.rglob("*.fna"))
        log_ok(f"  datasets · downloaded {len(fastas)} genome(s)")
    return fastas


def _build_index(
    fastas: List[Path], combined: Path, rpt_dir: Path, force: bool
) -> None:
    idx = Path(str(combined) + ".bwt.2bit.64")
    if idx.exists() and not force:
        log_info("  [bwa-mem2] index already exists — skipping")
        return

    log_info(f"  [bwa-mem2] building combined reference ({len(fastas)} genome(s))...")
    with open(combined, "w") as out:
        for fna in fastas:
            n = 0
            with open(fna) as f:
                for line in f:
                    out.write(f">{fna.stem}_{line[1:]}" if line.startswith(">") else line)
                    n += line.startswith(">")
            log_info(f"    + {fna.name}  ({n} seq(s), {human_size(fna)})")

    run_silent(["bwa-mem2", "index", str(combined)],
               log_file=rpt_dir / "bwa_index.log")
    run_silent(["samtools", "faidx", str(combined)],
               log_file=rpt_dir / "bwa_index.log")
    log_ok("  [bwa-mem2] index ready")


def _align_bwa(
    sample_id: str, r1: Path, r2: Path,
    combined: Path, tmp_dir: Path, log_f: Path, threads: int,
) -> Path:
    """
    Align reads with bwa-mem2 → unsorted BAM.

    sort -n REMOVED: samtools fastq with -N handles paired reads from an
    unsorted BAM without loss of correctness, avoiding the disk-bound
    O(n log n) sort step (Li et al. 2009, Bioinformatics 25:2078).
    The -M flag marks supplementary alignments as secondary for
    compatibility with downstream tools.
    """
    bam = tmp_dir / f"{sample_id}.bam"
    run_silent(
        f"bwa-mem2 mem -t {threads} -M {combined} {r1} {r2} "
        f"| samtools view -bS -@ 4 -o {bam}",
        log_file=log_f, shell=True,
    )
    return bam


def _extract_phage_bwa(
    bam: Path,
    r1_out: Path, r2_out: Path, singleton_out: Path,
    log_f: Path,
) -> dict:
    """
    Extract unmapped reads (phage) from BAM → FASTQ.gz.

    Paired phage reads (-f 12 -F 256):
        FLAG 4  : read unmapped
        FLAG 8  : mate unmapped
        Both flags set → neither read aligned to host → canonical phage pair.

    Singleton phage reads (-f 4 -F 8 -F 256):
        FLAG 4  : read unmapped
        FLAG 8  NOT set : mate DID map to host
        These are phage fragments whose pair happened to align to host,
        typically at DTR/ITR repeat boundaries. Retaining them improves
        terminal genome coverage, which CheckV uses to classify genomes as
        'Complete' (Nayfach et al. 2021, Nat Biotechnol 39:578).

    samtools fastq flags:
        -N  : append /1 /2 suffixes (Li et al. 2009); not available in 0.1.x
        --singleton : write singleton reads to a separate file

    References
    ----------
    Li et al. (2009) Bioinformatics 25:2078 — samtools.
    Nayfach et al. (2021) Nat Biotechnol 39:578 — DTR/ITR detection in CheckV.
    """
    # Paired unmapped reads → R1/R2 gzipped FASTQs
    run_silent(
        f"samtools view -bS -f 12 -F 256 {bam} "
        f"| samtools fastq -@ 4 -N "
        f"-1 {r1_out} -2 {r2_out} "
        f"--singleton {singleton_out}",
        log_file=log_f, shell=True,
    )

    # Count read categories for metrics
    def _c(flags: list) -> int:
        r = subprocess.run(
            ["samtools", "view", "-c"] + flags + [str(bam)],
            capture_output=True, text=True,
        )
        return int(r.stdout.strip() or 0)

    total      = _c([])
    host_r     = _c(["-F", "4"])                     # at least one read mapped
    phage_pair = _c(["-f", "12", "-F", "256"])        # both reads unmapped
    singletons = _c(["-f", "4", "-F", "8", "-F", "256"])  # only this read unmapped

    bam.unlink(missing_ok=True)

    pct_host      = f"{host_r      / total * 100:.2f}%" if total else "N/A"
    pct_phage     = f"{phage_pair  / total * 100:.2f}%" if total else "N/A"
    pct_singleton = f"{singletons  / total * 100:.2f}%" if total else "N/A"

    n_singletons = 0
    try:
        if singleton_out.exists():
            result = subprocess.run(
                ["bash", "-c", f"zcat {singleton_out} | wc -l"],
                capture_output=True, text=True,
            )
            n_singletons = int(result.stdout.strip() or 0) // 4
    except Exception:
        pass

    if n_singletons:
        log_ok(
            f"  [singletons] {n_singletons:,} singleton reads retained "
            f"({pct_singleton} of input) → pass to SPAdes --s1 "
            f"for DTR/ITR boundary coverage (Nayfach et al. 2021)"
        )

    return {
        "reads_in":       str(total),
        "reads_phage":    str(phage_pair),
        "reads_singleton":str(n_singletons),
        "pct_host":       pct_host,
        "pct_phage":      pct_phage,
        "pct_singleton":  pct_singleton,
    }


# ── Kraken2 pipeline ──────────────────────────────────────────────────────────

def _classify_kraken2(
    sample_id: str, r1: Path, r2: Path,
    db: Path, tmp_dir: Path, rpt_dir: Path, log_f: Path, threads: int,
) -> tuple:
    """
    Classify reads with Kraken2 using parameters optimised for phage retention.

    --confidence 0.5 (was 0.2):
        Higher specificity — reduces misclassification of phage reads as host.
        In purified_phage mode, retaining phage reads is the priority over
        exhaustive host removal (Wood et al. 2019, Genome Biology 20:257).

    --minimum-hit-groups 3 (was 2):
        Requires 3 distinct k-mer hit groups before calling a classification;
        reduces false positives at the cost of slightly lower host sensitivity,
        acceptable for purified preparations (Wood et al. 2019).
    """
    report = rpt_dir / f"{sample_id}_k2.report"
    k2_out = tmp_dir / f"{sample_id}_k2.tsv"

    run_silent([
        "kraken2",
        "--db",     str(db),
        "--report", str(report),
        "--output", str(k2_out),
        "--unclassified-out", f"{tmp_dir}/{sample_id}_unclassified#.fastq",
        "--paired",
        "--confidence",         "0.5",   # ← was 0.2; higher specificity (Wood 2019)
        "--minimum-hit-groups", "3",     # ← was 2; fewer false positives
        "--threads", str(threads),
        str(r1), str(r2),
    ], log_file=log_f, check=False)

    return report, k2_out


def _extract_phage_kraken2(
    sample_id: str, r1: Path, r2: Path,
    r1_out: Path, r2_out: Path,
    report: Path, k2_out: Path,
    tmp_dir: Path, log_f: Path,
) -> dict:
    """Extract unclassified + viral reads, merge to output FASTQ.gz."""
    un_r1     = tmp_dir / f"{sample_id}_unclassified_1.fastq"
    un_r2     = tmp_dir / f"{sample_id}_unclassified_2.fastq"
    viral_ids = tmp_dir / f"{sample_id}_viral_ids.txt"

    viral_taxids  = _parse_viral_taxids(report)
    n_viral_reads = _extract_viral_ids(k2_out, viral_taxids, viral_ids)
    log_info(
        f"  kraken2 · {len(viral_taxids)} viral taxid(s) | "
        f"{n_viral_reads:,} reads classified as Viruses"
    )
    _merge_reads(un_r1, un_r2, r1, r2, r1_out, r2_out,
                 viral_ids, n_viral_reads, log_f)

    # Parse stats from report
    n_total = n_unclass = n_viral = 0
    if report.exists():
        with open(report) as f:
            for line in f:
                p = line.strip().split("\t")
                if len(p) < 6:
                    continue
                try:
                    count = int(p[1])
                    name  = p[5].strip()
                    if name == "unclassified": n_unclass = count
                    if name == "root":         n_total   = count + n_unclass
                    if name == "Viruses":      n_viral   = count
                except ValueError:
                    pass

    n_kept   = n_unclass + n_viral
    pct_kept = f"{n_kept / n_total * 100:.2f}%" if n_total else "N/A"
    pct_host = f"{(n_total - n_kept) / n_total * 100:.2f}%" if n_total else "N/A"

    for p in (k2_out, viral_ids, un_r1, un_r2):
        if p and Path(p).exists():
            Path(p).unlink()

    return {
        "reads_in":        str(n_total),
        "reads_phage":     str(n_kept),
        "reads_singleton": "0",
        "pct_host":        pct_host,
        "pct_phage":       pct_kept,
        "pct_singleton":   "N/A",
        "k2_unclass":      str(n_unclass),
        "k2_viral":        str(n_viral),
    }


def _parse_viral_taxids(report: Path) -> set:
    taxids = set()
    if not report.exists():
        return taxids
    in_viral = False; viral_depth = None
    with open(report) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            try:
                taxid    = int(p[4].strip())
                name_raw = p[5]
                depth    = len(name_raw) - len(name_raw.lstrip())
                name     = name_raw.strip()
            except (ValueError, IndexError):
                continue
            if name == "Viruses" or taxid == TAXID_VIRUSES:
                in_viral = True; viral_depth = depth; taxids.add(taxid); continue
            if in_viral:
                if depth > viral_depth:
                    taxids.add(taxid)
                else:
                    in_viral = False; viral_depth = None
    return taxids


def _extract_viral_ids(k2_out: Path, viral_taxids: set, out_file: Path) -> int:
    n = 0
    if not k2_out.exists() or not viral_taxids:
        return n
    with open(k2_out) as f_in, open(out_file, "w") as f_out:
        for line in f_in:
            parts = line.split("\t")
            if len(parts) < 3 or parts[0].strip() != "C":
                continue
            try:
                if int(parts[2].strip()) in viral_taxids:
                    f_out.write(parts[1].strip() + "\n")
                    n += 1
            except ValueError:
                pass
    return n


def _merge_reads(
    un_r1: Path, un_r2: Path,
    src_r1: Path, src_r2: Path,
    r1_out: Path, r2_out: Path,
    viral_ids: Path, n_viral: int, log_f: Path,
) -> None:
    has_viral = viral_ids.exists() and n_viral > 0
    for un, src, out in [(un_r1, src_r1, r1_out), (un_r2, src_r2, r2_out)]:
        if un.exists() and has_viral:
            run_silent(
                f"{{ gzip -c {un}; seqtk subseq {src} {viral_ids} | gzip; }} > {out}",
                log_file=log_f, shell=True, check=False,
            )
        elif un.exists():
            run_silent(
                f"gzip -c {un} > {out}",
                log_file=log_f, shell=True, check=False,
            )
        elif has_viral:
            run_silent(
                f"seqtk subseq {src} {viral_ids} | gzip > {out}",
                log_file=log_f, shell=True, check=False,
            )


# ── Validation ────────────────────────────────────────────────────────────────

def _validate_output(
    sample_id: str,
    r1_out: Path, r2_out: Path,
    singleton_out: Path,
    stats: dict,
) -> None:
    """
    Validate output files and phage read counts.

    _MIN_READS = 50 000 is the empirical minimum for 50× coverage of a 150 kb
    phage genome at PE150, required for DTR/ITR detection in CheckV
    (Nayfach et al. 2021; Bankevich et al. 2012).
    """
    for label, path in [("R1", r1_out), ("R2", r2_out)]:
        if not path.exists() or path.stat().st_size == 0:
            log_warn(f"  [validate] {label} output missing or empty: {path}")

    try:
        n = int(stats.get("reads_phage", 0))
        if n < _MIN_READS:
            log_warn(
                f"  validate · only {n:,} phage reads retained — "
                f"< {_MIN_READS:,} reads may be insufficient for 50× coverage "
                f"of a 150 kb phage genome; Complete/HQ CheckV classification "
                f"may be compromised (Nayfach et al. 2021; Bankevich et al. 2012)."
            )
        else:
            log_ok(f"  validate · {n:,} phage reads retained — OK")
    except (ValueError, TypeError):
        pass


# ── Rich display helpers ──────────────────────────────────────────────────────

def _color_rate(val_str: str, good: float, warn: float) -> str:
    """Return Rich-colored string based on numeric thresholds."""
    try:
        v = float(val_str.rstrip("%"))
        if v >= good: return f"[bold green]{val_str}[/bold green]"
        if v >= warn: return f"[bold yellow]{val_str}[/bold yellow]"
        return f"[bold red]{val_str}[/bold red]"
    except (ValueError, AttributeError):
        return str(val_str)


def _fmt(val: str, unit: str = "") -> str:
    """Format metric value; dim 'N/A' entries."""
    if val in ("N/A", "", None):
        return "[dim]N/A[/dim]"
    return f"{val}{unit}"


def _print_input_table(sample_id: str, r1: Path, r2: Path) -> None:
    log_info(f"  R1 : {r1}  ({human_size(r1)})")
    log_info(f"  R2 : {r2}  ({human_size(r2)})")


def _print_summary_table(sample_id: str, stats: dict, mode: str) -> None:
    """Print host-removal metrics as compact log lines."""
    n_in  = stats.get("reads_in",        "N/A")
    n_out = stats.get("reads_phage",     "N/A")
    n_sin = stats.get("reads_singleton", "N/A")
    ph    = stats.get("pct_phage",       "N/A")
    ho    = stats.get("pct_host",        "N/A")

    # Kraken2 counts read pairs (not individual reads)
    reads_unit = "read pairs" if mode == "kraken2" else "reads"

    log_ok(
        f"  {reads_unit.capitalize()}: {n_in} → {n_out}  "
        f"({_color_rate(ph, _PHAGE_GOOD, _PHAGE_WARN)} phage retained)"
    )
    log_ok(
        f"  Host   : {_fmt(ho)} removed  |  "
        f"singletons={_fmt(n_sin)}  |  mode: {mode}"
    )
    if mode == "kraken2":
        log_info(
            f"  Note: Kraken2 reports {reads_unit}; "
            f"individual read count ≈ 2 × {n_in}"
        )
        unc = stats.get("k2_unclass", "N/A")
        vir = stats.get("k2_viral",   "N/A")
        log_ok(f"  kraken2: unclassified={_fmt(unc)}  viral={_fmt(vir)}")


def _check_warnings(stats: dict, mode: str) -> None:
    """Emit contextual warnings based on host-removal metrics."""
    try:
        ph = float(stats.get("pct_phage", "0").rstrip("%"))
        if ph < _PHAGE_CRITICAL:
            log_warn(
                f"  CRITICAL: phage fraction {stats.get('pct_phage')} < "
                f"{_PHAGE_CRITICAL}% — assembly failure highly probable. "
                "Verify host reference genome and sample purity."
            )
        elif ph < _PHAGE_WARN:
            log_warn(
                f"  Phage fraction {stats.get('pct_phage')} < {_PHAGE_WARN}% — "
                "possible wrong host reference, residual contamination, or "
                "low-purity preparation. Complete/HQ classification unlikely "
                "(Roux et al. 2019, eLife)."
            )
    except (ValueError, AttributeError):
        pass

    if mode == "kraken2":
        try:
            n_total   = int(stats.get("reads_in",   0))
            n_unclass = int(stats.get("k2_unclass", 0))
            if n_total and n_unclass / n_total < _K2_UNCLASS_WARN:
                log_warn(
                    f"  kraken2 · < {_K2_UNCLASS_WARN * 100:.0f}% unclassified — "
                    "phage may be partially represented in the Kraken2 database, "
                    "causing phage reads to be classified as host and discarded. "
                    "Consider switching to bwa-mem2 mode with --host-file."
                )
        except (ValueError, TypeError):
            pass


def _print_completion_panel(
    sample_id, r1_out, r2_out, singleton_out, rpt_dir, stats,
) -> None:
    n_in  = stats.get("reads_in",    "?")
    n_out = stats.get("reads_phage", "?")
    n_sin = stats.get("reads_singleton", "0")
    ph    = stats.get("pct_phage",   "?")
    ho    = stats.get("pct_host",    "?")

    text = Text()
    text.append("✓ ", style="bold green")
    text.append(f"{n_in} → ", style="dim white")
    text.append(f"{n_out}", style="bold green")
    text.append(f" retained  (phage: {ph}  |  host: {ho})\n", style="cyan")
    text.append(
        f"  singletons: {n_sin} retained for DTR/ITR coverage\n\n",
        style="dim white",
    )
    text.append("Output R1  : ", style="dim white")
    text.append(str(r1_out) + "\n", style="white")
    text.append("Output R2  : ", style="dim white")
    text.append(str(r2_out) + "\n", style="white")
    text.append("Singletons : ", style="dim white")
    text.append(str(singleton_out) + "\n", style="white")
    text.append("           : ", style="dim white")
    text.append("pass to SPAdes as --s1 (Nayfach et al. 2021)\n", style="dim white")
    text.append("Summary    : ", style="dim white")
    text.append(str(rpt_dir / "host_removal_summary.tsv"), style="white")

    console.print(Panel(
        text,
        title=f"[bold cyan]Host removal complete — {sample_id}[/bold cyan]",
        border_style="cyan",
        padding=(0, 2),
        width=90,
    ))


# ── TSV summary ───────────────────────────────────────────────────────────────

_TSV_HEADERS = [
    "sample", "mode",
    "reads_in", "reads_phage", "reads_singleton",
    "pct_host", "pct_phage", "pct_singleton",
    "k2_unclass", "k2_viral",
]


def _save_tsv(row: dict, path: Path) -> None:
    """
    Write / update the host-removal summary TSV.

    Existing rows are preserved and migrated to the current schema
    (missing columns filled with empty strings).
    """
    rows: dict[str, dict] = {}

    if path.exists():
        with open(path) as f:
            old_hdrs = f.readline().rstrip("\n").split("\t")
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if cols and cols[0]:
                    rows[cols[0]] = dict(zip(old_hdrs, cols))

    rows[row["sample"]] = {h: str(row.get(h, "")) for h in _TSV_HEADERS}

    with open(path, "w") as f:
        f.write("\t".join(_TSV_HEADERS) + "\n")
        for r in rows.values():
            f.write("\t".join(r.get(h, "") for h in _TSV_HEADERS) + "\n")

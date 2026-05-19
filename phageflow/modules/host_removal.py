"""PhageFlow Module 02 — Host read removal.

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
     Parameters: --confidence 0.2 --minimum-hit-groups 2
     (Wood et al. 2019, Genome Biology; aBIOTECH 2024)
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
_PHAGE_GOOD      = 70.0   # phage fraction: green
_PHAGE_WARN      = 30.0   # phage fraction: yellow; below → red + warning
_K2_UNCLASS_WARN = 0.30   # < 30% unclassified → phage may be in database
_MIN_READS       = 10_000 # fewer surviving reads → assembly-failure risk


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

    Returns
    -------
    (r1_out, r2_out) : paths to host-filtered reads
    """
    r1 = Path(r1); r2 = Path(r2)
    out_dir  = cfg.results(STEP)
    rpt_dir  = cfg.reports(STEP)
    host_dir = cfg.workdir / "data" / "host_genomes"
    tmp_dir  = out_dir / "tmp"
    mkdirs(out_dir, rpt_dir, host_dir, tmp_dir)

    r1_out = out_dir / f"{sample_id}_R1.fastq.gz"
    r2_out = out_dir / f"{sample_id}_R2.fastq.gz"
    log_f  = rpt_dir / f"{sample_id}_host_removal.log"

    log_step(f"Module 02 — host removal  [{sample_id}]")
    _print_input_table(sample_id, r1, r2)

    if r1_out.exists() and r1_out.stat().st_size > 0 and not force:
        log_info("  Already processed — skipping  (--force to re-run)")
        return r1_out, r2_out

    # ── Mode selection ─────────────────────────────────────────────────────────
    use_bwa = bool(host_file or accessions or accessions_file)

    if use_bwa:
        require_tools(*TOOLS_BWA)
        combined = host_dir / "combined_hosts.fasta"
        log_info("  Mode     : bwa-mem2 alignment  (bwa-mem2 + samtools)")
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
            return r1_out, r2_out
        log_info("  Mode     : kraken2 filter  (confidence=0.2 · minimum-hit-groups=2)")
        log_info(f"  Keeping  : unclassified + Viruses (taxid {TAXID_VIRUSES})")
        log_info("  Wood et al. 2019 (Genome Biology)")
        n_steps = 2

    # ── Pipeline steps ────────────────────────────────────────────────────────
    stats: dict = {}

    with Progress(
        SpinnerColumn(spinner_name="dots2"),
        TextColumn("[bold cyan]{task.description:<48}"),
        BarColumn(bar_width=20, complete_style="cyan", finished_style="green"),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        console=console,
        transient=True,
    ) as progress:
        task = progress.add_task("Initializing...", total=n_steps)

        if use_bwa:
            # 1/3 — resolve FASTAs + build index
            progress.update(task, description="[1/3] bwa-mem2 — preparing reference")
            fastas = _resolve_fastas(host_file, accessions, accessions_file,
                                     host_dir, rpt_dir)
            _build_index(fastas, combined, rpt_dir, force)
            progress.advance(task)

            # 2/3 — align reads
            progress.update(task, description="[2/3] bwa-mem2 — aligning reads")
            bam = _align_bwa(sample_id, r1, r2, combined, tmp_dir, log_f,
                             cfg.threads)
            progress.advance(task)

            # 3/3 — extract phage reads
            progress.update(task, description="[3/3] samtools — extracting phage reads")
            stats = _extract_phage_bwa(bam, r1_out, r2_out, log_f)
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

    # ── Metrics, display, save ─────────────────────────────────────────────────
    _validate_output(sample_id, r1_out, r2_out, stats)
    _save_tsv({**stats, "sample": sample_id, "mode": mode},
              rpt_dir / "host_removal_summary.tsv")
    _print_summary_table(sample_id, stats, mode)
    _check_warnings(stats, mode)
    _print_completion_panel(sample_id, r1_out, rpt_dir, stats)

    log_step(f"Module 02 completed ✓  [{sample_id}]")
    log_info(
        f"  Next: phageflow assembly --sample-id {sample_id} "
        f"--r1 {r1_out} --r2 {r2_out}"
    )
    return r1_out, r2_out


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
    if not fastas:
        log_error("  No host FASTAs found — check --host-file or --accessions.")
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
    run_silent(["samtools", "faidx",  str(combined)],
               log_file=rpt_dir / "bwa_index.log")
    log_ok("  [bwa-mem2] index ready")


def _align_bwa(
    sample_id: str, r1: Path, r2: Path,
    combined: Path, tmp_dir: Path, log_f: Path, threads: int,
) -> Path:
    """Align reads with bwa-mem2, sort by name → BAM."""
    bam = tmp_dir / f"{sample_id}.bam"
    run_silent(
        f"bwa-mem2 mem -t {threads} -M {combined} {r1} {r2} "
        f"| samtools view -bS -@ 4 "
        f"| samtools sort -n -@ 4 -T {tmp_dir}/{sample_id}_sort -o {bam}",
        log_file=log_f, shell=True,
    )
    return bam


def _extract_phage_bwa(
    bam: Path, r1_out: Path, r2_out: Path, log_f: Path,
) -> dict:
    """Extract unmapped read pairs (phage) from BAM → FASTQ.gz."""
    run_silent(
        f"samtools view -bS -f 12 -F 256 {bam} "
        f"| samtools fastq -@ 4 -1 {r1_out} -2 {r2_out} "
        f"-0 /dev/null -s /dev/null -n",
        log_file=log_f, shell=True,
    )

    def _c(flags):
        r = subprocess.run(["samtools", "view", "-c"] + flags + [str(bam)],
                           capture_output=True, text=True)
        return int(r.stdout.strip() or 0)

    total  = _c([])
    host_r = _c(["-F", "4"])
    phage  = _c(["-f", "12", "-F", "256"])
    bam.unlink(missing_ok=True)

    pct_host  = f"{host_r / total * 100:.2f}%" if total else "N/A"
    pct_phage = f"{phage  / total * 100:.2f}%" if total else "N/A"

    return {
        "reads_in":    str(total),
        "reads_phage": str(phage),
        "pct_host":    pct_host,
        "pct_phage":   pct_phage,
    }


# ── Kraken2 pipeline ──────────────────────────────────────────────────────────

def _classify_kraken2(
    sample_id: str, r1: Path, r2: Path,
    db: Path, tmp_dir: Path, rpt_dir: Path, log_f: Path, threads: int,
) -> tuple:
    """Classify reads with Kraken2. Returns (report_path, output_path)."""
    report = rpt_dir / f"{sample_id}_k2.report"
    k2_out = tmp_dir / f"{sample_id}_k2.tsv"

    run_silent([
        "kraken2",
        "--db",     str(db),
        "--report", str(report),
        "--output", str(k2_out),
        "--unclassified-out", f"{tmp_dir}/{sample_id}_unclassified#.fastq",
        "--paired",
        "--confidence",         "0.2",
        "--minimum-hit-groups", "2",
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

    # Cleanup temp files
    for p in (k2_out, viral_ids, un_r1, un_r2):
        if p and Path(p).exists():
            Path(p).unlink()

    return {
        "reads_in":    str(n_total),
        "reads_phage": str(n_kept),
        "pct_host":    pct_host,
        "pct_phage":   pct_kept,
        "k2_unclass":  str(n_unclass),
        "k2_viral":    str(n_viral),
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
    sample_id: str, r1_out: Path, r2_out: Path, stats: dict
) -> None:
    for label, path in [("R1", r1_out), ("R2", r2_out)]:
        if not path.exists() or path.stat().st_size == 0:
            log_warn(f"  [validate] {label} output missing or empty: {path}")
    try:
        n = int(stats.get("reads_phage", 0))
        if n < _MIN_READS:
            log_warn(
                f"  validate · only {n:,} phage reads retained — "
                f"assembly may fail with < {_MIN_READS:,} reads."
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
    n_in  = stats.get("reads_in",    "N/A")
    n_out = stats.get("reads_phage", "N/A")
    ph    = stats.get("pct_phage",   "N/A")
    ho    = stats.get("pct_host",    "N/A")

    log_ok(
        f"  Reads  : {n_in} → {n_out}  "
        f"({_color_rate(ph, _PHAGE_GOOD, _PHAGE_WARN)} phage retained)"
    )
    log_ok(f"  Host   : {_fmt(ho)} removed  |  mode: {mode}")

    if mode == "kraken2":
        unc = stats.get("k2_unclass", "N/A")
        vir = stats.get("k2_viral",   "N/A")
        log_ok(f"  kraken2: unclassified={_fmt(unc)}  viral={_fmt(vir)}")


def _check_warnings(stats: dict, mode: str) -> None:
    """Emit contextual warnings based on host-removal metrics."""
    try:
        ph = float(stats.get("pct_phage", "0").rstrip("%"))
        if ph < _PHAGE_WARN:
            log_warn(
                f"  Phage fraction {stats.get('pct_phage')} < {_PHAGE_WARN}% — "
                "possible contamination or wrong host reference."
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
                    "phage may be partially in the database."
                )
        except (ValueError, TypeError):
            pass


def _print_completion_panel(
    sample_id: str, r1_out: Path, rpt_dir: Path, stats: dict
) -> None:
    n_in  = stats.get("reads_in",    "?")
    n_out = stats.get("reads_phage", "?")
    ph    = stats.get("pct_phage",   "?")
    ho    = stats.get("pct_host",    "?")

    text = Text()
    text.append("✓ ", style="bold green")
    text.append(f"{n_in} reads  →  ", style="dim white")
    text.append(f"{n_out}", style="bold green")
    text.append(f"  (phage: {ph}  |  host: {ho})\n\n", style="cyan")
    text.append("Output R1 : ", style="dim white")
    text.append(str(r1_out) + "\n", style="white")
    text.append("Summary   : ", style="dim white")
    text.append(str(rpt_dir / "host_removal_summary.tsv"), style="white")

    console.print(Panel(
        text,
        title=f"[bold cyan]Host removal complete — {sample_id}[/bold cyan]",
        border_style="cyan",
        padding=(0, 2),
        width=72,
    ))


# ── TSV summary ───────────────────────────────────────────────────────────────

_TSV_HEADERS = [
    "sample", "mode",
    "reads_in", "reads_phage", "pct_host", "pct_phage",
    "k2_unclass", "k2_viral",
]


def _save_tsv(row: dict, path: Path) -> None:
    """
    Write / update the host-removal summary TSV.

    Existing rows are preserved and migrated to the current schema
    (missing columns filled with empty strings), so re-running with
    an older TSV does not lose previous samples.
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

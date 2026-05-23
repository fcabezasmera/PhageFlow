"""PhageFlow Module 02 — Host read removal.

Goal: maximise retention of phage reads — including reads at DTR/ITR boundaries —
while eliminating host contamination. Complete/HQ CheckV classification requires
uniform phage coverage across the entire genome, including terminal regions
(Nayfach et al. 2021, Nat Biotechnol 39:578).

Modes (priority order):

  1. --host-file  path
     Local FASTA (.fasta/.fna/.fa), folder of FASTAs, or text file of paths.
     Combined with DEFAULT_HOST_ACCESSIONS (e.g. E. coli K12) → bwa-mem2.

  2. --accessions GCF_1,GCF_2,...
     NCBI accessions + DEFAULT_HOST_ACCESSIONS → bwa-mem2.

  3. --accessions-file hosts.txt
     Text file with one GCF/GCA accession per line + DEFAULT_HOST_ACCESSIONS → bwa-mem2.

  4. (default — no reference provided)
     Kraken2 filter (confidence=0.5, min-hit-groups=3) → retain unclassified +
     Viruses (taxid 10239) → post-filter with bwa-mem2 vs DEFAULT_HOST_ACCESSIONS.

bwa-mem2 streaming pipeline
----------------------------
The alignment and extraction run in a single pipe — no BAM is written to disk:

  bwa-mem2 mem | samtools view | samtools collate | samtools fastq

  Filtering in samtools fastq (-f 4 -F 256):
    -f 4  : keep reads where THIS read is unmapped
    -F 256: exclude supplementary alignments
    Result:
      Both mates unmapped → -1 and -2 (phage read pairs)
      This read unmapped, mate mapped → -s (DTR/ITR boundary reads)

samtools version requirement
----------------------------
  samtools ≥ 1.15 required for -N flag and -s in samtools fastq.
  Fix: mamba install -n phageflow "samtools>=1.15"

Thresholds
----------
  _PHAGE_GOOD     = 90% : expected for well-prepared purified phage (Roux et al. 2019)
  _PHAGE_WARN     = 50% : real problem — wrong reference or contamination
  _PHAGE_CRITICAL = 20% : assembly failure highly probable
"""

from __future__ import annotations
import subprocess
import zipfile
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

# Accessions siempre incluidos como referencia de host genérica.
# Se suman a cualquier --accessions o --host-file que provea el usuario.
# En modo kraken2 se usan como post-filter bwa-mem2.
DEFAULT_HOST_ACCESSIONS: list[str] = [
    "GCF_000005845.2",   # Escherichia coli K-12 MG1655
]

# ── Thresholds ────────────────────────────────────────────────────────────────
_PHAGE_GOOD      = 90.0
_PHAGE_WARN      = 50.0
_PHAGE_CRITICAL  = 20.0
_K2_UNCLASS_WARN = 0.30
_MIN_READS       = 50_000

# samtools I/O concurrency caps (I/O-bound — no gain beyond these).
_ST_VIEW_THREADS    = 8
_ST_COLLATE_THREADS = 8
_ST_FASTQ_THREADS   = 4

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
    """Remove host reads for a single sample.

    Returns
    -------
    (r1_out, r2_out, singleton_out) : paths to host-filtered reads
    """
    r1 = Path(r1)
    r2 = Path(r2)
    out_dir  = cfg.results(STEP)
    rpt_dir  = cfg.reports(STEP) / sample_id
    host_dir = cfg.workdir / "data" / "host_genomes"
    tmp_dir  = out_dir / "tmp"
    mkdirs(out_dir, rpt_dir, host_dir, tmp_dir)

    r1_out        = out_dir / f"{sample_id}_R1.fastq.gz"
    r2_out        = out_dir / f"{sample_id}_R2.fastq.gz"
    singleton_out = out_dir / f"{sample_id}_singletons.fastq.gz"
    log_f         = rpt_dir / f"{sample_id}_host_removal.log"

    log_step(f"Module 02 — host removal  [{sample_id}]")
    log_info(f"  R1 : {r1}  ({human_size(r1)})")
    log_info(f"  R2 : {r2}  ({human_size(r2)})")

    if r1_out.exists() and r1_out.stat().st_size > 0 and not force:
        log_info("  Already processed — skipping  (--force to re-run)")
        return r1_out, r2_out, singleton_out

    use_bwa = bool(host_file or accessions or accessions_file)
    active_warnings: list[str] = []

    if use_bwa:
        require_tools(*TOOLS_BWA)
        _check_samtools_version()
        combined = host_dir / "combined_hosts.fasta"
        log_info("  Mode     : bwa-mem2 streaming  (no BAM intermediate)")
        log_info("  Pipeline : bwa-mem2 | view | collate | fastq  (single pipe)")
        log_info("  Singletons retained → DTR/ITR boundary coverage (Nayfach et al. 2021)")
        n_steps = 2
    else:
        require_tools(*TOOLS_K2)
        require_tools(*TOOLS_BWA)
        db = kraken_db or cfg.databases.kraken2
        if not db or not Path(db).exists():
            log_error(
                "  Kraken2 database not found. "
                "Provide --host-file, --accessions, or set databases.kraken2."
            )
            return r1_out, r2_out, singleton_out
        log_info("  Mode     : kraken2 + bwa-mem2 post-filter")
        log_info(f"  Keeping  : unclassified + Viruses (taxid {TAXID_VIRUSES})")
        log_info("  Ref      : Wood et al. 2019, Genome Biology 20:257")
        log_info(f"  Post-filter defaults : {', '.join(DEFAULT_HOST_ACCESSIONS)}")
        n_steps = 3

    stats: dict = {}

    with Progress(
        SpinnerColumn(spinner_name="dots2"),
        TextColumn("[bold cyan]{task.description:<56}"),
        BarColumn(bar_width=20, complete_style="cyan", finished_style="green"),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        console=console,
        transient=True,
    ) as progress:
        task = progress.add_task("Initializing...", total=n_steps)

        if use_bwa:
            # ── Modo bwa-mem2: host conocido ──────────────────────────────────
            progress.update(task, description="[1/2] bwa-mem2 — preparing reference + index")
            fastas = _resolve_fastas(host_file, accessions, accessions_file,
                                     host_dir, rpt_dir)
            if not fastas:
                log_error("  No host FASTAs resolved — check --host-file or --accessions.")
                raise FileNotFoundError("No host reference FASTAs found.")
            _build_index(fastas, combined, rpt_dir, force)
            progress.advance(task)

            progress.update(
                task,
                description="[2/2] bwa-mem2|collate|fastq — streaming alignment & extraction",
            )
            stats = _run_bwa_pipeline(
                sample_id, r1, r2,
                combined, r1_out, r2_out, singleton_out,
                tmp_dir, log_f, cfg.threads,
            )
            progress.advance(task)

        else:
            # ── Modo kraken2 + bwa-mem2 post-filter: host desconocido ─────────

            # Paso 1: clasificar con kraken2
            progress.update(task, description="[1/3] kraken2  — classifying reads")
            report, k2_out = _classify_kraken2(
                sample_id, r1, r2, Path(db), tmp_dir, rpt_dir, log_f, cfg.threads
            )
            progress.advance(task)

            # Paso 2: extraer reads de fago → paths intermedios en tmp/
            r1_k2 = tmp_dir / f"{sample_id}_k2_R1.fastq.gz"
            r2_k2 = tmp_dir / f"{sample_id}_k2_R2.fastq.gz"
            progress.update(task, description="[2/3] seqtk    — extracting phage reads")
            stats_k2 = _extract_phage_kraken2(
                sample_id, r1, r2, r1_k2, r2_k2,
                report, k2_out, tmp_dir, log_f,
            )
            progress.advance(task)

            # Paso 3: post-filter bwa-mem2 vs DEFAULT_HOST_ACCESSIONS
            progress.update(
                task,
                description="[3/3] bwa-mem2  — post-filter vs default hosts",
            )
            stats_bwa = _run_bwa_postfilter(
                sample_id, r1_k2, r2_k2,
                r1_out, r2_out, singleton_out,
                tmp_dir, log_f, host_dir, rpt_dir,
                cfg.threads, force,
            )
            progress.advance(task)

            # Merge: conserva reads_in y conteos k2, reemplaza con resultado final
            stats = {
                **stats_k2,
                "reads_phage":       stats_bwa["reads_phage"],
                "reads_singleton":   stats_bwa["reads_singleton"],
                "pct_phage":         stats_bwa["pct_phage"],
                "pct_host":          stats_bwa["pct_host"],
                "pct_singleton":     stats_bwa["pct_singleton"],
                "k2_pre_postfilter": stats_k2["reads_phage"],
            }

    import shutil as _sh
    _sh.rmtree(tmp_dir, ignore_errors=True)

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
    """Warn if legacy samtools 0.1.x is active."""
    try:
        r = subprocess.run(
            ["samtools", "--version"],
            capture_output=True, text=True, timeout=10,
        )
        parts = r.stdout.strip().split("\n")[0].split()
        if len(parts) >= 2:
            ver_str = parts[1].split("-")[0]
            major, minor = (int(x) for x in ver_str.split(".")[:2])
            if (major, minor) < _SAMTOOLS_MIN_VERSION:
                log_warn(
                    f"  CRITICAL: samtools {ver_str} is the legacy version. "
                    "Requires ≥1.15 for -N and -s. "
                    "Fix: mamba install -n phageflow 'samtools>=1.15'"
                )
    except Exception:
        log_warn("  Could not determine samtools version — ensure ≥1.15 is active.")


# ── bwa-mem2 streaming pipeline ───────────────────────────────────────────────

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
        valid = [
            Path(line.strip()) for line in open(p)
            if line.strip() and not line.startswith("#")
            and Path(line.strip()).suffix.lower() in FASTA_EXTS
            and Path(line.strip()).exists()
        ]
        log_info(f"  Reference : FASTA list  ({len(valid)} files)")
        return valid

    if accessions:
        combined_accs = list(dict.fromkeys(DEFAULT_HOST_ACCESSIONS + accessions))
        log_info(
            f"  Reference : NCBI accessions  ({len(accessions)} user + "
            f"{len(DEFAULT_HOST_ACCESSIONS)} default = {len(combined_accs)} genome(s))"
        )
        return _download_accessions(combined_accs, host_dir, rpt_dir)

    if accessions_file and Path(accessions_file).exists():
        accs = [
            line.strip() for line in open(accessions_file)
            if line.strip() and not line.startswith("#")
        ]
        combined_accs = list(dict.fromkeys(DEFAULT_HOST_ACCESSIONS + accs))
        log_info(
            f"  Reference : accessions file  ({len(accs)} user + "
            f"{len(DEFAULT_HOST_ACCESSIONS)} default = {len(combined_accs)} genome(s))"
        )
        return _download_accessions(combined_accs, host_dir, rpt_dir)

    return []


def _download_accessions(accs: List[str], host_dir: Path, rpt_dir: Path) -> List[Path]:
    require_tools("datasets")
    dl_dir   = host_dir / "ncbi_dataset"
    zip_out  = host_dir / "ncbi_download.zip"
    acc_file = host_dir / "accessions.txt"

    # Limpiar descargas previas para evitar acumulación entre muestras
    import shutil as _sh
    if dl_dir.exists():
        _sh.rmtree(dl_dir)

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
                    out.write(
                        f">{fna.stem}_{line[1:]}" if line.startswith(">") else line
                    )
                    n += line.startswith(">")
            log_info(f"    + {fna.name}  ({n} seq(s), {human_size(fna)})")
    run_silent(["bwa-mem2", "index", str(combined)], log_file=rpt_dir / "bwa_index.log")
    run_silent(["samtools", "faidx", str(combined)], log_file=rpt_dir / "bwa_index.log")
    log_ok("  [bwa-mem2] index ready")


def _run_bwa_pipeline(
    sample_id:    str,
    r1: Path, r2: Path,
    combined:     Path,
    r1_out:       Path,
    r2_out:       Path,
    singleton_out: Path,
    tmp_dir:      Path,
    log_f:        Path,
    threads:      int,
) -> dict:
    """Align + extract phage reads in a single streaming pipeline (no BAM on disk)."""
    collate_prefix = tmp_dir / f"{sample_id}.collate"
    st_view    = min(threads, _ST_VIEW_THREADS)
    st_collate = min(threads, _ST_COLLATE_THREADS)
    st_fastq   = min(threads, _ST_FASTQ_THREADS)

    cmd = (
        f"bwa-mem2 mem -t {threads} -M {combined} {r1} {r2} "
        f"| samtools view -@ {st_view} -bF 256 -u "
        f"| samtools collate -@ {st_collate} -O -u -T {collate_prefix} - "
        f"| samtools fastq -@ {st_fastq} -N "
        f"  -f 4 -F 256 "
        f"  -1 {r1_out} -2 {r2_out} "
        f"  -s {singleton_out} "
        f"  -0 /dev/null"
    )
    run_silent(cmd, log_file=log_f, shell=True)
    log_ok("  [bwa-mem2] streaming alignment + extraction complete")

    n_in_pairs  = _count_reads_fastq(r1,           threads)
    n_out_pairs = _count_reads_fastq(r1_out,        threads)
    n_sin       = _count_reads_fastq(singleton_out, threads)

    total_reads = n_in_pairs * 2
    phage_reads = n_out_pairs * 2 + n_sin
    host_reads  = max(0, total_reads - phage_reads)

    pct_phage = f"{phage_reads / total_reads * 100:.2f}%" if total_reads else "N/A"
    pct_host  = f"{host_reads  / total_reads * 100:.2f}%" if total_reads else "N/A"
    pct_sin   = f"{n_sin / total_reads * 100:.2f}%"       if total_reads else "N/A"

    if n_sin:
        log_ok(
            f"  [singletons] {n_sin:,} retained ({pct_sin} of input) "
            "— pass to assembly --s1 for DTR/ITR coverage (Nayfach et al. 2021)"
        )
    return {
        "reads_in":        str(total_reads),
        "reads_phage":     str(phage_reads),
        "reads_singleton": str(n_sin),
        "pct_host":        pct_host,
        "pct_phage":       pct_phage,
        "pct_singleton":   pct_sin,
    }


def _run_bwa_postfilter(
    sample_id:    str,
    r1_in:        Path,
    r2_in:        Path,
    r1_out:       Path,
    r2_out:       Path,
    singleton_out: Path,
    tmp_dir:      Path,
    log_f:        Path,
    host_dir:     Path,
    rpt_dir:      Path,
    threads:      int,
    force:        bool,
) -> dict:
    """Post-filter reads de kraken2 con bwa-mem2 vs DEFAULT_HOST_ACCESSIONS.

    Elimina contaminación residual de hosts genéricos (ej. E. coli K12) que
    Kraken2 puede dejar pasar por baja confianza o ausencia en su DB.

    Flujo:
      kraken2_R1/R2 (r1_in/r2_in, intermedios en tmp/)
        → bwa-mem2 | collate | fastq (-f 4 -F 256)
        → r1_out / r2_out / singleton_out (destino final)

    Los intermedios de kraken2 se eliminan al terminar.
    """
    log_info(
        f"  [post-filter] bwa-mem2 vs {len(DEFAULT_HOST_ACCESSIONS)} default host(s): "
        + ", ".join(DEFAULT_HOST_ACCESSIONS)
    )

    # Índice dedicado para defaults (separado del índice de hosts del usuario)
    combined = host_dir / "combined_defaults.fasta"
    fastas   = _download_accessions(DEFAULT_HOST_ACCESSIONS, host_dir, rpt_dir)
    _build_index(fastas, combined, rpt_dir, force)

    # Contar reads de entrada antes de alinear (para % sobre base kraken2)
    n_k2_pairs = _count_reads_fastq(r1_in, threads)

    collate_prefix = tmp_dir / f"{sample_id}.postfilter.collate"
    st_view    = min(threads, _ST_VIEW_THREADS)
    st_collate = min(threads, _ST_COLLATE_THREADS)
    st_fastq   = min(threads, _ST_FASTQ_THREADS)

    cmd = (
        f"bwa-mem2 mem -t {threads} -M {combined} {r1_in} {r2_in} "
        f"| samtools view -@ {st_view} -bF 256 -u "
        f"| samtools collate -@ {st_collate} -O -u -T {collate_prefix} - "
        f"| samtools fastq -@ {st_fastq} -N "
        f"  -f 4 -F 256 "
        f"  -1 {r1_out} -2 {r2_out} "
        f"  -s {singleton_out} "
        f"  -0 /dev/null"
    )
    run_silent(cmd, log_file=log_f, shell=True)
    log_ok("  [post-filter] bwa-mem2 post-filter complete")

    # Limpiar intermedios de kraken2
    for p in (r1_in, r2_in):
        if p.exists():
            p.unlink()

    n_out_pairs = _count_reads_fastq(r1_out,        threads)
    n_sin       = _count_reads_fastq(singleton_out, threads)

    total_k2    = n_k2_pairs * 2
    phage_reads = n_out_pairs * 2 + n_sin
    host_removed = max(0, total_k2 - phage_reads)

    pct_phage = f"{phage_reads  / total_k2 * 100:.2f}%" if total_k2 else "N/A"
    pct_host  = f"{host_removed / total_k2 * 100:.2f}%" if total_k2 else "N/A"
    pct_sin   = f"{n_sin        / total_k2 * 100:.2f}%" if total_k2 else "N/A"

    log_info(
        f"  [post-filter] {host_removed:,} reads eliminados "
        f"({pct_host} de reads post-kraken2)"
    )
    if n_sin:
        log_ok(
            f"  [singletons] {n_sin:,} retained ({pct_sin}) "
            "— pass to assembly --s1 (Nayfach et al. 2021)"
        )
    return {
        "reads_phage":     str(phage_reads),
        "reads_singleton": str(n_sin),
        "pct_phage":       pct_phage,
        "pct_host":        pct_host,
        "pct_singleton":   pct_sin,
    }


def _count_reads_fastq(path: Path, threads: int = 4) -> int:
    """Count reads in a gzipped FASTQ (lines ÷ 4). Uses pigz when available."""
    if not path.exists() or path.stat().st_size == 0:
        return 0
    try:
        import shutil as _sh
        decomp = (f"pigz -dc -p {threads}" if _sh.which("pigz") else "gzip -dc")
        r = subprocess.run(
            ["bash", "-c", f"{decomp} {path} | wc -l"],
            capture_output=True, text=True, timeout=600,
        )
        return int(r.stdout.strip()) // 4
    except Exception:
        return 0


# ── Kraken2 pipeline ──────────────────────────────────────────────────────────

def _classify_kraken2(
    sample_id: str, r1: Path, r2: Path,
    db: Path, tmp_dir: Path, rpt_dir: Path, log_f: Path, threads: int,
) -> tuple:
    report = rpt_dir / f"{sample_id}_k2.report"
    k2_out = tmp_dir  / f"{sample_id}_k2.tsv"
    run_silent([
        "kraken2",
        "--db",     str(db),
        "--report", str(report),
        "--output", str(k2_out),
        "--unclassified-out", f"{tmp_dir}/{sample_id}_unclassified#.fastq",
        "--paired",
        "--confidence",         "0.5",
        "--minimum-hit-groups", "3",
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
    taxids: set = set()
    if not report.exists():
        return taxids
    in_viral = False
    viral_depth: Optional[int] = None
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
                in_viral = True; viral_depth = depth; taxids.add(taxid)
                continue
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
    viral_taxids_str = {str(t) for t in viral_taxids}
    with (
        open(k2_out,   "r", buffering=1 << 20) as f_in,
        open(out_file, "w", buffering=1 << 20) as f_out,
    ):
        for line in f_in:
            if not line.startswith("C\t"):
                continue
            parts = line.split("\t")
            if len(parts) >= 3 and parts[2].strip() in viral_taxids_str:
                f_out.write(parts[1] + "\n")
                n += 1
    return n


def _merge_reads(
    un_r1: Path, un_r2: Path,
    src_r1: Path, src_r2: Path,
    r1_out: Path, r2_out: Path,
    viral_ids: Path, n_viral: int, log_f: Path,
) -> None:
    import shutil as _sh
    has_pigz  = bool(_sh.which("pigz"))
    zip_c     = "pigz -c"  if has_pigz else "gzip -c"
    zip_pipe  = "pigz"     if has_pigz else "gzip"
    has_viral = viral_ids.exists() and n_viral > 0

    cmds: list[str] = []
    for un, src, out in [(un_r1, src_r1, r1_out), (un_r2, src_r2, r2_out)]:
        if un.exists() and has_viral:
            cmds.append(
                f"{{ {zip_c} {un}; seqtk subseq {src} {viral_ids} | {zip_pipe}; }} > {out}"
            )
        elif un.exists():
            cmds.append(f"{zip_c} {un} > {out}")
        elif has_viral:
            cmds.append(f"seqtk subseq {src} {viral_ids} | {zip_pipe} > {out}")

    if cmds:
        run_silent(
            " & ".join(cmds) + " & wait",
            log_file=log_f, shell=True, check=False,
        )


# ── Validation ────────────────────────────────────────────────────────────────

def _validate_output(
    r1_out: Path, r2_out: Path, singleton_out: Path,
    stats: dict, active_warnings: list[str],
) -> None:
    for label, path in [("R1", r1_out), ("R2", r2_out)]:
        if not path.exists() or path.stat().st_size == 0:
            msg = f"{label} output missing or empty: {path}"
            log_warn(f"  [validate] {msg}")
            active_warnings.append(msg)
    try:
        n = int(stats.get("reads_phage", 0))
        if n < _MIN_READS:
            msg = (
                f"Only {n:,} phage reads retained — "
                f"< {_MIN_READS:,} may be insufficient for 50× on a 150 kb genome "
                "(Nayfach et al. 2021)."
            )
            log_warn(f"  [validate] {msg}")
            active_warnings.append(msg)
        else:
            log_ok(f"  [validate] {n:,} phage reads retained — OK")
    except (ValueError, TypeError):
        pass


# ── Warnings ──────────────────────────────────────────────────────────────────

def _check_warnings(stats: dict, mode: str, active_warnings: list[str]) -> None:
    try:
        ph = float(stats.get("pct_phage", "0").rstrip("%"))
        if ph < _PHAGE_CRITICAL:
            msg = (
                f"CRITICAL: phage fraction {stats.get('pct_phage')} < "
                f"{_PHAGE_CRITICAL}% — assembly failure highly probable."
            )
            log_warn(f"  {msg}")
            active_warnings.append(msg)
        elif ph < _PHAGE_WARN:
            msg = (
                f"Phage fraction {stats.get('pct_phage')} < {_PHAGE_WARN}% — "
                "possible wrong reference or contamination (Roux et al. 2019)."
            )
            log_warn(f"  {msg}")
            active_warnings.append(msg)
    except (ValueError, AttributeError):
        pass

    if mode == "kraken2":
        try:
            n_total   = int(stats.get("reads_in",   0))
            n_unclass = int(stats.get("k2_unclass", 0))
            if n_total and n_unclass / n_total < _K2_UNCLASS_WARN:
                msg = (
                    f"< {_K2_UNCLASS_WARN * 100:.0f}% unclassified — phage reads "
                    "may be in Kraken2 DB causing loss. "
                    "Consider bwa-mem2 mode with --host-file."
                )
                log_warn(f"  [kraken2] {msg}")
                active_warnings.append(msg)
        except (ValueError, TypeError):
            pass


# ── Completion panel ──────────────────────────────────────────────────────────

def _print_completion_panel(
    sample_id: str,
    r1_out: Path, r2_out: Path, singleton_out: Path,
    rpt_dir: Path, stats: dict, mode: str,
    active_warnings: list[str],
) -> None:

    def _c(val: str, good: float, warn: float) -> str:
        try:
            v = float(str(val).rstrip("%"))
            if v >= good: return f"[bold green]{val}[/bold green]"
            if v >= warn: return f"[bold yellow]{val}[/bold yellow]"
            return f"[bold red]{val}[/bold red]"
        except (ValueError, AttributeError):
            return str(val)

    def _f(val) -> str:
        return "[dim]N/A[/dim]" if str(val) in ("N/A", "", "None") else str(val)

    n_in  = _f(stats.get("reads_in",    "N/A"))
    n_out = _f(stats.get("reads_phage", "N/A"))
    n_sin = _f(stats.get("reads_singleton", "0"))
    ph    = stats.get("pct_phage", "N/A")
    ho    = stats.get("pct_host",  "N/A")
    reads_unit = "read pairs" if mode == "kraken2" else "reads"

    lines: list[str] = []
    lines.append("[bold]Key metrics[/bold]")
    lines.append(
        f"  [cyan]{reads_unit.capitalize():12s}:[/cyan]"
        f" {n_in} → {n_out}"
        f"  phage={_c(ph, _PHAGE_GOOD, _PHAGE_WARN)}"
        f"  host={_f(ho)}"
    )
    lines.append(
        f"  [cyan]Singletons  :[/cyan] {n_sin} retained"
        f"  [dim](pass to assembly --s1)[/dim]"
    )
    lines.append(f"  [cyan]Mode        :[/cyan] {mode}")
    if mode == "kraken2":
        k2_pre = _f(stats.get("k2_pre_postfilter", "N/A"))
        lines.append(
            f"  [dim]kraken2: unclassified={_f(stats.get('k2_unclass','N/A'))}"
            f"  viral={_f(stats.get('k2_viral','N/A'))}"
            f"  pre-postfilter={k2_pre}[/dim]"
        )

    lines.append("")
    lines.append("[bold]Output files[/bold]")
    lines.append(f"  Phage R1   : {r1_out}")
    lines.append(f"  Phage R2   : {r2_out}")
    lines.append(f"  Singletons : {singleton_out}")
    lines.append(f"  Summary    : {rpt_dir / 'host_removal_summary.tsv'}")

    if active_warnings:
        lines.append("")
        lines.append("[bold yellow]Warnings[/bold yellow]")
        for w in active_warnings:
            lines.append(f"  [yellow]⚠ {w}[/yellow]")

    console.print(Panel(
        "\n".join(lines),
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
    "k2_unclass", "k2_viral", "k2_pre_postfilter",
]


def _save_tsv(row: dict, path: Path) -> None:
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

"""PhageFlow command-line interface."""

from __future__ import annotations
import math
import multiprocessing
import re
import subprocess
import sys
from pathlib import Path
from typing import List, Optional

import click

from phageflow import __version__
from phageflow.utils.config import load_config, _default_threads
from phageflow.utils.logger import (
    log_header, log_info, log_ok, log_warn, log_error, log_step,
)

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

_AUTO_THREADS = _default_threads()

# ---------------------------------------------------------------------------
# Sample-ID inference
# ---------------------------------------------------------------------------

# Suffixes stripped from R1 filenames to produce a clean sample ID.
# Applied in order; first match wins.
_R1_PATTERNS = [
    r"[._\-]R1_001$",   # Illumina NextSeq:  sample_R1_001.fastq.gz
    r"[._\-]R1$",       # Standard:          sample_R1.fastq.gz
    r"[._\-]1_001$",    # Alternative:       sample_1_001.fastq.gz
    r"[._\-]1$",        # Short form:        sample_1.fastq.gz
    r"_read1$",         # Verbose:           sample_read1.fastq.gz
]

# FASTQ extensions to strip before applying R1 patterns.
_FASTQ_EXTS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


def _infer_sample_id(r1: str | Path) -> str:
    """Derive a sample ID from the R1 filename.

    Examples
    --------
    sampleA_R1_001.fastq.gz  →  sampleA
    sampleA_R1.fastq.gz      →  sampleA
    sampleA_1.fastq.gz       →  sampleA
    sampleA_1.fq.gz          →  sampleA

    Falls back to the full filename stem (minus extension) if no pattern
    matches — never crashes.
    """
    name = Path(r1).name.lower()
    stem = Path(r1).name

    # Strip FASTQ extension (case-insensitive)
    for ext in _FASTQ_EXTS:
        if name.endswith(ext):
            stem = stem[: -len(ext)]
            break

    # Strip R1 suffix
    for pattern in _R1_PATTERNS:
        new = re.sub(pattern, "", stem, flags=re.IGNORECASE)
        if new != stem:
            return new.strip("._-") or stem

    # No pattern matched — return stem as-is
    return stem


def _find_read_pairs(raw_dir: Path) -> list[tuple[Path, Path]]:
    """Scan *raw_dir* for paired FASTQ files.

    Matching strategy (in priority order):
      1. Files containing ``_R1`` / ``_R2`` (Illumina standard)
      2. Files ending in ``_1`` / ``_2`` (alternative convention)

    Returns a sorted list of ``(r1, r2)`` tuples.  Warns about R1 files
    without a matching R2 but never crashes.
    """
    pairs: list[tuple[Path, Path]] = []

    # Collect all FASTQ files
    fastqs = sorted(
        f for f in raw_dir.iterdir()
        if f.is_file() and any(
            f.name.lower().endswith(ext) for ext in _FASTQ_EXTS
        )
    )

    # Identify R1 files by pattern
    r1_files = [
        f for f in fastqs
        if re.search(r"[._\-]R1[._\-]|[._\-]R1\.", f.name, re.IGNORECASE)
        or re.search(r"[._\-]1[._\-]|[._\-]1\.", f.name, re.IGNORECASE)
    ]

    for r1 in r1_files:
        # Build the corresponding R2 name
        r2_name = re.sub(r"R1", "R2", r1.name, flags=re.IGNORECASE)
        if r2_name == r1.name:
            # No R1→R2 substitution worked; try _1→_2
            r2_name = re.sub(r"([._\-])1([._\-])", r"\g<1>2\2", r1.name)
        r2 = raw_dir / r2_name
        if r2.exists():
            pairs.append((r1, r2))
        else:
            log_warn(f"  R1 found but no matching R2: {r1.name}")

    return pairs


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _load(
    config:      str,
    workdir:     Optional[str] = None,
    output_dir:  Optional[str] = None,
    reports_dir: Optional[str] = None,
    threads:     Optional[int] = None,
    project:     Optional[str] = None,
) -> object:
    """Load config.yaml and apply CLI overrides."""
    
    if output_dir is None:
        log_error("--output-dir is required. Provide -o /path/to/output")
        sys.exit(1)
        
    if project:
        p = Path(project)
        config  = str(p / "config" / "config.yaml")
        workdir = str(p)
    cfg_path = Path(config)
    _using_bundled = False
    if not cfg_path.exists():
        from importlib.resources import files as _res_files
        _default = Path(str(_res_files("phageflow.config").joinpath("default_config.yaml")))
        if _default.exists():
            cfg_path = _default
            _using_bundled = True
        else:
            log_error(f"Config file not found: {cfg_path}")
            sys.exit(1)

    wd  = Path(workdir) if workdir else (Path(".").resolve() if _using_bundled else cfg_path.parent.parent)
    cfg = load_config(cfg_path, workdir=wd)

    if output_dir:
        base = Path(output_dir)
        cfg.set_output_dir(base / "results")
        if not reports_dir:
            cfg.set_reports_dir(base / "reports")

    if reports_dir:
        cfg.set_reports_dir(Path(reports_dir))

    if threads:
        cfg.threads = threads

    _log_threads(cfg)
    return cfg


def _log_threads(cfg) -> None:
    total = multiprocessing.cpu_count()
    if cfg.threads == _AUTO_THREADS:
        log_info(
            f"  Threads : {cfg.threads}/{total} logical CPUs "
            f"(auto 90% — override with -t N)"
        )
    else:
        log_info(f"  Threads : {cfg.threads}/{total} logical CPUs (manual override)")


# ---------------------------------------------------------------------------
# Shared option decorators
# ---------------------------------------------------------------------------

def common_options(f):
    """Options shared by every pipeline command."""
    f = click.option(
        "--project", default=None, type=click.Path(),
        help="Project directory (from phageflow init). Overrides --config and --workdir.",
    )(f)
    f = click.option(
        "-c", "--config",
        default="config/config.yaml", show_default=True,
        help="Path to config.yaml.",
    )(f)
    f = click.option(
        "-w", "--workdir", default=None,
        help="Pipeline working directory (must contain config/ and databases/).",
    )(f)
    f = click.option(
        "-o", "--output-dir", "output_dir",
        required=True,
        type=click.Path(),
        help="Base output directory.  results/ and reports/ are created inside it.",
    )(f)
    f = click.option(
        "--reports-dir", default=None, type=click.Path(),
        help="Override the reports directory independently of -o.",
    )(f)
    f = click.option(
        "-t", "--threads", default=None, type=int,
        help=(
            f"CPU threads  [default: auto {_AUTO_THREADS} = 90% of "
            f"{multiprocessing.cpu_count()} logical CPUs]."
        ),
    )(f)
    f = click.option(
        "--force", is_flag=True, default=False,
        help="Force re-run even if outputs already exist.",
    )(f)
    return f


def reads_options(f):
    """R1 / R2 / optional sample-id for modules that process paired reads."""
    f = click.option(
        "--sample-id", default=None,
        help=(
            "Sample identifier used for all output filenames.  "
            "Inferred from --r1 filename when omitted "
            "(e.g. sampleA_R1.fastq.gz → sampleA)."
        ),
    )(f)
    f = click.option(
        "--r1", "r1_path", required=True, type=click.Path(exists=True),
        help="R1 reads (FASTQ or FASTQ.gz).",
    )(f)
    f = click.option(
        "--r2", "r2_path", required=True, type=click.Path(exists=True),
        help="R2 reads (FASTQ or FASTQ.gz).",
    )(f)
    return f


def _resolve_sample_id(sample_id: Optional[str], r1_path: str) -> str:
    """Return *sample_id* if given, otherwise infer it from *r1_path*."""
    if sample_id:
        return sample_id
    inferred = _infer_sample_id(r1_path)
    log_info(f"  sample-id : {inferred}  (inferred from {Path(r1_path).name})")
    return inferred


# ---------------------------------------------------------------------------
# CLI group
# ---------------------------------------------------------------------------

@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(__version__, "-v", "--version", prog_name="PhageFlow")
def cli():
    """
    PhageFlow — bacteriophage genomics pipeline.

    \b
    Run modules in order for each sample:
      qc              Quality control and trimming (fastp + FastQC + MultiQC)
      host-removal    Remove host reads (bwa-mem2 or Kraken2)
      assembly        De novo assembly (SPAdes + MEGAHIT + cd-hit-est)
      viral-id        Viral identification (geNomad)
      quality         Genome quality and selection (CheckV)
      annotate        Structural annotation (Pharokka → Phold → Phynteny)
      safety          Biosafety screening (CARD + VFDB)
      report          Final HTML report per candidate genome

    \b
    Or run everything at once:
      run             Execute all modules for samples in a directory

    \b
      config          Copy default config.yaml template
      init            Initialise a new project directory
    """


# ---------------------------------------------------------------------------
# qc
# ---------------------------------------------------------------------------

@cli.command("qc")
@common_options
@reads_options
def cmd_qc(config, workdir, output_dir, reports_dir, threads, force, project,
           sample_id, r1_path, r2_path):
    """Quality control and trimming for a single sample.

    \b
    --sample-id is optional — derived from the R1 filename when omitted.

    \b
    Example (explicit):
      phageflow qc --r1 raw/sampleA_R1.fastq.gz --r2 raw/sampleA_R2.fastq.gz

    Example (explicit sample-id):
      phageflow qc --sample-id sampleA \\
        --r1 raw/sampleA_R1.fastq.gz --r2 raw/sampleA_R2.fastq.gz
    """
    cfg = _load(config, workdir, output_dir, reports_dir, threads, project)
    sid = _resolve_sample_id(sample_id, r1_path)
    from phageflow.modules.qc import run
    run(cfg, sample_id=sid, r1=Path(r1_path), r2=Path(r2_path), force=force)


# ---------------------------------------------------------------------------
# host-removal
# ---------------------------------------------------------------------------

@cli.command("host-removal")
@common_options
@reads_options
@click.option(
    "--host-file", default=None, type=click.Path(),
    help="Local host reference: FASTA file, folder of FASTAs, or text file of paths.",
)
@click.option(
    "--accessions", default=None,
    help="Comma-separated NCBI accessions to download (GCF/GCA).",
)
@click.option(
    "--accessions-file", default=None, type=click.Path(),
    help="Text file with one GCF/GCA accession per line.",
)
@click.option(
    "--kraken-db", default=None, type=click.Path(),
    help="Kraken2 database (used when no reference is provided).",
)
def cmd_host_removal(config, workdir, output_dir, reports_dir, threads, force, project,
                     sample_id, r1_path, r2_path,
                     host_file, accessions, accessions_file, kraken_db):
    """Remove host reads for a single sample.

    \b
    --sample-id is optional — derived from the R1 filename when omitted.

    \b
    Example:
      phageflow host-removal \\
        --r1 results/01_qc/sampleA_R1.fastq.gz \\
        --r2 results/01_qc/sampleA_R2.fastq.gz \\
        --accessions GCF_000005845.2
    """
    cfg  = _load(config, workdir, output_dir, reports_dir, threads, project)
    sid  = _resolve_sample_id(sample_id, r1_path)
    accs = [a.strip() for a in accessions.split(",")] if accessions else None
    from phageflow.modules.host_removal import run
    run(
        cfg, sample_id=sid,
        r1=Path(r1_path), r2=Path(r2_path),
        force=force,
        host_file=Path(host_file)             if host_file          else None,
        accessions=accs,
        accessions_file=Path(accessions_file) if accessions_file    else None,
        kraken_db=Path(kraken_db)             if kraken_db          else None,
    )


# ---------------------------------------------------------------------------
# assembly
# ---------------------------------------------------------------------------

@cli.command("assembly")
@common_options
@reads_options
@click.option(
    "--s1", "s1_path", default=None, type=click.Path(),
    help="Singleton reads from host-removal step (improves DTR/ITR coverage).",
)
def cmd_assembly(config, workdir, output_dir, reports_dir, threads, force, project,
                 sample_id, r1_path, r2_path, s1_path):
    """De novo assembly for a single sample.

    \b
    --sample-id is optional — derived from the R1 filename when omitted.
    """
    cfg = _load(config, workdir, output_dir, reports_dir, threads, project)
    sid = _resolve_sample_id(sample_id, r1_path)
    from phageflow.modules.assembly import run
    run(cfg, sample_id=sid,
        r1=Path(r1_path), r2=Path(r2_path),
        s1=Path(s1_path) if s1_path else None,
        force=force)


# ---------------------------------------------------------------------------
# viral-id
# ---------------------------------------------------------------------------

@cli.command("viral-id")
@common_options
@click.option("--sample-id", default=None,
              help="Sample identifier (derived from --contigs filename when omitted).")
@click.option(
    "--contigs", required=True, type=click.Path(exists=True),
    help="NR contigs FASTA (output of assembly step).",
)
def cmd_viral_id(config, workdir, output_dir, reports_dir, threads, force, project,
                 sample_id, contigs):
    """Viral identification with geNomad."""
    cfg = _load(config, workdir, output_dir, reports_dir, threads, project)
    sid = sample_id or _infer_sample_id(contigs)
    if not sample_id:
        log_info(f"  sample-id : {sid}  (inferred from {Path(contigs).name})")
    from phageflow.modules.viral_id import run
    run(cfg, sample_id=sid, contigs=Path(contigs), force=force)


# ---------------------------------------------------------------------------
# quality
# ---------------------------------------------------------------------------

@cli.command("quality")
@common_options
@click.option("--sample-id", default=None,
              help="Sample identifier (derived from --virus-fna filename when omitted).")
@click.option(
    "--virus-fna", required=True, type=click.Path(exists=True),
    help="Viral contigs FASTA (output of viral-id step).",
)
def cmd_quality(config, workdir, output_dir, reports_dir, threads, force, project,
                sample_id, virus_fna):
    """Genome quality assessment and selection with CheckV."""
    cfg = _load(config, workdir, output_dir, reports_dir, threads, project)
    sid = sample_id or Path(virus_fna).stem.removesuffix("_virus")
    if not sample_id:
        log_info(f"  sample-id : {sid}  (inferred from {Path(virus_fna).name})")
    from phageflow.modules.quality import run
    run(cfg, sample_id=sid, virus_fna=Path(virus_fna), force=force)


# ---------------------------------------------------------------------------
# annotate
# ---------------------------------------------------------------------------

@cli.command("annotate")
@common_options
@click.option("--sample-id", default=None,
              help="Sample identifier (derived from --genome filename when omitted).")
@click.option(
    "--genome", required=True, type=click.Path(exists=True),
    help="Genome FASTA from quality step (annotation_ready/).",
)
def cmd_annotate(config, workdir, output_dir, reports_dir, threads, force, project,
                 sample_id, genome):
    """Structural and functional annotation (Pharokka → Phold → Phynteny)."""
    cfg = _load(config, workdir, output_dir, reports_dir, threads, project)
    # genome path: results/05_quality/{sample_id}/annotation_ready/phages/{candidate}.fasta
    # Try to recover sample_id from path structure before falling back to stem.
    sid = sample_id or _sample_id_from_genome_path(genome)
    if not sample_id:
        log_info(f"  sample-id : {sid}  (inferred from path)")
    from phageflow.modules.annotate import run
    run(cfg, sample_id=sid, genome=Path(genome), force=force)


def _sample_id_from_genome_path(genome: str) -> str:
    """Walk up the genome path looking for the annotation_ready/phages structure."""
    p = Path(genome).resolve()
    parts = p.parts
    try:
        idx = parts.index("annotation_ready")
        # structure: .../{sample_id}/annotation_ready/phages/{candidate}.fasta
        return parts[idx - 1]
    except ValueError:
        return p.stem


# ---------------------------------------------------------------------------
# safety
# ---------------------------------------------------------------------------

@cli.command("safety")
@common_options
@click.option("--sample-id", default=None, help="Sample identifier.")
@click.option("--genome", required=True, type=click.Path(exists=True))
def cmd_safety(config, workdir, output_dir, reports_dir, threads, force, project,
               sample_id, genome):
    """Biosafety screening: CARD + VFDB."""
    cfg = _load(config, workdir, output_dir, reports_dir, threads, project)
    sid = sample_id or _sample_id_from_genome_path(genome)
    from phageflow.modules.safety import run
    run(cfg, sample_id=sid, genome=Path(genome), force=force)


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------

@cli.command("report")
@common_options
@click.option("--sample-id", default=None, help="Sample identifier.")
@click.option("--candidate-id", required=True,
              help="Candidate genome ID (e.g. Podoviridae_candidate_001).")
def cmd_report(config, workdir, output_dir, reports_dir, threads, force, project,
               sample_id, candidate_id):
    """Generate final HTML report for one annotated candidate genome."""
    cfg = _load(config, workdir, output_dir, reports_dir, threads, project)
    if not sample_id:
        log_error("--sample-id is required for the report command.")
        sys.exit(1)
    from phageflow.modules.report import run
    run(cfg, sample_id=sample_id, candidate_id=candidate_id, force=force)


# ---------------------------------------------------------------------------
# run  (master runner — directory-based, no samples.tsv)
# ---------------------------------------------------------------------------

@cli.command("run")
@common_options
@click.option(
    "--raw-dir", "raw_dir",
    default=None, type=click.Path(),
    help=(
        "Directory containing paired FASTQ files.  "
        "Pairs are auto-detected by *_R1* / *_R2* naming.  "
        "[default: <project>/raw/  or  raw/]"
    ),
)
@click.option(
    "--r1", "r1_list", multiple=True, type=click.Path(exists=True),
    help="R1 reads for a specific sample (repeat for multiple samples).",
)
@click.option(
    "--r2", "r2_list", multiple=True, type=click.Path(exists=True),
    help="R2 reads — must match --r1 in order.",
)
@click.option(
    "--host-file", default=None, type=click.Path(),
    help="Host reference FASTA / folder / path-list.",
)
@click.option(
    "--accessions", default=None,
    help="Comma-separated NCBI accessions for host download.",
)
@click.option(
    "--from-module",
    default=None,
    type=click.Choice(
        ["qc", "host-removal", "assembly", "viral-id",
         "quality", "annotate", "safety", "report"],
        case_sensitive=False,
    ),
    help="Resume pipeline from this module (skip earlier steps).",
)
def cmd_run(config, workdir, output_dir, reports_dir, threads, force, project,
            raw_dir, r1_list, r2_list, host_file, accessions, from_module):
    """Run the full pipeline for all samples.

    \b
    Sample discovery (priority order):
      1. --r1 / --r2 pairs supplied directly
      2. --raw-dir DIR  (scans DIR for *_R1*.fastq.gz + matching R2)
      3. <project>/raw/ or raw/  (default when using --project or init)

    \b
    Examples:
      # Auto-detect samples from raw/ directory
      phageflow run --project /data/myproject --accessions GCF_000005845.2

      # Explicit pairs (no directory scan)
      phageflow run \\
        --r1 raw/s1_R1.fastq.gz --r2 raw/s1_R2.fastq.gz \\
        --r1 raw/s2_R1.fastq.gz --r2 raw/s2_R2.fastq.gz \\
        --host-file /path/to/host.fasta

      # Resume from quality step
      phageflow run --project /data/myproject --from-module quality
    """
    cfg  = _load(config, workdir, output_dir, reports_dir, threads, project)
    accs = [a.strip() for a in accessions.split(",")] if accessions else None

    # ── Resolve pairs ─────────────────────────────────────────────────────────
    pairs: list[tuple[Path, Path]] = []

    if r1_list and r2_list:
        if len(r1_list) != len(r2_list):
            log_error("Number of --r1 and --r2 arguments must match.")
            sys.exit(1)
        pairs = [(Path(r1), Path(r2)) for r1, r2 in zip(r1_list, r2_list)]
        log_info(f"  Samples : {len(pairs)} pair(s) from explicit --r1/--r2 arguments")

    elif r1_list and not r2_list:
        log_error("--r1 provided without matching --r2.  Provide both.")
        sys.exit(1)

    else:
        # Auto-discover from a directory
        search_dirs: list[Path] = []
        if raw_dir:
            search_dirs = [Path(raw_dir)]
        elif project:
            search_dirs = [Path(project) / "raw"]
        else:
            # Try raw/ relative to workdir or cwd
            wd = Path(workdir) if workdir else Path(".")
            search_dirs = [wd / "raw", wd]

        for d in search_dirs:
            if d.is_dir():
                pairs = _find_read_pairs(d)
                if pairs:
                    log_info(f"  Samples : {len(pairs)} pair(s) discovered in {d}")
                    break
                log_info(f"  Scanned {d} — no FASTQ pairs found")

    if not pairs:
        log_error(
            "No sample pairs found.  Provide --r1/--r2, --raw-dir, "
            "or place FASTQ files in <project>/raw/."
        )
        sys.exit(1)

    # ── Run pipeline ──────────────────────────────────────────────────────────
    from phageflow.modules.runner import run
    run(
        cfg,
        pairs=pairs,
        force=force,
        host_file=Path(host_file) if host_file else None,
        accessions=accs,
        from_module=from_module,
    )

# ---------------------------------------------------------------------------
# config
# ---------------------------------------------------------------------------

@cli.command("config")
@click.option("-o", "--output", default="config/config.yaml", show_default=True)
def cmd_config(output):
    """Copy the default config.yaml template to OUTPUT."""
    import shutil
    from importlib.resources import files
    dst = Path(output)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(
        str(files("phageflow.config").joinpath("default_config.yaml")),
        dst,
    )
    log_ok(f"Config written to {dst}")


# ---------------------------------------------------------------------------
# init
# ---------------------------------------------------------------------------

@cli.command("init")
@click.argument("project_dir", type=click.Path())
@click.option("--force", is_flag=True, default=False,
              help="Re-initialize even if project already exists.")
def cmd_init(project_dir, force):
    """Initialize a PhageFlow project directory.

    \b
    Creates:
      PROJECT/
        config/config.yaml      ← copy of default config
        raw/                    ← place your FASTQ files here (auto-detected)
        results/                ← pipeline outputs
        reports/                ← logs, TSVs, HTML reports

    \b
    No samples manifest is needed — PhageFlow discovers FASTQ pairs in raw/.

    \b
    Next steps:
      1. Copy FASTQ files to PROJECT/raw/
      2. Edit PROJECT/config/config.yaml  (databases, threads, hosts …)
      3. phageflow run --project PROJECT --accessions GCF_XXXXXXXX.X
    """
    import shutil as _sh
    from importlib.resources import files as _files

    p = Path(project_dir)

    if (p / "config" / "config.yaml").exists() and not force:
        log_warn(
            f"  {p} already contains a config. "
            "Use --force to re-initialize."
        )
        return

    for sub in ("config", "raw", "results", "reports"):
        (p / sub).mkdir(parents=True, exist_ok=True)

    cfg_dst = p / "config" / "config.yaml"
    _sh.copy(
        str(_files("phageflow.config").joinpath("default_config.yaml")),
        cfg_dst,
    )

    (p / ".phageflow").write_text(f"project_dir: {p.resolve()}\n")

    log_ok(f"Project initialized at {p.resolve()}")
    log_info(f"  Config   : {cfg_dst}")
    log_info(f"  Raw reads: place FASTQ files in {p / 'raw'}/")
    log_info("")
    log_info("  Next steps:")
    log_info(f"    1. Copy FASTQ files to {p / 'raw'}/")
    log_info(f"    2. Edit {cfg_dst}")
    log_info(f"    3. phageflow run --project {p} --accessions GCF_XXXXXXXX.X")


# ---------------------------------------------------------------------------
# status
# ---------------------------------------------------------------------------

@cli.command("status")
@click.option("--project", default=None, type=click.Path())
@click.option("--sample-id", default=None)
def cmd_status(project, sample_id):
    """Show pipeline status for all samples in a project."""
    from rich.console import Console
    from rich.table import Table
    from phageflow.utils import status as _status

    reports_dir = Path(project) / "reports" if project else Path("reports")
    rows = _status.get_all(reports_dir)

    if not rows:
        log_info(f"No pipeline status found in {reports_dir}")
        return

    if sample_id:
        rows = [r for r in rows if r.get("sample_id") == sample_id]

    console = Console()
    table = Table(title="PhageFlow Pipeline Status", show_lines=True)
    table.add_column("Sample",   style="cyan",  min_width=12)
    table.add_column("Module",   style="white", min_width=14)
    table.add_column("Status",   style="bold",  min_width=10)
    table.add_column("Started",  style="dim",   min_width=17)
    table.add_column("Finished", style="dim",   min_width=17)
    table.add_column("Notes",    style="yellow")

    status_colors = {
        "completed": "[bold green]completed[/bold green]",
        "running":   "[bold cyan]running[/bold cyan]",
        "failed":    "[bold red]failed[/bold red]",
        "pending":   "[dim]pending[/dim]",
        "skipped":   "[dim]skipped[/dim]",
    }

    for row in rows:
        status_str = status_colors.get(row.get("status", ""), row.get("status", ""))
        table.add_row(
            row.get("sample_id", ""),
            row.get("module", ""),
            status_str,
            row.get("started_at", ""),
            row.get("finished_at", ""),
            row.get("notes", ""),
        )

    console.print(table)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    cli()

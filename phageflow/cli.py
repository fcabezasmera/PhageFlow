"""PhageFlow command-line interface."""

from __future__ import annotations
import math
import multiprocessing
import subprocess
import sys
from pathlib import Path
from typing import Optional

import click

from phageflow import __version__
from phageflow.utils.config import load_config, _default_threads
from phageflow.utils.logger import (
    log_header, log_info, log_ok, log_warn, log_error, log_step,
)

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

# Computed once at import time so the help text shows the real default.
_AUTO_THREADS = _default_threads()


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _load(
    config:      str,
    workdir:     Optional[str] = None,
    output_dir:  Optional[str] = None,
    reports_dir: Optional[str] = None,
    threads:     Optional[int] = None,
) -> object:
    """Load config.yaml and apply CLI overrides.

    Priority (highest → lowest):
        CLI flags  >  config.yaml values  >  auto-detected defaults
    """
    cfg_path = Path(config)
    if not cfg_path.exists():
        log_error(f"Config file not found: {cfg_path}")
        sys.exit(1)

    wd  = Path(workdir) if workdir else cfg_path.parent.parent
    cfg = load_config(cfg_path, workdir=wd)

    if output_dir:
        cfg.set_output_dir(Path(output_dir))
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
        "-c", "--config",
        default="config/config.yaml", show_default=True,
        help="Path to config.yaml.",
    )(f)
    f = click.option(
        "-w", "--workdir",
        default=None,
        help="Pipeline working directory (must contain config/ and databases/).",
    )(f)
    f = click.option(
        "-o", "--output-dir", "output_dir",
        default=None, type=click.Path(),
        help=(
            "Base output directory for results/  "
            "[default: workdir/results]. "
            "Overrides dirs.results in config.yaml."
        ),
    )(f)
    f = click.option(
        "--reports-dir",
        default=None, type=click.Path(),
        help=(
            "Base output directory for reports/  "
            "[default: workdir/reports]. "
            "Overrides dirs.reports in config.yaml."
        ),
    )(f)
    f = click.option(
        "-t", "--threads",
        default=None, type=int,
        help=(
            f"CPU threads  [default: auto {_AUTO_THREADS} = 90% of "
            f"{multiprocessing.cpu_count()} logical CPUs]. "
            "Overrides config.yaml and auto-detection."
        ),
    )(f)
    f = click.option(
        "--force", is_flag=True, default=False,
        help="Force re-run even if output already exists.",
    )(f)
    return f


def reads_options(f):
    """Options for modules that accept paired-end reads."""
    f = click.option(
        "--sample-id", required=True,
        help="Sample identifier (used for all output filenames).",
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
      safety          Biosafety screening (CARD + VFDB via Pharokka)
      report          Final HTML report per candidate genome

    \b
    Or run everything at once:
      run             Execute all modules for all samples in samples.tsv

    \b
    Utilities:
      check-tools     Verify all required tools are installed
      config          Copy default config.yaml template
      samples         Copy default samples.tsv template
    """


# ---------------------------------------------------------------------------
# qc
# ---------------------------------------------------------------------------

@cli.command("qc")
@common_options
@reads_options
def cmd_qc(config, workdir, output_dir, reports_dir, threads, force,
           sample_id, r1_path, r2_path):
    """Quality control and trimming for a single sample.

    \b
    Tools   : fastp · FastQC · MultiQC
    Params  : Q≥20 per-base · mean-Q≥25 · ≤10% low-qual · len≥75 bp
              poly-X filter · PE overlap correction · adapter auto-detect
    Ref     : Chen et al. 2018, Genome Biology 19:274

    \b
    Example:
      phageflow qc --sample-id s1 \\
        --r1 data/s1_R1.fastq.gz --r2 data/s1_R2.fastq.gz
    """
    cfg = _load(config, workdir, output_dir, reports_dir, threads)
    from phageflow.modules.qc import run
    run(cfg, sample_id=sample_id,
        r1=Path(r1_path), r2=Path(r2_path), force=force)


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
def cmd_host_removal(config, workdir, output_dir, reports_dir, threads, force,
                     sample_id, r1_path, r2_path,
                     host_file, accessions, accessions_file, kraken_db):
    """Remove host reads for a single sample.

    \b
    Modes (priority order):
      --host-file path        Local FASTA / folder / path-list
      --accessions GCF,...    Download specific NCBI genomes
      --accessions-file f.txt Accessions from file
      (no flag)               Auto-detect hosts with Kraken2

    Singleton reads (unmapped mates) are written to
    {sample}_singletons.fastq.gz — pass them to assembly --s1
    to improve DTR/ITR boundary coverage (Nayfach et al. 2021).

    \b
    Refs: Vasimuddin et al. 2019 (bwa-mem2)
          Li et al. 2009 (samtools)
          Wood et al. 2019 (Kraken2)

    \b
    Example:
      phageflow host-removal --sample-id s1 \\
        --r1 results/01_qc/s1_R1.fastq.gz \\
        --r2 results/01_qc/s1_R2.fastq.gz \\
        --host-file /path/to/host.fasta
    """
    cfg  = _load(config, workdir, output_dir, reports_dir, threads)
    accs = [a.strip() for a in accessions.split(",")] if accessions else None
    from phageflow.modules.host_removal import run
    run(
        cfg, sample_id=sample_id,
        r1=Path(r1_path), r2=Path(r2_path),
        force=force,
        host_file=Path(host_file)          if host_file          else None,
        accessions=accs,
        accessions_file=Path(accessions_file) if accessions_file else None,
        kraken_db=Path(kraken_db)          if kraken_db          else None,
    )


# ---------------------------------------------------------------------------
# assembly
# ---------------------------------------------------------------------------

@cli.command("assembly")
@common_options
@reads_options
@click.option(
    "--s1", "s1_path", default=None, type=click.Path(),
    help=(
        "Singleton reads from host-removal "
        "(results/02_host_removal/{sample}_singletons.fastq.gz). "
        "Passed to SPAdes --s1; improves DTR/ITR coverage "
        "(Nayfach et al. 2021)."
    ),
)
def cmd_assembly(config, workdir, output_dir, reports_dir, threads, force,
                 sample_id, r1_path, r2_path, s1_path):
    """De novo assembly for a single sample.

    \b
    Tools : SPAdes --isolate --only-assembler + MEGAHIT --no-mercy + cd-hit-est
    Refs  : Bankevich et al. 2012 · Li et al. 2015 · Fu et al. 2012

    \b
    Example:
      phageflow assembly --sample-id s1 \\
        --r1 results/02_host_removal/s1_R1.fastq.gz \\
        --r2 results/02_host_removal/s1_R2.fastq.gz \\
        --s1 results/02_host_removal/s1_singletons.fastq.gz
    """
    cfg = _load(config, workdir, output_dir, reports_dir, threads)
    from phageflow.modules.assembly import run
    run(cfg, sample_id=sample_id,
        r1=Path(r1_path), r2=Path(r2_path),
        s1=Path(s1_path) if s1_path else None,
        force=force)


# ---------------------------------------------------------------------------
# viral-id
# ---------------------------------------------------------------------------

@cli.command("viral-id")
@common_options
@click.option("--sample-id", required=True, help="Sample identifier.")
@click.option(
    "--contigs", required=True, type=click.Path(exists=True),
    help="NR contigs FASTA (output of assembly step).",
)
def cmd_viral_id(config, workdir, output_dir, reports_dir, threads, force,
                 sample_id, contigs):
    """Viral identification with geNomad.

    \b
    Key options (set in config.yaml):
      genomad.enable_score_calibration: true   (recommended for phage)
      genomad.composition: virome
      genomad.lenient_taxonomy: true           (genus-level co-binning)
      genomad.disable_find_proviruses: true    (purified phage)
      genomad.sensitivity: 4.2                 (increase to 5.7 for distant relatives)

    \b
    Ref: Camargo et al. 2023, Nat Biotechnol

    \b
    Example:
      phageflow viral-id --sample-id s1 \\
        --contigs results/03_assembly/combined/s1_contigs_nr.fasta
    """
    cfg = _load(config, workdir, output_dir, reports_dir, threads)
    from phageflow.modules.viral_id import run
    run(cfg, sample_id=sample_id, contigs=Path(contigs), force=force)


# ---------------------------------------------------------------------------
# quality
# ---------------------------------------------------------------------------

@cli.command("quality")
@common_options
@click.option("--sample-id", required=True, help="Sample identifier.")
@click.option(
    "--virus-fna", required=True, type=click.Path(exists=True),
    help="Viral contigs FASTA (output of viral-id step).",
)
def cmd_quality(config, workdir, output_dir, reports_dir, threads, force,
                sample_id, virus_fna):
    """Genome quality assessment and selection with CheckV.

    Selects Complete/HQ genomes for annotation.
    Rescues large ND contigs (myovirus rule) and draft co-bins.
    Concatemers (kmer_freq > threshold) are excluded automatically.

    \b
    Ref: Nayfach et al. 2021, Nat Biotechnol 39:578

    \b
    Example:
      phageflow quality --sample-id s1 \\
        --virus-fna results/04_viral_id/s1_virus.fna
    """
    cfg = _load(config, workdir, output_dir, reports_dir, threads)
    from phageflow.modules.quality import run
    run(cfg, sample_id=sample_id, virus_fna=Path(virus_fna), force=force)


# ---------------------------------------------------------------------------
# annotate
# ---------------------------------------------------------------------------

@cli.command("annotate")
@common_options
@click.option("--sample-id", required=True, help="Sample identifier.")
@click.option(
    "--genome", required=True, type=click.Path(exists=True),
    help="Genome FASTA from quality step (annotation_ready/).",
)
def cmd_annotate(config, workdir, output_dir, reports_dir, threads, force,
                 sample_id, genome):
    """Structural and functional annotation.

    \b
    Pipeline:
      Pharokka  → base GBK (PHANOTATE + PHROGs + tRNA/CRISPR)
      Phold     → improved /product via ProstT5 + Foldseek
      Phynteny  → /phynteny_* qualifiers via transformer + ESM2
      Merge     → canonical GBK (CDS from Phynteny · non-CDS from Pharokka)

    Genetic code is read from geNomad output per candidate.
    dnaapler reorientation applied only to DTR (circular) genomes.

    \b
    Refs: Bouras et al. 2023 (Pharokka) · Bouras et al. 2025 (Phold)
          Grigson et al. 2025 (Phynteny) · Bouras et al. 2024 (dnaapler)

    \b
    Example:
      phageflow annotate --sample-id s1 \\
        --genome results/05_quality/s1/annotation_ready/phages/Podoviridae_candidate_001.fasta
    """
    cfg = _load(config, workdir, output_dir, reports_dir, threads)
    from phageflow.modules.annotate import run
    run(cfg, sample_id=sample_id, genome=Path(genome), force=force)


# ---------------------------------------------------------------------------
# safety
# ---------------------------------------------------------------------------

@cli.command("safety")
@common_options
@click.option("--sample-id", required=True, help="Sample identifier.")
@click.option(
    "--genome", required=True, type=click.Path(exists=True),
    help="Genome FASTA.",
)
def cmd_safety(config, workdir, output_dir, reports_dir, threads, force,
               sample_id, genome):
    """Biosafety screening: CARD + VFDB via Pharokka output.

    \b
    Example:
      phageflow safety --sample-id s1 --genome candidate.fasta
    """
    cfg = _load(config, workdir, output_dir, reports_dir, threads)
    from phageflow.modules.safety import run
    run(cfg, sample_id=sample_id, genome=Path(genome), force=force)


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------

@cli.command("report")
@common_options
@click.option("--sample-id", required=True, help="Sample identifier.")
@click.option(
    "--candidate-id", required=True,
    help="Candidate genome ID (e.g. Podoviridae_candidate_001).",
)
def cmd_report(config, workdir, output_dir, reports_dir, threads, force,
               sample_id, candidate_id):
    """Generate final HTML report for one annotated candidate genome.

    Aggregates metrics from all previous modules into a single
    self-contained HTML file per candidate.

    \b
    Example:
      phageflow report --sample-id s1 \\
        --candidate-id Podoviridae_candidate_001
    """
    cfg = _load(config, workdir, output_dir, reports_dir, threads)
    from phageflow.modules.report import run
    run(cfg, sample_id=sample_id, candidate_id=candidate_id, force=force)


# ---------------------------------------------------------------------------
# run  (master runner — executes all modules for all samples)
# ---------------------------------------------------------------------------

@cli.command("run")
@common_options
@click.option(
    "--host-file", default=None, type=click.Path(),
    help="Host reference FASTA / folder / path-list (passed to host-removal).",
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
def cmd_run(config, workdir, output_dir, reports_dir, threads, force,
            host_file, accessions, from_module):
    """Run the full pipeline for all samples in samples.tsv.

    Executes in order: qc → host-removal → assembly → viral-id → quality
    → [for each candidate: annotate → safety → report]

    Failures in one sample are logged and the pipeline continues with
    the next sample. Use --from-module to resume an interrupted run.

    \b
    Example:
      phageflow run --host-file /path/to/host.fasta

      phageflow run --host-file /path/to/host.fasta \\
        --from-module quality -o /data/results
    """
    cfg  = _load(config, workdir, output_dir, reports_dir, threads)
    accs = [a.strip() for a in accessions.split(",")] if accessions else None
    from phageflow.modules.runner import run
    run(
        cfg, force=force,
        host_file=Path(host_file) if host_file else None,
        accessions=accs,
        from_module=from_module,
    )


# ---------------------------------------------------------------------------
# check-tools
# ---------------------------------------------------------------------------

@cli.command("check-tools")
@click.option("-c", "--config", default="config/config.yaml", show_default=True)
def cmd_check_tools(config):
    """Verify all required tools are installed and report versions."""
    from phageflow.utils.tools import check_tool, get_tool_version
    log_header(__version__)
    log_step("Checking required tools")

    # samtools legacy check FIRST — critical for host-removal.
    _check_samtools_version()

    TOOLS: dict[str, list[tuple[str, Optional[str], str]]] = {
        # category: [(tool, version_flag, install_hint)]
        "QC": [
            ("fastp",   "--version", ""),
            ("fastqc",  "--version", ""),
            ("multiqc", "--version", ""),
        ],
        "Host removal": [
            ("bwa-mem2",  "version",    ""),
            ("samtools",  "--version",  "must be ≥1.15"),
            ("seqtk",     None,         ""),
            ("pigz",      "--version",  "optional — parallel gzip"),
        ],
        "Assembly": [
            ("spades.py",  "--version", ""),
            ("megahit",    "--version", ""),
            ("cd-hit-est", None,        ""),
        ],
        "Viral ID": [
            ("genomad", "--version", ""),
        ],
        "Quality": [
            ("checkv", None,        ""),
            ("mash",   "--version", "required for circular-safe dereplication"),
            ("seqkit", "--version", ""),
        ],
        "Annotation": [
            ("pharokka.py",          "--version", ""),
            ("phold",                "--version", ""),
            ("phynteny_transformer", "--version", ""),
            ("dnaapler",             "--version", ""),
        ],
        "Safety": [
            ("abricate", "--version", ""),
        ],
        "Databases": [
            ("datasets", "version",    "NCBI datasets CLI"),
            ("kraken2",  "--version",  "optional — auto host detection"),
        ],
    }

    all_ok = True
    for category, tool_list in TOOLS.items():
        for tool, flag, hint in tool_list:
            found = check_tool(tool)
            if found:
                ver = get_tool_version(tool, flag) if flag else "(installed)"
                log_ok(f"  [{category:14}] {tool:28s} {ver}")
            else:
                suffix = f"  ← {hint}" if hint else ""
                log_warn(f"  [{category:14}] {tool:28s} NOT FOUND{suffix}")
                if tool not in ("pigz", "kraken2", "datasets"):
                    all_ok = False

    n_total   = multiprocessing.cpu_count()
    n_default = _default_threads()
    print()
    log_info(
        f"System  : {n_total} logical CPU(s) | "
        f"default threads = {n_default} (90% auto-detected)"
    )
    print()
    if all_ok:
        log_ok("All required tools found.")
    else:
        log_warn(
            "Some required tools are missing. "
            "Activate the phageflow conda environment: conda activate phageflow"
        )


def _check_samtools_version() -> None:
    """Warn if legacy samtools 0.1.x is active.

    abricate / blast-legacy pull in samtools=0.1.19 as a conda dependency.
    This version does NOT implement -N or --singleton, which are required
    by the bwa-mem2 host-removal pipeline (Module 02).

    Fix:
        mamba install -n phageflow "samtools>=1.15"
    """
    try:
        r = subprocess.run(
            ["samtools", "--version"],
            capture_output=True, text=True, timeout=10,
        )
        first = r.stdout.strip().split("\n")[0]
        parts = first.split()
        if len(parts) >= 2:
            ver_str = parts[1].split("-")[0]
            major, minor = (int(x) for x in ver_str.split(".")[:2])
            if (major, minor) < (1, 15):
                log_warn(
                    f"\n  *** CRITICAL: samtools {ver_str} on PATH is the LEGACY version ***\n"
                    "  The host-removal pipeline requires samtools ≥1.15 for\n"
                    "  the -N and --singleton flags (Li et al. 2009).\n"
                    "  Fix: mamba install -n phageflow 'samtools>=1.15'\n"
                )
    except Exception:
        log_warn("  Could not determine samtools version — ensure ≥1.15 is active.")


# ---------------------------------------------------------------------------
# config / samples templates
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


@cli.command("samples")
@click.option("-o", "--output", default="config/samples.tsv", show_default=True)
def cmd_samples(output):
    """Copy the default samples.tsv template to OUTPUT."""
    import shutil
    from importlib.resources import files
    dst = Path(output)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(
        str(files("phageflow.config").joinpath("default_samples.tsv")),
        dst,
    )
    log_ok(f"Samples template written to {dst}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    cli()

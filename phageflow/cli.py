"""PhageFlow command-line interface."""

from __future__ import annotations
import sys
from pathlib import Path
from typing import Optional

import click

from phageflow import __version__
from phageflow.utils.logger import (
    log_header, log_info, log_ok, log_warn, log_error, log_step
)
from phageflow.utils.config import load_config, load_samples

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


def _load(config: str, workdir: Optional[str]) -> tuple:
    cfg_path = Path(config)
    if not cfg_path.exists():
        log_error(f"Config file not found: {cfg_path}")
        sys.exit(1)
    wd  = Path(workdir) if workdir else cfg_path.parent.parent
    cfg = load_config(cfg_path, workdir=wd)
    return cfg


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(__version__, "-v", "--version", prog_name="PhageFlow")
def cli():
    """
    PhageFlow — modular bacteriophage genomics pipeline.

    \b
    Modules (run in order, one sample at a time):
      qc              Quality control and trimming
      host-removal    Remove host reads (local FASTA, NCBI accessions, or Kraken2)
      assembly        De novo assembly (SPAdes + MEGAHIT + cd-hit-est)
      viral-id        Viral identification (geNomad)
      quality         Genome quality and selection (CheckV)
      annotate        Structural annotation (Pharokka + Phold)
      safety          Biosafety screening (CARD + VFDB)
      lifecycle       Lifecycle prediction (BACPHLIP)
    """
    pass


def common_options(f):
    f = click.option("-c", "--config", default="config/config.yaml",
                     show_default=True, help="Path to config.yaml")(f)
    f = click.option("-w", "--workdir", default=None,
                     help="Pipeline working directory")(f)
    f = click.option("-t", "--threads", default=None, type=int,
                     help="Threads (overrides config)")(f)
    f = click.option("--force", is_flag=True, default=False,
                     help="Force re-run even if output exists")(f)
    return f


def reads_options(f):
    """Options for modules that accept raw reads input."""
    f = click.option("--sample-id", required=True,
                     help="Sample identifier (used for output filenames)")(f)
    f = click.option("--r1", "r1_path", required=True,
                     type=click.Path(exists=True),
                     help="R1 reads (FASTQ or FASTQ.gz)")(f)
    f = click.option("--r2", "r2_path", required=True,
                     type=click.Path(exists=True),
                     help="R2 reads (FASTQ or FASTQ.gz)")(f)
    return f


# ---------------------------------------------------------------------------
# qc
# ---------------------------------------------------------------------------

@cli.command("qc")
@common_options
@reads_options
def cmd_qc(config, workdir, threads, force, sample_id, r1_path, r2_path):
    """Quality control and trimming for a single sample.

    \b
    Tools: fastp + FastQC + MultiQC
    Parameters optimized for purified phage PE150 sequencing.
    Reference: Chen et al. 2018, Genome Biology.

    \b
    Example:
      phageflow qc --sample-id s1 \\
        --r1 reads_R1.fastq.gz --r2 reads_R2.fastq.gz
    """
    cfg = _load(config, workdir)
    if threads:
        cfg.threads = threads
    from phageflow.modules.qc import run
    run(cfg, sample_id=sample_id,
        r1=Path(r1_path), r2=Path(r2_path), force=force)


# ---------------------------------------------------------------------------
# host-removal
# ---------------------------------------------------------------------------

@cli.command("host-removal")
@common_options
@reads_options
@click.option("--host-file", default=None, type=click.Path(),
              help="Local host reference: single FASTA, folder of FASTAs, "
                   "or text file with one FASTA path per line.")
@click.option("--accessions", default=None,
              help="Comma-separated NCBI accessions (GCF/GCA). "
                   "Example: GCF_000013465.1,GCF_000007785.1")
@click.option("--accessions-file", default=None, type=click.Path(),
              help="Text file with one GCF/GCA accession per line")
@click.option("--kraken-db", default=None, type=click.Path(),
              help="Kraken2 database (used when no other reference provided)")
def cmd_host_removal(config, workdir, threads, force,
                     sample_id, r1_path, r2_path,
                     host_file, accessions, accessions_file,
                     kraken_db):
    """Remove host reads using bwa-mem2 for a single sample.

    \b
    Four modes (priority order):
      --host-file  path       Local FASTA, folder, or path list
      --accessions GCF_1,...  Download specific NCBI genomes
      --accessions-file f.txt Read accessions from file
      (no flag)               Auto-detect hosts with Kraken2

    \b
    Examples:
      phageflow host-removal --sample-id s1 --r1 R1.fq --r2 R2.fq \\
        --host-file /path/to/host.fasta
      phageflow host-removal --sample-id s1 --r1 R1.fq --r2 R2.fq \\
        --accessions GCF_000013465.1,GCF_000007785.1
      phageflow host-removal --sample-id s1 --r1 R1.fq --r2 R2.fq \\
        --kraken-db /path/to/db --kraken-mode filter
    """
    cfg = _load(config, workdir)
    if threads:
        cfg.threads = threads
    from phageflow.modules.host_removal import run
    accs = [a.strip() for a in accessions.split(",")] if accessions else None
    run(cfg, sample_id=sample_id,
        r1=Path(r1_path), r2=Path(r2_path),
        force=force,
        host_file=Path(host_file) if host_file else None,
        accessions=accs,
        accessions_file=Path(accessions_file) if accessions_file else None,
        kraken_db=Path(kraken_db) if kraken_db else None)


# ---------------------------------------------------------------------------
# assembly
# ---------------------------------------------------------------------------

@cli.command("assembly")
@common_options
@reads_options
def cmd_assembly(config, workdir, threads, force, sample_id, r1_path, r2_path):
    """De novo assembly for a single sample.

    \b
    Tools: metaSPAdes --only-assembler + MEGAHIT --no-mercy + cd-hit-est 100%
    References: Roux et al. 2019 (eLife); Li et al. 2015 (Bioinformatics).

    \b
    Example:
      phageflow assembly --sample-id s1 \\
        --r1 results/02_host_removal/s1_R1.fastq.gz \\
        --r2 results/02_host_removal/s1_R2.fastq.gz
    """
    cfg = _load(config, workdir)
    if threads:
        cfg.threads = threads
    from phageflow.modules.assembly import run
    run(cfg, sample_id=sample_id,
        r1=Path(r1_path), r2=Path(r2_path), force=force)


# ---------------------------------------------------------------------------
# viral-id
# ---------------------------------------------------------------------------

@cli.command("viral-id")
@common_options
@click.option("--sample-id", required=True, help="Sample identifier")
@click.option("--contigs", required=True, type=click.Path(exists=True),
              help="NR contigs FASTA (output of assembly step)")
def cmd_viral_id(config, workdir, threads, force, sample_id, contigs):
    """Viral identification with geNomad.

    \b
    Example:
      phageflow viral-id --sample-id s1 \\
        --contigs results/03_assembly/combined/s1_contigs_nr.fasta
    """
    cfg = _load(config, workdir)
    if threads:
        cfg.threads = threads
    from phageflow.modules.viral_id import run
    run(cfg, sample_id=sample_id, contigs=Path(contigs), force=force)


# ---------------------------------------------------------------------------
# quality
# ---------------------------------------------------------------------------

@cli.command("quality")
@common_options
@click.option("--sample-id", required=True, help="Sample identifier")
@click.option("--virus-fna", required=True, type=click.Path(exists=True),
              help="Viral contigs FASTA (output of viral-id step)")
def cmd_quality(config, workdir, threads, force, sample_id, virus_fna):
    """Genome quality assessment and selection with CheckV.

    \b
    Example:
      phageflow quality --sample-id s1 \\
        --virus-fna results/04_viral_id/s1_virus.fna
    """
    cfg = _load(config, workdir)
    if threads:
        cfg.threads = threads
    from phageflow.modules.quality import run
    run(cfg, sample_id=sample_id, virus_fna=Path(virus_fna), force=force)


# ---------------------------------------------------------------------------
# annotate
# ---------------------------------------------------------------------------

@cli.command("annotate")
@common_options
@click.option("--sample-id", required=True, help="Sample identifier")
@click.option("--genome", required=True, type=click.Path(exists=True),
              help="Genome FASTA (HQ genome from quality step)")
def cmd_annotate(config, workdir, threads, force, sample_id, genome):
    """Structural annotation with Pharokka + Phold.

    \b
    Example:
      phageflow annotate --sample-id s1 \\
        --genome results/05_quality/annotation_ready/s1_HQ.fasta
    """
    cfg = _load(config, workdir)
    if threads:
        cfg.threads = threads
    from phageflow.modules.annotate import run
    run(cfg, sample_id=sample_id, genome=Path(genome), force=force)


# ---------------------------------------------------------------------------
# safety
# ---------------------------------------------------------------------------

@cli.command("safety")
@common_options
@click.option("--sample-id", required=True, help="Sample identifier")
@click.option("--genome", required=True, type=click.Path(exists=True),
              help="Genome FASTA")
def cmd_safety(config, workdir, threads, force, sample_id, genome):
    """Biosafety screening: CARD + VFDB + integrase detection.

    \b
    Example:
      phageflow safety --sample-id s1 \\
        --genome results/05_quality/annotation_ready/s1_HQ.fasta
    """
    cfg = _load(config, workdir)
    if threads:
        cfg.threads = threads
    from phageflow.modules.safety import run
    run(cfg, sample_id=sample_id, genome=Path(genome), force=force)


# ---------------------------------------------------------------------------
# lifecycle
# ---------------------------------------------------------------------------

@cli.command("lifecycle")
@common_options
@click.option("--sample-id", required=True, help="Sample identifier")
@click.option("--genome", required=True, type=click.Path(exists=True),
              help="Genome FASTA")
def cmd_lifecycle(config, workdir, threads, force, sample_id, genome):
    """Lifecycle prediction with BACPHLIP (activate bacphlip_env first).

    \b
    Example:
      conda activate bacphlip_env
      phageflow lifecycle --sample-id s1 \\
        --genome results/05_quality/annotation_ready/s1_HQ.fasta
    """
    cfg = _load(config, workdir)
    if threads:
        cfg.threads = threads
    from phageflow.modules.lifecycle import run
    run(cfg, sample_id=sample_id, genome=Path(genome), force=force)


# ---------------------------------------------------------------------------
# check-tools
# ---------------------------------------------------------------------------

@cli.command("check-tools")
@click.option("-c", "--config", default="config/config.yaml")
def cmd_check_tools(config):
    """Check all required external tools are available in PATH."""
    from phageflow.utils.tools import check_tool, get_tool_version

    log_header(__version__)
    log_step("Checking required tools")

    TOOLS = {
        "QC":         [("fastp","--version"),("fastqc","--version"),("multiqc","--version")],
        "Alignment":  [("bwa-mem2","version"),("samtools",None),("seqtk",None)],
        "Assembly":   [("spades.py","--version"),("megahit","--version"),("cd-hit-est",None)],
        "Viral ID":   [("genomad","--version")],
        "Quality":    [("checkv",None)],
        "Annotation": [("pharokka.py","--version"),("phold","--version")],
        "Safety":     [("abricate","--version")],
        "Lifecycle":  [("bacphlip",None)],
        "Datasets":   [("datasets","version"),("kraken2","--version")],
    }

    all_ok = True
    for category, tool_list in TOOLS.items():
        for tool, flag in tool_list:
            found = check_tool(tool)
            if found:
                ver = get_tool_version(tool, flag) if flag else "(installed)"
                log_ok(f"  [{category:12}] {tool:20s} {ver}")
            else:
                log_warn(f"  [{category:12}] {tool:20s} NOT FOUND")
                all_ok = False

    print()
    if all_ok:
        log_ok("All tools found.")
    else:
        log_warn("Some tools missing — activate the correct conda environment.")


# ---------------------------------------------------------------------------
# config / samples templates
# ---------------------------------------------------------------------------

@cli.command("config")
@click.option("-o", "--output", default="config/config.yaml")
def cmd_config(output):
    """Copy the default config.yaml template."""
    import shutil
    from importlib.resources import files
    dst = Path(output)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(str(files("phageflow.config").joinpath("default_config.yaml")), dst)
    log_ok(f"Config written to {dst}")


@cli.command("samples")
@click.option("-o", "--output", default="config/samples.tsv")
def cmd_samples(output):
    """Copy the default samples.tsv template."""
    import shutil
    from importlib.resources import files
    dst = Path(output)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(str(files("phageflow.config").joinpath("default_samples.tsv")), dst)
    log_ok(f"Samples template written to {dst}")


if __name__ == "__main__":
    cli()

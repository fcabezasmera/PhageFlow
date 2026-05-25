"""PhageFlow Module 03b — Pre-assembly coverage profiling and validation.

Goal: map post-host-removal reads to the assembled contigs to compute
per-contig coverage metrics before viral identification. Provides early
warning of assembly quality issues and identifies contigs with suspiciously
high or low coverage.

This module runs AFTER assembly (03) and BEFORE viral identification (04).
It does not filter or remove any contigs — it only reports metrics and
warns when coverage patterns suggest problems.

Rationale
---------
Coverage profiling before viral ID provides three benefits:

1. Assembly quality assessment
   Fragmented assemblies (N50 < 5 kb) benefit from coverage context:
   a contig with high coverage but short length may be a real viral
   element, while a long contig with near-zero coverage is likely an
   assembly artifact. Roux et al. 2017 (Nature 537:689).

2. Ultra-high coverage detection
   In purified phage preparations, coverage > 1000x can degrade SPAdes
   assembly quality due to accumulation of systematic errors in de Bruijn
   graph traversal. Wick et al. 2021 (Genome Biol 22:241).
   Contigs with coverage > max_coverage_warn are flagged for review.

3. Multi-genome mixture detection
   Discrete coverage peaks suggest multiple phages co-purified.
   A bimodal coverage distribution (e.g. 500x and 2000x) indicates
   at least two distinct genomes in the sample.

Tool: CoverM contig (Woodcroft et al. 2023, github.com/wwood/CoverM)
  --coupled       : paired-end reads (bwa-mem2 internal alignment)
  --reference     : assembled contigs FASTA
  Methods computed:
    mean           : average depth across all positions
    trimmed_mean   : mean after removing top/bottom 5% positions
                     (robust to terminal artifacts in phage genomes)
    covered_bases  : positions covered at >= 1x (for breadth calculation)
    variance       : depth variance (high CV = uneven = possible chimera)

  Trimmed mean rationale: phage terminal repeats (DTR/ITR) are assembled
  twice and may appear at double depth at contig ends, inflating mean
  coverage. Trimmed mean (5-95%%) is more representative of true genome
  coverage. Bankevich et al. 2012 (J Comput Biol 19:455).

Coverage thresholds
-------------------
  min_coverage_warn   : contigs with mean < 5x flagged (assembly artifact)
  max_coverage_warn   : contigs with mean > 1000x flagged (ultra-high coverage)
  min_breadth_warn    : contigs with breadth < 80% flagged (uncovered regions)
  cv_warn             : CV > 2.0 flagged (highly uneven — possible chimera)
                        CV = std / mean

Output
------
  results/03b_coverage/{sample_id}/
    {sample_id}_coverage.tsv   : per-contig coverage metrics
  reports/03b_coverage/{sample_id}/
    coverage_summary.tsv       : module-level statistics
    {sample_id}_coverage.log   : CoverM log
"""

from __future__ import annotations
import math
import subprocess
import shutil
from pathlib import Path
from typing import Optional

from rich.panel import Panel
from rich.progress import (
    Progress, SpinnerColumn, TextColumn,
    BarColumn, MofNCompleteColumn, TimeElapsedColumn,
)

from phageflow.utils.config import Config
from phageflow.utils.logger import (
    log_step, log_info, log_ok, log_warn, log_error, console,
)
from phageflow.utils.tools import require_tools, run_silent, mkdirs

STEP  = "03b_coverage"
TOOLS = ["coverm"]

# Coverage thresholds
_MIN_COV_WARN    =     5.0   # < 5x → possible artifact
_MAX_COV_WARN    =  1000.0   # > 1000x → ultra-high coverage
_MIN_BREADTH     =     0.80  # < 80% covered → uncovered regions
_CV_WARN         =     2.0   # CV > 2.0 → highly uneven depth

# CoverM alignment filters
_MIN_READ_ID     =    95.0   # minimum read percent identity (%)
_MIN_READ_ALIGNED =   50.0   # minimum read aligned percent (%)


def run(
    cfg:       Config,
    sample_id: str,
    r1:        Path,
    r2:        Path,
    contigs:   Path,
    force:     bool = False,
) -> Path:
    """Compute per-contig coverage metrics for assembled contigs.

    Parameters
    ----------
    r1, r2  : post-host-removal paired reads (from Module 02)
    contigs : NR assembled contigs (from Module 03)

    Returns
    -------
    Path to per-contig coverage TSV
    """
    r1      = Path(r1)
    r2      = Path(r2)
    contigs = Path(contigs)

    out_dir = cfg.results(STEP) / sample_id
    rpt_dir = cfg.reports(STEP) / sample_id
    mkdirs(out_dir, rpt_dir)

    coverage_tsv = out_dir / f"{sample_id}_coverage.tsv"
    log_f        = rpt_dir / f"{sample_id}_coverage.log"

    log_step(f"Module 03b — coverage profiling  [{sample_id}]")
    log_info(f"  Contigs : {contigs}  ({_count_fasta(contigs):,} sequences)")
    log_info(f"  R1      : {r1}")
    log_info(f"  R2      : {r2}")

    if coverage_tsv.exists() and coverage_tsv.stat().st_size > 0 and not force:
        log_info("  Already computed — skipping  (--force to re-run)")
        return coverage_tsv

    if not contigs.exists() or contigs.stat().st_size == 0:
        log_warn("  Contigs file empty — skipping coverage profiling")
        return coverage_tsv

    require_tools(*TOOLS)
    active_warnings: list[str] = []

    with Progress(
        SpinnerColumn(spinner_name="dots2"),
        TextColumn("[bold cyan]{task.description:<60}"),
        BarColumn(bar_width=20, complete_style="cyan", finished_style="green"),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        console=console,
        transient=True,
    ) as progress:
        task = progress.add_task("Initializing...", total=3)

        # ── Step 1: CoverM ────────────────────────────────────────────────────
        progress.update(task, description="[1/3] CoverM  — mapping + coverage metrics")
        _run_coverm(r1, r2, contigs, coverage_tsv, log_f, cfg.threads)
        progress.advance(task)

        # ── Step 2: Parse + compute derived metrics ───────────────────────────
        progress.update(task, description="[2/3] Parsing coverage metrics")
        metrics = _parse_coverage(coverage_tsv)
        progress.advance(task)

        # ── Step 3: Warn + summarize ──────────────────────────────────────────
        progress.update(task, description="[3/3] Evaluating coverage distribution")
        stats = _evaluate_coverage(metrics, active_warnings)
        _save_summary_tsv(
            {**stats, "sample": sample_id},
            rpt_dir / "coverage_summary.tsv",
        )
        progress.advance(task)

    _print_completion_panel(
        sample_id, coverage_tsv, rpt_dir, stats, active_warnings,
    )
    log_step(f"Module 03b completed ✓  [{sample_id}]")
    return coverage_tsv


# ── CoverM runner ─────────────────────────────────────────────────────────────

def _run_coverm(
    r1:          Path,
    r2:          Path,
    contigs:     Path,
    output_tsv:  Path,
    log_f:       Path,
    threads:     int,
) -> None:
    """Run CoverM contig with multiple coverage methods.

    Methods:
      mean          : average depth (includes terminal repeat artifacts)
      trimmed_mean  : robust mean (5-95% trimmed — recommended for phages)
      covered_bases : positions covered ≥ 1x (for breadth computation)
      variance      : depth variance (for CV computation)

    Alignment filters:
      --min-read-percent-identity 95   : high-identity alignments only
      --min-read-aligned-percent 50    : at least 50% of read aligned
    """
    cmd = [
        "coverm", "contig",
        "--coupled", str(r1), str(r2),
        "--reference", str(contigs),
        "--methods", "mean", "trimmed_mean", "covered_bases", "variance",
        "--min-read-percent-identity", str(_MIN_READ_ID),
        "--min-read-aligned-percent",  str(_MIN_READ_ALIGNED),
        "--trim-min", "5",
        "--trim-max", "95",
        "--threads", str(min(threads, 16)),
        "--output-file", str(output_tsv),
    ]
    try:
        run_silent(cmd, log_file=log_f)
        log_ok("  [CoverM] coverage metrics computed")
    except Exception as e:
        log_warn(f"  [CoverM] failed: {e}")


# ── Coverage parsing ──────────────────────────────────────────────────────────

def _parse_coverage(coverage_tsv: Path) -> list[dict]:
    """Parse CoverM output and compute derived metrics.

    CoverM output columns (--methods mean trimmed_mean covered_bases variance):
      Contig  Mean  Trimmed Mean  Covered Bases  Variance

    Derived metrics:
      breadth : covered_bases / contig_length  (requires length from FASTA)
      cv      : sqrt(variance) / mean  (coefficient of variation)
    """
    if not coverage_tsv.exists() or coverage_tsv.stat().st_size == 0:
        return []

    metrics = []
    with open(coverage_tsv, newline="") as f:
        header = f.readline().rstrip("\n").split("\t")

        # CoverM uses full path as column prefix:
        # "uce01_contigs.fasta/uce01_R1.fastq.gz Mean"
        # Normalize headers by extracting just the metric suffix.
        def _metric(col: str) -> str:
            """Strip file path prefix from CoverM column header."""
            for suffix in ("Mean", "Trimmed Mean", "Covered Bases", "Variance"):
                if col.endswith(suffix):
                    return suffix
            return col
        norm_header = [_metric(h) for h in header]

        for line in f:
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0]:
                continue
            row = dict(zip(norm_header, parts))

            contig = parts[0]
            mean_cov  = _safe_float(row.get("Mean",          "0"))
            trimmed   = _safe_float(row.get("Trimmed Mean",  "0"))
            cov_bases = _safe_float(row.get("Covered Bases", "0"))
            variance  = _safe_float(row.get("Variance",      "0"))

            cv = math.sqrt(variance) / mean_cov if mean_cov > 0 else 0.0

            metrics.append({
                "contig":        contig,
                "mean_cov":      mean_cov,
                "trimmed_mean":  trimmed,
                "covered_bases": int(cov_bases),
                "variance":      variance,
                "cv":            cv,
            })

    return metrics


def _safe_float(val: str) -> float:
    try:
        return float(val)
    except (ValueError, TypeError):
        return 0.0


# ── Coverage evaluation ───────────────────────────────────────────────────────

def _evaluate_coverage(
    metrics:         list[dict],
    active_warnings: list[str],
) -> dict:
    """Evaluate coverage distribution and emit warnings.

    Returns summary statistics dict.
    """
    if not metrics:
        return {
            "n_contigs": "0", "n_low_cov": "0", "n_high_cov": "0",
            "n_low_breadth": "0", "n_high_cv": "0",
            "median_cov": "0", "max_cov": "0", "min_cov": "0",
        }

    coverages   = [m["mean_cov"]     for m in metrics]
    sorted_covs = sorted(coverages)
    n           = len(coverages)
    median_cov  = sorted_covs[n // 2]
    max_cov     = max(coverages)
    min_cov     = min(coverages)

    n_low_cov    = 0
    n_high_cov   = 0
    n_low_breadth = 0
    n_high_cv    = 0

    for m in metrics:
        mean    = m["mean_cov"]
        cv      = m["cv"]
        covered = m["covered_bases"]

        if mean < _MIN_COV_WARN and mean > 0:
            n_low_cov += 1

        if mean > _MAX_COV_WARN:
            n_high_cov += 1
            log_info(
                f"  [cov] {m['contig']}: {mean:.0f}x — ultra-high coverage. "
                "Possible concatemer or sequencing artifact."
            )

        if cv > _CV_WARN and mean > 10:
            n_high_cv += 1
            log_warn(
                f"  [cov] {m['contig']}: CV={cv:.2f} — highly uneven depth. "
                "Possible chimeric assembly or terminal repeat artifact."
            )

    # Summary warnings
    if n_low_cov > 0:
        pct = n_low_cov / n * 100
        msg = (
            f"{n_low_cov}/{n} contigs ({pct:.1f}%) have mean coverage "
            f"< {_MIN_COV_WARN}x — likely assembly artifacts or very low abundance."
        )
        log_warn(f"  [cov] {msg}")
        active_warnings.append(msg)

    if n_high_cov > 0:
        msg = (
            f"{n_high_cov}/{n} contigs have mean coverage "
            f"> {_MAX_COV_WARN:.0f}x — ultra-high coverage. "
            "Consider checking for concatemers in quality module."
        )
        log_warn(f"  [cov] {msg}")
        active_warnings.append(msg)

    if n_high_cv > 0:
        msg = (
            f"{n_high_cv}/{n} contigs have CV > {_CV_WARN} — "
            "uneven depth. Review assembly graph for chimeric contigs."
        )
        log_warn(f"  [cov] {msg}")
        active_warnings.append(msg)

    log_ok(
        f"  [cov] median={median_cov:.0f}x  "
        f"min={min_cov:.0f}x  max={max_cov:.0f}x  "
        f"n={n:,}"
    )

    return {
        "n_contigs":    str(n),
        "n_low_cov":    str(n_low_cov),
        "n_high_cov":   str(n_high_cov),
        "n_low_breadth": str(n_low_breadth),
        "n_high_cv":    str(n_high_cv),
        "median_cov":   f"{median_cov:.1f}",
        "max_cov":      f"{max_cov:.1f}",
        "min_cov":      f"{min_cov:.1f}",
    }


# ── Helpers ───────────────────────────────────────────────────────────────────

def _count_fasta(path: Path) -> int:
    if not path.exists():
        return 0
    n = 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n


# ── Completion panel ──────────────────────────────────────────────────────────

def _print_completion_panel(
    sample_id:       str,
    coverage_tsv:    Path,
    rpt_dir:         Path,
    stats:           dict,
    active_warnings: list[str],
) -> None:
    lines = [
        "[bold]Coverage summary[/bold]",
        f"  [cyan]Contigs profiled :[/cyan] {stats.get('n_contigs', '0')}",
        f"  [cyan]Median coverage  :[/cyan] {stats.get('median_cov', '?')}x",
        f"  [cyan]Range            :[/cyan] {stats.get('min_cov', '?')}x – {stats.get('max_cov', '?')}x",
        f"  [cyan]Low cov (<{_MIN_COV_WARN:.0f}x)  :[/cyan] {stats.get('n_low_cov', '0')} contigs",
        f"  [cyan]High cov (>{_MAX_COV_WARN:.0f}x):[/cyan] {stats.get('n_high_cov', '0')} contigs",
        f"  [cyan]High CV (>{_CV_WARN})  :[/cyan] {stats.get('n_high_cv', '0')} contigs",
        "",
        "[bold]Output[/bold]",
        f"  Coverage TSV : {coverage_tsv}",
        f"  Summary      : {rpt_dir / 'coverage_summary.tsv'}",
    ]
    if active_warnings:
        lines += ["", "[bold yellow]Warnings[/bold yellow]"]
        lines += [f"  [yellow]⚠ {w}[/yellow]" for w in active_warnings]

    console.print(Panel(
        "\n".join(lines),
        title=f"[bold cyan]Coverage profiling — {sample_id}[/bold cyan]",
        border_style="cyan", padding=(0, 2), width=120,
    ))


# ── TSV summary ───────────────────────────────────────────────────────────────

_TSV_HEADERS = [
    "sample", "n_contigs", "n_low_cov", "n_high_cov",
    "n_low_breadth", "n_high_cv",
    "median_cov", "max_cov", "min_cov",
]


def _save_summary_tsv(row: dict, path: Path) -> None:
    rows: dict[str, dict] = {}
    if path.exists():
        with open(path) as f:
            old_h = f.readline().rstrip("\n").split("\t")
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if cols and cols[0]:
                    rows[cols[0]] = dict(zip(old_h, cols))
    rows[row["sample"]] = {h: str(row.get(h, "")) for h in _TSV_HEADERS}
    with open(path, "w") as f:
        f.write("\t".join(_TSV_HEADERS) + "\n")
        for r in rows.values():
            f.write("\t".join(r.get(h, "") for h in _TSV_HEADERS) + "\n")

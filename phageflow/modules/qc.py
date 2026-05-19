"""PhageFlow Module 01 — Quality control and trimming.

Designed for purified phage Illumina PE sequencing.

Tools:
    fastp   : adapter trimming, sliding-window quality trimming,
              low-complexity filtering (Chen et al. 2018, Genome Biology)
    FastQC  : per-read quality metrics (Andrews 2010)
    MultiQC : aggregate QC report (Ewels et al. 2016, Bioinformatics)

fastp parameters (literature-based for phage sequencing):
    --cut_right / --cut_right_window_size 4 / --cut_right_mean_quality 20
        Sliding-window 3' trimming (Bolger et al. 2014; Chen et al. 2018)
    --qualified_quality_phred 20
        Q20 threshold: 99% base call accuracy (Illumina quality guidelines)
    --unqualified_percent_limit 20
        Max 20% low-quality bases per read (Chen et al. 2018)
    --length_required 75
        Minimum 75 bp for PE150 data (Bankevich et al. 2012)
    --low_complexity_filter --complexity_threshold 30
        Removes homopolymer/repetitive reads (Roux et al. 2019, eLife)
    --n_base_limit 5
        Remove reads with >5 N calls (standard practice)
    --detect_adapter_for_pe
        Automatic adapter detection for paired-end data (Chen et al. 2018)
"""

from __future__ import annotations
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from rich.table import Table
from rich.progress import (
    Progress, SpinnerColumn, TextColumn,
    BarColumn, TaskProgressColumn, TimeElapsedColumn,
)

from phageflow.utils.config import Config
from phageflow.utils.logger import log_step, log_info, log_ok, log_warn, console
from phageflow.utils.tools import require_tools, run_silent, human_size, mkdirs

STEP  = "01_qc"
TOOLS = ["fastp", "fastqc", "multiqc"]

# Thresholds for coloring pass rate in the summary table
_PASS_GOOD = 90.0   # green
_PASS_WARN = 80.0   # yellow (below → red)


def run(
    cfg:       Config,
    sample_id: str,
    r1:        Path,
    r2:        Path,
    force:     bool = False,
) -> tuple[Path, Path]:
    """
    Run QC + trimming for a single sample.

    Returns
    -------
    (r1_out, r2_out) : paths to trimmed reads
    """
    require_tools(*TOOLS)

    out_dir = cfg.results(STEP)
    rpt_dir = cfg.reports(STEP)
    mkdirs(out_dir, rpt_dir)

    r1 = Path(r1); r2 = Path(r2)

    # ── Validate inputs early ─────────────────────────────────────────────────
    for label, path in [("R1", r1), ("R2", r2)]:
        if not path.exists():
            raise FileNotFoundError(f"[{sample_id}] {label} not found: {path}")
        if path.stat().st_size == 0:
            raise ValueError(f"[{sample_id}] {label} is empty (0 bytes): {path}")

    r1_out = out_dir / f"{sample_id}_R1.fastq.gz"
    r2_out = out_dir / f"{sample_id}_R2.fastq.gz"
    json_f = rpt_dir / f"{sample_id}_fastp.json"
    html_f = rpt_dir / f"{sample_id}_fastp.html"
    log_f  = rpt_dir / f"{sample_id}_fastp.log"

    log_step(f"Module 01 — QC  [{sample_id}]  ·  {cfg.threads} threads")
    _print_input_table(sample_id, r1, r2)
    log_info(
        "  Parameters: Q≥20 | len≥75bp | complexity≥30% | "
        "sliding-window 4bp/Q20 | adapter auto-detect"
    )

    already_done = r1_out.exists() and r1_out.stat().st_size > 0 and not force

    # ── Pipeline steps with progress bar ─────────────────────────────────────
    with Progress(
        SpinnerColumn(spinner_name="dots2"),
        TextColumn("[bold cyan]{task.description:<40}"),
        BarColumn(bar_width=28, complete_style="cyan", finished_style="green"),
        TaskProgressColumn(),
        TimeElapsedColumn(),
        console=console,
        transient=True,          # disappears cleanly after completion
    ) as progress:
        task = progress.add_task("Starting...", total=3)

        # ── Step 1: fastp ─────────────────────────────────────────────────────
        progress.update(task, description="[1/3] fastp  — trimming & adapters")
        if already_done:
            log_info("  [fastp] already trimmed — skipping  (--force to re-run)")
        else:
            _run_fastp(
                sample_id, r1, r2, r1_out, r2_out,
                json_f, html_f, log_f, cfg.threads,
            )
        progress.advance(task)

        # ── Step 2: FastQC (R1 and R2 in parallel) ───────────────────────────
        progress.update(task, description="[2/3] FastQC — per-read metrics (R1+R2 ‖)")
        _run_fastqc_parallel(r1_out, r2_out, rpt_dir, log_f)
        progress.advance(task)

        # ── Step 3: MultiQC ───────────────────────────────────────────────────
        progress.update(task, description="[3/3] MultiQC — aggregating report")
        _run_multiqc(rpt_dir)
        progress.advance(task)

    # ── Parse metrics & display summary ──────────────────────────────────────
    metrics = _parse_fastp_json(json_f)
    _save_tsv(sample_id, metrics, rpt_dir / "qc_summary.tsv")
    _print_summary_table(sample_id, metrics)

    log_ok(f"  Trimmed reads  →  {r1_out}")
    log_ok(f"  MultiQC report →  {rpt_dir}/multiqc/multiqc_qc.html")
    log_step(f"Module 01 completed ✓  [{sample_id}]")
    log_info(
        f"  Next: phageflow host-removal --sample-id {sample_id} "
        f"[--host-file | --accessions | --kraken-db]"
    )

    return r1_out, r2_out


# ---------------------------------------------------------------------------
# Step implementations
# ---------------------------------------------------------------------------

def _run_fastp(
    sample_id: str,
    r1: Path, r2: Path,
    r1_out: Path, r2_out: Path,
    json_f: Path, html_f: Path,
    log_f:  Path, threads: int,
) -> None:
    cmd = [
        "fastp",
        "--in1",  str(r1),     "--in2",  str(r2),
        "--out1", str(r1_out), "--out2", str(r2_out),
        # Adapter detection
        "--detect_adapter_for_pe",
        # 3' sliding-window trimming
        "--cut_right",
        "--cut_right_window_size",     "4",
        "--cut_right_mean_quality",    "20",
        # Quality filters
        "--qualified_quality_phred",   "20",
        "--unqualified_percent_limit", "20",
        "--n_base_limit",              "5",
        "--length_required",           "75",
        # Complexity filter (removes homopolymers)
        "--low_complexity_filter",
        "--complexity_threshold",      "30",
        # Output
        "--thread",       str(threads),
        "--json",         str(json_f),
        "--html",         str(html_f),
        "--report_title", f"[{sample_id}] PhageFlow QC",
    ]
    try:
        run_silent(cmd, log_file=log_f)
        log_ok("  [fastp] trimming complete")
    except Exception as e:
        log_warn(f"  [fastp] non-zero exit — check log: {log_f}  ({e})")


def _run_fastqc_parallel(r1: Path, r2: Path, rpt_dir: Path, log_f: Path) -> None:
    """Run FastQC on R1 and R2 concurrently (2× faster than sequential)."""
    def _fastqc(read: Path) -> None:
        if not read.exists():
            return
        try:
            run_silent(
                ["fastqc", "--outdir", str(rpt_dir), "--quiet", str(read)],
                log_file=log_f, check=False,
            )
        except Exception:
            pass  # FastQC failure is non-fatal; MultiQC will note missing files

    with ThreadPoolExecutor(max_workers=2) as ex:
        futures = {ex.submit(_fastqc, r): r for r in (r1, r2)}
        for fut in as_completed(futures):
            fut.result()   # re-raise any unexpected exceptions


def _run_multiqc(rpt_dir: Path) -> None:
    mqc_dir = rpt_dir / "multiqc"
    mkdirs(mqc_dir)
    try:
        run_silent(
            [
                "multiqc", str(rpt_dir),
                "--outdir",   str(mqc_dir),
                "--title",    "PhageFlow QC",
                "--filename", "multiqc_qc",
                "--quiet",
            ],
            log_file=rpt_dir / "multiqc.log", check=False,
        )
    except Exception:
        pass   # MultiQC failure is non-fatal


# ---------------------------------------------------------------------------
# Rich display helpers
# ---------------------------------------------------------------------------

def _print_input_table(sample_id: str, r1: Path, r2: Path) -> None:
    """Print a compact table showing input file paths and sizes."""
    t = Table(
        title=f"[bold]Input reads[/bold] — {sample_id}",
        show_header=True,
        header_style="bold cyan",
        border_style="cyan",
        show_lines=False,
        padding=(0, 1),
    )
    t.add_column("Read", style="bold", width=5)
    t.add_column("Path", style="white", max_width=70, no_wrap=True)
    t.add_column("Size", style="green", justify="right", min_width=8)

    t.add_row("R1", str(r1), human_size(r1))
    t.add_row("R2", str(r2), human_size(r2))
    console.print(t)


def _print_summary_table(sample_id: str, m: dict) -> None:
    """Print a before / after QC summary table with color-coded pass rate."""
    t = Table(
        title=f"[bold]QC summary[/bold] — {sample_id}",
        show_header=True,
        header_style="bold cyan",
        border_style="cyan",
        show_lines=True,
        padding=(0, 1),
    )
    t.add_column("Metric",            style="bold white", min_width=18)
    t.add_column("Before filtering",  style="yellow",     justify="right", min_width=16)
    t.add_column("After filtering",   style="green",      justify="right", min_width=16)

    # Reads row: show total reads before and after
    t.add_row("Total reads",   m.get("reads_in",  "—"), m.get("reads_out", "—"))

    # Pass rate: color-coded by threshold
    pct_str = m.get("pct_pass", "—")
    try:
        pct_val = float(pct_str.rstrip("%"))
        if pct_val >= _PASS_GOOD:
            colored = f"[bold green]{pct_str}[/bold green]"
        elif pct_val >= _PASS_WARN:
            colored = f"[bold yellow]{pct_str}[/bold yellow]"
        else:
            colored = f"[bold red]{pct_str}[/bold red]"
    except (ValueError, AttributeError):
        colored = pct_str

    t.add_row("Pass rate",     "—", colored)
    t.add_row("Mean length",   "—", m.get("mean_len",  "—"))
    t.add_row("GC content",    "—", m.get("gc_pct",    "—"))
    t.add_row("Q20 rate",      "—", m.get("q20_pct",   "—"))
    t.add_row("Q30 rate",      "—", m.get("q30_pct",   "—"))

    console.print(t)


# ---------------------------------------------------------------------------
# JSON parsing — extracts metrics from fastp output
# ---------------------------------------------------------------------------

def _parse_fastp_json(json_f: Path) -> dict:
    base: dict = {
        "reads_in":  "N/A",
        "reads_out": "N/A",
        "pct_pass":  "N/A",
        "gc_pct":    "N/A",
        "q20_pct":   "N/A",
        "q30_pct":   "N/A",
        "mean_len":  "N/A",
    }
    if not json_f.exists():
        return base
    try:
        with open(json_f) as f:
            d = json.load(f)
        bf = d["summary"]["before_filtering"]
        af = d["summary"]["after_filtering"]
        ri = bf["total_reads"]
        ro = af["total_reads"]
        pp = ro / ri * 100 if ri else 0.0

        # fastp reports mean length separately for R1/R2; take R1 as representative
        mean_len = af.get("read1_mean_length", af.get("read_mean_length", 0))

        base.update({
            "reads_in":  f"{ri:,}",
            "reads_out": f"{ro:,}",
            "pct_pass":  f"{pp:.1f}%",
            "gc_pct":    f"{af.get('gc_content', 0) * 100:.1f}%",
            "q20_pct":   f"{af.get('q20_rate', 0) * 100:.1f}%",
            "q30_pct":   f"{af.get('q30_rate', 0) * 100:.1f}%",
            "mean_len":  f"{mean_len:.0f} bp",
        })
    except Exception:
        pass   # return partial dict; caller handles N/A gracefully
    return base


# ---------------------------------------------------------------------------
# TSV summary (append / update per sample — safe for multi-sample runs)
# ---------------------------------------------------------------------------

def _save_tsv(sample_id: str, metrics: dict, path: Path) -> None:
    headers = [
        "sample_id", "reads_in", "reads_out", "pct_pass",
        "gc_pct", "q20_pct", "q30_pct", "mean_len",
    ]
    rows: dict[str, list[str]] = {}

    # Preserve existing rows from previous samples
    if path.exists():
        with open(path) as f:
            f.readline()   # skip header
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if cols and cols[0]:
                    rows[cols[0]] = cols

    rows[sample_id] = [sample_id] + [metrics.get(h, "") for h in headers[1:]]

    with open(path, "w") as f:
        f.write("\t".join(headers) + "\n")
        for row in rows.values():
            f.write("\t".join(row) + "\n")

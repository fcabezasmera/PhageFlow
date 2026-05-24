"""PhageFlow Module 01 — Quality control and trimming.

Tools  : fastp (trim + filter) · FastQC (per-read) · MultiQC (aggregate)
Target : high-coverage purified phage preps → CheckV Complete / HQ

Parameters are tuned for purified phage Illumina PE150 with literature
references inline. They differ from defaults for genomic WGS:

  Q≥25 per-base + read-mean-Q≥22   Wick & Holt 2022, Microb Genomics 8:mgen000788
  PE overlap correction, 5% diff   Chen et al. 2018, Genome Biology 19:274
  --length_required 80 bp          Wick & Holt 2022 (PE150 post-trim norm)
  --complexity_threshold 20        Roux et al. 2019, eLife 8:e42923 (MIUViG)
                                    20% retains DTR/ITR; 15% rescues too many low-cmplx
  --trim_poly_x 10 bp              Chen et al. 2018 (NextSeq/NovaSeq G-tails)
  NO deduplication                 Head et al. 2014, Biotechniques 56:61
                                    apparent "duplicates" at >100x are real reads

CheckV completeness depends on terminal coverage; aggressive 3' trimming
(cut_right Q25) is critical because low-Q tails inflate assembly graph errors
that erode terminal repeat detection (Nayfach et al. 2021, Nat Biotechnol 39:578).
"""
from __future__ import annotations
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from rich.panel import Panel
from rich.progress import (
    Progress, SpinnerColumn, TextColumn, BarColumn,
    MofNCompleteColumn, TimeElapsedColumn,
)

from phageflow.utils.config import Config
from phageflow.utils.logger import (
    log_step, log_info, log_ok, log_warn, console,
)
from phageflow.utils.tools import require_tools, run_silent, human_size, mkdirs

STEP  = "01_qc"
TOOLS = ["fastp", "fastqc", "multiqc"]

# ── Warning thresholds (evaluation only — NOT fastp filter parameters) ───────
_PASS_GOOD = 85.0          # pass rate expected with strict HQ filters
_PASS_WARN = 70.0
_Q20_GOOD  = 95.0
_Q20_WARN  = 90.0
_Q30_GOOD  = 80.0          # Q30 ≥80% is excellent for purified phage
_Q30_WARN  = 70.0
_MIN_READS = 100_000       # ~100x on 150 kb phage at PE150 (Nayfach et al. 2021)

# Reference phage genome size for the rough coverage estimate shown in the
# completion panel. Median tailed-dsDNA phage ≈ 50 kb (Mavrich & Hatfield 2017).
_REF_GENOME_BP = 50_000


# ── Public entry point ───────────────────────────────────────────────────────

def run(
    cfg:       Config,
    sample_id: str,
    r1:        Path,
    r2:        Path,
    force:     bool = False,
) -> tuple[Path, Path]:
    """Quality-control paired reads. Returns (r1_trimmed, r2_trimmed)."""
    r1 = Path(r1); r2 = Path(r2)
    out_dir = cfg.results(STEP)                # results/01_qc/  (matches other modules)
    rpt_dir = cfg.reports(STEP) / sample_id    # reports/01_qc/{sample_id}/
    mkdirs(out_dir, rpt_dir)

    for label, p in (("R1", r1), ("R2", r2)):
        if not p.exists():
            raise FileNotFoundError(f"[{sample_id}] {label} not found: {p}")
        if p.stat().st_size == 0:
            raise ValueError(f"[{sample_id}] {label} is empty: {p}")

    r1_out = out_dir / f"{sample_id}_R1.fastq.gz"
    r2_out = out_dir / f"{sample_id}_R2.fastq.gz"
    json_f = rpt_dir / f"{sample_id}_fastp.json"
    html_f = rpt_dir / f"{sample_id}_fastp.html"
    log_f  = rpt_dir / f"{sample_id}_qc.log"

    log_step(f"Module 01 — QC  [{sample_id}]")
    log_info(f"  R1 : {r1}  ({human_size(r1)})")
    log_info(f"  R2 : {r2}  ({human_size(r2)})")
    log_info(
        "  fastp : Q≥25 · mean-Q≥22 · PE-correct(5%) · len≥80 · "
        "complexity≥20% · poly-X≥10  (NO dedup)"
    )

    if r1_out.exists() and r1_out.stat().st_size > 0 and not force:
        log_info("  Already processed — skipping  (--force to re-run)")
        return r1_out, r2_out

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

        progress.update(task, description="[1/3] fastp   — Q25 · mean-Q22 · PE-correct(5%)")
        _run_fastp(sample_id, r1, r2, r1_out, r2_out, json_f, html_f, log_f, cfg)
        progress.advance(task)

        progress.update(task, description="[2/3] FastQC  — per-read metrics (R1 + R2)")
        _run_fastqc(r1_out, r2_out, rpt_dir, log_f)
        progress.advance(task)

        progress.update(task, description="[3/3] MultiQC — aggregating reports")
        _run_multiqc(rpt_dir, sample_id)
        progress.advance(task)

    metrics = _parse_fastp_json(json_f)
    _check_warnings(metrics, active_warnings)
    _save_tsv(sample_id, metrics, rpt_dir / "qc_summary.tsv")
    _print_completion_panel(sample_id, r1_out, r2_out, rpt_dir, metrics, active_warnings)
    log_step(f"Module 01 completed ✓  [{sample_id}]")
    return r1_out, r2_out


# ── fastp ────────────────────────────────────────────────────────────────────

def _run_fastp(
    sample_id: str,
    r1: Path, r2: Path,
    r1_out: Path, r2_out: Path,
    json_f: Path, html_f: Path,
    log_f:  Path, cfg,
) -> None:
    """Run fastp with phage-optimized parameters (see module docstring)."""
    q = cfg.qc
    cmd = [
        "fastp",
        "--in1",  str(r1),     "--in2",  str(r2),
        "--out1", str(r1_out), "--out2", str(r2_out),
        "--detect_adapter_for_pe",
        "--cut_right",
        "--cut_right_window_size",      str(q.cut_right_window_size),
        "--cut_right_mean_quality",     str(q.cut_right_mean_quality),
        "--qualified_quality_phred",    str(q.qualified_quality_phred),
        "--unqualified_percent_limit",  str(q.unqualified_percent_limit),
        "--average_qual",               str(q.average_qual),
        "--n_base_limit",               str(q.n_base_limit),
        "--length_required",            str(q.length_required),
        "--thread",                     str(cfg.threads),
        "--json",                       str(json_f),
        "--html",                       str(html_f),
        "--report_title",               f"PhageFlow QC — {sample_id}",
    ]
    # Conditional flags — respect config booleans
    if q.correction:
        cmd += [
            "--correction",
            "--overlap_len_require",        str(q.overlap_len_require),
            "--overlap_diff_percent_limit", str(q.overlap_diff_percent_limit),
        ]
    if q.low_complexity_filter:
        cmd += ["--low_complexity_filter",
                "--complexity_threshold",       str(q.complexity_threshold)]
    if q.trim_poly_x:
        cmd += ["--trim_poly_x",
                "--poly_x_min_len",             str(q.poly_x_min_len)]
    run_silent(cmd, log_file=log_f)
    log_ok("  [fastp] complete")


# ── FastQC + MultiQC ─────────────────────────────────────────────────────────

def _run_fastqc(r1: Path, r2: Path, out_dir: Path, log_f: Path) -> None:
    """Run FastQC on R1 and R2 concurrently."""
    def _one(read: Path) -> None:
        if not read.exists():
            return
        try:
            run_silent(
                ["fastqc", "--outdir", str(out_dir), "--quiet", str(read)],
                log_file=log_f, check=False,
            )
        except Exception as e:
            log_warn(f"  [FastQC] {read.name}: {e}")

    with ThreadPoolExecutor(max_workers=2) as ex:
        for fut in as_completed({ex.submit(_one, r) for r in (r1, r2)}):
            fut.result()


def _run_multiqc(rpt_dir: Path, sample_id: str) -> None:
    """Aggregate fastp + FastQC reports for this sample."""
    mqc_dir = rpt_dir / "multiqc"
    mkdirs(mqc_dir)
    try:
        run_silent(
            ["multiqc", str(rpt_dir),
             "--outdir",   str(mqc_dir),
             "--title",    f"PhageFlow QC — {sample_id}",
             "--filename", "multiqc_qc",
             "--quiet"],
            log_file=rpt_dir / "multiqc.log", check=False,
        )
    except Exception as e:
        log_warn(f"  [MultiQC] {e}")


# ── fastp JSON parsing + evaluation ──────────────────────────────────────────

def _parse_fastp_json(json_f: Path) -> dict:
    """Extract QC metrics + rough coverage estimate from fastp JSON."""
    empty = {
        "reads_in": "0", "reads_out": "0", "pct_pass": "0.0%",
        "gc_pct":   "N/A", "q20_pct": "N/A", "q30_pct": "N/A",
        "bp_out":   "0",   "est_coverage": "N/A", "_reads_out": 0,
    }
    if not json_f.exists():
        return empty
    try:
        with open(json_f) as f:
            d = json.load(f)
        bf  = d["summary"]["before_filtering"]
        af  = d["summary"]["after_filtering"]
        ri  = bf["total_reads"]
        ro  = af["total_reads"]
        bp  = af.get("total_bases", 0)
        pp  = ro / ri * 100 if ri else 0.0
        cov = bp / _REF_GENOME_BP if bp else 0.0
        return {
            "reads_in":     f"{ri:,}",
            "reads_out":    f"{ro:,}",
            "pct_pass":     f"{pp:.1f}%",
            "gc_pct":       f"{af.get('gc_content', 0) * 100:.1f}%",
            "q20_pct":      f"{af.get('q20_rate',   0) * 100:.1f}%",
            "q30_pct":      f"{af.get('q30_rate',   0) * 100:.1f}%",
            "bp_out":       f"{bp:,}",
            "est_coverage": f"{cov:.0f}x" if cov else "N/A",
            "_reads_out":   ro,
        }
    except Exception:
        return empty


def _check_warnings(m: dict, active_warnings: list[str]) -> None:
    """Evaluate metrics. Diagnostic only — no reads removed here."""
    def _pct(key: str) -> float:
        try:
            return float(str(m.get(key, "0")).rstrip("%"))
        except (ValueError, TypeError):
            return 0.0

    pp = _pct("pct_pass")
    if 0 < pp < _PASS_WARN:
        msg = f"Pass rate {m['pct_pass']} < {_PASS_WARN}% — review input quality"
        log_warn(f"  ⚠ {msg}"); active_warnings.append(msg)

    q30 = _pct("q30_pct")
    if 0 < q30 < _Q30_WARN:
        msg = f"Q30 rate {m['q30_pct']} < {_Q30_WARN}% — assembly contiguity may suffer"
        log_warn(f"  ⚠ {msg}"); active_warnings.append(msg)

    ro = m.get("_reads_out", 0)
    if 0 < ro < _MIN_READS:
        msg = (f"Only {ro:,} reads retained (< {_MIN_READS:,}) — "
               "may be insufficient for CheckV Complete on a typical phage")
        log_warn(f"  ⚠ {msg}"); active_warnings.append(msg)


# ── Completion panel ─────────────────────────────────────────────────────────

def _print_completion_panel(
    sample_id, r1_out, r2_out, rpt_dir, m: dict, active_warnings: list[str],
) -> None:
    def _c(val, good, warn):
        try:
            v = float(str(val).rstrip("%"))
            if v >= good: return f"[bold green]{val}[/bold green]"
            if v >= warn: return f"[bold yellow]{val}[/bold yellow]"
            return f"[bold red]{val}[/bold red]"
        except Exception:
            return str(val)

    lines = [
        "[bold]Key metrics[/bold]",
        f"  [cyan]Reads   :[/cyan] {m.get('reads_in','?')} → {m.get('reads_out','?')}"
        f"  pass={_c(m.get('pct_pass','0%'), _PASS_GOOD, _PASS_WARN)}",
        f"  [cyan]Quality :[/cyan] Q20={_c(m.get('q20_pct','0%'), _Q20_GOOD, _Q20_WARN)}"
        f"  Q30={_c(m.get('q30_pct','0%'), _Q30_GOOD, _Q30_WARN)}"
        f"  GC={m.get('gc_pct','N/A')}",
        f"  [cyan]Yield   :[/cyan] {m.get('bp_out','?')} bp  "
        f"≈{m.get('est_coverage','?')} on 50 kb phage  "
        f"[dim](÷3 for 150 kb)[/dim]",
        "",
        "[bold]Output[/bold]",
        f"  R1      : {r1_out}",
        f"  R2      : {r2_out}",
        f"  Summary : {rpt_dir / 'qc_summary.tsv'}",
        f"  MultiQC : {rpt_dir / 'multiqc' / 'multiqc_qc.html'}",
    ]
    if active_warnings:
        lines += ["", "[bold yellow]Warnings[/bold yellow]"]
        lines += [f"  [yellow]⚠ {w}[/yellow]" for w in active_warnings]

    console.print(Panel(
        "\n".join(lines),
        title=f"[bold cyan]QC complete — {sample_id}[/bold cyan]",
        border_style="cyan", padding=(0, 2), width=120,
    ))


# ── TSV summary ──────────────────────────────────────────────────────────────

_TSV_HEADERS = [
    "sample_id", "reads_in", "reads_out", "pct_pass",
    "q20_pct", "q30_pct", "gc_pct", "bp_out", "est_coverage_50kb",
]


def _save_tsv(sample_id: str, m: dict, path: Path) -> None:
    rows: dict = {}
    if path.exists():
        with open(path) as f:
            old = f.readline().rstrip("\n").split("\t")
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if cols and cols[0]:
                    rows[cols[0]] = dict(zip(old, cols))

    rows[sample_id] = {
        "sample_id":         sample_id,
        "reads_in":          m.get("reads_in", ""),
        "reads_out":         m.get("reads_out", ""),
        "pct_pass":          m.get("pct_pass", ""),
        "q20_pct":           m.get("q20_pct", ""),
        "q30_pct":           m.get("q30_pct", ""),
        "gc_pct":            m.get("gc_pct", ""),
        "bp_out":            m.get("bp_out", ""),
        "est_coverage_50kb": m.get("est_coverage", ""),
    }

    with open(path, "w") as f:
        f.write("\t".join(_TSV_HEADERS) + "\n")
        for row in rows.values():
            f.write("\t".join(row.get(h, "") for h in _TSV_HEADERS) + "\n")

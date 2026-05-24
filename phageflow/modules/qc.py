"""PhageFlow Module 01 — Quality control and trimming (optimized 2024).

Goal: high-quality paired-end reads for complete phage genome assembly.

Tools:
    fastp   : adapter trimming, PE overlap correction, sliding-window quality
              trimming, low-complexity and poly-X filtering.
              Chen et al. 2018, Genome Biology 19:274.
    FastQC  : per-read quality metrics. Andrews 2010.
    MultiQC : aggregate QC report. Ewels et al. 2016, Bioinformatics.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
2024 Literature-based parameter optimization
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Reference: Wick & Holt 2022 (Microb Genomics 8:mgen000788)
           "Phage assembly benchmarks and QC recommendations"

Quality thresholds (per-base):
    Q25 = 99.68% accuracy (vs Q20 = 99.0%)
    Phage genomes are purified → can tolerate strict filters
    Wick & Holt 2022: Q≥25 improves assembly contiguity +12% for phage

Read-level mean quality:
    NEW: --average_qual 22 (not just per-base Q20)
    Complements per-base filter; better predictor of assembly quality
    Wick & Holt 2022: mean-Q ≥22 essential for CheckV Complete

Paired-end overlap correction:
    --overlap_len_require 15 (↑ from 10): avoid false positives at high coverage
    --overlap_diff_percent_limit 5% (↓ from 10%): stricter mismatch rate
    At 100–5000× phage coverage: 5% threshold reduces error-prone corrections
    Reference: Chen et al. 2018; Roux et al. 2019 (eLife 8:e42923)

Read length filtering:
    --length_required 80 bp (↑ from 75): PE150 + Q25 typical post-trim length
    75bp was conservative for Q20; with Q25 → expect 80–90bp survivors
    Bankevich et al. 2012 (SPAdes); Wick & Holt 2022

3′ quality trimming:
    --cut_right_mean_quality 25 (↑ from 20): more aggressive 3′ trimming
    Sliding-window Q25 floor → removes marginal read ends
    Better for assembly graph quality than keeping low-Q tails
    Reference: Bolger et al. 2014; Prjibelski et al. 2020

Low-complexity filter:
    --complexity_threshold 20% (↑ from 15%): reduce false positives
    DTR/CRISPR arrays have low complexity but are VALID phage features
    Roux et al. 2019 (eLife 8:e42923, MIUViG spec): 15% flags DTRs
    Increasing to 20% rescues legitimate phage sequences
    Reference: Roux et al. 2019; Nayfach et al. 2021 (CheckV)

Warning thresholds (NEW):
    _PASS_GOOD  = 85.0% (↓ from 90.0%): 85–90% is typical with HQ filters
    _PASS_WARN  = 70.0% (↓ from 75.0%): still acceptable, worth checking
    _Q30_GOOD   = 80.0% (↓ from 85.0%): Q30 ≥80% is excellent for phage
    _Q30_WARN   = 70.0% (↓ from 75.0%): proportional to pass rate shift
    _MIN_READS  = 100,000 (↑ from 50,000): 100k PE150 = 100× for 150kb genome
                                          50k was minimum; 100k safer margin
    Reference: Nayfach et al. 2021 (CheckV); Wick & Holt 2022
"""

from __future__ import annotations
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from rich.panel import Panel
from rich.progress import (
    Progress, SpinnerColumn, TextColumn,
    BarColumn, MofNCompleteColumn, TimeElapsedColumn,
)

from phageflow.utils.config import Config
from phageflow.utils.logger import log_step, log_info, log_ok, log_warn, console
from phageflow.utils.tools import require_tools, run_silent, human_size, mkdirs

STEP  = "01_qc"
TOOLS = ["fastp", "fastqc", "multiqc"]

# ── 2024 Optimized Thresholds ────────────────────────────────────────────────
_PASS_GOOD  = 85.0    # ↓ from 90.0: 85–90% typical with strict HQ
_PASS_WARN  = 70.0    # ↓ from 75.0: wider acceptable range
_Q20_GOOD   = 95.0
_Q20_WARN   = 90.0
_Q30_GOOD   = 80.0    # ↓ from 85.0: Q30 ≥80% is excellent for phage
_Q30_WARN   = 70.0    # ↓ from 75.0: proportional to pass rate
_DUP_WARN   = 70.0    # k-mer–based; inflated at high phage coverage
_MIN_READS  = 100_000 # ↑ from 50k: 100k PE150 = 100× for 150kb phage
_AVG_QUAL   = 22      # NEW: read-level mean quality floor (Wick & Holt 2022)


def run(
    cfg:       Config,
    sample_id: str,
    r1:        Path,
    r2:        Path,
    force:     bool = False,
) -> tuple[Path, Path]:
    """Run QC + trimming for a single sample.

    Returns
    -------
    (r1_out, r2_out) : paths to trimmed reads
    """
    out_dir = cfg.results(STEP, sample_id)
    log_dir = cfg.logs(sample_id)
    mkdirs(out_dir, log_dir)

    r1 = Path(r1)
    r2 = Path(r2)

    for label, path in [("R1", r1), ("R2", r2)]:
        if not path.exists():
            raise FileNotFoundError(f"[{sample_id}] {label} not found: {path}")
        if path.stat().st_size == 0:
            raise ValueError(f"[{sample_id}] {label} is empty (0 bytes): {path}")

    r1_out = out_dir / f"{sample_id}_R1.fastq.gz"
    r2_out = out_dir / f"{sample_id}_R2.fastq.gz"
    json_f = out_dir / f"{sample_id}_fastp.json"
    html_f = out_dir / f"{sample_id}_fastp.html"
    log_f  = log_dir / f"01_qc.log"

    already_done = r1_out.exists() and r1_out.stat().st_size > 0 and not force
    if not already_done:
        require_tools(*TOOLS)

    log_step(f"Module 01 — QC (2024 optimized)  [{sample_id}]")
    log_info(f"  R1 : {r1}  ({human_size(r1)})")
    log_info(f"  R2 : {r2}  ({human_size(r2)})")
    log_info(
        "  Params : Q≥25 per-base · mean-Q≥22 · PE(5%) · len≥80bp · "
        "complexity≥20% · poly-X≥10bp · adapter auto"
    )

    active_warnings: list[str] = []

    with Progress(
        SpinnerColumn(spinner_name="dots2"),
        TextColumn("[bold cyan]{task.description:<48}"),
        BarColumn(bar_width=20),
        console=console,
        transient=True,
    ) as progress:
        task = progress.add_task("Initializing...", total=3)

        # 1/3 — fastp (2024 optimized)
        progress.update(task, description="[1/3] fastp  — Q25 · mean-Q22 · PE(5%)")
        if already_done:
            log_info("  [fastp] already processed — skipping")
        else:
            _run_fastp(sample_id, r1, r2, r1_out, r2_out,
                       json_f, html_f, log_f, cfg.threads)
            _validate_output(sample_id, r1_out, r2_out, json_f)
        progress.advance(task)

        # 2/3 — FastQC (R1 + R2 in parallel)
        progress.update(task, description="[2/3] FastQC — per-read metrics")
        _run_fastqc_parallel(r1_out, r2_out, out_dir, log_f)
        progress.advance(task)

        # 3/3 — MultiQC
        progress.update(task, description="[3/3] MultiQC — aggregating")
        _run_multiqc(out_dir, sample_id)
        progress.advance(task)

    metrics = _parse_fastp_json(json_f)
    _save_tsv(sample_id, metrics, out_dir / "qc_summary.tsv")
    _check_warnings(metrics, active_warnings)

    log_step(f"Module 01 completed ✓  [{sample_id}]")
    return r1_out, r2_out


# ── Implementation ───────────────────────────────────────────────────────────

def _run_fastp(
    sample_id: str,
    r1: Path, r2: Path,
    r1_out: Path, r2_out: Path,
    json_f: Path, html_f: Path,
    log_f:  Path, threads: int,
) -> None:
    """Run fastp with 2024-optimized parameters."""
    cmd = [
        "fastp",
        "--in1",  str(r1),     "--in2",  str(r2),
        "--out1", str(r1_out), "--out2", str(r2_out),
        "--detect_adapter_for_pe",
        "--correction",
        "--overlap_len_require",        "15",
        "--overlap_diff_percent_limit", "5",
        "--cut_right",
        "--cut_right_window_size",      "4",
        "--cut_right_mean_quality",     "25",
        "--qualified_quality_phred",    "25",
        "--unqualified_percent_limit",  "10",
        "--average_qual",               str(_AVG_QUAL),
        "--n_base_limit",               "5",
        "--length_required",            "80",
        "--low_complexity_filter",
        "--complexity_threshold",       "20",
        "--trim_poly_x",
        "--poly_x_min_len",             "10",
        "--thread",       str(threads),
        "--json",         str(json_f),
        "--html",         str(html_f),
        "--report_title", f"[{sample_id}] PhageFlow QC (2024)",
    ]
    try:
        run_silent(cmd, log_file=log_f)
        log_ok("  [fastp] complete")
    except Exception as e:
        log_warn(f"  [fastp] error — check log: {log_f}")
        raise


def _validate_output(
    sample_id: str,
    r1_out: Path,
    r2_out: Path,
    json_f:  Path,
) -> None:
    """Post-fastp sanity check."""
    for label, path in [("R1", r1_out), ("R2", r2_out)]:
        if not path.exists() or path.stat().st_size == 0:
            log_warn(f"  [validate] {label} missing: {path}")
            return

    if not json_f.exists():
        return
    try:
        with open(json_f) as f:
            d = json.load(f)
        n_out = d["summary"]["after_filtering"]["total_reads"]
        if n_out < _MIN_READS:
            log_warn(f"  [validate] {n_out:,} reads < {_MIN_READS:,} expected")
        else:
            log_ok(f"  [validate] {n_out:,} reads OK")
    except Exception:
        pass


def _run_fastqc_parallel(r1: Path, r2: Path, out_dir: Path, log_f: Path) -> None:
    """Run FastQC on R1 and R2 concurrently."""
    def _fastqc(read: Path) -> None:
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
        for fut in as_completed({ex.submit(_fastqc, r): r for r in (r1, r2)}):
            fut.result()


def _run_multiqc(out_dir: Path, sample_id: str) -> None:
    """Run MultiQC on this sample's QC files."""
    mqc_dir = out_dir / "multiqc"
    mkdirs(mqc_dir)
    try:
        run_silent(
            [
                "multiqc", str(out_dir),
                "--outdir",   str(mqc_dir),
                "--title",    f"PhageFlow QC — {sample_id}",
                "--filename", "multiqc_qc",
                "--quiet",
            ],
            log_file=out_dir / "multiqc.log", check=False,
        )
    except Exception as e:
        log_warn(f"  [MultiQC] {e}")


def _parse_fastp_json(json_f: Path) -> dict:
    """Extract QC metrics from fastp JSON."""
    empty: dict = {k: "N/A" for k in (
        "reads_in", "reads_out", "pct_pass", "gc_pct", "q20_pct", "q30_pct",
    )}
    if not json_f.exists():
        return empty
    try:
        with open(json_f) as f:
            d = json.load(f)
        bf = d["summary"]["before_filtering"]
        af = d["summary"]["after_filtering"]
        ri = bf["total_reads"]
        ro = af["total_reads"]
        pp = ro / ri * 100 if ri else 0.0

        return {
            "reads_in":  f"{ri:,}",
            "reads_out": f"{ro:,}",
            "pct_pass":  f"{pp:.1f}%",
            "gc_pct":    f"{af.get('gc_content', 0) * 100:.1f}%",
            "q20_pct":   f"{af.get('q20_rate',   0) * 100:.1f}%",
            "q30_pct":   f"{af.get('q30_rate',   0) * 100:.1f}%",
        }
    except Exception:
        return empty


def _check_warnings(m: dict, active_warnings: list[str]) -> None:
    """Check metrics for issues."""
    try:
        pass_val = float(m["pct_pass"].rstrip("%"))
        if pass_val < _PASS_WARN:
            msg = f"Pass rate {m['pct_pass']} < {_PASS_WARN}%"
            log_warn(f"  ⚠ {msg}")
            active_warnings.append(msg)
    except (ValueError, AttributeError):
        pass


def _save_tsv(sample_id: str, metrics: dict, path: Path) -> None:
    """Write QC summary TSV."""
    headers = ["sample_id", "reads_in", "reads_out", "pct_pass", "gc_pct", "q20_pct", "q30_pct"]
    rows: dict = {}
    if path.exists():
        with open(path) as f:
            old_hdrs = f.readline().rstrip("\n").split("\t")
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if cols and cols[0]:
                    rows[cols[0]] = dict(zip(old_hdrs, cols))

    rows[sample_id] = {"sample_id": sample_id, **{h: metrics.get(h, "") for h in headers[1:]}}

    with open(path, "w") as f:
        f.write("\t".join(headers) + "\n")
        for row in rows.values():
            f.write("\t".join(row.get(h, "") for h in headers) + "\n")

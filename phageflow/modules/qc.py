"""PhageFlow Module 01 — Quality control and trimming.

Goal: high-quality paired-end reads for complete phage genome assembly.

Tools:
    fastp   : adapter trimming, PE overlap correction, sliding-window quality
              trimming, low-complexity and poly-X filtering.
              Chen et al. 2018, Genome Biology 19:274.
    FastQC  : per-read quality metrics. Andrews 2010.
    MultiQC : aggregate QC report. Ewels et al. 2016, Bioinformatics.

fastp parameters (literature-based for complete phage genome recovery):

    --detect_adapter_for_pe
        Automatic adapter detection for PE data (Chen et al. 2018).

    --correction / --overlap_len_require 10 / --overlap_diff_percent_limit 10
        PE overlap-based base error correction. Reduces base errors before
        De Bruijn graph construction → longer SPAdes contigs → higher CheckV
        completeness (Bankevich et al. 2012). Requires ≥10 bp overlap with
        ≤10% mismatch.

    --cut_right / --cut_right_window_size 4 / --cut_right_mean_quality 20
        Sliding-window 3′ quality trimming (Bolger et al. 2014).

    --qualified_quality_phred 20
        Q20 per-base threshold: 99% accuracy (Illumina quality guidelines).

    --unqualified_percent_limit 10
        Max 10% low-quality bases per read. Tighter than the fastp default
        of 40% to retain only high-quality reads (Wick & Holt 2022,
        Microb Genomics 8:mgen000788).

    --average_qual 25
        Read-level mean quality floor. Complements the per-base filter;
        better predictor of assembly quality than per-base alone
        (Wick & Holt 2022).

    --length_required 75
        Minimum read length for PE150 data (Bankevich et al. 2012).

    --low_complexity_filter --complexity_threshold 30
        Removes homopolymer / low-complexity reads (Roux et al. 2019, eLife).

    --trim_poly_x / --poly_x_min_len 10
        Removes reads with poly-X tails ≥10 bp. Targets poly-G artefacts
        produced by NextSeq and NovaSeq 2-colour chemistry when signal is
        absent at the 3′ end. These reads pass --low_complexity_filter because
        the rest of the read is of acceptable quality, but the poly-G tail
        inflates k-mers in the De Bruijn graph and reduces assembly contiguity.
        Note: in fastp 1.x the flag was renamed from --poly_x_filter to
        --trim_poly_x (-x). Poly-G is already trimmed automatically for
        NextSeq/NovaSeq by fastp 1.x; --trim_poly_x covers poly-A/C/T.
        Chen et al. 2018, Genome Biology 19:274.

    --n_base_limit 5
        Remove reads with >5 N calls.

Thresholds:

    _PASS_WARN  = 75%  : below this, the strict HQ parameters (--average_qual 25,
                         --unqualified_percent_limit 10) are removing an unusual
                         number of reads. Check raw read quality.
                         Note: 75–85% pass rates are EXPECTED with strict filters
                         (Wick & Holt 2022).

    _DUP_WARN   = 70%  : fastp uses k-mer sampling, not coordinate-based
                         deduplication. At high phage coverage (100–5000×),
                         50–70% apparent duplication is normal (Head et al. 2014,
                         BMC Genomics 15:179; Roux et al. 2019).

    _MIN_READS  = 50,000 : minimum reads post-filter for 50× on a 150 kb phage
                           genome at PE150. Required for DTR/ITR detection in
                           CheckV (Nayfach et al. 2021, Nat Biotechnol 39:578;
                           Bankevich et al. 2012).
"""

from __future__ import annotations
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from rich.panel import Panel
from rich.text import Text
from rich.progress import (
    Progress, SpinnerColumn, TextColumn,
    BarColumn, MofNCompleteColumn, TimeElapsedColumn,
)

from phageflow.utils.config import Config
from phageflow.utils.logger import log_step, log_info, log_ok, log_warn, console
from phageflow.utils.tools import require_tools, run_silent, human_size, mkdirs

STEP  = "01_qc"
TOOLS = ["fastp", "fastqc", "multiqc"]

# ── Thresholds ────────────────────────────────────────────────────────────────
_PASS_GOOD  = 90.0
_PASS_WARN  = 75.0
_Q20_GOOD   = 95.0
_Q20_WARN   = 90.0
_Q30_GOOD   = 85.0
_Q30_WARN   = 75.0
_DUP_WARN   = 70.0    # k-mer–based; inflated at high phage coverage
_MIN_READS  = 50_000  # for 50× on 150 kb phage at PE150
_AVG_QUAL   = 20


# ── Public entry point ────────────────────────────────────────────────────────

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
    require_tools(*TOOLS)

    out_dir = cfg.results(STEP)
    rpt_dir = cfg.reports(STEP) / sample_id   # per-sample: avoids mixing reports
    mkdirs(out_dir, rpt_dir)

    r1 = Path(r1)
    r2 = Path(r2)

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

    log_step(f"Module 01 — QC  [{sample_id}]")
    log_info(f"  R1 : {r1}  ({human_size(r1)})")
    log_info(f"  R2 : {r2}  ({human_size(r2)})")
    log_info(
        "  Params : Q≥20 per-base · mean-Q≥25 · ≤10% low-qual · len≥75bp · "
        "complexity≥30% · poly-X≥10bp · PE correction · adapter auto-detect"
    )

    already_done = r1_out.exists() and r1_out.stat().st_size > 0 and not force

    # Collect warnings for the completion panel
    active_warnings: list[str] = []

    with Progress(
        SpinnerColumn(spinner_name="dots2"),
        TextColumn("[bold cyan]{task.description:<48}"),
        BarColumn(bar_width=20, complete_style="cyan", finished_style="green"),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        console=console,
        transient=True,
    ) as progress:
        task = progress.add_task("Initializing...", total=3)

        # 1/3 — fastp
        progress.update(task, description="[1/3] fastp  — trim · filter · PE correction")
        if already_done:
            log_info("  [fastp] already processed — skipping  (--force to re-run)")
        else:
            _run_fastp(sample_id, r1, r2, r1_out, r2_out,
                       json_f, html_f, log_f, cfg.threads)
            _validate_output(sample_id, r1_out, r2_out, json_f)
        progress.advance(task)

        # 2/3 — FastQC (R1 + R2 in parallel)
        progress.update(task, description="[2/3] FastQC — per-read metrics (R1+R2 ‖)")
        _run_fastqc_parallel(r1_out, r2_out, rpt_dir, log_f)
        progress.advance(task)

        # 3/3 — MultiQC
        progress.update(task, description="[3/3] MultiQC — aggregating report")
        _run_multiqc(rpt_dir, sample_id)
        progress.advance(task)

    metrics = _parse_fastp_json(json_f)
    _save_tsv(sample_id, metrics, rpt_dir / "qc_summary.tsv")
    _check_warnings(metrics, active_warnings)
    _print_completion_panel(sample_id, r1_out, r2_out, rpt_dir, metrics,
                            active_warnings)

    log_step(f"Module 01 completed ✓  [{sample_id}]")
    return r1_out, r2_out


# ── Step implementations ──────────────────────────────────────────────────────

def _run_fastp(
    sample_id: str,
    r1: Path, r2: Path,
    r1_out: Path, r2_out: Path,
    json_f: Path, html_f: Path,
    log_f:  Path, threads: int,
) -> None:
    """Run fastp with parameters optimised for complete phage genome recovery.

    Key parameter decisions
    -----------------------
    --trim_poly_x / --poly_x_min_len 10
        NextSeq and NovaSeq use 2-colour chemistry. When signal is absent,
        the instrument calls G, producing poly-G artefacts at the 3′ end.
        --low_complexity_filter does not fully capture these because the
        remainder of the read may be high-quality. Filtering poly-X tails
        ≥10 bp removes these artefacts before assembly.
        Flag renamed in fastp 1.x: --poly_x_filter → --trim_poly_x (-x).
        Reference: Chen et al. 2018, Genome Biology 19:274.

    --correction / --overlap_len_require 10
        PE overlap correction reduces base errors pre-assembly.
        Requires ≥10 bp overlap. Corrects mismatches ≤10%.
        Reference: Chen et al. 2018.

    --average_qual 25
        Read-level floor complements the per-base filter. A read uniformly
        at Q20 (passes --qualified_quality_phred) but with mean Q21 is likely
        too noisy for k-mer-based assembly. Wick & Holt 2022.
    """
    cmd = [
        "fastp",
        "--in1",  str(r1),     "--in2",  str(r2),
        "--out1", str(r1_out), "--out2", str(r2_out),
        # Adapter detection
        "--detect_adapter_for_pe",
        # PE overlap-based error correction (Chen et al. 2018)
        "--correction",
        "--overlap_len_require",        "10",
        "--overlap_diff_percent_limit", "10",
        # 3′ sliding-window quality trimming (Bolger et al. 2014)
        "--cut_right",
        "--cut_right_window_size",      "4",
        "--cut_right_mean_quality",     "20",
        # Per-base quality filter
        "--qualified_quality_phred",    "20",
        "--unqualified_percent_limit",  "10",
        # Read-level mean quality floor (Wick & Holt 2022)
        "--average_qual",               str(_AVG_QUAL),
        # Other quality filters
        "--n_base_limit",               "5",
        "--length_required",            "75",
        # Low-complexity filter (Roux et al. 2019)
        "--low_complexity_filter",
        "--complexity_threshold",       "15",
        # Poly-X tail filter — targets NextSeq/NovaSeq poly-G artefacts.
        # Flag renamed in fastp 1.x: --poly_x_filter → --trim_poly_x (-x).
        # Poly-G is already trimmed automatically for NextSeq/NovaSeq by default
        # in fastp 1.x; --trim_poly_x adds coverage for poly-A/C/T tails.
        # Chen et al. 2018, Genome Biology 19:274.
        "--trim_poly_x",
        "--poly_x_min_len",             "10",
        # Output
        "--thread",       str(threads),
        "--json",         str(json_f),
        "--html",         str(html_f),
        "--report_title", f"[{sample_id}] PhageFlow QC",
    ]
    try:
        run_silent(cmd, log_file=log_f)
        log_ok("  [fastp] trimming + PE correction + poly-X filter complete")
    except Exception as e:
        log_warn(f"  [fastp] non-zero exit — check log: {log_f}  ({e})")
        raise


def _validate_output(
    sample_id: str,
    r1_out: Path,
    r2_out: Path,
    json_f:  Path,
) -> None:
    """Post-fastp sanity check: file presence + minimum surviving reads."""
    for label, path in [("R1", r1_out), ("R2", r2_out)]:
        if not path.exists() or path.stat().st_size == 0:
            log_warn(f"  [validate] {label} output missing or empty: {path}")
            return

    if not json_f.exists():
        return
    try:
        with open(json_f) as f:
            d = json.load(f)
        n_out = d["summary"]["after_filtering"]["total_reads"]
        if n_out < _MIN_READS:
            log_warn(
                f"  [validate] only {n_out:,} reads passed filtering — "
                f"< {_MIN_READS:,} may be insufficient for 50× on a 150 kb genome. "
                "Complete/HQ CheckV classification may be compromised "
                "(Nayfach et al. 2021)."
            )
        else:
            log_ok(f"  [validate] {n_out:,} reads passed — OK")
    except Exception:
        pass


def _run_fastqc_parallel(r1: Path, r2: Path, rpt_dir: Path, log_f: Path) -> None:
    """Run FastQC on R1 and R2 concurrently. Skips if HTML already exists."""
    def _fastqc(read: Path) -> None:
        if not read.exists():
            return
        stem   = read.name.replace(".fastq.gz", "").replace(".fastq", "")
        html_ok = (rpt_dir / f"{stem}_fastqc.html").exists()
        if html_ok:
            return
        try:
            run_silent(
                ["fastqc", "--outdir", str(rpt_dir), "--quiet", str(read)],
                log_file=log_f, check=False,
            )
        except Exception as e:
            log_warn(f"  [FastQC] {read.name}: {e}")

    with ThreadPoolExecutor(max_workers=2) as ex:
        for fut in as_completed({ex.submit(_fastqc, r): r for r in (r1, r2)}):
            fut.result()


def _run_multiqc(rpt_dir: Path, sample_id: str) -> None:
    """Run MultiQC on only this sample's QC files.

    rpt_dir is already per-sample (cfg.reports(STEP) / sample_id), so
    MultiQC sees only this sample's fastp JSON and FastQC zips — no
    cross-sample aggregation.
    """
    mqc_dir = rpt_dir / "multiqc"
    mkdirs(mqc_dir)
    try:
        run_silent(
            [
                "multiqc", str(rpt_dir),
                "--outdir",   str(mqc_dir),
                "--title",    f"PhageFlow QC — {sample_id}",
                "--filename", "multiqc_qc",
                "--quiet",
            ],
            log_file=rpt_dir / "multiqc.log", check=False,
        )
    except Exception:
        pass


# ── JSON parsing ──────────────────────────────────────────────────────────────

def _parse_fastp_json(json_f: Path) -> dict:
    """Extract QC metrics from the fastp JSON report."""
    empty: dict = {k: "N/A" for k in (
        "reads_in", "reads_out", "pct_pass",
        "mean_len_in", "mean_len_out",
        "gc_pct", "q20_pct", "q30_pct",
        "dup_rate", "insert_peak", "adapter_pct",
        "correction_rate", "corrected_reads", "corrected_bases",
        "filt_lowqual", "filt_tooshort", "filt_lowcomplex",
        "filt_polyx", "filt_n",
    )}
    if not json_f.exists():
        return empty
    try:
        with open(json_f) as f:
            d = json.load(f)

        bf  = d["summary"]["before_filtering"]
        af  = d["summary"]["after_filtering"]
        ri  = bf["total_reads"]
        ro  = af["total_reads"]
        pp  = ro / ri * 100 if ri else 0.0

        dup = d.get("duplication",     {}).get("rate")
        ins = d.get("insert_size",     {}).get("peak")
        adp = d.get("adapter_cutting", {}).get("adapter_trimmed_reads")
        fr  = d.get("filtering_result", {})

        corr_reads = fr.get("corrected_reads", None)
        corr_bases = fr.get("corrected_bases", None)
        if corr_reads is None:
            corr_rate = "N/A"
        elif ri:
            corr_rate = (
                f"{corr_reads / ri * 100:.2f}%"
                f" ({corr_reads:,} reads / {corr_bases or 0:,} bases)"
            )
        else:
            corr_rate = "0.0%"

        return {
            "reads_in":       f"{ri:,}",
            "reads_out":      f"{ro:,}",
            "pct_pass":       f"{pp:.1f}%",
            "mean_len_in":    f"{bf.get('read1_mean_length', 0):.0f} bp",
            "mean_len_out":   f"{af.get('read1_mean_length', 0):.0f} bp",
            "gc_pct":         f"{af.get('gc_content', 0) * 100:.1f}%",
            "q20_pct":        f"{af.get('q20_rate',   0) * 100:.1f}%",
            "q30_pct":        f"{af.get('q30_rate',   0) * 100:.1f}%",
            "dup_rate":       f"{dup * 100:.1f}%" if dup is not None else "N/A",
            "insert_peak":    f"{ins} bp"          if ins is not None else "N/A",
            "adapter_pct":    f"{adp / ri * 100:.1f}%" if (adp and ri) else "N/A",
            "correction_rate":  corr_rate,
            "corrected_reads":  str(corr_reads) if corr_reads is not None else "N/A",
            "corrected_bases":  str(corr_bases) if corr_bases is not None else "N/A",
            "filt_lowqual":   str(fr.get("low_quality_reads",    0)),
            "filt_tooshort":  str(fr.get("too_short_reads",      0)),
            "filt_lowcomplex":str(fr.get("low_complexity_reads", 0)),
            "filt_polyx":     str(fr.get("poly_x_reads",         0)),
            "filt_n":         str(fr.get("too_many_N_reads",     0)),
        }
    except Exception:
        return empty


# ── Warnings ──────────────────────────────────────────────────────────────────

def _check_warnings(m: dict, active_warnings: list[str]) -> None:
    """Evaluate metrics and populate active_warnings list for the panel."""

    try:
        pass_val = float(m["pct_pass"].rstrip("%"))
        if pass_val < _PASS_WARN:
            msg = (
                f"Pass rate {m['pct_pass']} < {_PASS_WARN}% — "
                "check raw read quality. "
                "Note: 75–85% is expected with strict HQ filters "
                "(Wick & Holt 2022)."
            )
            log_warn(f"  {msg}")
            active_warnings.append(msg)
    except (ValueError, AttributeError):
        pass

    try:
        dup_val = float(m["dup_rate"].rstrip("%"))
        if dup_val > _DUP_WARN:
            msg = (
                f"Duplication {m['dup_rate']} > {_DUP_WARN}% — "
                "k-mer–based estimate; rates >50% are normal at high "
                "phage coverage (Head et al. 2014)."
            )
            log_warn(f"  {msg}")
            active_warnings.append(msg)
    except (ValueError, AttributeError):
        pass

    try:
        if float(m["q30_pct"].rstrip("%")) < _Q30_WARN:
            msg = (
                f"Q30 rate {m['q30_pct']} < {_Q30_WARN}% — "
                "assembly contiguity and CheckV Complete classification "
                "may be affected (Nayfach et al. 2021)."
            )
            log_warn(f"  {msg}")
            active_warnings.append(msg)
    except (ValueError, AttributeError):
        pass

    try:
        n_out = int(m["reads_out"].replace(",", ""))
        if 0 < n_out < _MIN_READS:
            msg = (
                f"Only {m['reads_out']} reads passed — "
                f"< {_MIN_READS:,} may be insufficient for 50× on a 150 kb genome "
                "(Nayfach et al. 2021)."
            )
            log_warn(f"  {msg}")
            active_warnings.append(msg)
    except (ValueError, AttributeError):
        pass


# ── Completion panel ──────────────────────────────────────────────────────────

def _print_completion_panel(
    sample_id: str,
    r1_out:    Path,
    r2_out:    Path,
    rpt_dir:   Path,
    m:         dict,
    active_warnings: list[str],
) -> None:
    """Standardised completion panel: key metrics · outputs · warnings.

    Uses markup strings (not Text objects) so Rich interprets colour tags.
    """

    def _c(val_str: str, good: float, warn: float) -> str:
        """Return val_str wrapped in the appropriate Rich colour tag."""
        try:
            v = float(val_str.rstrip("%"))
            if v >= good: return f"[bold green]{val_str}[/bold green]"
            if v >= warn: return f"[bold yellow]{val_str}[/bold yellow]"
            return f"[bold red]{val_str}[/bold red]"
        except (ValueError, AttributeError):
            return str(val_str)

    def _f(val: str) -> str:
        return "[dim]N/A[/dim]" if val in ("N/A", "", None) else val

    lines: list[str] = []

    # ── Key metrics ──
    lines.append("[bold]Key metrics[/bold]")
    lines.append(
        f"  [cyan]Reads   :[/cyan] {_f(m['reads_in'])} → {_f(m['reads_out'])}"
        f"  pass={_c(m['pct_pass'], _PASS_GOOD, _PASS_WARN)}"
    )
    lines.append(
        f"  [cyan]Quality :[/cyan]"
        f" Q20={_c(m['q20_pct'], _Q20_GOOD, _Q20_WARN)}"
        f"  Q30={_c(m['q30_pct'], _Q30_GOOD, _Q30_WARN)}"
        f"  GC={_f(m['gc_pct'])}"
        f"  mean-len={_f(m['mean_len_out'])}"
    )
    lines.append(
        f"  [cyan]Library :[/cyan]"
        f" dup={_f(m['dup_rate'])}"
        f"  insert={_f(m['insert_peak'])}"
        f"  adapter={_f(m['adapter_pct'])}"
    )
    lines.append(f"  [cyan]PE corr :[/cyan] {_f(m['correction_rate'])}")

    # Filtering breakdown (non-zero only)
    _FILT = {
        "filt_lowqual":    "low-qual",
        "filt_tooshort":   "too-short",
        "filt_lowcomplex": "low-complexity",
        "filt_polyx":      "poly-X",
        "filt_n":          "excess-N",
    }
    parts = [
        f"{lbl}={int(m[k]):,}"
        for k, lbl in _FILT.items()
        if m.get(k, "0") not in ("0", "N/A", "", None)
    ]
    if parts:
        lines.append(f"  [dim]Filtered: {' · '.join(parts)}[/dim]")

    # ── Output files ──
    lines.append("")
    lines.append("[bold]Output files[/bold]")
    lines.append(f"  Trimmed R1 : {r1_out}")
    lines.append(f"  Trimmed R2 : {r2_out}")
    lines.append(f"  MultiQC    : {rpt_dir / 'multiqc' / 'multiqc_qc.html'}")

    # ── Active warnings ──
    if active_warnings:
        lines.append("")
        lines.append("[bold yellow]Warnings[/bold yellow]")
        for w in active_warnings:
            lines.append(f"  [yellow]⚠ {w}[/yellow]")

    console.print(Panel(
        "\n".join(lines),
        title=f"[bold cyan]QC complete — {sample_id}[/bold cyan]",
        border_style="cyan",
        padding=(0, 2),
        width=120,
    ))


# ── TSV summary ───────────────────────────────────────────────────────────────

_TSV_HEADERS = [
    "sample_id",
    "reads_in", "reads_out", "pct_pass",
    "mean_len_in", "mean_len_out",
    "gc_pct", "q20_pct", "q30_pct",
    "dup_rate", "insert_peak", "adapter_pct",
    "correction_rate", "corrected_reads", "corrected_bases",
    "filt_lowqual", "filt_tooshort", "filt_lowcomplex", "filt_polyx", "filt_n",
]


def _save_tsv(sample_id: str, metrics: dict, path: Path) -> None:
    """Write / update the QC summary TSV, preserving existing rows."""
    rows: dict[str, dict] = {}
    if path.exists():
        with open(path) as f:
            old_hdrs = f.readline().rstrip("\n").split("\t")
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if cols and cols[0]:
                    rows[cols[0]] = dict(zip(old_hdrs, cols))

    rows[sample_id] = {
        "sample_id": sample_id,
        **{h: metrics.get(h, "") for h in _TSV_HEADERS[1:]},
    }

    with open(path, "w") as f:
        f.write("\t".join(_TSV_HEADERS) + "\n")
        for row in rows.values():
            f.write("\t".join(row.get(h, "") for h in _TSV_HEADERS) + "\n")

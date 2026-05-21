"""PhageFlow Module 01 — Quality control and trimming.

Optimized for purified phage Illumina PE150 sequencing.
Goal: complete and High-quality genomes via CheckV (Nayfach et al. 2021).

Tools:
    fastp   : adapter trimming, PE overlap correction, sliding-window quality
              trimming, low-complexity filtering
              (Chen et al. 2018, Genome Biology 19:274)
    FastQC  : per-read quality metrics (Andrews 2010)
    MultiQC : aggregate QC report (Ewels et al. 2016, Bioinformatics)

fastp parameters (literature-based for phage complete genome recovery):

    --correction / --overlap_len_require 10 / --overlap_diff_percent_limit 10
        PE overlap-based base error correction (Chen et al. 2018).
        Reduces base errors before De Bruijn graph construction → longer SPAdes
        contigs → higher CheckV completeness scores (Bankevich et al. 2012).
        Corrects mismatches in overlapping PE read pairs; requires ≥10 bp overlap
        and tolerates ≤10% mismatch in the overlap region.

    --cut_right / --cut_right_window_size 4 / --cut_right_mean_quality 20
        Sliding-window 3' trimming (Bolger et al. 2014; Chen et al. 2018).

    --qualified_quality_phred 20
        Q20 threshold: 99% base call accuracy (Illumina quality guidelines).

    --unqualified_percent_limit 10
        Max 10% low-quality bases per read — tighter than the fastp default of
        40% to retain only high-quality reads for complete genome assembly
        (Wick & Holt 2022, Microb Genomics 8:mgen000788).

    --average_qual 25
        Discard reads with mean Phred quality < Q25. Per-base thresholds can
        pass reads that are uniformly mediocre (many Q20-Q24 bases, none below
        Q20); a mean-quality floor complements the per-base filter and is a
        better predictor of assembly quality (Wick & Holt 2022).

    --length_required 75
        Minimum 75 bp for PE150 data (Bankevich et al. 2012).

    --low_complexity_filter --complexity_threshold 30
        Removes homopolymer/repetitive reads (Roux et al. 2019, eLife).

    --n_base_limit 5
        Remove reads with >5 N calls (standard practice).

    --detect_adapter_for_pe
        Automatic adapter detection for paired-end data (Chen et al. 2018).

Thresholds (mode: purified_phage):

    _DUP_WARN  = 70% : fastp estimates duplication via k-mer sampling, NOT
        coordinate-based (unlike Picard). At the coverage typical of purified
        phage preparations (100–5000×), many independent fragments share
        identical k-mers from the phage genome, inflating the apparent
        duplication rate. Rates of 50–70% are expected at high coverage without
        PCR artefacts. The threshold is set to 70% to avoid false positives
        (Head et al. 2014, BMC Genomics 15:179; Roux et al. 2019).

    _MIN_READS = 50 000 : minimum reads after filtering to trigger the low-read
        warning. Rationale — to assemble a complete 150 kb phage genome
        (e.g. Herelleviridae, Ackermannviridae) at ≥50× coverage with PE150:
            150 000 bp × 50× / 150 bp × 2 reads ≈ 50 000 read pairs.
        50× is the empirical minimum for SPAdes to resolve DTR/ITR terminal
        repeats that CheckV uses to classify genomes as "Complete"
        (Nayfach et al. 2021, Nat Biotechnol 39:578; Bankevich et al. 2012).
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
_PASS_GOOD  = 90.0   # pass rate: green
_PASS_WARN  = 75.0   # pass rate: yellow; below → red + warning message
                     # Lowered from 80% to 75%: --unqualified_percent_limit 10
                     # and --average_qual 25 intentionally reject more reads
                     # than the fastp defaults. A 75–85% pass rate is expected
                     # and desirable with strict HQ parameters
                     # (Wick & Holt 2022, Microb Genomics 8:mgen000788).
_Q20_GOOD   = 95.0
_Q20_WARN   = 90.0
_Q30_GOOD   = 85.0
_Q30_WARN   = 75.0
_DUP_WARN   = 70.0   # fastp k-mer-based; inflated at high phage coverage
                     # (Head et al. 2014, BMC Genomics; Roux et al. 2019)
_MIN_READS  = 50_000 # minimum for 50× on 150 kb phage genome at PE150
                     # required for DTR/ITR detection (Nayfach et al. 2021)
_AVG_QUAL   = 25     # mean Phred per read; complements per-base filter
                     # (Wick & Holt 2022, Microb Genomics)


# ── Public entry point ────────────────────────────────────────────────────────

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

    # ── Early input validation ─────────────────────────────────────────────────
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
    _print_input_table(sample_id, r1, r2)
    log_info(
        "  Parameters: Q≥20 (per-base) | mean-Q≥25 | ≤10% low-qual bases | "
        "len≥75 bp | complexity≥30% | PE correction | sliding-window 4/Q20 | "
        "adapter auto-detect"
    )
    log_info(
        "  PE correction: --correction --overlap_len_require 10 "
        "--overlap_diff_percent_limit 10  (Chen et al. 2018)"
    )

    already_done = r1_out.exists() and r1_out.stat().st_size > 0 and not force

    # ── Pipeline steps ────────────────────────────────────────────────────────
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
        progress.update(task, description="[1/3] fastp  — trimming, correction & filtering")
        if already_done:
            log_info("  [fastp] already trimmed — skipping  (--force to re-run)")
        else:
            _run_fastp(sample_id, r1, r2, r1_out, r2_out, json_f, html_f, log_f, cfg.threads)
            _validate_output(sample_id, r1_out, r2_out, json_f)
        progress.advance(task)

        # 2/3 — FastQC (R1 + R2 in parallel, with output caching)
        progress.update(task, description="[2/3] FastQC — per-read metrics (R1+R2 ‖)")
        _run_fastqc_parallel(r1_out, r2_out, rpt_dir, log_f)
        progress.advance(task)

        # 3/3 — MultiQC
        progress.update(task, description="[3/3] MultiQC — aggregating report")
        _run_multiqc(rpt_dir)
        progress.advance(task)

    # ── Metrics, display, save ─────────────────────────────────────────────────
    metrics = _parse_fastp_json(json_f)
    _save_tsv(sample_id, metrics, rpt_dir / "qc_summary.tsv")
    _print_summary_table(sample_id, metrics)
    _check_warnings(metrics)
    _print_completion_panel(sample_id, r1_out, r2_out, rpt_dir, metrics)

    log_step(f"Module 01 completed ✓  [{sample_id}]")
    log_info(
        f"  Next: phageflow host-removal --sample-id {sample_id} "
        f"[--host-file | --accessions | --kraken-db]"
    )
    return r1_out, r2_out


# ── Step implementations ──────────────────────────────────────────────────────

def _run_fastp(
    sample_id: str,
    r1: Path, r2: Path,
    r1_out: Path, r2_out: Path,
    json_f: Path, html_f: Path,
    log_f:  Path, threads: int,
) -> None:
    """
    Run fastp with parameters optimised for complete phage genome recovery.

    Key additions vs. generic QC:
        --correction               : PE overlap-based error correction
        --overlap_len_require 10   : min overlap bp to attempt correction
        --overlap_diff_percent_limit 10 : max mismatch % in overlap
        --average_qual 25          : read-level mean quality floor
        --unqualified_percent_limit 10  : tighter than default 40%

    References
    ----------
    Chen et al. (2018) Genome Biology 19:274 — fastp design and correction.
    Wick & Holt (2022) Microb Genomics 8:mgen000788 — per-read quality filters.
    Bankevich et al. (2012) J Comput Biol — k-mer quality and SPAdes assembly.
    """
    cmd = [
        "fastp",
        "--in1",  str(r1),     "--in2",  str(r2),
        "--out1", str(r1_out), "--out2", str(r2_out),
        # Adapter detection
        "--detect_adapter_for_pe",
        # PE overlap-based error correction (Chen et al. 2018)
        "--correction",
        "--overlap_len_require",         "10",
        "--overlap_diff_percent_limit",  "10",
        # 3' sliding-window quality trimming (Bolger et al. 2014)
        "--cut_right",
        "--cut_right_window_size",       "4",
        "--cut_right_mean_quality",      "20",
        # Per-base quality filter
        "--qualified_quality_phred",     "20",
        "--unqualified_percent_limit",   "10",   # tighter for HQ assembly
        # Read-level mean quality floor (Wick & Holt 2022)
        "--average_qual",                str(_AVG_QUAL),
        # Other filters
        "--n_base_limit",                "5",
        "--length_required",             "75",
        "--low_complexity_filter",
        "--complexity_threshold",        "30",
        # Output
        "--thread",       str(threads),
        "--json",         str(json_f),
        "--html",         str(html_f),
        "--report_title", f"[{sample_id}] PhageFlow QC",
    ]
    try:
        run_silent(cmd, log_file=log_f)
        log_ok("  fastp · trimming + PE correction complete")
    except Exception as e:
        log_warn(f"  fastp · non-zero exit — check log: {log_f}  ({e})")
        raise


def _validate_output(
    sample_id: str,
    r1_out: Path, r2_out: Path,
    json_f:  Path,
) -> None:
    """
    Post-fastp sanity checks: file presence + minimum surviving reads.

    _MIN_READS = 50 000 is the empirical minimum for 50× coverage of a
    150 kb phage genome (PE150), required for DTR/ITR detection in CheckV
    (Nayfach et al. 2021). Genomes below this coverage threshold are unlikely
    to receive a 'Complete' or 'High-quality' classification.
    """
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
                f"  validate · only {n_out:,} reads passed filtering — "
                f"< {_MIN_READS:,} reads may be insufficient for 50× coverage "
                f"of a 150 kb phage genome; Complete/HQ CheckV classification "
                f"may be compromised (Nayfach et al. 2021)."
            )
        else:
            log_ok(f"  validate · {n_out:,} reads passed — OK")
    except Exception:
        pass


def _run_fastqc_parallel(r1: Path, r2: Path, rpt_dir: Path, log_f: Path) -> None:
    """Run FastQC on R1 and R2 concurrently. Skips if HTML output already exists."""
    def _fastqc(read: Path) -> None:
        if not read.exists():
            return
        # Cache check: FastQC names output as <stem>_fastqc.html
        stem    = read.name.replace(".fastq.gz", "").replace(".fastq", "")
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
        pass


# ── JSON parsing ──────────────────────────────────────────────────────────────

def _parse_fastp_json(json_f: Path) -> dict:
    """
    Extract QC metrics from fastp JSON output.

    Fields returned
    ---------------
    reads_in / reads_out / pct_pass
    mean_len_in / mean_len_out
    gc_pct / q20_pct / q30_pct
    dup_rate          — duplication rate (k-mer-based; inflated at high coverage)
    insert_peak       — insert size peak (bp)
    adapter_pct       — % reads with adapter trimmed
    correction_rate   — % read pairs corrected by PE overlap correction (new)
    filt_lowqual      — reads removed: low quality
    filt_tooshort     — reads removed: too short
    filt_lowcomplex   — reads removed: low complexity
    filt_n            — reads removed: excess N bases
    """
    empty: dict = {k: "N/A" for k in (
        "reads_in", "reads_out", "pct_pass",
        "mean_len_in", "mean_len_out",
        "gc_pct", "q20_pct", "q30_pct",
        "dup_rate", "insert_peak", "adapter_pct",
        "correction_rate",
        "filt_lowqual", "filt_tooshort", "filt_lowcomplex", "filt_n",
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

        ml_in  = bf.get("read1_mean_length", 0)
        ml_out = af.get("read1_mean_length", 0)

        dup    = d.get("duplication",    {}).get("rate")
        ins    = d.get("insert_size",    {}).get("peak")
        adp_r  = d.get("adapter_cutting", {}).get("adapter_trimmed_reads")
        fr     = d.get("filtering_result", {})

        # PE correction stats (present only when --correction is active).
        # corrected_reads == 0 is a valid result (no overlapping pairs found);
        # must show "0.0%" not "N/A" so the user knows correction ran.
        corr_data  = d.get("correction", {})
        corr_pairs = corr_data.get("corrected_reads", None)
        if corr_pairs is None:
            # --correction was not active or fastp version does not emit this key
            corr_rate = "N/A"
        elif ri:
            corr_rate = f"{corr_pairs / ri * 100:.1f}%"
        else:
            corr_rate = "0.0%"

        return {
            "reads_in":        f"{ri:,}",
            "reads_out":       f"{ro:,}",
            "pct_pass":        f"{pp:.1f}%",
            "mean_len_in":     f"{ml_in:.0f} bp",
            "mean_len_out":    f"{ml_out:.0f} bp",
            "gc_pct":          f"{af.get('gc_content', 0) * 100:.1f}%",
            "q20_pct":         f"{af.get('q20_rate',   0) * 100:.1f}%",
            "q30_pct":         f"{af.get('q30_rate',   0) * 100:.1f}%",
            "dup_rate":        f"{dup * 100:.1f}%" if dup  is not None else "N/A",
            "insert_peak":     f"{ins} bp"          if ins  is not None else "N/A",
            "adapter_pct":     f"{adp_r / ri * 100:.1f}%" if (adp_r and ri) else "N/A",
            "correction_rate": corr_rate,
            "filt_lowqual":    str(fr.get("low_quality_reads",    0)),
            "filt_tooshort":   str(fr.get("too_short_reads",      0)),
            "filt_lowcomplex": str(fr.get("low_complexity_reads", 0)),
            "filt_n":          str(fr.get("too_many_N_reads",     0)),
        }
    except Exception:
        return empty


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


def _print_summary_table(sample_id: str, m: dict) -> None:
    """Print QC metrics as compact log lines."""
    log_ok(
        f"  Reads    : {m['reads_in']} → {m['reads_out']}  "
        f"({_color_rate(m['pct_pass'], _PASS_GOOD, _PASS_WARN)} pass)"
    )
    log_ok(
        f"  Quality  : Q20={m['q20_pct']}  Q30={m['q30_pct']}  "
        f"GC={m['gc_pct']}  len={m['mean_len_out']}"
    )
    log_ok(
        f"  Library  : dup={_fmt(m['dup_rate'])}  "
        f"insert={_fmt(m['insert_peak'])}  adapter={_fmt(m['adapter_pct'])}"
    )
    log_ok(
        f"  Correction: PE overlap corrected={_fmt(m['correction_rate'])}  "
        f"(Chen et al. 2018)"
    )

    # Filtering breakdown — only non-zero entries
    _FILTER_LABELS = {
        "filt_lowqual":    "low-quality",
        "filt_tooshort":   "too-short",
        "filt_lowcomplex": "low-complexity",
        "filt_n":          "excess-N",
    }
    parts = [
        f"{label}={int(m[key]):,}"
        for key, label in _FILTER_LABELS.items()
        if m.get(key, "0") not in ("0", "N/A", "", None)
    ]
    if parts:
        log_info(f"  Filtered : {' | '.join(parts)}")


def _check_warnings(m: dict) -> None:
    """Emit contextual warnings based on QC metrics."""
    try:
        pass_val = float(m["pct_pass"].rstrip("%"))
        if pass_val < _PASS_WARN:
            log_warn(
                f"  Pass rate {m['pct_pass']} < {_PASS_WARN}% — "
                "check raw read quality or residual host contamination. "
                "Note: --unqualified_percent_limit 10 and --average_qual 25 "
                "are intentionally strict; 75–85% pass rates are expected "
                "and indicate effective removal of low-quality reads "
                "(Wick & Holt 2022)."
            )
        elif pass_val < _PASS_GOOD:
            log_info(
                f"  Pass rate {m['pct_pass']}: acceptable with strict HQ filters "
                f"(--unqualified_percent_limit 10 / --average_qual 25). "
                f"Q30={m.get('q30_pct','?')} confirms retained reads are high quality."
            )
    except (ValueError, AttributeError):
        pass

    try:
        dup_val = float(m["dup_rate"].rstrip("%"))
        if dup_val > _DUP_WARN:
            log_warn(
                f"  Duplication {m['dup_rate']} > {_DUP_WARN}% — "
                "fastp uses k-mer sampling; rates >50% are expected at high "
                "phage coverage (100–5000×) without PCR artefacts. "
                "Consider Picard MarkDuplicates for coordinate-based estimation "
                "if PCR bias is a concern (Head et al. 2014, BMC Genomics)."
            )
    except (ValueError, AttributeError):
        pass

    try:
        if float(m["q30_pct"].rstrip("%")) < _Q30_WARN:
            log_warn(
                f"  Q30 rate {m['q30_pct']} < {_Q30_WARN}% — "
                "base call quality is low; SPAdes assembly contiguity will be "
                "reduced and Complete/HQ genome classification in CheckV "
                "may be compromised (Nayfach et al. 2021)."
            )
    except (ValueError, AttributeError):
        pass

    try:
        n_out_str = m.get("reads_out", "0").replace(",", "")
        n_out = int(n_out_str) if n_out_str not in ("N/A", "") else 0
        if 0 < n_out < _MIN_READS:
            log_warn(
                f"  Only {m['reads_out']} reads passed filtering — "
                f"< {_MIN_READS:,} reads may be insufficient for 50× coverage "
                f"of a 150 kb phage genome; Complete/HQ CheckV classification "
                f"requires adequate DTR/ITR coverage (Nayfach et al. 2021; "
                f"Bankevich et al. 2012)."
            )
    except (ValueError, AttributeError):
        pass


def _print_completion_panel(sample_id, r1_out, r2_out, rpt_dir, m: dict) -> None:
    """Clean panel summarising output paths and key stats."""
    text = Text()
    text.append("✓ ", style="bold green")
    text.append(f"{m.get('reads_in', '?')} reads  →  ", style="dim white")
    text.append(f"{m.get('reads_out', '?')}", style="bold green")
    text.append(f"  (pass rate: {m.get('pct_pass', '?')})\n", style="cyan")
    text.append(
        f"  PE correction: {m.get('correction_rate', 'N/A')} reads corrected\n\n",
        style="dim white",
    )
    text.append("Trimmed R1 : ", style="dim white")
    text.append(str(r1_out) + "\n", style="white")
    text.append("Trimmed R2 : ", style="dim white")
    text.append(str(r2_out) + "\n", style="white")
    text.append("MultiQC    : ", style="dim white")
    text.append(str(rpt_dir / "multiqc" / "multiqc_qc.html"), style="white")

    console.print(Panel(
        text,
        title=f"[bold cyan]QC complete — {sample_id}[/bold cyan]",
        border_style="cyan",
        padding=(0, 2),
        width=72,
    ))


# ── TSV summary ───────────────────────────────────────────────────────────────

_TSV_HEADERS = [
    "sample_id",
    "reads_in", "reads_out", "pct_pass",
    "mean_len_in", "mean_len_out",
    "gc_pct", "q20_pct", "q30_pct",
    "dup_rate", "insert_peak", "adapter_pct",
    "correction_rate",
    "filt_lowqual", "filt_tooshort", "filt_lowcomplex", "filt_n",
]


def _save_tsv(sample_id: str, metrics: dict, path: Path) -> None:
    """
    Write / update the QC summary TSV.

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

    rows[sample_id] = {
        "sample_id": sample_id,
        **{h: metrics.get(h, "") for h in _TSV_HEADERS[1:]},
    }

    with open(path, "w") as f:
        f.write("\t".join(_TSV_HEADERS) + "\n")
        for row in rows.values():
            f.write("\t".join(row.get(h, "") for h in _TSV_HEADERS) + "\n")

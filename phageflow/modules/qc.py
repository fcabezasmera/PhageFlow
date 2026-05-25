"""PhageFlow Module 01 — Read quality control (fastp + FastQC + MultiQC).

Goal: trim adapters, filter low-quality reads, and remove low-complexity
sequences from paired-end Illumina reads prior to host removal and assembly.
Inputs are assumed to be raw reads of unknown composition (metagenome,
virome, or purified preparation). No phage-specific assumptions are made here.

Tool: fastp (Chen et al. 2018, Bioinformatics 34:i884)
  Adapter trimming  : auto-detect PE adapters (--detect_adapter_for_pe)
  Quality trimming  : sliding-window 3\'→5\' (--cut_right)
  PE correction     : overlap-based base correction for paired reads
  Low-complexity    : homopolymer / dust filter (Roux et al. 2019 MIUViG)
  Poly-X removal    : poly-G/C/A/T tail trimming (NovaSeq 2-color artefacts)

Parameter rationale
  average_qual 20       Q20 = 99% per-base accuracy; MIUViG minimum standard
                        (Roux et al. 2019, Nat Biotechnol 37:505).
                        Stricter thresholds (Q22+) discard reads from
                        divergent viral regions without improving assembly.
  length_required 50    Reads ≥50 bp assemble reliably in SPAdes/MEGAHIT
                        (Holtgrewe et al. 2013, PLoS ONE 8:e61458).
                        Retains reads from terminal repeats (DTR/ITR) and
                        small viral genomes (Microviridae ~3-6 kb,
                        Inoviridae ~6-9 kb).
  low_complexity 20%    Removes homopolymers and dust reads that inflate
                        assembly graphs without contributing sequence
                        information (Roux et al. 2019 MIUViG §2.1).
  correction true       PE overlap correction reduces SNP-like errors at
                        read ends, improving k-mer graph accuracy.
  trim_poly_x true      NovaSeq and NextSeq (2-color SBS) append poly-G
                        artefacts to reads shorter than the flow cell cycle
                        count (Illumina tech note 2017).

Duplication note
  High duplication (>50%) is expected and normal in purified phage
  preparations at high coverage. Optical/PCR duplicates are NOT removed
  here — deduplication would discard real reads in high-coverage samples
  and is not recommended for viral genome assembly
  (Head et al. 2014, Biotechniques 56:61).

Coverage estimates use genome size reference ranges:
  Small (Microviridae/Inoviridae) : ~5 kb
  Typical tailed phage            : ~50 kb
  Large (Jumbo phage)             : ~200-540 kb
  (Al-Shayeb et al. 2020, Nature 578:425)
"""
from __future__ import annotations
import gzip
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from rich.panel import Panel
from rich.progress import (
    Progress, SpinnerColumn, TextColumn, BarColumn,
    MofNCompleteColumn, TimeElapsedColumn,
)

from phageflow.utils.config import Config
from phageflow.utils.logger import log_step, log_info, log_ok, log_warn, console
from phageflow.utils.tools import require_tools, run_silent, human_size, mkdirs

STEP  = "01_qc"
TOOLS = ["fastp", "fastqc", "multiqc"]

_PASS_GOOD = 85.0
_PASS_WARN = 75.0
_Q20_GOOD  = 95.0
_Q20_WARN  = 90.0
_Q30_GOOD  = 80.0
_Q30_WARN  = 70.0
_MIN_READS          = 100_000
# Phage genome size references for coverage estimation display
# Al-Shayeb et al. 2020 (Nature 578:425); Ackermann 2007 (Arch Virol 152:227)
_REF_GENOME_SMALL   =   5_000   # Microviridae / Inoviridae
_REF_GENOME_TYPICAL =  50_000   # typical tailed phage
_REF_GENOME_LARGE   = 200_000   # jumbo phage (conservative)

def _detect_read_length(r1: Path, n_reads: int = 1000) -> int:
    """Auto-detect read length from first n_reads of R1 FASTQ.

    Samples the first n_reads sequences and returns the median read length.
    Used to adjust QC and assembly parameters for PE150/PE250/PE300.

    Samples 1000 reads by default — sufficient for stable median estimate
    while keeping I/O minimal (Schirmer et al. 2015, Nucleic Acids Res 43:e83).

    Returns
    -------
    int : detected read length bucket (150, 250, or 300)
          bucketed to nearest standard Illumina PE length.
    """
    lengths: list[int] = []
    try:
        opener = gzip.open if str(r1).endswith(".gz") else open
        with opener(r1, "rt", errors="replace") as f:
            i = 0
            for line in f:
                if i % 4 == 1:  # FASTQ sequence line
                    lengths.append(len(line.rstrip()))
                    if len(lengths) >= n_reads:
                        break
                i += 1
    except Exception:
        return 150  # safe default

    if not lengths:
        return 150

    median = sorted(lengths)[len(lengths) // 2]

    # Bucket to nearest standard Illumina PE length
    if median >= 275:
        return 300
    elif median >= 200:
        return 250
    else:
        return 150


def _read_length_params(read_length: int) -> dict:
    """Return QC parameter overrides for a given read length.

    PE150 (default): standard Illumina short-read parameters.
    PE250: MiSeq v2/v3 — longer reads, more overlap, larger windows.
    PE300: MiSeq v3 600-cycle — maximum read length, extensive overlap.

    Parameter rationale
    -------------------
    length_required:
        Minimum post-trim read length. Set to ~1/3 of read length to
        retain reads that are still long enough for reliable k-mer
        assembly. Too short → spurious k-mer matches; too long → discards
        valid trimmed reads from low-quality 3' ends.
        Holtgrewe et al. 2013 (PLoS ONE 8:e61458).

    overlap_len_require:
        PE250/PE300 produce genuine read overlaps of 50-150 bp in
        amplicon-like preparations. Requiring longer overlaps for
        correction reduces false corrections from chance overlaps.
        fastp documentation; Chen et al. 2018 (Bioinformatics 34:i884).

    cut_right_window_size:
        Sliding window for 3' quality trimming. Larger windows smooth
        out local quality dips in longer reads, avoiding over-trimming.
        Bolger et al. 2014 (Bioinformatics 30:2114) — Trimmomatic.
    """
    if read_length >= 275:  # PE300
        return {
            "length_required":   100,
            "overlap_len_require": 50,
            "cut_right_window_size": 8,
        }
    elif read_length >= 200:  # PE250
        return {
            "length_required":   75,
            "overlap_len_require": 30,
            "cut_right_window_size": 6,
        }
    else:  # PE150 — use config defaults unchanged
        return {}


def run(
    cfg:       Config,
    sample_id: str,
    r1:        Path,
    r2:        Path,
    force:     bool = False,
) -> tuple[Path, Path]:
    """Quality-control paired reads. Returns (r1_trimmed, r2_trimmed)."""
    r1 = Path(r1); r2 = Path(r2)
    out_dir = cfg.results(STEP)
    rpt_dir = cfg.reports(STEP) / sample_id
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
    log_info(f"  {sample_id}  ·  R1: {human_size(r1)}  ·  R2: {human_size(r2)}")

    if r1_out.exists() and r1_out.stat().st_size > 0 and not force:
        log_info("  Already processed — skipping  (--force to re-run)")
        return r1_out, r2_out

    require_tools(*TOOLS)
    active_warnings: list[str] = []

    # Auto-detect read length and adjust parameters accordingly
    detected_rl = _detect_read_length(r1)
    rl_overrides = _read_length_params(detected_rl)
    if rl_overrides:
        log_info(
            f"  Read length : {detected_rl} bp detected (PE{detected_rl}) — "
            f"adjusting: length_required={rl_overrides['length_required']} "
            f"overlap={rl_overrides['overlap_len_require']} "
            f"window={rl_overrides['cut_right_window_size']}"
        )
    else:
        log_info(f"  Read length : {detected_rl} bp detected (PE150) — using default parameters")

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

        progress.update(task, description="[1/3] fastp   — trimming + quality filter")
        _run_fastp(sample_id, r1, r2, r1_out, r2_out, json_f, html_f, log_f, cfg, rl_overrides)
        progress.advance(task)
        log_ok("  [1/3] fastp   — complete")

        progress.update(task, description="[2/3] FastQC  — per-read QC")
        _run_fastqc(r1_out, r2_out, rpt_dir, log_f)
        progress.advance(task)
        log_ok("  [2/3] FastQC  — reports written")

        progress.update(task, description="[3/3] MultiQC — aggregate report")
        _run_multiqc(rpt_dir, sample_id)
        progress.advance(task)
        log_ok("  [3/3] MultiQC — report ready")

    metrics = _parse_fastp_json(json_f)
    _check_warnings(metrics, active_warnings)
    _save_tsv(sample_id, metrics, rpt_dir / "qc_summary.tsv")
    _print_completion_panel(sample_id, r1_out, r2_out, rpt_dir, metrics, active_warnings)
    log_step(f"Module 01 completed ✓  [{sample_id}]")
    return r1_out, r2_out

def _run_fastp(
    sample_id: str,
    r1: Path, r2: Path,
    r1_out: Path, r2_out: Path,
    json_f: Path, html_f: Path,
    log_f:  Path, cfg,
    rl_overrides: dict | None = None,
) -> None:
    q = cfg.qc
    # Apply read-length-specific overrides (PE250/PE300 adjustments)
    length_required      = (rl_overrides or {}).get("length_required",      q.length_required)
    overlap_len_require  = (rl_overrides or {}).get("overlap_len_require",   q.overlap_len_require)
    cut_right_window_size = (rl_overrides or {}).get("cut_right_window_size", q.cut_right_window_size)
    cmd = [
        "fastp",
        "--in1",  str(r1),     "--in2",  str(r2),
        "--out1", str(r1_out), "--out2", str(r2_out),
        "--detect_adapter_for_pe",
        "--cut_right",
        "--cut_right_window_size",      str(cut_right_window_size),
        "--cut_right_mean_quality",     str(q.cut_right_mean_quality),
        "--qualified_quality_phred",    str(q.qualified_quality_phred),
        "--unqualified_percent_limit",  str(q.unqualified_percent_limit),
        "--average_qual",               str(q.average_qual),
        "--n_base_limit",               str(q.n_base_limit),
        "--length_required",            str(length_required),
        "--thread",                     str(cfg.threads),
        "--json",                       str(json_f),
        "--html",                       str(html_f),
        "--report_title",               f"PhageFlow QC — {sample_id}",
    ]
    if q.correction:
        cmd += [
            "--correction",
            "--overlap_len_require",        str(overlap_len_require),
            "--overlap_diff_percent_limit", str(q.overlap_diff_percent_limit),
        ]
    if q.low_complexity_filter:
        cmd += [
            "--low_complexity_filter",
            "--complexity_threshold", str(q.complexity_threshold),
        ]
    if q.trim_poly_x:
        cmd += [
            "--trim_poly_x",
            "--poly_x_min_len", str(q.poly_x_min_len),
        ]
    run_silent(cmd, log_file=log_f)

def _run_fastqc(r1: Path, r2: Path, out_dir: Path, log_f: Path) -> None:
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

def _parse_fastp_json(json_f: Path) -> dict:
    empty = {
        "reads_in": "0", "reads_out": "0", "pct_pass": "0.0%",
        "gc_pct":   "N/A", "q20_pct": "N/A", "q30_pct": "N/A",
        "bp_out":   "0",   "dup_pct": "N/A",
        "est_coverage_small": "N/A", "est_coverage_typical": "N/A", "est_coverage_large": "N/A",
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
        cov_small   = bp / _REF_GENOME_SMALL   if bp else 0.0
        cov_typical = bp / _REF_GENOME_TYPICAL if bp else 0.0
        cov_large   = bp / _REF_GENOME_LARGE   if bp else 0.0
        dup = d.get("duplication", {}).get("rate", 0) * 100
        return {
            "reads_in":     f"{ri:,}",
            "reads_out":    f"{ro:,}",
            "pct_pass":     f"{pp:.1f}%",
            "gc_pct":       f"{af.get('gc_content', 0) * 100:.1f}%",
            "q20_pct":      f"{af.get('q20_rate',   0) * 100:.1f}%",
            "q30_pct":      f"{af.get('q30_rate',   0) * 100:.1f}%",
            "bp_out":       f"{bp:,}",
            "dup_pct":      f"{dup:.1f}%",
            "est_coverage_small":   f"{cov_small:.0f}x"   if cov_small   else "N/A",
            "est_coverage_typical": f"{cov_typical:.0f}x" if cov_typical else "N/A",
            "est_coverage_large":   f"{cov_large:.0f}x"   if cov_large   else "N/A",
        }
    except Exception:
        return empty

def _check_warnings(m: dict, active_warnings: list[str]) -> None:
    def _pct(key: str) -> float:
        try:
            return float(str(m.get(key, "0")).rstrip("%"))
        except (ValueError, TypeError):
            return 0.0

    # Skip all warnings if fastp produced no output (empty JSON fallback)
    try:
        ri = int(str(m.get("reads_in", "0")).replace(",", ""))
    except (ValueError, TypeError):
        ri = 0
    if ri == 0:
        return

    pp = _pct("pct_pass")
    if pp < _PASS_WARN:
        msg = f"Pass rate {m['pct_pass']} < {_PASS_WARN}% — review input quality"
        log_warn(f"  ⚠ {msg}"); active_warnings.append(msg)

    q30 = _pct("q30_pct")
    if q30 < _Q30_WARN:
        msg = f"Q30 rate {m['q30_pct']} < {_Q30_WARN}% — assembly contiguity may suffer"
        log_warn(f"  ⚠ {msg}"); active_warnings.append(msg)

    dup = _pct("dup_pct")
    if dup > 50.0:
        log_info(
            f"  ℹ Duplication {m.get('dup_pct','?')} — elevated duplication is "
            "normal in high-coverage samples and purified viral preparations. "
            "Reads are NOT deduplicated (Head et al. 2014, Biotechniques 56:61)."
        )

    try:
        ro = int(str(m.get("reads_out", "0")).replace(",", ""))
    except (ValueError, TypeError):
        ro = 0
    if ro < _MIN_READS:
        msg = (
            f"Only {ro:,} reads retained (< {_MIN_READS:,}) — "
            "may be insufficient for CheckV Complete on a typical phage"
        )
        log_warn(f"  ⚠ {msg}"); active_warnings.append(msg)

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
        f"  GC={m.get('gc_pct','N/A')}  dup={m.get('dup_pct','N/A')}",
        f"  [cyan]Yield   :[/cyan] {m.get('bp_out','?')} bp",
        f"  [cyan]Coverage:[/cyan] [dim]~5 kb:[/dim] {m.get('est_coverage_small','?')}"
        f"  [dim]~50 kb:[/dim] {m.get('est_coverage_typical','?')}"
        f"  [dim]~200 kb:[/dim] {m.get('est_coverage_large','?')}",
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

_TSV_HEADERS = [
    "sample_id", "reads_in", "reads_out", "pct_pass",
    "q20_pct", "q30_pct", "gc_pct", "dup_pct", "bp_out", "est_coverage_50kb",
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
        "dup_pct":           m.get("dup_pct", ""),
        "bp_out":            m.get("bp_out", ""),
        "est_coverage_50kb": m.get("est_coverage", ""),
    }

    with open(path, "w") as f:
        f.write("\t".join(_TSV_HEADERS) + "\n")
        for row in rows.values():
            f.write("\t".join(row.get(h, "") for h in _TSV_HEADERS) + "\n")

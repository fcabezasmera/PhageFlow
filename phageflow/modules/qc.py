"""PhageFlow Module 01 — Quality control and trimming.

Designed for purified phage Illumina PE sequencing.

Tools:
    fastp   : adapter trimming, sliding-window quality trimming,
              low-complexity filtering (Chen et al. 2018, Genome Biology)
    FastQC  : per-read quality metrics (Andrews 2010)
    MultiQC : aggregate QC report (Ewels et al. 2016, Bioinformatics)

fastp parameters (literature-based for phage sequencing):
    --cut_right / --cut_right_window_size 4 / --cut_right_mean_quality 20
        Sliding-window 3' trimming — equivalent to Trimmomatic SLIDINGWINDOW
        (Bolger et al. 2014, Bioinformatics; recommended by Chen et al. 2018)
    --qualified_quality_phred 20
        Q20 threshold: 99% base call accuracy (Illumina quality guidelines)
    --unqualified_percent_limit 20
        Max 20% low-quality bases per read — more stringent than default 40%
        (Chen et al. 2018 recommend 20% for high-quality datasets)
    --length_required 75
        Minimum 75 bp for PE150 data ensures reliable k-mer coverage
        for phage genome assembly (Bankevich et al. 2012, J Comp Biol)
    --low_complexity_filter --complexity_threshold 30
        Removes homopolymer/repetitive reads enriched in phage tail fibers
        (Roux et al. 2019, eLife; Li et al. 2015, Bioinformatics)
    --n_base_limit 5
        Remove reads with >5 N calls (standard practice)
    --detect_adapter_for_pe
        Automatic adapter detection for paired-end data (Chen et al. 2018)
"""

from __future__ import annotations
import json
from pathlib import Path
from typing import Optional

from phageflow.utils.config import Config
from phageflow.utils.logger import log_step, log_info, log_ok, log_warn
from phageflow.utils.tools import require_tools, run_silent, human_size, mkdirs

STEP  = "01_qc"
TOOLS = ["fastp", "fastqc", "multiqc"]


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
    for label, path in [("r1", r1), ("r2", r2)]:
        if not path.exists():
            raise FileNotFoundError(f"[{sample_id}] {label} not found: {path}")

    r1_out = out_dir / f"{sample_id}_R1.fastq.gz"
    r2_out = out_dir / f"{sample_id}_R2.fastq.gz"
    json_f = rpt_dir / f"{sample_id}_fastp.json"
    html_f = rpt_dir / f"{sample_id}_fastp.html"
    log_f  = rpt_dir / f"{sample_id}_fastp.log"

    log_step(f"Module 01 — QC [{sample_id}] · {cfg.threads} threads")
    log_info(f"  R1 : {human_size(r1)}")
    log_info(f"  R2 : {human_size(r2)}")
    log_info(f"  Params: Q≥20 | len≥75bp | complexity≥30% | sliding-window 4bp")

    if r1_out.exists() and r1_out.stat().st_size > 0 and not force:
        log_info("  Already processed — skipping (use --force to re-run)")
    else:
        _run_fastp(sample_id, r1, r2, r1_out, r2_out, json_f, html_f, log_f, cfg.threads)
        _run_fastqc(r1_out, r2_out, rpt_dir, log_f)

    metrics = _parse_fastp_json(json_f)
    _save_tsv(sample_id, metrics, rpt_dir / "qc_summary.tsv")

    log_ok(
        f"  [{sample_id}] "
        f"{metrics['reads_in']} → {metrics['reads_out']} reads "
        f"({metrics['pct_pass']} pass | "
        f"GC={metrics['gc_pct']} Q20={metrics['q20_pct']} Q30={metrics['q30_pct']})"
    )

    # MultiQC — regenerate after each sample
    _run_multiqc(rpt_dir)

    log_ok(f"  Trimmed reads : {r1_out}")
    log_ok(f"  Report HTML   : {rpt_dir}/multiqc/multiqc_qc.html")
    log_step(f"Module 01 completed ✓  [{sample_id}]")
    log_info(f"Next: phageflow host-removal --sample-id {sample_id} "
             f"[--host-file|--accessions|--kraken-db]")

    return r1_out, r2_out


def _run_fastp(
    sample_id: str,
    r1: Path, r2: Path,
    r1_out: Path, r2_out: Path,
    json_f: Path, html_f: Path,
    log_f:  Path, threads: int,
) -> None:
    cmd = [
        "fastp",
        "--in1",  str(r1), "--in2",  str(r2),
        "--out1", str(r1_out), "--out2", str(r2_out),
        "--detect_adapter_for_pe",
        "--cut_right",
        "--cut_right_window_size",  "4",
        "--cut_right_mean_quality", "20",
        "--qualified_quality_phred",   "20",
        "--unqualified_percent_limit", "20",
        "--n_base_limit",    "5",
        "--length_required", "75",
        "--low_complexity_filter",
        "--complexity_threshold", "30",
        "--thread", str(threads),
        "--json",  str(json_f),
        "--html",  str(html_f),
        "--report_title", f"[{sample_id}] PhageFlow",
    ]
    try:
        run_silent(cmd, log_file=log_f)
        log_ok("  fastp OK")
    except Exception as e:
        log_warn(f"  fastp warning: {e}")


def _run_fastqc(r1: Path, r2: Path, rpt_dir: Path, log_f: Path) -> None:
    try:
        run_silent(
            ["fastqc", "--threads", "2", "--outdir", str(rpt_dir),
             "--quiet", str(r1), str(r2)],
            log_file=log_f, check=False
        )
    except Exception:
        pass


def _run_multiqc(rpt_dir: Path) -> None:
    mqc_dir = rpt_dir / "multiqc"
    mkdirs(mqc_dir)
    try:
        run_silent(
            ["multiqc", str(rpt_dir), "--outdir", str(mqc_dir),
             "--title", "PhageFlow QC", "--filename", "multiqc_qc", "--quiet"],
            log_file=rpt_dir / "multiqc.log", check=False
        )
    except Exception:
        pass


def _parse_fastp_json(json_f: Path) -> dict:
    base = {"reads_in":"N/A","reads_out":"N/A","pct_pass":"N/A",
            "gc_pct":"N/A","q20_pct":"N/A","q30_pct":"N/A"}
    if not json_f.exists():
        return base
    try:
        d  = json.load(open(json_f))
        bf = d["summary"]["before_filtering"]
        af = d["summary"]["after_filtering"]
        ri = bf["total_reads"]; ro = af["total_reads"]
        pp = ro / ri * 100 if ri else 0
        base.update({
            "reads_in":  str(ri), "reads_out": str(ro),
            "pct_pass":  f"{pp:.1f}%",
            "gc_pct":    f"{af.get('gc_content',0)*100:.1f}%",
            "q20_pct":   f"{af.get('q20_rate',0)*100:.1f}%",
            "q30_pct":   f"{af.get('q30_rate',0)*100:.1f}%",
        })
    except Exception:
        pass
    return base


def _save_tsv(sample_id: str, metrics: dict, path: Path) -> None:
    """Append or create the QC summary TSV."""
    headers = ["sample_id","reads_in","reads_out","pct_pass",
               "gc_pct","q20_pct","q30_pct"]
    # Load existing rows
    rows = {}
    if path.exists():
        with open(path) as f:
            f.readline()  # header
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if cols:
                    rows[cols[0]] = cols
    rows[sample_id] = [sample_id] + [metrics.get(h,"") for h in headers[1:]]
    with open(path,"w") as f:
        f.write("\t".join(headers)+"\n")
        for row in rows.values():
            f.write("\t".join(row)+"\n")

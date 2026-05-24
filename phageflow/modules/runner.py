"""PhageFlow master runner — executes all modules for a list of (r1, r2) pairs.

Called by ``phageflow run``.  Receives pairs auto-discovered from a directory
(or supplied explicitly via --r1/--r2) — no samples.tsv required.
"""

from __future__ import annotations
import traceback
from pathlib import Path
from typing import Optional

from phageflow.utils.config import Config
from phageflow.utils.logger import log_step, log_info, log_ok, log_warn, log_error
from phageflow.cli import _infer_sample_id

# Module order for --from-module resume
_MODULE_ORDER = [
    "qc",
    "host-removal",
    "assembly",
    "viral-id",
    "quality",
    "annotate",
    "safety",
    "report",
]


def run(
    cfg:         Config,
    pairs:       list[tuple[Path, Path]],
    force:       bool = False,
    host_file:   Optional[Path] = None,
    accessions:  Optional[list[str]] = None,
    from_module: Optional[str] = None,
) -> None:
    """Execute the full pipeline for every (r1, r2) pair.

    Parameters
    ----------
    cfg         : loaded Config object
    pairs       : list of (r1, r2) Path tuples; sample_id inferred from r1
    force       : pass --force to every module
    host_file   : passed to host-removal
    accessions  : passed to host-removal
    from_module : skip all modules that come before this one
    """
    start_idx = 0
    if from_module:
        try:
            start_idx = _MODULE_ORDER.index(from_module.lower())
        except ValueError:
            log_warn(f"  Unknown --from-module '{from_module}' — starting from qc")

    log_step(f"PhageFlow run — {len(pairs)} sample(s)")
    if from_module:
        log_info(f"  Resuming from module : {from_module}")

    results: dict[str, str] = {}  # sample_id → "ok" | "failed"

    for r1, r2 in pairs:
        sample_id = _infer_sample_id(r1)
        log_step(f"Sample : {sample_id}")
        log_info(f"  R1 : {r1}")
        log_info(f"  R2 : {r2}")

        try:
            _run_sample(
                cfg, sample_id, r1, r2,
                force, host_file, accessions, start_idx,
            )
            results[sample_id] = "ok"
            log_ok(f"  Sample {sample_id} completed successfully")
        except Exception as exc:
            log_error(f"  Sample {sample_id} FAILED: {exc}")
            traceback.print_exc()
            results[sample_id] = f"failed: {exc}"

    # Summary
    log_step("Run summary")
    for sid, status in results.items():
        if status == "ok":
            log_ok(f"  {sid:30s} OK")
        else:
            log_warn(f"  {sid:30s} {status}")


def _run_sample(
    cfg:         Config,
    sample_id:   str,
    r1:          Path,
    r2:          Path,
    force:       bool,
    host_file:   Optional[Path],
    accessions:  Optional[list[str]],
    start_idx:   int,
) -> None:
    """Run all pipeline modules for a single sample."""

    def _skip(module: str) -> bool:
        return _MODULE_ORDER.index(module) < start_idx

    # ── 01 QC ─────────────────────────────────────────────────────────────────
    if not _skip("qc"):
        from phageflow.modules.qc import run as qc_run
        r1_qc, r2_qc = qc_run(cfg, sample_id=sample_id, r1=r1, r2=r2, force=force)
    else:
        r1_qc = cfg.results("01_qc") / f"{sample_id}_R1.fastq.gz"
        r2_qc = cfg.results("01_qc") / f"{sample_id}_R2.fastq.gz"

    # ── 02 Host removal ───────────────────────────────────────────────────────
    if not _skip("host-removal"):
        from phageflow.modules.host_removal import run as hr_run
        r1_hr, r2_hr, s1 = hr_run(
            cfg, sample_id=sample_id,
            r1=r1_qc, r2=r2_qc,
            force=force,
            host_file=host_file,
            accessions=accessions,
        )
    else:
        r1_hr = cfg.results("02_host_removal") / f"{sample_id}_R1.fastq.gz"
        r2_hr = cfg.results("02_host_removal") / f"{sample_id}_R2.fastq.gz"
        s1    = cfg.results("02_host_removal") / f"{sample_id}_singletons.fastq.gz"

    # ── 03 Assembly ───────────────────────────────────────────────────────────
    if not _skip("assembly"):
        from phageflow.modules.assembly import run as asm_run
        contigs = asm_run(
            cfg, sample_id=sample_id,
            r1=r1_hr, r2=r2_hr,
            s1=s1 if s1.exists() and s1.stat().st_size > 0 else None,
            force=force,
        )
    else:
        contigs = cfg.results("03_assembly") / f"{sample_id}_contigs.fasta"

    # ── 04 Viral ID ───────────────────────────────────────────────────────────
    if not _skip("viral-id"):
        from phageflow.modules.viral_id import run as vid_run
        virus_fna, _ = vid_run(cfg, sample_id=sample_id, contigs=contigs, force=force)
    else:
        virus_fna = cfg.results("04_viral_id") / f"{sample_id}_virus.fna"

    # ── 05 Quality ────────────────────────────────────────────────────────────
    if not _skip("quality"):
        from phageflow.modules.quality import run as qual_run
        ann_ready, _ = qual_run(cfg, sample_id=sample_id, virus_fna=virus_fna, force=force)
    else:
        ann_ready = (
            cfg.results("05_quality") / sample_id / "annotation_ready" / "phages"
        )

    # ── 06 Annotate (per candidate) ───────────────────────────────────────────
    if not _skip("annotate"):
        candidates = sorted(ann_ready.glob("*.fasta")) if ann_ready.exists() else []
        if not candidates:
            log_warn(f"  [{sample_id}] No annotation_ready candidates found — skipping annotate")
        for genome in candidates:
            from phageflow.modules.annotate import run as ann_run
            ann_run(cfg, sample_id=sample_id, genome=genome, force=force)

    # ── 07 Safety (per candidate) ─────────────────────────────────────────────
    if not _skip("safety"):
        candidates = sorted(ann_ready.glob("*.fasta")) if ann_ready.exists() else []
        for genome in candidates:
            try:
                from phageflow.modules.safety import run as safe_run
                safe_run(cfg, sample_id=sample_id, genome=genome, force=force)
            except Exception as e:
                log_warn(f"  [{sample_id}] safety skipped for {genome.stem}: {e}")

    # ── 08 Report (per candidate) ─────────────────────────────────────────────
    if not _skip("report"):
        candidates = sorted(ann_ready.glob("*.fasta")) if ann_ready.exists() else []
        for genome in candidates:
            try:
                from phageflow.modules.report import run as rpt_run
                rpt_run(cfg, sample_id=sample_id,
                        candidate_id=genome.stem, force=force)
            except Exception as e:
                log_warn(f"  [{sample_id}] report skipped for {genome.stem}: {e}")

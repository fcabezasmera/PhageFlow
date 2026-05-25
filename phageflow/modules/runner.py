"""PhageFlow master runner — executes all modules for a list of (r1, r2) pairs.

Called by ``phageflow run``. Receives pairs auto-discovered from a directory
(or supplied explicitly via --r1/--r2) — no samples.tsv required.

Pipeline order
--------------
  01  qc              : adapter trimming, quality filter, poly-X removal
  02  host-removal    : remove host reads (bwa-mem2 or Kraken2)
  03  assembly        : SPAdes + MEGAHIT + cd-hit-est NR
  03b coverage        : per-contig coverage profiling (CoverM)
  04  viral-id        : geNomad viral classification
  05  quality         : CheckV quality tiers + dereplication
  06  annotate        : Pharokka + Phold per candidate
  07  resistance      : AMR + virulence + ACR + defense screening

Resume support
--------------
  --from-module NAME : skip all modules before NAME.
  Intermediate outputs are inferred from standard paths.
  Example: resume from assembly after fixing host-removal parameters:
    phageflow run --from-module assembly --r1 ... --r2 ... -o ...

Error handling
--------------
  Per-sample errors are caught and logged without stopping other samples.
  Per-candidate errors (annotate, resistance) are caught individually.
  Final summary shows OK / FAILED status per sample.
"""

from __future__ import annotations
import traceback
from pathlib import Path
from typing import Optional

from phageflow.utils.config import Config
from phageflow.utils.logger import (
    log_step, log_info, log_ok, log_warn, log_error,
)

# Module order for --from-module resume
_MODULE_ORDER = [
    "qc",
    "host-removal",
    "assembly",
    "coverage",
    "viral-id",
    "quality",
    "annotate",
    "resistance",
]


def run(
    cfg:         Config,
    pairs:       list[tuple[Path, Path]],
    force:       bool                = False,
    host_file:   Optional[Path]      = None,
    accessions:  Optional[list[str]] = None,
    from_module: Optional[str]       = None,
) -> dict[str, str]:
    """Execute the full pipeline for every (r1, r2) pair.

    Parameters
    ----------
    cfg         : loaded Config object
    pairs       : list of (r1, r2) Path tuples; sample_id inferred from r1
    force       : pass --force to every module
    host_file   : passed to host-removal (local FASTA reference)
    accessions  : passed to host-removal (NCBI accessions to download)
    from_module : resume from this module, skipping earlier ones

    Returns
    -------
    dict mapping sample_id → "ok" | "failed: <reason>"
    """
    from phageflow.cli import _infer_sample_id

    start_idx = 0
    if from_module:
        key = from_module.lower()
        if key in _MODULE_ORDER:
            start_idx = _MODULE_ORDER.index(key)
        else:
            log_warn(
                f"  Unknown --from-module '{from_module}'. "
                f"Valid options: {', '.join(_MODULE_ORDER)}"
            )

    log_step(f"PhageFlow run — {len(pairs)} sample(s)")
    if from_module and start_idx > 0:
        log_info(f"  Resuming from  : {from_module} (skipping {_MODULE_ORDER[:start_idx]})")
    if host_file:
        log_info(f"  Host reference : {host_file}")
    if accessions:
        log_info(f"  Accessions     : {', '.join(accessions)}")

    results: dict[str, str] = {}

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
            log_ok(f"  Sample {sample_id} — pipeline complete")
        except Exception as exc:
            log_error(f"  Sample {sample_id} FAILED: {exc}")
            traceback.print_exc()
            results[sample_id] = f"failed: {exc}"

    # ── Final summary ─────────────────────────────────────────────────────────
    log_step("PhageFlow run — summary")
    n_ok   = sum(1 for v in results.values() if v == "ok")
    n_fail = len(results) - n_ok
    for sid, status in results.items():
        if status == "ok":
            log_ok(f"  {sid:40s} ✓ OK")
        else:
            log_warn(f"  {sid:40s} ✗ {status}")
    log_info(f"  {n_ok}/{len(results)} samples completed successfully")
    if n_fail:
        log_warn(f"  {n_fail} sample(s) failed — check logs above")

    return results


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
    """Execute all pipeline modules for a single sample."""

    def _skip(module: str) -> bool:
        """Return True if this module should be skipped (resume mode)."""
        try:
            return _MODULE_ORDER.index(module) < start_idx
        except ValueError:
            return False

    # ── 01 QC ─────────────────────────────────────────────────────────────────
    if not _skip("qc"):
        from phageflow.modules.qc import run as qc_run
        r1_qc, r2_qc = qc_run(
            cfg, sample_id=sample_id, r1=r1, r2=r2, force=force,
        )
    else:
        log_info("  [skip] qc")
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
        log_info("  [skip] host-removal")
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
        log_info("  [skip] assembly")
        contigs = cfg.results("03_assembly") / f"{sample_id}_contigs.fasta"

    # ── 03b Coverage ──────────────────────────────────────────────────────────
    if not _skip("coverage"):
        try:
            from phageflow.modules.coverage import run as cov_run
            cov_run(
                cfg, sample_id=sample_id,
                r1=r1_hr, r2=r2_hr,
                contigs=contigs,
                force=force,
            )
        except Exception as e:
            log_warn(f"  [coverage] non-fatal error: {e} — continuing")
    else:
        log_info("  [skip] coverage")

    # ── 04 Viral ID ───────────────────────────────────────────────────────────
    if not _skip("viral-id"):
        from phageflow.modules.viral_id import run as vid_run
        virus_fna, _ = vid_run(
            cfg, sample_id=sample_id, contigs=contigs, force=force,
        )
    else:
        log_info("  [skip] viral-id")
        virus_fna = cfg.results("04_viral_id") / f"{sample_id}_virus.fna"

    # ── 05 Quality ────────────────────────────────────────────────────────────
    if not _skip("quality"):
        from phageflow.modules.quality import run as qual_run
        ann_ready, _ = qual_run(
            cfg, sample_id=sample_id, virus_fna=virus_fna, force=force,
        )
    else:
        log_info("  [skip] quality")
        ann_ready = (
            cfg.results("05_quality") / sample_id / "annotation_ready" / "phages"
        )

    # Candidates for per-genome modules
    candidates = sorted(ann_ready.glob("*.fasta")) if ann_ready.exists() else []
    if not candidates:
        log_warn(f"  [{sample_id}] No annotation_ready candidates — skipping annotate + resistance")
        return

    log_info(f"  [{sample_id}] {len(candidates)} candidate(s) → annotate + resistance")

    # ── 06 Annotate (per candidate) ───────────────────────────────────────────
    if not _skip("annotate"):
        from phageflow.modules.annotate import run as ann_run
        for genome in candidates:
            try:
                ann_run(cfg, sample_id=sample_id, genome=genome, force=force)
            except Exception as e:
                log_warn(f"  [annotate] {genome.stem}: {e} — continuing")
    else:
        log_info("  [skip] annotate")

    # ── 07 Resistance (per candidate) ─────────────────────────────────────────
    if not _skip("resistance"):
        from phageflow.modules.resistance import run as res_run
        for genome in candidates:
            try:
                res_run(
                    cfg, sample_id=sample_id,
                    candidate_id=genome.stem,
                    force=force,
                )
            except Exception as e:
                log_warn(f"  [resistance] {genome.stem}: {e} — continuing")
    else:
        log_info("  [skip] resistance")

"""PhageFlow Module 08 — Lifecycle prediction.

Tool:
    BACPHLIP v0.9+ (Hockenberry & Wilke 2021, PeerJ):
        Random forest classifier trained on 62 lifestyle-associated
        protein families (HMM profiles derived from over 1,000 manually
        curated phage genomes).

        Input  : phage genome FASTA (single or multi-contig)
        Output : {genome}.bacphlip — TSV with Virulent / Temperate scores

Interpretation:
    Virulent  score ≥ 0.9 : confident lytic phage (preferred for therapy)
    Temperate score ≥ 0.9 : confident lysogen (requires expert review)
    Neither   score ≥ 0.9 : ambiguous — insufficient lifestyle markers
                            (common in highly novel or fragmented genomes)

Biological context:
    BACPHLIP detects lifestyle through lysogeny-associated HMMs
    (integrases, CI repressors, Cro proteins, excisionases) and
    lysis-associated HMMs (holins, endolysins, spanins).
    Hockenberry & Wilke 2021: 98% accuracy on complete genomes;
    performance drops for fragmented or highly divergent sequences.

    Complement with:
      - CheckV completeness (Module 05): low-completeness genomes → lower confidence
      - Safety module (Module 07): integrase presence → temperate signal
      - Pharokka annotation (Module 06): manual check of CI repressor / Cro genes

Activation note:
    BACPHLIP requires its own conda environment (bacphlip_env).
    Activate before running: conda activate bacphlip_env
"""

from __future__ import annotations
import shutil
from pathlib import Path

from phageflow.utils.config import Config
from phageflow.utils.logger import log_step, log_info, log_ok, log_warn, log_error
from phageflow.utils.tools import require_tools, run_silent, human_size, mkdirs, fasta_stats

STEP  = "08_lifecycle"
TOOLS = ["bacphlip"]

_CONFIDENCE_THRESHOLD = 0.9


def run(
    cfg:       Config,
    sample_id: str,
    genome:    Path,
    force:     bool = False,
) -> dict:
    """
    Predict lifecycle (virulent / temperate) for a phage genome.

    Parameters
    ----------
    genome : phage genome FASTA (ideally HQ / complete from quality step)

    Returns
    -------
    result : dict with 'virulent_score', 'temperate_score', 'prediction', 'confidence'
    """
    require_tools(*TOOLS)

    genome  = Path(genome)
    out_dir = cfg.results(STEP)
    rpt_dir = cfg.reports(STEP)
    mkdirs(out_dir, rpt_dir)

    log_step(f"Module 08 — lifecycle prediction [{sample_id}] (BACPHLIP)")

    if not genome.exists():
        log_error(f"  Genome not found: {genome}")
        raise FileNotFoundError(genome)

    s = fasta_stats(genome)
    log_info(f"  Genome : {human_size(genome)} | {s['n']} contig(s) | "
             f"total={s['total_bp']}bp")
    log_info("  Tool: BACPHLIP RF classifier (Hockenberry & Wilke 2021, PeerJ)")
    log_info(f"  Confidence threshold: ≥{_CONFIDENCE_THRESHOLD}")

    # BACPHLIP writes output next to the input FASTA by default
    # We copy the genome to the output directory to keep results organised
    genome_copy = out_dir / f"{sample_id}.fasta"
    bacphlip_out = out_dir / f"{sample_id}.fasta.bacphlip"

    if bacphlip_out.exists() and bacphlip_out.stat().st_size > 0 and not force:
        log_info("  Already processed — skipping (use --force to re-run)")
    else:
        shutil.copy(genome, genome_copy)
        log_info("  [BACPHLIP] running classifier...")
        try:
            run_silent(
                ["bacphlip", "-i", str(genome_copy), "--multi_fasta"],
                log_file=rpt_dir / f"{sample_id}_bacphlip.log",
            )
            log_ok("  [BACPHLIP] OK")
        except Exception as e:
            log_warn(f"  [BACPHLIP] warning: {e}")

    result = _parse_bacphlip(sample_id, bacphlip_out, s)
    _print_result(sample_id, result, s)

    # Cross-check with safety module integrase output
    safety_tsv = cfg.reports("07_safety") / "safety_summary.tsv"
    _crosscheck_integrase(sample_id, result, safety_tsv)

    _save_tsv(result, rpt_dir / "lifecycle_summary.tsv")

    log_step(f"Module 08 completed ✓  [{sample_id}]")
    log_ok(f"  Output: {bacphlip_out}")
    log_info("Pipeline complete. Check reports/ for summaries.")

    return result


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def _parse_bacphlip(sample_id: str, bacphlip_out: Path, fasta_s: dict) -> dict:
    """
    Parse BACPHLIP output.

    BACPHLIP output format (tab-separated):
        contig_id  Virulent  Temperate
        (one row per contig in multi_fasta mode)

    For multi-contig inputs, aggregate by weighted average (by contig length)
    — longer contigs are more likely to carry lifestyle markers.
    """
    result = {
        "sample_id":       sample_id,
        "n_contigs":       fasta_s.get("n", 0),
        "genome_bp":       fasta_s.get("total_bp", 0),
        "virulent_score":  None,
        "temperate_score": None,
        "prediction":      "unknown",
        "confidence":      "low",
        "note":            "",
    }

    if not bacphlip_out.exists():
        result["note"] = "BACPHLIP output not found"
        return result

    rows = []
    try:
        with open(bacphlip_out) as f:
            header = f.readline().strip().split("\t")
            col = {h: i for i, h in enumerate(header)}
            for line in f:
                if not line.strip():
                    continue
                parts = line.strip().split("\t")
                try:
                    vir  = float(parts[col.get("Virulent",  1)])
                    temp = float(parts[col.get("Temperate", 2)])
                    rows.append((vir, temp))
                except (ValueError, IndexError):
                    pass
    except Exception as e:
        result["note"] = f"Parse error: {e}"
        return result

    if not rows:
        result["note"] = "No predictions in BACPHLIP output"
        return result

    # Average across contigs (unweighted for simplicity; adequate for complete genomes)
    avg_vir  = sum(r[0] for r in rows) / len(rows)
    avg_temp = sum(r[1] for r in rows) / len(rows)

    result["virulent_score"]  = round(avg_vir,  4)
    result["temperate_score"] = round(avg_temp, 4)

    if avg_vir >= _CONFIDENCE_THRESHOLD:
        result["prediction"]  = "Virulent"
        result["confidence"]  = "high"
    elif avg_temp >= _CONFIDENCE_THRESHOLD:
        result["prediction"]  = "Temperate"
        result["confidence"]  = "high"
    elif avg_vir > avg_temp:
        result["prediction"]  = "Virulent (ambiguous)"
        result["confidence"]  = "low"
    else:
        result["prediction"]  = "Temperate (ambiguous)"
        result["confidence"]  = "low"

    if len(rows) > 1:
        result["note"] = f"Averaged over {len(rows)} contigs"

    return result


def _print_result(sample_id: str, result: dict, fasta_s: dict) -> None:
    pred = result["prediction"]
    vir  = result.get("virulent_score")
    temp = result.get("temperate_score")
    conf = result.get("confidence", "low")

    score_str = ""
    if vir is not None:
        score_str = f"Virulent={vir:.3f} | Temperate={temp:.3f}"

    if "Virulent" in pred and conf == "high":
        log_ok(f"  [{sample_id}] ✓ {pred} | {score_str}")
    elif "Temperate" in pred and conf == "high":
        log_warn(f"  [{sample_id}] ⚠  {pred} | {score_str}")
    else:
        log_warn(f"  [{sample_id}] ? {pred} | {score_str} (ambiguous)")

    if fasta_s.get("total_bp", 0) < 10_000:
        log_warn(
            "  Genome <10 kb — BACPHLIP confidence reduced for short contigs. "
            "Increase completeness threshold in CheckV step."
        )


def _crosscheck_integrase(sample_id: str, result: dict, safety_tsv: Path) -> None:
    """Warn if lifecycle prediction conflicts with integrase detection."""
    if not safety_tsv.exists():
        return
    try:
        with open(safety_tsv) as f:
            header = f.readline().strip().split("\t")
            col = {h: i for i, h in enumerate(header)}
            for line in f:
                parts = line.strip().split("\t")
                if not parts or parts[0] != sample_id:
                    continue
                n_integ = int(parts[col.get("integrase_genes", 3)] or 0)
                if n_integ > 0 and result.get("prediction", "").startswith("Virulent"):
                    log_warn(
                        f"  ⚠  Conflict: BACPHLIP predicts Virulent but "
                        f"{n_integ} integrase/lysogeny gene(s) detected "
                        f"(Safety module). Manual review recommended."
                    )
    except Exception:
        pass


# ---------------------------------------------------------------------------
# Summary TSV
# ---------------------------------------------------------------------------

def _save_tsv(result: dict, path: Path) -> None:
    headers = [
        "sample_id", "n_contigs", "genome_bp",
        "virulent_score", "temperate_score",
        "prediction", "confidence", "note",
    ]
    rows = {}
    if path.exists():
        with open(path) as f:
            f.readline()
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if cols:
                    rows[cols[0]] = cols

    rows[result["sample_id"]] = [str(result.get(h, "")) for h in headers]

    with open(path, "w") as f:
        f.write("\t".join(headers) + "\n")
        for row in rows.values():
            f.write("\t".join(row) + "\n")

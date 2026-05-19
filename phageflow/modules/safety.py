"""PhageFlow Module 07 — Biosafety screening.

Screens phage genomes for:

  1. Antimicrobial Resistance Genes (ARGs) — CARD database
     ABRicate + CARD (Alcock et al. 2023, Nucleic Acids Research):
       Identifies AMR genes that could be horizontally transferred to
       bacterial hosts during phage therapy (phage-mediated ARG transfer;
       Enault et al. 2017, ISME J; Colavecchio et al. 2017, Front Microbiol).
       Threshold: --minid 80 --mincov 80 (standard for ARG detection;
       McArthur et al. 2013, Antimicrobial Agents and Chemotherapy).

  2. Virulence Factors (VFs) — VFDB database
     ABRicate + VFDB (Liu et al. 2022, Nucleic Acids Research):
       Detects virulence factor genes that may enhance pathogen fitness.
       Lysogenic phages can transduce virulence factors (Brüssow et al. 2004,
       Microbiol Mol Biol Rev; Boyd & Brussow 2002, Trends Microbiol).
       Threshold: --minid 80 --mincov 80.

  3. Integrase / Lysogeny Genes
     Detected from Pharokka annotation (PHROG categories):
       - Integration/excision (PHROG category: "integration and excision")
       - Lysogeny module (PHROG category: "lysogeny")
       Presence suggests potential for chromosomal integration and
       lysogenic lifestyle (phage therapy risk factor;
       Reardon 2014, Nature; Chan et al. 2013, Science Trans Med).

Safety verdict:
  PASS   : No ARG, no VF, no integrase detected
  CAUTION: Integrase or VF detected (warrants expert review)
  FAIL   : ARG detected (biosafety concern; exclude from therapy use)
"""

from __future__ import annotations
import csv
from pathlib import Path
from typing import List

from phageflow.utils.config import Config
from phageflow.utils.logger import log_step, log_info, log_ok, log_warn, log_error
from phageflow.utils.tools import require_tools, run_silent, human_size, mkdirs

STEP  = "07_safety"
TOOLS = ["abricate"]

_ABRICATE_MINID  = "80"
_ABRICATE_MINCOV = "80"

# PHROG categories that indicate lysogeny / integration potential
_LYSOGENY_PHROG_CATEGORIES = {
    "integration and excision",
    "lysogeny",
    "transcription regulation",   # CI repressor, Cro
}


def run(
    cfg:       Config,
    sample_id: str,
    genome:    Path,
    force:     bool = False,
) -> dict:
    """
    Screen a phage genome for ARGs, virulence factors, and lysogeny genes.

    Returns
    -------
    result : dict with keys 'arg_hits', 'vf_hits', 'integrase', 'verdict'
    """
    require_tools(*TOOLS)

    genome  = Path(genome)
    out_dir = cfg.results(STEP)
    rpt_dir = cfg.reports(STEP)
    mkdirs(out_dir, rpt_dir)

    log_step(f"Module 07 — biosafety screening [{sample_id}]")

    if not genome.exists():
        log_error(f"  Genome not found: {genome}")
        raise FileNotFoundError(genome)

    log_info(f"  Genome : {human_size(genome)}")
    log_info(f"  Thresholds: --minid {_ABRICATE_MINID} --mincov {_ABRICATE_MINCOV}")
    log_info("  Databases: CARD (ARG), VFDB (virulence)")

    # ── CARD (AMR) ────────────────────────────────────────────────────────────
    card_tsv = out_dir / f"{sample_id}_CARD.tsv"

    if card_tsv.exists() and card_tsv.stat().st_size > 0 and not force:
        log_info("  [CARD] already exists — skipping")
    else:
        log_info("  [CARD] scanning for AMR genes...")
        try:
            run_silent(
                f"abricate --db card --minid {_ABRICATE_MINID} "
                f"--mincov {_ABRICATE_MINCOV} {genome} > {card_tsv}",
                log_file=rpt_dir / f"{sample_id}_abricate.log",
                shell=True,
            )
        except Exception as e:
            log_warn(f"  [CARD] warning: {e}")

    card_hits = _parse_abricate(card_tsv)
    if card_hits:
        log_warn(f"  [CARD] ⚠  {len(card_hits)} ARG hit(s) detected:")
        for h in card_hits:
            log_warn(f"     · {h['gene']} | {h['pct_id']}%id | "
                     f"{h['pct_cov']}%cov | {h['product']}")
    else:
        log_ok("  [CARD] No ARGs detected ✓")

    # ── VFDB (virulence factors) ───────────────────────────────────────────────
    vfdb_tsv = out_dir / f"{sample_id}_VFDB.tsv"

    if vfdb_tsv.exists() and vfdb_tsv.stat().st_size > 0 and not force:
        log_info("  [VFDB] already exists — skipping")
    else:
        log_info("  [VFDB] scanning for virulence factors...")
        try:
            run_silent(
                f"abricate --db vfdb --minid {_ABRICATE_MINID} "
                f"--mincov {_ABRICATE_MINCOV} {genome} > {vfdb_tsv}",
                log_file=rpt_dir / f"{sample_id}_abricate.log",
                shell=True,
            )
        except Exception as e:
            log_warn(f"  [VFDB] warning: {e}")

    vf_hits = _parse_abricate(vfdb_tsv)
    if vf_hits:
        log_warn(f"  [VFDB] ⚠  {len(vf_hits)} virulence factor(s) detected:")
        for h in vf_hits:
            log_warn(f"     · {h['gene']} | {h['pct_id']}%id | "
                     f"{h['pct_cov']}%cov | {h['product']}")
    else:
        log_ok("  [VFDB] No virulence factors detected ✓")

    # ── Integrase / Lysogeny detection (from Pharokka output) ─────────────────
    pharokka_tsv = (cfg.results("06_annotation") / sample_id
                    / "pharokka" / f"{sample_id}_cds_functions.tsv")

    integrase_genes: List[dict] = []
    if pharokka_tsv.exists():
        integrase_genes = _parse_integrase(pharokka_tsv)
        if integrase_genes:
            log_warn(f"  [Integrase] ⚠  {len(integrase_genes)} lysogeny-related gene(s):")
            for g in integrase_genes:
                log_warn(f"     · {g['gene_id']} | {g['category']} | {g['product']}")
        else:
            log_ok("  [Integrase] No integrase/lysogeny genes detected ✓")
    else:
        log_info(
            "  [Integrase] Pharokka TSV not found — "
            "run 'phageflow annotate' first for integrase detection"
        )

    # ── Safety verdict ────────────────────────────────────────────────────────
    verdict = _compute_verdict(card_hits, vf_hits, integrase_genes)

    if verdict == "PASS":
        log_ok(f"  Safety verdict: ✓ PASS")
    elif verdict == "CAUTION":
        log_warn(f"  Safety verdict: ⚠  CAUTION — review before therapeutic use")
    else:
        log_warn(f"  Safety verdict: ✗ FAIL — ARG detected, exclude from therapy")

    result = {
        "sample_id":        sample_id,
        "arg_hits":         len(card_hits),
        "vf_hits":          len(vf_hits),
        "integrase_genes":  len(integrase_genes),
        "verdict":          verdict,
    }

    _save_tsv(result, card_hits, vf_hits, integrase_genes,
              rpt_dir / "safety_summary.tsv")

    log_step(f"Module 07 completed ✓  [{sample_id}]")
    log_info("Next: phageflow lifecycle --sample-id <id> --genome <genome.fasta>")

    return result


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def _parse_abricate(tsv: Path) -> List[dict]:
    """Parse abricate output TSV. Returns list of hit dicts."""
    hits = []
    if not tsv.exists():
        return hits
    try:
        with open(tsv) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if row.get("#FILE", "").startswith("#"):
                    continue
                hits.append({
                    "contig":   row.get("SEQUENCE", ""),
                    "gene":     row.get("GENE", ""),
                    "pct_id":   row.get("%IDENTITY", ""),
                    "pct_cov":  row.get("%COVERAGE", ""),
                    "product":  row.get("PRODUCT", ""),
                    "database": row.get("DATABASE", ""),
                    "resistance": row.get("RESISTANCE", ""),
                })
    except Exception:
        pass
    return hits


def _parse_integrase(tsv: Path) -> List[dict]:
    """Detect lysogeny-related genes from Pharokka CDS functions TSV."""
    genes = []
    try:
        with open(tsv) as f:
            header_line = f.readline().strip().split("\t")
            col = {h: i for i, h in enumerate(header_line)}
            for line in f:
                if not line.strip():
                    continue
                row = line.strip().split("\t")
                cat = (row[col["category"]].lower()
                       if "category" in col and len(row) > col["category"] else "")
                if any(kw in cat for kw in _LYSOGENY_PHROG_CATEGORIES):
                    genes.append({
                        "gene_id": row[col.get("gene_id", 0)] if row else "",
                        "category": cat,
                        "product": (row[col["product"]]
                                    if "product" in col and len(row) > col["product"] else ""),
                    })
    except Exception:
        pass
    return genes


def _compute_verdict(card_hits, vf_hits, integrase_genes) -> str:
    if card_hits:
        return "FAIL"
    if vf_hits or integrase_genes:
        return "CAUTION"
    return "PASS"


# ---------------------------------------------------------------------------
# Summary TSV
# ---------------------------------------------------------------------------

def _save_tsv(result, card_hits, vf_hits, integrase_genes, path):
    """Write/update the safety summary TSV."""
    # Main summary
    headers = [
        "sample_id", "arg_hits", "vf_hits",
        "integrase_genes", "verdict",
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

    # Detailed hit files (overwrite per sample)
    sid = result["sample_id"]
    detail_path = path.parent / f"{sid}_safety_details.tsv"
    detail_headers = ["sample_id", "type", "gene", "contig",
                      "pct_id", "pct_cov", "product"]
    with open(detail_path, "w") as f:
        f.write("\t".join(detail_headers) + "\n")
        for h in card_hits:
            f.write("\t".join([
                sid, "ARG", h["gene"], h["contig"],
                h["pct_id"], h["pct_cov"], h["product"],
            ]) + "\n")
        for h in vf_hits:
            f.write("\t".join([
                sid, "VF", h["gene"], h["contig"],
                h["pct_id"], h["pct_cov"], h["product"],
            ]) + "\n")
        for g in integrase_genes:
            f.write("\t".join([
                sid, "Integrase/Lysogeny", g["gene_id"],
                "", "", "", g["product"],
            ]) + "\n")

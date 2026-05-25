"""PhageFlow Module 07 — Resistance, virulence and defense element screening.

Goal: parse outputs from Pharokka (sequence-based) and Phold (structure-based)
to identify AMR genes, virulence factors, anti-CRISPR proteins, defense systems,
and toxin-antitoxin systems in annotation-ready phage genomes.

No additional tools are required — all data is parsed from existing annotation
outputs produced by Module 06 (annotate.py).

Data sources
------------
Pharokka (sequence-based, MMseqs2):
  _cds_final_merged_output.tsv  : per-CDS CARD and VFDB hits
    CARD columns : CARD_hit, CARD_alnScore, CARD_seqIdentity, CARD_eVal,
                   ARO_Accession, CARD_short_name, AMR_Gene_Family,
                   Drug_Class, Resistance_Mechanism
    VFDB columns : vfdb_hit, vfdb_alnScore, vfdb_seqIdentity, vfdb_eVal,
                   vfdb_short_name, vfdb_description, vfdb_species

Phold (structure-based, Foldseek):
  sub_db_tophits/card_cds_predictions.tsv         : AMR (structural)
  sub_db_tophits/vfdb_cds_predictions.tsv         : virulence (structural)
  sub_db_tophits/acr_cds_predictions.tsv          : anti-CRISPR proteins
  sub_db_tophits/defensefinder_cds_predictions.tsv: defense systems
  sub_db_tophits/netflax_cds_predictions.tsv      : toxin-antitoxin systems

Dual-evidence approach
----------------------
A finding is reported at two confidence levels:
  HIGH   : detected by BOTH Pharokka (sequence) AND Phold (structure)
  MEDIUM : detected by Pharokka sequence search only
  LOW    : detected by Phold structure search only (below sequence threshold)

This mirrors the annotation tier logic from annotate.py and reduces both
false positives (structural hits without sequence support) and false negatives
(divergent sequences missed by MMseqs2 but detected via 3Di structure).

AMR in phages
-------------
Phages can carry and mobilise AMR genes via generalised transduction
(Colomer-Lluch et al. 2011, Antimicrob Agents Chemother 55:4247;
 Enright et al. 2002, Proc Natl Acad Sci 99:7687).
AMR detection in phage genomes is a biosafety screening step — the presence
of an AMR gene is a red flag for therapeutic use.

Anti-CRISPR proteins
--------------------
Anti-CRISPR (Acr) proteins inhibit CRISPR-Cas immune systems in bacterial
hosts. They are relevant for phage therapy (phage resistance to host immunity)
and for understanding phage-host dynamics.
Bondy-Denomy et al. 2013 (Nature 493:429); Pawluk et al. 2016 (Science 351:aad8405).

Defense systems
---------------
DefenseFinder detects >150 prokaryotic anti-phage defense systems.
Reported here as context for phage-host interaction interpretation.
Abby et al. 2022 (Nucleic Acids Res 50:W245).

Toxin-antitoxin (NetFlax)
--------------------------
TA systems in phage genomes can contribute to addiction modules and
host manipulation. NetFlax detects phage TA systems via structural homology.
"""

from __future__ import annotations
import csv
from pathlib import Path
from typing import Optional

from rich.panel import Panel
from rich.progress import (
    Progress, SpinnerColumn, TextColumn,
    BarColumn, MofNCompleteColumn, TimeElapsedColumn,
)

from phageflow.utils.config import Config
from phageflow.utils.logger import (
    log_step, log_info, log_ok, log_warn, console,
)
from phageflow.utils.tools import mkdirs

STEP = "07_resistance"

# Confidence thresholds for reporting
_EVALUE_THRESH    = 1e-5   # sequence-based (Pharokka)
_STRUCT_EVALUE    = 1e-3   # structure-based (Phold) — more permissive
_FIDENT_THRESH    = 0.30   # minimum fraction identity for structure hits
_QCOV_THRESH      = 0.40   # minimum query coverage for structure hits

# ACR-specific thresholds — more permissive than general structure thresholds.
# Anti-CRISPR proteins are among the most sequence-divergent phage proteins;
# structural homology at low sequence identity is still meaningful.
# Pawluk et al. 2016 (Science 351:aad8405); Stanley & Maxwell 2018 (Cell 175:1481).
_ACR_FIDENT_THRESH = 0.20   # lower identity threshold for ACR detection
_ACR_QCOV_THRESH   = 0.35   # lower coverage threshold for ACR detection


def run(
    cfg:          Config,
    sample_id:    str,
    candidate_id: str,
    force:        bool = False,
) -> Path:
    """Screen a single annotation-ready candidate for resistance and virulence.

    Parses Pharokka and Phold outputs from Module 06.
    No external tools required.

    Parameters
    ----------
    candidate_id : stem of the candidate FASTA (e.g. Ackermannviridae_candidate_001)

    Returns
    -------
    Path to resistance_summary.tsv
    """
    out_dir = cfg.results(STEP) / candidate_id
    rpt_dir = cfg.reports(STEP) / candidate_id
    mkdirs(out_dir, rpt_dir)

    summary_tsv = rpt_dir / f"{candidate_id}_resistance_summary.tsv"

    log_step(f"Module 07 — resistance screening  [{candidate_id}]")

    if summary_tsv.exists() and not force:
        log_info("  Already screened — skipping  (--force to re-run)")
        return summary_tsv

    # Locate annotation outputs
    ann_dir     = cfg.results("06_annotation") / candidate_id
    pharokka_dir = ann_dir / "pharokka"
    phold_dir    = ann_dir / "phold"

    if not ann_dir.exists():
        log_warn(
            f"  Annotation directory not found: {ann_dir}. "
            "Run annotate module first."
        )
        return summary_tsv

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
        task = progress.add_task("Initializing...", total=5)

        # ── Step 1: AMR (CARD) ────────────────────────────────────────────────
        progress.update(task, description="[1/5] AMR screening (CARD — seq + struct)")
        amr_hits = _screen_amr(pharokka_dir, phold_dir, candidate_id)
        progress.advance(task)

        # ── Step 2: Virulence (VFDB) ──────────────────────────────────────────
        progress.update(task, description="[2/5] Virulence factors (VFDB — seq + struct)")
        vf_hits = _screen_virulence(pharokka_dir, phold_dir, candidate_id)
        progress.advance(task)

        # ── Step 3: Anti-CRISPR ───────────────────────────────────────────────
        progress.update(task, description="[3/5] Anti-CRISPR proteins (Phold ACR)")
        acr_hits = _screen_acr(phold_dir)
        progress.advance(task)

        # ── Step 4: Defense systems ───────────────────────────────────────────
        progress.update(task, description="[4/5] Defense systems (DefenseFinder)")
        defense_hits = _screen_defense(phold_dir)
        progress.advance(task)

        # ── Step 5: Toxin-antitoxin ───────────────────────────────────────────
        progress.update(task, description="[5/5] Toxin-antitoxin systems (NetFlax)")
        ta_hits = _screen_netflax(phold_dir)
        progress.advance(task)

    # ── Log findings ──────────────────────────────────────────────────────────
    _log_findings("AMR (CARD)",      amr_hits,     active_warnings, flag=True)
    _log_findings("Virulence (VFDB)", vf_hits,     active_warnings, flag=True)
    _log_findings("Anti-CRISPR",     acr_hits,     active_warnings, flag=False)
    _log_findings("Defense systems", defense_hits, active_warnings, flag=False)
    _log_findings("Toxin-antitoxin", ta_hits,      active_warnings, flag=False)

    # ── Write outputs ─────────────────────────────────────────────────────────
    _write_detail_tsvs(
        out_dir, candidate_id,
        amr_hits, vf_hits, acr_hits, defense_hits, ta_hits,
    )
    _write_summary_tsv(
        summary_tsv, candidate_id, sample_id,
        amr_hits, vf_hits, acr_hits, defense_hits, ta_hits,
    )
    _update_aggregate_tsv(
        cfg.reports(STEP) / "resistance_aggregate.tsv",
        candidate_id, sample_id,
        amr_hits, vf_hits, acr_hits, defense_hits, ta_hits,
    )

    _print_completion_panel(
        candidate_id, out_dir, rpt_dir,
        amr_hits, vf_hits, acr_hits, defense_hits, ta_hits,
        active_warnings,
    )

    log_step(f"Module 07 completed ✓  [{candidate_id}]")
    return summary_tsv


# ── AMR screening ─────────────────────────────────────────────────────────────

def _screen_amr(pharokka_dir: Path, phold_dir: Path, candidate_id: str) -> list[dict]:
    """Parse CARD hits from Pharokka (sequence) and Phold (structure).

    Returns merged list with confidence level per hit.
    """
    hits: dict[str, dict] = {}

    # Pharokka sequence-based CARD hits
    merged = pharokka_dir / f"{candidate_id}_cds_final_merged_output.tsv"
    if merged.exists():
        with open(merged, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                card_hit = (row.get("CARD_hit") or "").strip()
                if not card_hit or card_hit.lower() in ("none", ""):
                    continue
                try:
                    evalue = float(row.get("CARD_eVal") or 1)
                except ValueError:
                    evalue = 1.0
                if evalue > _EVALUE_THRESH:
                    continue
                cds_id = (row.get("gene") or row.get("locus_tag") or "").strip()
                hits[cds_id] = {
                    "cds_id":     cds_id,
                    "hit":        card_hit,
                    "short_name": (row.get("CARD_short_name") or "").strip(),
                    "aro":        (row.get("ARO_Accession") or "").strip(),
                    "gene_family": (row.get("AMR_Gene_Family") or "").strip(),
                    "drug_class": (row.get("Drug_Class") or "").strip(),
                    "mechanism":  (row.get("Resistance_Mechanism") or "").strip(),
                    "evalue":     str(evalue),
                    "pident":     (row.get("CARD_seqIdentity") or "").strip(),
                    "source":     "pharokka_seq",
                    "confidence": "medium",
                    "database":   "CARD",
                }

    # Phold structure-based CARD hits
    phold_card = phold_dir / "sub_db_tophits" / "card_cds_predictions.tsv"
    if phold_card.exists() and phold_card.stat().st_size > 0:
        with open(phold_card, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                cds_id = (row.get("cds_id") or "").strip()
                hit    = (row.get("tophit_protein") or row.get("annotation") or "").strip()
                if not hit:
                    continue
                try:
                    evalue = float(row.get("evalue") or 1)
                    fident = float(row.get("fident") or 0)
                    qcov   = float(row.get("qCov") or 0)
                except ValueError:
                    continue
                if evalue > _STRUCT_EVALUE or fident < _FIDENT_THRESH or qcov < _QCOV_THRESH:
                    continue
                if cds_id in hits:
                    hits[cds_id]["confidence"] = "high"
                    hits[cds_id]["source"]     = "pharokka_seq+phold_struct"
                else:
                    hits[cds_id] = {
                        "cds_id":     cds_id,
                        "hit":        hit,
                        "short_name": hit,
                        "aro":        "",
                        "gene_family": "",
                        "drug_class": "",
                        "mechanism":  "",
                        "evalue":     str(evalue),
                        "pident":     f"{fident:.3f}",
                        "source":     "phold_struct",
                        "confidence": "low",
                        "database":   "CARD",
                    }

    return list(hits.values())


# ── Virulence screening ───────────────────────────────────────────────────────

def _screen_virulence(pharokka_dir: Path, phold_dir: Path, candidate_id: str) -> list[dict]:
    """Parse VFDB hits from Pharokka (sequence) and Phold (structure)."""
    hits: dict[str, dict] = {}

    # Pharokka sequence-based VFDB hits
    merged = pharokka_dir / f"{candidate_id}_cds_final_merged_output.tsv"
    if merged.exists():
        with open(merged, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                vf_hit = (row.get("vfdb_hit") or "").strip()
                if not vf_hit or vf_hit.lower() in ("none", ""):
                    continue
                try:
                    evalue = float(row.get("vfdb_eVal") or 1)
                except ValueError:
                    evalue = 1.0
                if evalue > _EVALUE_THRESH:
                    continue
                cds_id = (row.get("gene") or row.get("locus_tag") or "").strip()
                hits[cds_id] = {
                    "cds_id":      cds_id,
                    "hit":         vf_hit,
                    "short_name":  (row.get("vfdb_short_name") or "").strip(),
                    "description": (row.get("vfdb_description") or "").strip(),
                    "species":     (row.get("vfdb_species") or "").strip(),
                    "evalue":      str(evalue),
                    "pident":      (row.get("vfdb_seqIdentity") or "").strip(),
                    "source":      "pharokka_seq",
                    "confidence":  "medium",
                    "database":    "VFDB",
                }

    # Phold structure-based VFDB hits
    phold_vfdb = phold_dir / "sub_db_tophits" / "vfdb_cds_predictions.tsv"
    if phold_vfdb.exists() and phold_vfdb.stat().st_size > 0:
        with open(phold_vfdb, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                cds_id = (row.get("cds_id") or "").strip()
                hit    = (row.get("tophit_protein") or row.get("annotation") or "").strip()
                if not hit:
                    continue
                try:
                    evalue = float(row.get("evalue") or 1)
                    fident = float(row.get("fident") or 0)
                    qcov   = float(row.get("qCov") or 0)
                except ValueError:
                    continue
                if evalue > _STRUCT_EVALUE or fident < _FIDENT_THRESH or qcov < _QCOV_THRESH:
                    continue
                if cds_id in hits:
                    hits[cds_id]["confidence"] = "high"
                    hits[cds_id]["source"]     = "pharokka_seq+phold_struct"
                else:
                    hits[cds_id] = {
                        "cds_id":      cds_id,
                        "hit":         hit,
                        "short_name":  hit,
                        "description": "",
                        "species":     "",
                        "evalue":      str(evalue),
                        "pident":      f"{fident:.3f}",
                        "source":      "phold_struct",
                        "confidence":  "low",
                        "database":    "VFDB",
                    }

    return list(hits.values())


# ── Anti-CRISPR screening ─────────────────────────────────────────────────────

def _screen_acr(phold_dir: Path) -> list[dict]:
    """Parse anti-CRISPR hits from Phold ACR database.

    Bondy-Denomy et al. 2013 (Nature 493:429).
    Pawluk et al. 2016 (Science 351:aad8405).
    """
    hits = []
    tsv = phold_dir / "sub_db_tophits" / "acr_cds_predictions.tsv"
    if not tsv.exists() or tsv.stat().st_size == 0:
        return hits

    with open(tsv, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            cds_id = (row.get("cds_id") or "").strip()
            acr_id = (row.get("anti_CRISPR_id") or "").strip()
            if not cds_id or not acr_id:
                continue
            try:
                evalue = float(row.get("evalue") or 1)
                fident = float(row.get("fident") or 0)
                qcov   = float(row.get("qCov") or 0)
            except ValueError:
                continue
            if evalue > _STRUCT_EVALUE or fident < _ACR_FIDENT_THRESH or qcov < _ACR_QCOV_THRESH:
                continue
            hits.append({
                "cds_id":     cds_id,
                "acr_id":     acr_id,
                "family":     (row.get("Family") or "").strip(),
                "anti_type":  (row.get("Anti_type") or "").strip(),
                "species":    (row.get("Species") or "").strip(),
                "evalue":     str(evalue),
                "pident":     f"{fident:.3f}",
                "qcov":       f"{qcov:.3f}",
                "confidence": (row.get("annotation_confidence") or "").strip(),
                "pubmed":     (row.get("Pubmed") or "").strip(),
                "database":   "ACR",
            })

    return hits


# ── Defense system screening ──────────────────────────────────────────────────

def _screen_defense(phold_dir: Path) -> list[dict]:
    """Parse defense system hits from Phold DefenseFinder database.

    Abby et al. 2022 (Nucleic Acids Res 50:W245).
    """
    hits = []
    tsv = phold_dir / "sub_db_tophits" / "defensefinder_cds_predictions.tsv"
    if not tsv.exists() or tsv.stat().st_size == 0:
        return hits

    with open(tsv, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            cds_id = (row.get("cds_id") or "").strip()
            if not cds_id:
                continue
            try:
                evalue = float(row.get("evalue") or 1)
                fident = float(row.get("fident") or 0)
                qcov   = float(row.get("qCov") or 0)
            except ValueError:
                continue
            if evalue > _STRUCT_EVALUE or fident < _FIDENT_THRESH or qcov < _QCOV_THRESH:
                continue
            hits.append({
                "cds_id":     cds_id,
                "hit":        (row.get("tophit_protein") or row.get("product") or "").strip(),
                "evalue":     str(evalue),
                "pident":     f"{fident:.3f}",
                "qcov":       f"{qcov:.3f}",
                "confidence": (row.get("annotation_confidence") or "").strip(),
                "database":   "DefenseFinder",
            })

    return hits


# ── Toxin-antitoxin screening ─────────────────────────────────────────────────

def _screen_netflax(phold_dir: Path) -> list[dict]:
    """Parse toxin-antitoxin hits from Phold NetFlax database."""
    hits = []
    tsv = phold_dir / "sub_db_tophits" / "netflax_cds_predictions.tsv"
    if not tsv.exists() or tsv.stat().st_size == 0:
        return hits

    with open(tsv, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            cds_id = (row.get("cds_id") or "").strip()
            if not cds_id:
                continue
            try:
                evalue = float(row.get("evalue") or 1)
                fident = float(row.get("fident") or 0)
                qcov   = float(row.get("qCov") or 0)
            except ValueError:
                continue
            if evalue > _STRUCT_EVALUE or fident < _FIDENT_THRESH or qcov < _QCOV_THRESH:
                continue
            hits.append({
                "cds_id":     cds_id,
                "hit":        (row.get("tophit_protein") or row.get("product") or "").strip(),
                "evalue":     str(evalue),
                "pident":     f"{fident:.3f}",
                "qcov":       f"{qcov:.3f}",
                "confidence": (row.get("annotation_confidence") or "").strip(),
                "database":   "NetFlax",
            })

    return hits


# ── Logging ───────────────────────────────────────────────────────────────────

def _log_findings(
    label:           str,
    hits:            list[dict],
    active_warnings: list[str],
    flag:            bool,
) -> None:
    if not hits:
        log_ok(f"  [{label}] no hits detected")
        return
    n = len(hits)
    if flag:
        msg = f"{n} {label} hit(s) detected — review before therapeutic use"
        log_warn(f"  ⚠ {msg}")
        active_warnings.append(msg)
        for h in hits:
            name = h.get("short_name") or h.get("hit") or h.get("acr_id") or "?"
            conf = h.get("confidence", "")
            log_warn(f"    {h.get('cds_id','')} → {name} [{conf}]")
    else:
        log_ok(f"  [{label}] {n} hit(s) detected")
        for h in hits:
            name = h.get("short_name") or h.get("hit") or h.get("acr_id") or "?"
            conf = h.get("confidence", "")
            log_info(f"    {h.get('cds_id','')} → {name} [{conf}]")


# ── Output writers ────────────────────────────────────────────────────────────

def _write_detail_tsvs(
    out_dir:      Path,
    candidate_id: str,
    amr_hits:     list[dict],
    vf_hits:      list[dict],
    acr_hits:     list[dict],
    defense_hits: list[dict],
    ta_hits:      list[dict],
) -> None:
    """Write per-category detail TSVs."""
    _write_tsv(
        out_dir / f"{candidate_id}_amr.tsv",
        amr_hits,
        ["cds_id", "hit", "short_name", "aro", "gene_family",
         "drug_class", "mechanism", "evalue", "pident", "confidence", "source", "database"],
    )
    _write_tsv(
        out_dir / f"{candidate_id}_virulence.tsv",
        vf_hits,
        ["cds_id", "hit", "short_name", "description", "species",
         "evalue", "pident", "confidence", "source", "database"],
    )
    _write_tsv(
        out_dir / f"{candidate_id}_acr.tsv",
        acr_hits,
        ["cds_id", "acr_id", "family", "anti_type", "species",
         "evalue", "pident", "qcov", "confidence", "pubmed", "database"],
    )
    _write_tsv(
        out_dir / f"{candidate_id}_defense.tsv",
        defense_hits,
        ["cds_id", "hit", "evalue", "pident", "qcov", "confidence", "database"],
    )
    _write_tsv(
        out_dir / f"{candidate_id}_toxin_antitoxin.tsv",
        ta_hits,
        ["cds_id", "hit", "evalue", "pident", "qcov", "confidence", "database"],
    )


def _write_tsv(path: Path, rows: list[dict], headers: list[str]) -> None:
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _write_summary_tsv(
    path:         Path,
    candidate_id: str,
    sample_id:    str,
    amr_hits:     list[dict],
    vf_hits:      list[dict],
    acr_hits:     list[dict],
    defense_hits: list[dict],
    ta_hits:      list[dict],
) -> None:
    n_amr_high = sum(1 for h in amr_hits if h.get("confidence") == "high")
    n_vf_high  = sum(1 for h in vf_hits  if h.get("confidence") == "high")
    headers = [
        "candidate_id", "sample_id",
        "n_amr", "n_amr_high_conf",
        "n_virulence", "n_virulence_high_conf",
        "n_acr", "n_defense", "n_toxin_antitoxin",
        "biosafety_flag",
    ]
    flag = "YES" if (amr_hits or vf_hits) else "NO"
    row = {
        "candidate_id":        candidate_id,
        "sample_id":           sample_id,
        "n_amr":               str(len(amr_hits)),
        "n_amr_high_conf":     str(n_amr_high),
        "n_virulence":         str(len(vf_hits)),
        "n_virulence_high_conf": str(n_vf_high),
        "n_acr":               str(len(acr_hits)),
        "n_defense":           str(len(defense_hits)),
        "n_toxin_antitoxin":   str(len(ta_hits)),
        "biosafety_flag":      flag,
    }
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


def _update_aggregate_tsv(
    path:         Path,
    candidate_id: str,
    sample_id:    str,
    amr_hits:     list[dict],
    vf_hits:      list[dict],
    acr_hits:     list[dict],
    defense_hits: list[dict],
    ta_hits:      list[dict],
) -> None:
    headers = [
        "candidate_id", "sample_id",
        "n_amr", "n_amr_high_conf",
        "n_virulence", "n_virulence_high_conf",
        "n_acr", "n_defense", "n_toxin_antitoxin",
        "biosafety_flag",
    ]
    rows: dict[str, dict] = {}
    if path.exists():
        with open(path) as f:
            old_h = f.readline().rstrip("\n").split("\t")
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if cols and cols[0]:
                    rows[cols[0]] = dict(zip(old_h, cols))

    n_amr_high = sum(1 for h in amr_hits if h.get("confidence") == "high")
    n_vf_high  = sum(1 for h in vf_hits  if h.get("confidence") == "high")
    flag = "YES" if (amr_hits or vf_hits) else "NO"
    rows[candidate_id] = {
        "candidate_id":          candidate_id,
        "sample_id":             sample_id,
        "n_amr":                 str(len(amr_hits)),
        "n_amr_high_conf":       str(n_amr_high),
        "n_virulence":           str(len(vf_hits)),
        "n_virulence_high_conf": str(n_vf_high),
        "n_acr":                 str(len(acr_hits)),
        "n_defense":             str(len(defense_hits)),
        "n_toxin_antitoxin":     str(len(ta_hits)),
        "biosafety_flag":        flag,
    }
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=headers, delimiter="\t")
        writer.writeheader()
        for r in rows.values():
            writer.writerow({h: r.get(h, "") for h in headers})


# ── Completion panel ──────────────────────────────────────────────────────────

def _print_completion_panel(
    candidate_id: str,
    out_dir:      Path,
    rpt_dir:      Path,
    amr_hits:     list[dict],
    vf_hits:      list[dict],
    acr_hits:     list[dict],
    defense_hits: list[dict],
    ta_hits:      list[dict],
    active_warnings: list[str],
) -> None:
    flag = "[bold red]⚠ YES[/bold red]" if (amr_hits or vf_hits) else "[bold green]✓ NO[/bold green]"

    def _hit_lines(hits: list[dict], name_key: str) -> list[str]:
        lines = []
        for h in hits[:5]:  # show max 5
            name = h.get(name_key) or h.get("hit") or h.get("acr_id") or "?"
            conf = h.get("confidence", "")
            lines.append(f"    [dim]{h.get('cds_id','')}[/dim] → {name} [{conf}]")
        if len(hits) > 5:
            lines.append(f"    [dim]... and {len(hits)-5} more[/dim]")
        return lines

    lines = [
        f"  [cyan]Biosafety flag    :[/cyan] {flag}",
        f"  [cyan]AMR genes (CARD)  :[/cyan] {len(amr_hits)}"
        + (f"  [dim](high={sum(1 for h in amr_hits if h.get('confidence')=='high')})[/dim]" if amr_hits else ""),
        *_hit_lines(amr_hits, "short_name"),
        f"  [cyan]Virulence (VFDB)  :[/cyan] {len(vf_hits)}"
        + (f"  [dim](high={sum(1 for h in vf_hits if h.get('confidence')=='high')})[/dim]" if vf_hits else ""),
        *_hit_lines(vf_hits, "short_name"),
        f"  [cyan]Anti-CRISPR (ACR) :[/cyan] {len(acr_hits)}",
        *_hit_lines(acr_hits, "acr_id"),
        f"  [cyan]Defense systems   :[/cyan] {len(defense_hits)}",
        f"  [cyan]Toxin-antitoxin   :[/cyan] {len(ta_hits)}",
        "",
        "[bold]Output files[/bold]",
        f"  AMR      : {out_dir}/{candidate_id}_amr.tsv",
        f"  Virulence: {out_dir}/{candidate_id}_virulence.tsv",
        f"  ACR      : {out_dir}/{candidate_id}_acr.tsv",
        f"  Summary  : {rpt_dir}/{candidate_id}_resistance_summary.tsv",
    ]
    if active_warnings:
        lines += ["", "[bold yellow]Warnings[/bold yellow]"]
        lines += [f"  [yellow]⚠ {w}[/yellow]" for w in active_warnings]

    console.print(Panel(
        "\n".join(lines),
        title=f"[bold cyan]Resistance screening — {candidate_id}[/bold cyan]",
        border_style="cyan", padding=(0, 2), width=120,
    ))

"""PhageFlow Module 04 — Viral identification with geNomad.

Tool:
    geNomad (Camargo et al. 2023, Nature Biotechnology):
        End-to-end classification of contigs as virus / plasmid / chromosome.
        Detects integrated proviruses and assigns viral taxonomy.

        --splits 8         : splits input for parallel processing.
        --min-score 0.7    : virus score threshold (default); configurable via config.yaml.
        --cleanup          : removes large intermediate files after run.

Input  : results/03_assembly/combined/{sample}_contigs_nr.fasta
Output : results/04_viral_id/{sample}_virus.fna  ← CheckV input
         reports/04_viral_id/genomad_summary.tsv

Activation note:
    geNomad requires its own conda environment.
    conda activate genomad
"""

from __future__ import annotations
import shutil
from pathlib import Path

from phageflow.utils.config import Config
from phageflow.utils.logger import log_step, log_info, log_ok, log_warn
from phageflow.utils.tools import require_tools, run_silent, mkdirs, fasta_stats

STEP  = "04_viral_id"
TOOLS = ["genomad"]


def run(
    cfg:       Config,
    sample_id: str,
    contigs:   Path,
    force:     bool = False,
) -> Path:
    """
    Identify viral contigs with geNomad for a single sample.

    Parameters
    ----------
    sample_id : sample identifier (used as output filename prefix)
    contigs   : NR contigs FASTA (output of assembly step)

    Returns
    -------
    virus_fna : Path to viral contigs FASTA (input for CheckV / quality step)
    """
    require_tools(*TOOLS)

    contigs = Path(contigs)
    out_dir = cfg.results(STEP)
    rpt_dir = cfg.reports(STEP)
    mkdirs(out_dir, rpt_dir)

    db = cfg.databases.genomad
    if not db.exists():
        log_warn(f"  geNomad database not found: {db}")

    log_step(f"Module 04 — geNomad viral identification [{sample_id}] · {cfg.threads} threads")

    if not contigs.exists():
        raise FileNotFoundError(f"[{sample_id}] contigs not found: {contigs}")

    s = fasta_stats(contigs)
    log_info(f"  Contigs : {s['n']} sequences | total={s['total_bp']}bp | "
             f"largest={s['largest_bp']}bp")
    log_info(f"  geNomad min-score: {cfg.genomad.min_score} (Camargo et al. 2023)")
    log_info("  Classifying: virus / plasmid / chromosome + provirus detection")

    # ── geNomad end-to-end ────────────────────────────────────────────────────
    sdir   = out_dir / sample_id            # working dir per sample
    prefix = contigs.stem                   # e.g. s1_contigs_nr

    # geNomad output paths (fixed by geNomad naming convention)
    v_fna = sdir / f"{prefix}_summary" / f"{prefix}_virus.fna"
    v_tsv = sdir / f"{prefix}_summary" / f"{prefix}_virus_summary.tsv"

    if v_fna.exists() and v_fna.stat().st_size > 0 and not force:
        log_info("  Already processed — skipping (use --force to re-run)")
    else:
        mkdirs(sdir)
        log_info("  [geNomad] running end-to-end classification...")
        try:
            run_silent([
                "genomad", "end-to-end",
                "--cleanup",
                "--splits",    "8",
                "--min-score", str(cfg.genomad.min_score),
                "--threads",   str(cfg.threads),
                str(contigs), str(sdir), str(db),
            ], log_file=rpt_dir / f"{sample_id}_genomad.log")
            log_ok("  [geNomad] OK")
        except Exception as e:
            log_warn(f"  [geNomad] warning: {e}")

    # ── Parse results ─────────────────────────────────────────────────────────
    if not v_fna.exists():
        log_warn(f"  [geNomad] virus FASTA not found: {v_fna}")
        log_warn("  Check the geNomad log in reports/04_viral_id/")
        virus_fna = out_dir / f"{sample_id}_virus.fna"
        return virus_fna   # return expected path even if empty

    row = _parse_genomad_output(sample_id, contigs, sdir, prefix, v_fna, v_tsv)

    # Copy virus FASTA to flat output dir (CheckV input convention)
    virus_fna = out_dir / f"{sample_id}_virus.fna"
    shutil.copy(v_fna, virus_fna)

    log_ok(
        f"  [{sample_id}] viral={row['viral_ctg']} contigs ({row['viral_bp']}bp) "
        f"| plasmid={row['plasmid']} | provirus={row['provirus']} "
        f"| best_family={row['best_family']}"
    )

    _save_tsv(row, rpt_dir / "genomad_summary.tsv")

    log_step(f"Module 04 completed ✓  [{sample_id}]")
    log_ok(f"  Virus FASTA : {virus_fna}")
    log_info("Next: phageflow quality --sample-id <id> --virus-fna <virus.fna>  "
             "(activate phageflow env)")

    return virus_fna


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def _parse_genomad_output(
    sample_id: str,
    fa:        Path,
    sdir:      Path,
    prefix:    str,
    v_fna:     Path,
    v_tsv:     Path,
) -> dict:
    n_input  = fasta_stats(fa)["n"]
    n_viral  = fasta_stats(v_fna)["n"]     if v_fna.exists()  else 0
    viral_bp = fasta_stats(v_fna)["total_bp"] if v_fna.exists() else 0

    # Plasmid count
    p_tsv  = sdir / f"{prefix}_summary" / f"{prefix}_plasmid_summary.tsv"
    n_plas = _count_tsv_rows(p_tsv)

    # Provirus count
    pr_tsv = sdir / f"{prefix}_find_proviruses" / f"{prefix}_provirus.tsv"
    n_prov = _count_tsv_rows(pr_tsv)

    # Best viral family (highest virus_score)
    best_fam = "unclassified"
    if v_tsv.exists():
        best_score = 0.0
        try:
            with open(v_tsv) as f:
                header = f.readline().strip().split("\t")
                col = {h: i for i, h in enumerate(header)}
                for line in f:
                    row = line.strip().split("\t")
                    if not row:
                        continue
                    try:
                        score = float(row[col.get("virus_score", 6)])
                        if score > best_score:
                            best_score = score
                            tax = row[col.get("taxonomy", 10)] if len(row) > col.get("taxonomy", 10) else ""
                            parts = [p.strip() for p in tax.split(";") if p.strip()]
                            if parts:
                                best_fam = parts[-1]
                    except (ValueError, IndexError):
                        pass
        except Exception:
            pass

    return {
        "sample":      sample_id,
        "input_ctg":   n_input,
        "viral_ctg":   n_viral,
        "viral_bp":    viral_bp,
        "plasmid":     n_plas,
        "provirus":    n_prov,
        "best_family": best_fam,
    }


def _count_tsv_rows(path: Path) -> int:
    """Count non-header, non-comment rows in a TSV."""
    if not path.exists():
        return 0
    try:
        count = sum(
            1 for line in open(path)
            if line.strip() and not line.startswith("#") and not line.startswith("seq_name")
        )
        return max(0, count - 1)   # subtract header
    except Exception:
        return 0


# ---------------------------------------------------------------------------
# Summary TSV (append / update per sample)
# ---------------------------------------------------------------------------

def _save_tsv(row: dict, path: Path) -> None:
    headers = ["sample", "input_ctg", "viral_ctg", "viral_bp",
               "plasmid", "provirus", "best_family"]
    rows = {}
    if path.exists():
        with open(path) as f:
            f.readline()
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if cols:
                    rows[cols[0]] = cols
    rows[row["sample"]] = [str(row.get(h, "")) for h in headers]
    with open(path, "w") as f:
        f.write("\t".join(headers) + "\n")
        for r in rows.values():
            f.write("\t".join(r) + "\n")

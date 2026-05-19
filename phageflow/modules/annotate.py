"""PhageFlow Module 06 — Structural and functional annotation.

Tools and parameters (literature-based):

Pharokka v1.6+ (Bouras et al. 2023, Bioinformatics):
    A fast phage annotation tool using PHANOTATE for gene calling and
    PHROG (Phage Orthologous Groups) for functional annotation.

    --single                 : single contig mode for isolated phage genomes.
    --meta                   : metagenomic mode for multi-contig inputs.
    --dnaapler all           : attempts reorientation to terminase large subunit
                               (standard convention for phage genome deposition;
                               Casjens & Thuman-Commike 2011, Virology).
    --gene_predict_hmm       : PHANOTATE + HMMER hybrid calling for better
                               sensitivity on divergent phages.

    Outputs: {sample}.gbk, {sample}.gff, {sample}_pharokka.tsv,
             {sample}_cds_functions.tsv

Phold v0.2+ (Bouras et al. 2024):
    Protein structure-based annotation using ProstT5 language model
    and Foldseek structural alignment to PHROGs.
    Complements Pharokka by assigning function to hypothetical proteins
    that lack sequence homology.

    --cpu                    : number of threads (structure predictions are
                               memory-intensive — use fewer threads than assembly).
    --skip_extra_annotations : skip tRNA/tmRNA (already done in Pharokka).

    Run after Pharokka — takes {sample}.gbk as input.
    Outputs updated {sample}_phold.gbk with enhanced annotations.

Integration strategy:
    Pharokka → base annotation (sequence homology to PHROGs)
    Phold    → structural function calls for hypothetical proteins
    Combined outputs are used for safety screening and lifecycle prediction.
"""

from __future__ import annotations
from pathlib import Path

from phageflow.utils.config import Config
from phageflow.utils.logger import log_step, log_info, log_ok, log_warn, log_error
from phageflow.utils.tools import require_tools, run_silent, human_size, mkdirs, fasta_stats

STEP  = "06_annotation"
TOOLS = ["pharokka.py", "phold"]


def run(
    cfg:       Config,
    sample_id: str,
    genome:    Path,
    force:     bool = False,
) -> Path:
    """
    Annotate a single phage genome with Pharokka + Phold.

    Parameters
    ----------
    genome : HQ genome FASTA (output of quality step)

    Returns
    -------
    gbk : Path to the final annotated GenBank file
    """
    require_tools(*TOOLS)

    genome  = Path(genome)
    out_dir = cfg.results(STEP)
    rpt_dir = cfg.reports(STEP)
    sdir    = out_dir / sample_id
    mkdirs(sdir, rpt_dir)

    log_step(f"Module 06 — annotation [{sample_id}] · {cfg.threads} threads")

    if not genome.exists():
        log_error(f"  Genome not found: {genome}")
        raise FileNotFoundError(genome)

    s = fasta_stats(genome)
    log_info(f"  Genome : {human_size(genome)} | {s['n']} contig(s) | "
             f"largest={s['largest_bp']}bp")
    log_info("  [1/2] Pharokka: PHANOTATE + PHROG (Bouras et al. 2023)")
    log_info("  [2/2] Phold: structure-based functional annotation (Bouras et al. 2024)")

    pharokka_db = cfg.databases.pharokka
    phold_db    = cfg.databases.phold

    if not pharokka_db.exists():
        log_warn(f"  Pharokka database not found: {pharokka_db}")
    if not phold_db.exists():
        log_warn(f"  Phold database not found: {phold_db}")

    # ── Pharokka ──────────────────────────────────────────────────────────────
    pharokka_dir = sdir / "pharokka"
    gbk_out      = pharokka_dir / f"{sample_id}.gbk"
    mkdirs(pharokka_dir)

    if gbk_out.exists() and gbk_out.stat().st_size > 0 and not force:
        log_info("  [Pharokka] already exists — skipping")
    else:
        n_ctg = s["n"]
        mode_flag = "--meta" if n_ctg > 1 else "--single"
        log_info(f"  [Pharokka] mode: {'meta' if n_ctg > 1 else 'single'} "
                 f"({n_ctg} contig{'s' if n_ctg > 1 else ''})")
        try:
            run_silent([
                "pharokka.py",
                "-i", str(genome),
                "-o", str(pharokka_dir),
                "-d", str(pharokka_db),
                "-p", sample_id,
                mode_flag,
                "--dnaapler", "all",
                "--gene_predict_hmm",
                "--threads", str(cfg.threads),
                "--force",
            ], log_file=rpt_dir / f"{sample_id}_pharokka.log")
            log_ok("  [Pharokka] OK")
        except Exception as e:
            log_warn(f"  [Pharokka] warning: {e}")

    pharokka_stats = _parse_pharokka_tsv(pharokka_dir / f"{sample_id}_cds_functions.tsv")
    if pharokka_stats:
        log_ok(
            f"  [Pharokka] {pharokka_stats['total_cds']} CDS | "
            f"annotated={pharokka_stats['annotated']} "
            f"({pharokka_stats['pct_annotated']}) | "
            f"hypothetical={pharokka_stats['hypothetical']}"
        )

    # ── Phold ─────────────────────────────────────────────────────────────────
    phold_dir = sdir / "phold"
    gbk_phold = phold_dir / f"{sample_id}_phold.gbk"
    mkdirs(phold_dir)

    if gbk_phold.exists() and gbk_phold.stat().st_size > 0 and not force:
        log_info("  [Phold] already exists — skipping")
    elif gbk_out.exists():
        log_info("  [Phold] structural annotation...")
        try:
            run_silent([
                "phold", "run",
                "-i", str(gbk_out),
                "-o", str(phold_dir),
                "-p", sample_id,
                "-d", str(phold_db),
                "--cpu", str(cfg.threads),
                "--skip_extra_annotations",
                "--force",
            ], log_file=rpt_dir / f"{sample_id}_phold.log")
            log_ok("  [Phold] OK")
        except Exception as e:
            log_warn(f"  [Phold] warning: {e}")
    else:
        log_warn("  [Phold] skipped — no Pharokka GBK found")

    phold_stats = _parse_phold_tsv(phold_dir / f"{sample_id}_phold.tsv")
    if phold_stats:
        log_ok(
            f"  [Phold] upgraded={phold_stats['upgraded']} hypothetical proteins "
            f"→ function assigned"
        )

    # ── Summary TSV ───────────────────────────────────────────────────────────
    final_gbk = gbk_phold if gbk_phold.exists() else gbk_out
    _save_tsv(
        sample_id, genome,
        pharokka_stats or {}, phold_stats or {},
        rpt_dir / "annotation_summary.tsv",
    )

    log_step(f"Module 06 completed ✓  [{sample_id}]")
    log_ok(f"  GBK (final) : {final_gbk}")
    log_ok(f"  GFF         : {pharokka_dir / (sample_id + '.gff')}")
    log_info("Next: phageflow safety --sample-id <id> --genome <genome.fasta>")

    return final_gbk


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def _parse_pharokka_tsv(tsv: Path) -> dict:
    """Parse pharokka _cds_functions.tsv to get annotation statistics."""
    if not tsv.exists():
        return {}
    total = annotated = 0
    try:
        with open(tsv) as f:
            header_line = f.readline().strip().split("\t")
            col = {h: i for i, h in enumerate(header_line)}
            for line in f:
                if not line.strip():
                    continue
                row = line.strip().split("\t")
                total += 1
                cat = row[col.get("category", 1)] if len(row) > 1 else ""
                if cat and cat.lower() not in ("unknown function", "hypothetical protein", ""):
                    annotated += 1
    except Exception:
        return {}
    hypo = total - annotated
    pct  = f"{annotated / total * 100:.1f}%" if total else "0.0%"
    return {
        "total_cds": total, "annotated": annotated,
        "hypothetical": hypo, "pct_annotated": pct,
    }


def _parse_phold_tsv(tsv: Path) -> dict:
    """Count how many hypothetical proteins phold upgraded."""
    if not tsv.exists():
        return {}
    upgraded = 0
    try:
        with open(tsv) as f:
            header_line = f.readline().strip().split("\t")
            col = {h: i for i, h in enumerate(header_line)}
            for line in f:
                if not line.strip():
                    continue
                row = line.strip().split("\t")
                phrog = row[col.get("phrog", 0)] if row else ""
                if phrog and phrog not in ("No_PHROG", "", "NA"):
                    upgraded += 1
    except Exception:
        return {}
    return {"upgraded": upgraded}


# ---------------------------------------------------------------------------
# Summary TSV
# ---------------------------------------------------------------------------

def _save_tsv(sample_id, genome, pharokka, phold, path):
    headers = [
        "sample_id", "genome_bp", "total_cds",
        "annotated", "pct_annotated", "hypothetical", "phold_upgraded",
    ]
    rows = {}
    if path.exists():
        with open(path) as f:
            f.readline()
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if cols:
                    rows[cols[0]] = cols

    s = fasta_stats(genome)
    rows[sample_id] = [
        sample_id,
        str(s.get("total_bp", 0)),
        str(pharokka.get("total_cds", "")),
        str(pharokka.get("annotated", "")),
        pharokka.get("pct_annotated", ""),
        str(pharokka.get("hypothetical", "")),
        str(phold.get("upgraded", "")),
    ]

    with open(path, "w") as f:
        f.write("\t".join(headers) + "\n")
        for row in rows.values():
            f.write("\t".join(row) + "\n")

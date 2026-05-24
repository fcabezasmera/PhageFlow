"""PhageFlow Module 06 — Structural and functional annotation.

Three-tier annotation cascade with per-CDS delta tracking from tool TSVs.

Pipeline:
    Pharokka → Phold → Phynteny Transformer → plot (phold or pyGenomeViz)

Output per candidate genome (candidate_id = genome filename stem):
    results/06_annotation/{candidate_id}/
        pharokka/   gene calling + PHROG (MMseqs2 + PyHMMER)
        phold/      structure-based upgrade via ProstT5 + Foldseek
        phynteny/   synteny + ESM2 upgrade (final GBK)
        plots/      genome map — circular (phold plot) or linear (pyGenomeViz)

    reports/06_annotation/
        annotation_summary.tsv          per-candidate aggregate stats (all runs)
        {candidate_id}/
            {candidate_id}_pharokka.log
            {candidate_id}_phold.log
            {candidate_id}_phynteny.log
            {candidate_id}_plot.log
            {candidate_id}_delta_report.tsv

─────────────────────────────────────────────────────────────────────────────
Plot strategy
─────────────────────────────────────────────────────────────────────────────
  Single-contig genome (n_ctg == 1):
      phold plot — circular map (PNG + SVG).
      Reads /product, annotation_confidence, and Foldseek qualifiers.
      Reference: phold.readthedocs.io/en/latest/tutorial/

  Multi-contig genome (n_ctg > 1):
      pyGenomeViz — single linear figure with one track per contig, drawn to
      scale, colored by PHROG category (same palette as phold).
      Rationale: phold plot generates one circular map per GBK record,
      producing n × 2 files with no inter-contig context. A single linear
      figure is more informative for fragmented assemblies.
      Reference: Shimoyama 2022, bioRxiv 2022.11.24.517870

  In both cases the input is the Phold GBK because:
      1. phold plot reads /product, annotation_confidence, and Foldseek
         qualifiers written by Phold — these drive colour coding and labels.
      2. Phynteny does NOT modify /product; its predictions live only in
         /phynteny_category, /phynteny_score, /phynteny_confidence.
      3. Phold explicitly preserves non-CDS features (tRNA/tmRNA/CRISPR)
         from the Pharokka input; Phynteny preservation is undocumented.

─────────────────────────────────────────────────────────────────────────────
Per-CDS delta tracking strategy
─────────────────────────────────────────────────────────────────────────────
Each tier exposes a dedicated per-CDS TSV:

  Tier 1 — Pharokka: {prefix}_cds_final_merged_output.tsv
      Columns of interest: locus_tag (or gene), annot, top_hit, phrog,
                           category, vfdb_hit, card_hit
      Reference: Bouras et al. 2023, Bioinformatics 39:btac776

  Tier 2 — Phold: {prefix}_per_cds_predictions.tsv
      Columns of interest: locus_tag, product, phrog_category,
                           annotation_confidence, foldseek_evalue,
                           foldseek_bitscore, foldseek_lddt, tm_score
      annotation_confidence: high / medium / low — primary quality indicator.
      Reference: Bouras et al. 2025, bioRxiv 2025.08.05.668817

  Tier 3 — Phynteny Transformer: phynteny_per_cds_funcions.tsv
      Columns of interest: locus_tag, phynteny_category,
                           phynteny_score, phynteny_confidence
      Note: 'funcions' is a confirmed typo in the tool output filename.
      Phynteny adds /phynteny_* qualifiers to the GBK but does NOT modify
      /product or /function — the prediction is purely additive.
      Reference: Grigson et al. 2025, bioRxiv 2025.07.28.667340

─────────────────────────────────────────────────────────────────────────────
Non-CDS features (tRNA / tmRNA / CRISPR)
─────────────────────────────────────────────────────────────────────────────
  Phold explicitly preserves non-CDS features from the Pharokka GBK.
  Phynteny preservation is NOT documented and must be validated at runtime.
  The delta report counts non-CDS features from the Phold GBK.
  A validation step warns if Phynteny altered the non-CDS count.
"""

from __future__ import annotations
import csv
import re
import shutil
from pathlib import Path
from typing import Optional

from rich.panel import Panel
from rich.text import Text
from rich.progress import (
    Progress, SpinnerColumn, TextColumn,
    BarColumn, MofNCompleteColumn, TimeElapsedColumn,
)

from phageflow.utils.config import Config
from phageflow.utils.logger import (
    log_step, log_info, log_ok, log_warn, log_error, console,
)
from phageflow.utils.tools import (
    require_tools, run_silent, human_size, mkdirs, fasta_stats,
)

STEP  = "06_annotation"
TOOLS = ["pharokka.py", "phold", "phynteny_transformer"]

# ── Phold parameters ──────────────────────────────────────────────────────────
_PHOLD_SENSITIVITY = "9.5"   # default in phold; explicit for reproducibility
_PHOLD_EVALUE      = "1e-3"  # phold default
_PHOLD_CARD_EVALUE = "1e-3"  # phold default

# ── phold plot parameters (single-contig only) ────────────────────────────────
_PLOT_DPI         = "600"
_PLOT_ANNOTATIONS = "1"
_PLOT_LABEL_SIZE  = "8"
_PLOT_INTERVAL    = "5000"

# ── Annotation category sets ──────────────────────────────────────────────────
_UNKNOWN_CATS = frozenset({
    "unknown function", "hypothetical protein", "", "na", "none",
    "unknown", "uncharacterized protein",
})

# ── Phold confidence tiers ────────────────────────────────────────────────────
_PHOLD_HIGH_CONF   = "high"
_PHOLD_MEDIUM_CONF = "medium"
_PHOLD_LOW_CONF    = "low"

# ── PHROG color palette (consistent with phold plot) ─────────────────────────
# Source: phold colour scheme / PHROG category colours
_PHROG_COLORS: dict[str, str] = {
    "head and packaging":                              "#F0A500",
    "integration and excision":                        "#A020F0",
    "lysis":                                           "#FF4500",
    "moron, auxiliary metabolic and host takeover":    "#228B22",
    "other":                                           "#A9A9A9",
    "tail":                                            "#1E90FF",
    "transcription regulation":                        "#FF69B4",
    "dna, rna and nucleotide metabolism":              "#8B0000",
    "connector":                                       "#20B2AA",
    "unknown function":                                "#D3D3D3",
}


# ── Public entry point ────────────────────────────────────────────────────────

def run(
    cfg:       Config,
    sample_id: str,
    genome:    Path,
    force:     bool = False,
) -> Path:
    """Annotate a single phage candidate genome (three tiers + plot)."""
    require_tools(*TOOLS)

    genome       = Path(genome)
    candidate_id = genome.stem

    out_dir  = cfg.results(STEP)
    rpt_base = cfg.reports(STEP)           # aggregate TSVs live here
    rpt_dir  = rpt_base / candidate_id     # per-candidate logs and delta report
    cdir     = out_dir / candidate_id
    mkdirs(cdir, rpt_base, rpt_dir)

    log_step(f"Module 06 — annotation [{candidate_id}]")

    if not genome.exists():
        log_error(f"  Genome not found: {genome}")
        raise FileNotFoundError(genome)

    s     = fasta_stats(genome)
    n_ctg = s["n"]
    mode  = "meta" if n_ctg > 1 else "single"

    use_gpu    = getattr(cfg.annotate, "phold_gpu",           True)
    batch_size = getattr(cfg.annotate, "phold_batch_size",    1)
    phynteny_c = getattr(cfg.annotate, "phynteny_confidence", 0.7)
    gpu_tag    = "--foldseek_gpu" if use_gpu else "--cpu"

    plot_strategy = (
        "phold plot (circular)"
        if n_ctg == 1
        else f"pyGenomeViz linear ({n_ctg} tracks)"
    )

    log_info(
        f"  Genome   : {human_size(genome)} | {n_ctg} contig(s) | "
        f"largest={s['largest_bp']:,} bp | mode={mode}"
    )
    log_info(f"  [1/4] Pharokka  : PHANOTATE/prodigal-gv + PHROGs MMseqs2 + PyHMMER")
    log_info(f"  [2/4] Phold     : ProstT5 + Foldseek --hyps {gpu_tag} (s={_PHOLD_SENSITIVITY})")
    log_info(f"  [3/4] Phynteny  : Transformer + ESM2 (confidence ≥{phynteny_c})")
    log_info(f"  [4/4] Plot      : {plot_strategy} from Phold GBK")

    pharokka_dir = cdir / "pharokka"
    phold_dir    = cdir / "phold"
    phynteny_dir = cdir / "phynteny"
    plots_dir    = cdir / "plots"
    mkdirs(pharokka_dir, phold_dir, phynteny_dir, plots_dir)

    gbk_pharokka = pharokka_dir / f"{candidate_id}.gbk"
    gbk_phold    = phold_dir    / f"{candidate_id}.gbk"
    gbk_phynteny: Optional[Path] = None

    tier1: dict = {}
    tier2: dict = {}
    tier3: dict = {}

    with Progress(
        SpinnerColumn(spinner_name="dots2"),
        TextColumn("[bold cyan]{task.description:<60}"),
        BarColumn(bar_width=20, complete_style="cyan", finished_style="green"),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        console=console,
        transient=True,
    ) as progress:
        task = progress.add_task("Initializing...", total=4)

        # 1/4 — Pharokka
        progress.update(task, description="[1/4] Pharokka  — gene calling + PHROG")
        _run_pharokka(cfg, candidate_id, genome, pharokka_dir, gbk_pharokka,
                      n_ctg, rpt_dir, force)
        tier1 = _parse_pharokka_per_cds(pharokka_dir, candidate_id)
        if tier1:
            log_ok(
                f"  [Tier 1 Pharokka] {tier1['total_cds']} CDS | "
                f"annotated={tier1['annotated']} ({tier1['pct_annotated']}) | "
                f"hypothetical={tier1['hypothetical']} | "
                f"tRNA={tier1['n_trna']} tmRNA={tier1['n_tmrna']} "
                f"CRISPR={tier1['n_crispr']}"
            )
        progress.advance(task)

        # 2/4 — Phold
        progress.update(task, description="[2/4] Phold     — structure-based upgrade")
        _run_phold(cfg, candidate_id, gbk_pharokka, phold_dir, gbk_phold,
                   rpt_dir, force)
        tier2 = _parse_phold_per_cds(phold_dir, candidate_id, tier1)
        if tier2:
            cum = tier1.get("annotated", 0) + tier2["upgraded"]
            tot = tier1.get("total_cds", 0)
            pct = f"{cum / tot * 100:.1f}%" if tot else "?"
            log_ok(
                f"  [Tier 2 Phold] Δ+{tier2['upgraded']} upgraded | "
                f"high={tier2['high_conf']} medium={tier2['medium_conf']} "
                f"low={tier2['low_conf']} | cumulative {cum}/{tot} ({pct})"
            )
            _validate_noncds_in_phold(gbk_pharokka, gbk_phold, tier1)
        progress.advance(task)

        # 3/4 — Phynteny Transformer
        progress.update(task, description="[3/4] Phynteny  — synteny + ESM2 upgrade")
        _run_phynteny(cfg, candidate_id, gbk_phold, phynteny_dir, rpt_dir, force)
        gbk_phynteny = _find_phynteny_gbk(phynteny_dir, candidate_id)
        tier3 = _parse_phynteny_per_cds(
            phynteny_dir, candidate_id, phynteny_c, tier2,
        )
        if tier3:
            p1, p2, p3 = (
                tier1.get("annotated", 0),
                tier2.get("upgraded", 0),
                tier3["upgraded"],
            )
            cum = p1 + p2 + p3
            tot = tier1.get("total_cds", 0)
            pct = f"{cum / tot * 100:.1f}%" if tot else "?"
            log_ok(
                f"  [Tier 3 Phynteny] Δ+{p3} predicted | "
                f"cumulative {cum}/{tot} ({pct})"
            )
            if gbk_phynteny:
                _validate_noncds_in_phynteny(gbk_phold, gbk_phynteny)
        progress.advance(task)

        # 4/4 — Plot (phold circular OR pyGenomeViz linear)
        progress.update(task, description=f"[4/4] plot      — {plot_strategy}")
        _run_plots(candidate_id, gbk_phold, plots_dir, rpt_dir, n_ctg, force)
        progress.advance(task)

    # ── Build and save per-CDS delta report ───────────────────────────────────
    delta_path = rpt_dir / f"{candidate_id}_delta_report.tsv"
    _save_delta_report(
        candidate_id, sample_id,
        tier1, tier2, tier3,
        phynteny_c,         # pass config value — not the module constant
        delta_path,
    )

    # ── Aggregate summary TSV (shared across all candidates) ──────────────────
    _save_summary_tsv(
        candidate_id, sample_id, genome,
        tier1, tier2, tier3,
        rpt_base / "annotation_summary.tsv",
    )

    _print_completion_panel(
        candidate_id, gbk_phold, gbk_phynteny,
        pharokka_dir, plots_dir,
        tier1, tier2, tier3,
        n_ctg,
    )

    log_step(f"Module 06 completed ✓  [{candidate_id}]")
    log_ok(f"  GBK (Phold)     : {gbk_phold}  ← use for downstream / phold plot")
    log_ok(
        f"  GBK (Phynteny)  : {gbk_phynteny or 'not produced'}"
        f"  ← /phynteny_* qualifiers only"
    )
    log_ok(f"  Delta report    : {delta_path}")
    log_ok(
        f"  GFF             : "
        f"{pharokka_dir / (candidate_id + '.gff')}"
    )
    log_info("Next: phageflow safety --sample-id <id> --genome <genome.fasta>")

    return gbk_phold


# ── Step 1: Pharokka ──────────────────────────────────────────────────────────

def _run_pharokka(cfg, candidate_id, genome, pharokka_dir, gbk_out,
                  n_ctg, rpt_dir, force):
    if gbk_out.exists() and gbk_out.stat().st_size > 0 and not force:
        log_info("  [Pharokka] already exists — skipping  (--force to re-run)")
        return

    pharokka_db = cfg.databases.pharokka
    if not pharokka_db.exists():
        log_warn(f"  [Pharokka] database not found: {pharokka_db}")

    cmd = [
        "pharokka.py",
        "-i",        str(genome),
        "-o",        str(pharokka_dir),
        "-d",        str(pharokka_db),
        "-p",        candidate_id,
        "--threads", str(cfg.threads),
        "--force",
    ]

    if n_ctg > 1:
        # --meta: prodigal-gv gene caller; --meta_hmm for PyHMMER in meta mode
        # --dnaapler is NOT added: it only works with --single (single contig)
        # and would fail or be silently ignored in meta mode
        cmd += ["-m", "--meta_hmm"]
        log_info(
            f"  [Pharokka] mode=meta ({n_ctg} contigs) | "
            f"gene caller=prodigal-gv | MMseqs2 + PyHMMER (--meta_hmm)"
        )
    else:
        # --dnaapler: reorient circular genome to canonical origin
        # Only valid in single-contig mode (Pharokka docs; Bouras et al. 2023)
        cmd.append("--dnaapler")
        log_info("  [Pharokka] mode=single | gene caller=PHANOTATE | MMseqs2 + PyHMMER")

    try:
        run_silent(cmd, log_file=rpt_dir / f"{candidate_id}_pharokka.log")
        log_ok("  [Pharokka] OK")
    except Exception as e:
        log_warn(f"  [Pharokka] warning: {e}")


# ── Step 2: Phold ─────────────────────────────────────────────────────────────

def _run_phold(cfg, candidate_id, input_gbk, phold_dir, gbk_out, rpt_dir, force):
    if gbk_out.exists() and gbk_out.stat().st_size > 0 and not force:
        log_info("  [Phold] already exists — skipping  (--force to re-run)")
        return

    if not input_gbk.exists():
        log_warn("  [Phold] skipped — Pharokka GBK not found")
        return

    phold_db   = cfg.databases.phold
    use_gpu    = getattr(cfg.annotate, "phold_gpu",        True)
    batch_size = getattr(cfg.annotate, "phold_batch_size", 1)
    # phold_finetune default matches AnnotateConfig and config.yaml (True)
    finetune   = getattr(cfg.annotate, "phold_finetune",   True)

    cmd = [
        "phold", "run",
        "-i", str(input_gbk),
        "-o", str(phold_dir),
        "-p", candidate_id,
        "-d", str(phold_db),
        "-t", str(cfg.threads),
        "--hyps",
        "-s", _PHOLD_SENSITIVITY,
        "-e", _PHOLD_EVALUE,
        "--card_vfdb_evalue", _PHOLD_CARD_EVALUE,
        "--batch_size", str(batch_size),
        "--force",
        "--foldseek_gpu" if use_gpu else "--cpu",
    ]

    if finetune:
        cmd.append("--finetune")
        log_info("  [Phold] --finetune enabled (phage-finetuned ProstT5)")

    log_info(
        f"  [Phold] {('--foldseek_gpu' if use_gpu else '--cpu')} | "
        f"--hyps | sensitivity={_PHOLD_SENSITIVITY} | "
        f"card_vfdb_evalue={_PHOLD_CARD_EVALUE}"
    )

    try:
        run_silent(cmd, log_file=rpt_dir / f"{candidate_id}_phold.log")
        log_ok("  [Phold] OK")
    except Exception as e:
        log_warn(f"  [Phold] warning: {e}")


# ── Step 3: Phynteny Transformer ──────────────────────────────────────────────

def _run_phynteny(cfg, candidate_id, input_gbk, phynteny_dir, rpt_dir, force):
    if not input_gbk.exists():
        log_warn("  [Phynteny] skipped — Phold GBK not found")
        return

    existing = list(phynteny_dir.glob("*.gbk"))
    if existing and not force:
        log_info("  [Phynteny] already exists — skipping  (--force to re-run)")
        return

    models_path = cfg.databases.phynteny / "models"
    if not models_path.exists():
        log_warn(f"  [Phynteny] models not found: {models_path}")

    confidence = getattr(cfg.annotate, "phynteny_confidence", 0.7)

    cmd = [
        "phynteny_transformer",
        str(input_gbk),
        "-o",                     str(phynteny_dir),
        "--prefix",               candidate_id,
        "-m",                     str(models_path),
        "--confidence-threshold", str(confidence),
        "-f",   # --force: overwrite existing output directory
    ]

    log_info(
        f"  [Phynteny] input=Phold GBK | confidence ≥{confidence} | "
        f"output: /phynteny_category, /phynteny_score, /phynteny_confidence"
    )
    log_info(
        "  [Phynteny] NOTE: /product and /function are NOT modified by Phynteny"
    )

    try:
        run_silent(cmd, log_file=rpt_dir / f"{candidate_id}_phynteny.log")
        log_ok("  [Phynteny] OK")
    except Exception as e:
        log_warn(f"  [Phynteny] warning: {e}")


def _find_phynteny_gbk(phynteny_dir: Path, candidate_id: str) -> Optional[Path]:
    """Locate Phynteny output GBK — checks known naming patterns before glob.

    phynteny_transformer outputs 'phynteny_transformer.gbk' by default
    when no --prefix is given. With --prefix, the name may vary.
    Candidates ordered from most to least specific.
    """
    candidates = [
        phynteny_dir / "phynteny_transformer.gbk",   # default tool output name
        phynteny_dir / "phynteny.gbk",
        phynteny_dir / f"{candidate_id}.gbk",
        phynteny_dir / f"{candidate_id}_phynteny_transformer.gbk",
    ]
    found = next((p for p in candidates if p.exists()), None)
    if found is None:
        gbks  = list(phynteny_dir.glob("*.gbk"))
        found = gbks[0] if gbks else None
    return found


# ── Step 4: Plot routing ──────────────────────────────────────────────────────

def _run_plots(
    candidate_id: str,
    gbk_phold:   Path,
    plots_dir:   Path,
    rpt_dir:     Path,
    n_ctg:       int,
    force:       bool,
) -> None:
    """
    Route to the appropriate genome visualisation tool.

    Single-contig (n_ctg == 1):
        phold plot — circular map (PNG + SVG).
        phold reads /product and Foldseek qualifiers that Phynteny does not
        expose, so the Phold GBK is always the correct input.

    Multi-contig (n_ctg > 1):
        pyGenomeViz — single linear figure with one track per contig drawn
        to scale, colored by PHROG category (same palette as phold plot).
        Rationale: phold --all on a multi-record GBK produces n × 2 files
        (PNG + SVG per contig) with no inter-contig context. One figure
        with n tracks is more informative for fragmented assemblies.
        Reference: Shimoyama 2022, bioRxiv 2022.11.24.517870
    """
    existing = list(plots_dir.glob("*.png")) + list(plots_dir.glob("*.svg"))
    if existing and not force:
        log_info(f"  [plot] already exists ({len(existing)} file(s)) — skipping")
        return

    if not gbk_phold.exists():
        log_warn("  [plot] skipped — Phold GBK not found")
        return

    if n_ctg == 1:
        _run_phold_plot_single(candidate_id, gbk_phold, plots_dir, rpt_dir)
    else:
        _run_pygenomeviz_multi(candidate_id, gbk_phold, plots_dir, rpt_dir, n_ctg)


def _run_phold_plot_single(
    candidate_id: str,
    gbk_phold:   Path,
    plots_dir:   Path,
    rpt_dir:     Path,
) -> None:
    """
    phold plot for single-contig complete/circular genomes.

    --all omitted: generates PNG + SVG (2 files) without the additional
    per-feature annotation files that --all adds — sufficient for a
    single-contig genome and avoids redundant output.
    """
    log_info(
        f"  [phold plot] single-contig circular map | "
        f"input: {gbk_phold.name}"
    )
    cmd = [
        "phold", "plot",
        "-i",            str(gbk_phold),
        "-o",            str(plots_dir),
        "-p",            candidate_id,
        "--dpi",         _PLOT_DPI,
        "--annotations", _PLOT_ANNOTATIONS,
        "--label_size",  _PLOT_LABEL_SIZE,
        "--interval",    _PLOT_INTERVAL,
        "-f",   # --force: overwrite existing output
    ]
    try:
        run_silent(cmd, log_file=rpt_dir / f"{candidate_id}_plot.log")
        pngs = list(plots_dir.glob("*.png"))
        svgs = list(plots_dir.glob("*.svg"))
        log_ok(
            f"  [phold plot] {len(pngs)} PNG + {len(svgs)} SVG → {plots_dir}"
        )
    except Exception as e:
        log_warn(f"  [phold plot] warning: {e}")


def _run_pygenomeviz_multi(
    candidate_id: str,
    gbk_phold:   Path,
    plots_dir:   Path,
    rpt_dir:     Path,
    n_ctg:       int,
) -> None:
    """
    pyGenomeViz linear figure for multi-contig (fragmented) genomes.

    One horizontal track per contig, drawn to scale. CDS colored by PHROG
    category (/product) using the same palette as phold plot. tRNA and tmRNA
    marked with a distinct color. A shared legend is added at the bottom.

    Reference: Shimoyama 2022, bioRxiv 2022.11.24.517870
    """
    log_info(
        f"  [pyGenomeViz] multi-contig linear map | {n_ctg} tracks | "
        f"input: {gbk_phold.name}"
    )
    out_png = plots_dir / f"{candidate_id}_linear.png"
    out_svg = plots_dir / f"{candidate_id}_linear.svg"

    try:
        from Bio import SeqIO
        from pygenomeviz import GenomeViz
        import matplotlib.patches as mpatches
        import matplotlib.pyplot as plt

        records = list(SeqIO.parse(str(gbk_phold), "genbank"))
        if not records:
            log_warn("  [pyGenomeViz] no records found in Phold GBK")
            return

        gv = GenomeViz(
            fig_width           = 16,
            fig_track_height    = 0.5,
            track_align_type    = "left",
            feature_track_ratio = 0.3,
        )

        for rec in records:
            track = gv.add_feature_track(
                name      = rec.id,
                segments  = len(rec.seq),
                labelsize = 8,
            )
            for feat in rec.features:
                if feat.type == "CDS":
                    product = feat.qualifiers.get(
                        "product", ["unknown function"]
                    )[0].lower()
                    color = next(
                        (v for k, v in _PHROG_COLORS.items() if k in product),
                        _PHROG_COLORS["unknown function"],
                    )
                    track.add_features(feat, color=color, linewidth=0.3)

                elif feat.type in ("tRNA", "tmRNA"):
                    track.add_features(
                        feat,
                        color     = "#000000",
                        linewidth = 0.8,
                    )

        # Build figure and add shared legend
        fig = gv.plotfig()
        legend_handles = [
            mpatches.Patch(color=color, label=cat.title())
            for cat, color in _PHROG_COLORS.items()
        ] + [
            mpatches.Patch(color="#000000", label="tRNA / tmRNA"),
        ]
        fig.legend(
            handles        = legend_handles,
            loc            = "lower center",
            ncol           = 3,
            fontsize       = 7,
            frameon        = False,
            bbox_to_anchor = (0.5, -0.02),
        )
        fig.savefig(str(out_png), dpi=int(_PLOT_DPI), bbox_inches="tight")
        fig.savefig(str(out_svg),                     bbox_inches="tight")
        plt.close(fig)
        log_ok(f"  [pyGenomeViz] 2 file(s) → {plots_dir}  (PNG + SVG)")

    except Exception as e:
        log_warn(
            f"  [pyGenomeViz] warning: {e} — "
            f"falling back to phold plot (individual contig maps)"
        )
        _run_phold_plot_single(candidate_id, gbk_phold, plots_dir, rpt_dir)


# ── Per-CDS TSV parsers ───────────────────────────────────────────────────────

def _parse_pharokka_per_cds(pharokka_dir: Path, candidate_id: str) -> dict:
    """
    Parse Pharokka per-CDS data from:
      {prefix}_cds_final_merged_output.tsv  — per-CDS MMseqs2/PyHMMER hits
      {prefix}_cds_functions.tsv            — per-category/feature counts

    The locus_tag column may be named 'locus_tag', 'gene', or 'ID'
    depending on the Pharokka version. All variants are tried in order.

    Reference: Bouras et al. 2023, Bioinformatics 39:btac776
    """
    merged_tsv  = pharokka_dir / f"{candidate_id}_cds_final_merged_output.tsv"
    summary_tsv = pharokka_dir / f"{candidate_id}_cds_functions.tsv"

    # Candidate column names for the locus identifier (version-dependent)
    _LOCUS_COLS = ("locus_tag", "gene", "ID", "id", "seq_id")

    cds_rows: list[dict] = []
    locus_col: str = "locus_tag"   # default; will be resolved below

    if merged_tsv.exists():
        with open(merged_tsv, newline="") as f:
            reader  = csv.DictReader(f, delimiter="\t")
            headers = reader.fieldnames or []
            # Resolve actual locus column name from file headers
            locus_col = next(
                (c for c in _LOCUS_COLS if c in headers),
                _LOCUS_COLS[0],
            )
            cds_rows = list(reader)

    # Normalise all rows so downstream code can always use "locus_tag"
    if locus_col != "locus_tag":
        for row in cds_rows:
            row["locus_tag"] = row.get(locus_col, "")

    total     = len(cds_rows)
    annotated = sum(
        1 for r in cds_rows
        if r.get("annot", "").strip().lower() not in _UNKNOWN_CATS
        and r.get("annot", "").strip() != ""
    )
    hypo = total - annotated
    pct  = f"{annotated / total * 100:.1f}%" if total else "0.0%"

    n_trna = n_tmrna = n_crispr = 0
    if summary_tsv.exists():
        with open(summary_tsv, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                try:
                    n_trna   += int(row.get("tRNAs",   0) or 0)
                    n_tmrna  += int(row.get("tmRNAs",  0) or 0)
                    n_crispr += int(row.get("CRISPRs", 0) or 0)
                except (ValueError, TypeError):
                    pass

    return {
        "total_cds":     total,
        "annotated":     annotated,
        "hypothetical":  hypo,
        "pct_annotated": pct,
        "n_trna":        n_trna,
        "n_tmrna":       n_tmrna,
        "n_crispr":      n_crispr,
        "cds_rows":      cds_rows,
    }


def _parse_phold_per_cds(
    phold_dir:   Path,
    candidate_id: str,
    tier1:       dict,
) -> dict:
    """
    Parse Phold per-CDS data from {prefix}_per_cds_predictions.tsv.

    annotation_confidence (high / medium / low) is the primary quality
    indicator per phold documentation.

    Delta logic (--hyps mode):
        Only CDS that were hypothetical in Pharokka (tier1) can be upgraded
        by Phold. Any such CDS that receives a non-hypothetical /product
        from Phold counts as a genuine upgrade.

    Reference: Bouras et al. 2025, bioRxiv 2025.08.05.668817
    """
    tsv = _find_phold_per_cds_tsv(phold_dir, candidate_id)
    if not tsv:
        return {}

    pharokka_hypo: set[str] = set()
    for r in tier1.get("cds_rows", []):
        annot = r.get("annot", "").strip().lower()
        if annot in _UNKNOWN_CATS or annot == "":
            pharokka_hypo.add(r.get("locus_tag", "").strip())

    phold_rows: list[dict] = []
    upgraded = high_conf = medium_conf = low_conf = 0

    with open(tsv, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            phold_rows.append(row)
            locus   = row.get("locus_tag", "").strip()
            product = row.get("product",   "").strip().lower()
            conf    = row.get("annotation_confidence", "").strip().lower()

            if locus in pharokka_hypo and product not in _UNKNOWN_CATS and product:
                upgraded += 1
                if conf == _PHOLD_HIGH_CONF:
                    high_conf += 1
                elif conf == _PHOLD_MEDIUM_CONF:
                    medium_conf += 1
                else:
                    low_conf += 1

    return {
        "upgraded":    upgraded,
        "high_conf":   high_conf,
        "medium_conf": medium_conf,
        "low_conf":    low_conf,
        "phold_rows":  phold_rows,
    }


def _parse_phynteny_per_cds(
    phynteny_dir:  Path,
    candidate_id:  str,
    confidence:    float,
    tier2:         dict,
) -> dict:
    """
    Parse Phynteny Transformer per-CDS data from phynteny_per_cds_funcions.tsv.

    NOTE: 'funcions' is a confirmed typo in the tool's output filename
    (verified against the susiegriggo/Phynteny_transformer repository README).

    Phynteny adds /phynteny_category, /phynteny_score, /phynteny_confidence
    qualifiers to the GBK but does NOT modify /product or /function.

    Delta logic:
        Count CDS that (a) were still hypothetical after Phold, AND
        (b) received a Phynteny prediction with confidence ≥ threshold.

    Reference: Grigson et al. 2025, bioRxiv 2025.07.28.667340
    """
    tsv_candidates = [
        phynteny_dir / "phynteny_per_cds_funcions.tsv",           # confirmed tool typo
        phynteny_dir / f"{candidate_id}_phynteny_per_cds_funcions.tsv",
        phynteny_dir / "phynteny_per_cds_functions.tsv",           # corrected spelling fallback
        phynteny_dir / f"{candidate_id}_phynteny_per_cds_functions.tsv",
    ]
    tsv = next((p for p in tsv_candidates if p.exists()), None)

    if not tsv:
        log_info(
            "  [Phynteny] per-CDS TSV not found — "
            "falling back to GBK qualifier parsing for delta count"
        )
        return _parse_phynteny_gbk_qualifiers_fallback(
            phynteny_dir, candidate_id, confidence, tier2,
        )

    phold_still_hypo: set[str] = set()
    for r in tier2.get("phold_rows", []):
        product = r.get("product", "").strip().lower()
        if product in _UNKNOWN_CATS or product == "":
            phold_still_hypo.add(r.get("locus_tag", "").strip())

    phynteny_rows: list[dict] = []
    upgraded = 0

    with open(tsv, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            phynteny_rows.append(row)
            locus = row.get("locus_tag", "").strip()
            cat   = row.get("phynteny_category", "").strip().lower()
            try:
                conf_val = float(row.get("phynteny_confidence", 0) or 0)
            except (ValueError, TypeError):
                conf_val = 0.0

            if (
                locus in phold_still_hypo
                and cat
                and cat not in _UNKNOWN_CATS
                and conf_val >= confidence
            ):
                upgraded += 1

    return {
        "upgraded":      upgraded,
        "phynteny_rows": phynteny_rows,
        "source":        "tsv",
    }


def _parse_phynteny_gbk_qualifiers_fallback(
    phynteny_dir:  Path,
    candidate_id:  str,
    confidence:    float,
    tier2:         dict,
) -> dict:
    """
    Fallback: parse /phynteny_* qualifiers directly from the Phynteny GBK
    when the per-CDS TSV is absent.

    Preferred only as a fallback — the TSV is more structured and less
    fragile than regex-based GBK parsing.
    """
    gbk_phynteny = _find_phynteny_gbk(phynteny_dir, candidate_id)
    if not gbk_phynteny:
        return {}

    phold_still_hypo: set[str] = set()
    for r in tier2.get("phold_rows", []):
        product = r.get("product", "").strip().lower()
        if product in _UNKNOWN_CATS or product == "":
            phold_still_hypo.add(r.get("locus_tag", "").strip())

    try:
        text = gbk_phynteny.read_text(errors="replace")
    except Exception:
        return {}

    upgraded = 0
    for block in text.split("     CDS ")[1:]:
        def _q(name: str) -> str:
            m = re.search(rf'/{name}="((?:[^"\\]|\\.|\n\s+)*)"', block)
            if not m:
                return ""
            return re.sub(r"\n\s+", " ", m.group(1)).strip()

        locus = _q("locus_tag")
        cat   = _q("phynteny_category").lower()
        try:
            conf_val = float(_q("phynteny_confidence") or 0)
        except (ValueError, TypeError):
            conf_val = 0.0

        if (
            locus in phold_still_hypo
            and cat
            and cat not in _UNKNOWN_CATS
            and conf_val >= confidence
        ):
            upgraded += 1

    return {
        "upgraded":      upgraded,
        "phynteny_rows": [],
        "source":        "gbk_fallback",
    }


# ── Non-CDS feature validation ────────────────────────────────────────────────

def _count_noncds_gbk(gbk_path: Path) -> dict[str, int]:
    """Count non-CDS feature lines in a GBK file by feature type keyword."""
    if not gbk_path or not gbk_path.exists():
        return {}
    counts: dict[str, int] = {}
    text = gbk_path.read_text(errors="replace")
    for feat_type in ("tRNA", "tmRNA", "repeat_region", "misc_feature"):
        counts[feat_type] = len(
            re.findall(rf"^\s+{feat_type}\s", text, re.MULTILINE)
        )
    return counts


def _validate_noncds_in_phold(
    gbk_pharokka: Path, gbk_phold: Path, tier1: dict,
) -> None:
    """Verify Phold preserved non-CDS features from the Pharokka GBK."""
    pk = _count_noncds_gbk(gbk_pharokka)
    ph = _count_noncds_gbk(gbk_phold)
    for feat, n_pk in pk.items():
        n_ph = ph.get(feat, 0)
        if n_pk > 0 and n_ph == 0:
            log_warn(
                f"  [non-CDS] {feat}: in Pharokka ({n_pk}) "
                f"but missing in Phold GBK — check Phold output"
            )
        elif n_pk != n_ph:
            log_warn(
                f"  [non-CDS] {feat}: Pharokka={n_pk} vs Phold={n_ph}"
            )


def _validate_noncds_in_phynteny(
    gbk_phold: Path, gbk_phynteny: Path,
) -> None:
    """
    Check whether Phynteny preserved non-CDS features from the Phold GBK.
    Phynteny does NOT document this behaviour — validation is mandatory.
    If features are lost, the Phold GBK remains the authoritative source.
    """
    ph = _count_noncds_gbk(gbk_phold)
    py = _count_noncds_gbk(gbk_phynteny)
    lost: list[str] = []
    for feat, n_ph in ph.items():
        n_py = py.get(feat, 0)
        if n_ph > 0 and n_py == 0:
            lost.append(f"{feat}({n_ph}→0)")
        elif n_ph != n_py:
            lost.append(f"{feat}({n_ph}→{n_py})")
    if lost:
        log_warn(
            f"  [non-CDS] Phynteny GBK lost features: {', '.join(lost)}. "
            f"Phold GBK is authoritative for plotting and downstream use."
        )
    else:
        log_ok("  [non-CDS] Phynteny preserved all non-CDS features from Phold GBK")


# ── Per-CDS delta report ──────────────────────────────────────────────────────

_DELTA_HEADERS = [
    "candidate_id",
    "locus_tag",
    # Tier 1 — Pharokka
    "t1_annot",
    "t1_phrog",
    "t1_category",
    "t1_vfdb_hit",
    "t1_card_hit",
    # Tier 2 — Phold
    "t2_product",
    "t2_phrog_category",
    "t2_annotation_confidence",
    "t2_foldseek_evalue",
    "t2_foldseek_lddt",
    "t2_tm_score",
    # Tier 3 — Phynteny
    "t3_phynteny_category",
    "t3_phynteny_confidence",
    "t3_phynteny_score",
    # Summary
    "improved_at_tier",   # "1", "2", "3", or "" (still hypothetical)
]


def _save_delta_report(
    candidate_id: str,
    sample_id:    str,
    tier1:        dict,
    tier2:        dict,
    tier3:        dict,
    phynteny_c:   float,    # from config — NOT the module constant
    path:         Path,
) -> None:
    """
    Write a per-CDS delta report joining all three tier TSVs on locus_tag.

    improved_at_tier:
        "1" — Pharokka gave a non-hypothetical annotation
        "2" — Phold upgraded a Pharokka hypothetical
        "3" — Phynteny predicted a category with confidence ≥ phynteny_c
        ""  — still hypothetical / unknown after all three tiers
    """
    p1_map: dict[str, dict] = {
        r.get("locus_tag", "").strip(): r
        for r in tier1.get("cds_rows", [])
    }
    p2_map: dict[str, dict] = {
        r.get("locus_tag", "").strip(): r
        for r in tier2.get("phold_rows", [])
    }
    p3_map: dict[str, dict] = {
        r.get("locus_tag", "").strip(): r
        for r in tier3.get("phynteny_rows", [])
    }

    all_loci = sorted(set(p1_map) | set(p2_map) | set(p3_map))

    rows: list[dict] = []
    for locus in all_loci:
        r1 = p1_map.get(locus, {})
        r2 = p2_map.get(locus, {})
        r3 = p3_map.get(locus, {})

        t1_annot = r1.get("annot", "").strip()
        t2_prod  = r2.get("product", "").strip()
        t3_cat   = r3.get("phynteny_category", "").strip()

        t1_annotated = t1_annot.lower() not in _UNKNOWN_CATS and t1_annot != ""
        t2_annotated = t2_prod.lower()  not in _UNKNOWN_CATS and t2_prod  != ""
        try:
            t3_conf = float(r3.get("phynteny_confidence", 0) or 0)
        except (ValueError, TypeError):
            t3_conf = 0.0
        # Use phynteny_c from config — not a hardcoded module constant
        t3_predicted = (
            t3_cat.lower() not in _UNKNOWN_CATS
            and t3_cat != ""
            and t3_conf >= phynteny_c
        )

        if t1_annotated:
            improved = "1"
        elif t2_annotated:
            improved = "2"
        elif t3_predicted:
            improved = "3"
        else:
            improved = ""

        rows.append({
            "candidate_id":             candidate_id,
            "locus_tag":                locus,
            "t1_annot":                 t1_annot,
            "t1_phrog":                 r1.get("phrog",    ""),
            "t1_category":              r1.get("category", ""),
            "t1_vfdb_hit":              r1.get("vfdb_hit", ""),
            "t1_card_hit":              r1.get("card_hit", ""),
            "t2_product":               t2_prod,
            "t2_phrog_category":        r2.get("phrog_category",       ""),
            "t2_annotation_confidence": r2.get("annotation_confidence",""),
            "t2_foldseek_evalue":       r2.get("foldseek_evalue",      ""),
            "t2_foldseek_lddt":         r2.get("foldseek_lddt",        ""),
            "t2_tm_score":              r2.get("tm_score",             ""),
            "t3_phynteny_category":     t3_cat,
            "t3_phynteny_confidence":   r3.get("phynteny_confidence",  ""),
            "t3_phynteny_score":        r3.get("phynteny_score",       ""),
            "improved_at_tier":         improved,
        })

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=_DELTA_HEADERS, delimiter="\t", extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)

    n1 = sum(1 for r in rows if r["improved_at_tier"] == "1")
    n2 = sum(1 for r in rows if r["improved_at_tier"] == "2")
    n3 = sum(1 for r in rows if r["improved_at_tier"] == "3")
    n0 = sum(1 for r in rows if r["improved_at_tier"] == "")
    log_ok(
        f"  [delta report] {len(rows)} CDS | "
        f"tier1={n1} tier2={n2} tier3={n3} still_unknown={n0} "
        f"→ {path.name}"
    )


# ── Helper: locate phold per-CDS TSV ─────────────────────────────────────────

def _find_phold_per_cds_tsv(
    phold_dir: Path, candidate_id: str,
) -> Optional[Path]:
    """Locate the phold per-CDS predictions TSV (named by phold prefix)."""
    candidates = [
        phold_dir / f"{candidate_id}_per_cds_predictions.tsv",
        phold_dir / "phold_per_cds_predictions.tsv",
    ]
    found = next((p for p in candidates if p.exists()), None)
    if found is None:
        tsvs  = list(phold_dir.glob("*per_cds_predictions*.tsv"))
        found = tsvs[0] if tsvs else None
    if found is None:
        log_warn(
            "  [Phold] per_cds_predictions.tsv not found — "
            "Tier 2 delta will be 0. Check Phold output directory."
        )
    return found


# ── Rich display helpers ──────────────────────────────────────────────────────

def _print_completion_panel(
    candidate_id:  str,
    gbk_phold:     Path,
    gbk_phynteny:  Optional[Path],
    pharokka_dir:  Path,
    plots_dir:     Path,
    tier1:         dict,
    tier2:         dict,
    tier3:         dict,
    n_ctg:         int,
) -> None:
    total = tier1.get("total_cds",  0)
    p1    = tier1.get("annotated",  0)
    p2    = tier2.get("upgraded",   0)
    p3    = tier3.get("upgraded",   0)
    final = p1 + p2 + p3
    hypo  = max(0, total - final)
    plots = list(plots_dir.glob("*.png")) + list(plots_dir.glob("*.svg"))

    plot_mode = (
        "circular (phold plot)"
        if n_ctg == 1
        else f"linear multi-track (pyGenomeViz, {n_ctg} contigs)"
    )

    text = Text()
    text.append("Annotation improvement (per-CDS TSVs):\n", style="bold")
    text.append(
        f"  Tier 1 Pharokka  : {p1:>4}/{total} "
        f"({tier1.get('pct_annotated','?')}) | "
        f"tRNA={tier1.get('n_trna',0)} "
        f"tmRNA={tier1.get('n_tmrna',0)} "
        f"CRISPR={tier1.get('n_crispr',0)}\n",
        style="cyan",
    )
    if total:
        pct2 = f"{(p1 + p2) / total * 100:.1f}%"
        text.append(
            f"  Tier 2 Phold     : +{p2:>3}  → {p1 + p2:>4}/{total} ({pct2}) | "
            f"high={tier2.get('high_conf',0)} "
            f"medium={tier2.get('medium_conf',0)} "
            f"low={tier2.get('low_conf',0)}\n",
            style="cyan",
        )
        pct3 = f"{final / total * 100:.1f}%"
        text.append(
            f"  Tier 3 Phynteny  : +{p3:>3}  → {final:>4}/{total} ({pct3})\n",
            style="bold green",
        )
    text.append(
        f"  Still unknown    : {hypo:>4}/{total}\n\n",
        style="dim white",
    )
    text.append("GBK (Phold)    : ", style="dim white")
    text.append(str(gbk_phold) + "  ← canonical output\n", style="white")
    text.append("GBK (Phynteny) : ", style="dim white")
    text.append(
        (str(gbk_phynteny) if gbk_phynteny else "not produced")
        + "  ← /phynteny_* only\n",
        style="white",
    )
    text.append("GFF            : ", style="dim white")
    text.append(
        str(pharokka_dir / f"{candidate_id}.gff") + "\n",
        style="white",
    )
    text.append("Plots          : ", style="dim white")
    text.append(
        f"{plots_dir}  ({len(plots)} file(s), {plot_mode})\n",
        style="white",
    )

    console.print(Panel(
        text,
        title=f"[bold cyan]Annotation complete — {candidate_id}[/bold cyan]",
        border_style="cyan",
        padding=(0, 2),
        width=120,
    ))


# ── Aggregate summary TSV ─────────────────────────────────────────────────────

_SUMMARY_HEADERS = [
    "candidate_id", "sample_id", "genome_bp",
    "total_cds",
    "n_trna", "n_tmrna", "n_crispr",
    "tier1_annotated", "tier1_pct",
    "tier2_upgraded", "tier2_high", "tier2_medium", "tier2_low",
    "tier3_upgraded",
    "final_annotated", "final_pct", "final_hypothetical",
    "tier2_delta_source",
    "tier3_delta_source",
]


def _save_summary_tsv(
    candidate_id: str,
    sample_id:    str,
    genome:       Path,
    tier1:        dict,
    tier2:        dict,
    tier3:        dict,
    path:         Path,
) -> None:
    """Append or update the aggregate annotation_summary.tsv."""
    rows: dict[str, dict] = {}
    if path.exists():
        with open(path) as f:
            old_hdrs = f.readline().rstrip("\n").split("\t")
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if cols and cols[0]:
                    rows[cols[0]] = dict(zip(old_hdrs, cols))

    s     = fasta_stats(genome)
    total = tier1.get("total_cds",  0)
    p1    = tier1.get("annotated",  0)
    p2    = tier2.get("upgraded",   0)
    p3    = tier3.get("upgraded",   0)
    final = p1 + p2 + p3
    hypo  = max(0, total - final)

    rows[candidate_id] = {
        "candidate_id":       candidate_id,
        "sample_id":          sample_id,
        "genome_bp":          str(s.get("total_bp", 0)),
        "total_cds":          str(total),
        "n_trna":             str(tier1.get("n_trna",    0)),
        "n_tmrna":            str(tier1.get("n_tmrna",   0)),
        "n_crispr":           str(tier1.get("n_crispr",  0)),
        "tier1_annotated":    str(p1),
        "tier1_pct":          tier1.get("pct_annotated", ""),
        "tier2_upgraded":     str(p2),
        "tier2_high":         str(tier2.get("high_conf",   0)),
        "tier2_medium":       str(tier2.get("medium_conf", 0)),
        "tier2_low":          str(tier2.get("low_conf",    0)),
        "tier3_upgraded":     str(p3),
        "final_annotated":    str(final),
        "final_pct":          f"{final / total * 100:.1f}%" if total else "",
        "final_hypothetical": str(hypo),
        "tier2_delta_source": "tsv",
        "tier3_delta_source": tier3.get("source", ""),
    }

    with open(path, "w") as f:
        f.write("\t".join(_SUMMARY_HEADERS) + "\n")
        for row in rows.values():
            f.write(
                "\t".join(row.get(h, "") for h in _SUMMARY_HEADERS) + "\n"
            )

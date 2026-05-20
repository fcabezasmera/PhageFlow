"""PhageFlow Module 06 — Structural and functional annotation.

Three-tier annotation cascade with per-step delta tracking.

Pipeline:
    Pharokka → Phold → Phynteny Transformer → phold plot (from Phynteny GBK)

Output per candidate genome (candidate_id = genome filename stem):
    results/06_annotation/{candidate_id}/
        pharokka/   gene calling + PHROG (MMseqs2 + PyHMMER)
        phold/      structure-based upgrade via ProstT5 + Foldseek
        phynteny/   synteny + ESM2 upgrade (final GBK)
        plots/      circular map from phynteny GBK (600 dpi)

Tier 1 — Pharokka v1.9+ (Bouras et al. 2023, Bioinformatics 39:btac776):
    Single-contig : PHANOTATE + MMseqs2 + PyHMMER + --dnaapler
    Multi-contig  : prodigal-gv + -m + --meta_hmm
    --dnaapler is a boolean flag in v1.4+ — no argument.

Tier 2 — Phold v1.2+ (Bouras et al. 2024, Bioinformatics):
    Input : Pharokka GBK
    --hyps      : targets hypothetical proteins only (Bouras 2024 Supp)
    --finetune  : phage-finetuned ProstT5 (better phage fold coverage)
    -s 9.5      : maximum Foldseek sensitivity (Steinegger 2017)
    --card_vfdb_evalue 1e-10 : strict AMR/VF threshold
    Output: {candidate_id}.gbk in phold dir

Tier 3 — Phynteny Transformer (Grigson et al. 2025, bioRxiv):
    Input : Phold GBK  ← confirmed as intended workflow in benchmarking
            (Zenodo v2: "Phold GBK used as input for Phynteny")
    Integrates positional encoding + bidirectional LSTM + circular attention
    transformer. Trained on 280 000+ phage genomes. AUC > 0.84 across
    9 PHROG categories. Average +14% annotation improvement.
    Output: phynteny_transformer.gbk (final annotated GBK)

Plot — phold plot (from Phynteny GBK = final annotations):
    --all : one circular map per contig (essential for multi-contig bins)
    --dpi 600 : publication-quality resolution
"""

from __future__ import annotations
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
from phageflow.utils.logger import log_step, log_info, log_ok, log_warn, log_error, console
from phageflow.utils.tools import require_tools, run_silent, human_size, mkdirs, fasta_stats

STEP  = "06_annotation"
TOOLS = ["pharokka.py", "phold", "phynteny_transformer"]

# ── Phold parameters ──────────────────────────────────────────────────────────
_PHOLD_SENSITIVITY = "9.5"    # max Foldseek sensitivity (Steinegger 2017)
_PHOLD_EVALUE      = "1e-3"
_PHOLD_CARD_EVALUE = "1e-10"  # strict AMR/VF threshold

# ── Phold plot parameters ─────────────────────────────────────────────────────
_PLOT_DPI         = "600"
_PLOT_ANNOTATIONS = "1"       # label all annotated functions
_PLOT_LABEL_SIZE  = "8"
_PLOT_INTERVAL    = "5000"    # axis tick every 5 kb

# ── Phynteny parameters ───────────────────────────────────────────────────────
_PHYNTENY_CONFIDENCE = 0.7


# ── Public entry point ────────────────────────────────────────────────────────

def run(
    cfg:       Config,
    sample_id: str,
    genome:    Path,
    force:     bool = False,
) -> Path:
    """
    Annotate a single phage candidate genome (three tiers + plot).

    Parameters
    ----------
    genome : FASTA from quality/annotation_ready/

    Returns
    -------
    final_gbk : Phynteny GBK if available, else Phold GBK, else Pharokka GBK
    """
    require_tools(*TOOLS)

    genome       = Path(genome)
    candidate_id = genome.stem

    out_dir = cfg.results(STEP)
    rpt_dir = cfg.reports(STEP)
    cdir    = out_dir / candidate_id
    mkdirs(cdir, rpt_dir)

    log_step(f"Module 06 — annotation [{candidate_id}]")

    if not genome.exists():
        log_error(f"  Genome not found: {genome}")
        raise FileNotFoundError(genome)

    s     = fasta_stats(genome)
    n_ctg = s["n"]
    mode  = "meta" if n_ctg > 1 else "single"

    use_gpu    = getattr(cfg.annotate, "phold_gpu",           True)
    batch_size = getattr(cfg.annotate, "phold_batch_size",    1)
    phynteny_c = getattr(cfg.annotate, "phynteny_confidence", _PHYNTENY_CONFIDENCE)
    gpu_tag    = "--foldseek_gpu" if use_gpu else "--cpu"

    log_info(
        f"  Genome   : {human_size(genome)} | {n_ctg} contig(s) | "
        f"largest={s['largest_bp']:,} bp | mode={mode}"
    )
    log_info(f"  [1/4] Pharokka  : PHANOTATE/prodigal-gv + PHROGs MMseqs2 + PyHMMER")
    log_info(f"  [2/4] Phold     : ProstT5 + Foldseek --hyps --finetune {gpu_tag}")
    log_info(f"  [3/4] Phynteny  : LSTM + transformer circular attention (confidence ≥{phynteny_c})")
    log_info(f"  [4/4] Phold plot: circular map from Phynteny GBK (--all --dpi {_PLOT_DPI})")

    pharokka_dir = cdir / "pharokka"
    phold_dir    = cdir / "phold"
    phynteny_dir = cdir / "phynteny"
    plots_dir    = cdir / "plots"
    mkdirs(pharokka_dir, phold_dir, phynteny_dir, plots_dir)

    # Output file paths
    gbk_pharokka = pharokka_dir / f"{candidate_id}.gbk"
    gbk_phold    = phold_dir    / f"{candidate_id}.gbk"   # phold run -p {id} → {id}.gbk

    pharokka_stats: dict = {}
    phold_stats:    dict = {}
    phynteny_stats: dict = {}

    # ── Pipeline steps ────────────────────────────────────────────────────────
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
        pharokka_stats = _parse_pharokka_tsv(
            pharokka_dir / f"{candidate_id}_cds_functions.tsv"
        )
        if pharokka_stats:
            log_ok(
                f"  [Tier 1 Pharokka] {pharokka_stats['total_cds']} CDS | "
                f"annotated={pharokka_stats['annotated']} ({pharokka_stats['pct_annotated']}) | "
                f"hypothetical={pharokka_stats['hypothetical']}"
            )
        progress.advance(task)

        # 2/4 — Phold
        progress.update(task, description="[2/4] Phold     — structure-based upgrade")
        _run_phold(cfg, candidate_id, gbk_pharokka, phold_dir, gbk_phold, rpt_dir, force)
        phold_stats = _parse_phold_per_cds(phold_dir, candidate_id)
        if phold_stats:
            cum = pharokka_stats.get("annotated", 0) + phold_stats["upgraded"]
            tot = pharokka_stats.get("total_cds", 0)
            pct = f"{cum / tot * 100:.1f}%" if tot else "?"
            log_ok(
                f"  [Tier 2 Phold] Δ+{phold_stats['upgraded']} upgraded → "
                f"cumulative {cum}/{tot} ({pct})"
            )
        progress.advance(task)

        # 3/4 — Phynteny Transformer
        progress.update(task, description="[3/4] Phynteny  — synteny + ESM2 upgrade")
        _run_phynteny(cfg, candidate_id, gbk_phold, phynteny_dir, rpt_dir, force)
        gbk_phynteny = _find_phynteny_gbk(phynteny_dir, candidate_id)
        phynteny_stats = _parse_phynteny_per_cds(phynteny_dir, candidate_id)
        if phynteny_stats:
            cum = (pharokka_stats.get("annotated", 0) +
                   phold_stats.get("upgraded", 0) +
                   phynteny_stats["upgraded"])
            tot = pharokka_stats.get("total_cds", 0)
            pct = f"{cum / tot * 100:.1f}%" if tot else "?"
            log_ok(
                f"  [Tier 3 Phynteny] Δ+{phynteny_stats['upgraded']} predicted → "
                f"cumulative {cum}/{tot} ({pct})"
            )
        progress.advance(task)

        # 4/4 — Phold plot (from Phynteny GBK = final annotations)
        progress.update(task, description="[4/4] Phold plot — circular genome map")
        plot_input = gbk_phynteny if gbk_phynteny else gbk_phold
        _run_phold_plot(candidate_id, plot_input, plots_dir, rpt_dir, force)
        progress.advance(task)

    # ── Final output ──────────────────────────────────────────────────────────
    final_gbk = gbk_phynteny or gbk_phold if gbk_phold.exists() else gbk_pharokka

    _save_tsv(
        candidate_id, sample_id, genome,
        pharokka_stats, phold_stats, phynteny_stats,
        rpt_dir / "annotation_summary.tsv",
    )

    _print_completion_panel(
        candidate_id, final_gbk, pharokka_dir, plots_dir, phynteny_dir,
        pharokka_stats, phold_stats, phynteny_stats,
    )

    log_step(f"Module 06 completed ✓  [{candidate_id}]")
    log_ok(f"  GBK (final)  : {final_gbk}")
    log_ok(f"  GFF          : {pharokka_dir / (candidate_id + '.gff')}")
    log_ok(f"  Plots        : {plots_dir}")
    log_info("Next: phageflow safety --sample-id <id> --genome <genome.fasta>")

    return final_gbk


# ── Step 1: Pharokka ──────────────────────────────────────────────────────────

def _run_pharokka(cfg, candidate_id, genome, pharokka_dir, gbk_out, n_ctg, rpt_dir, force):
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
        "--dnaapler",               # boolean flag in v1.4+ — reorients to terminase
        "--threads", str(cfg.threads),
        "--force",
    ]

    if n_ctg > 1:
        cmd += ["-m", "--meta_hmm"]
        log_info(
            f"  [Pharokka] mode=meta ({n_ctg} contigs) | "
            f"gene caller=prodigal-gv | MMseqs2 + PyHMMER (--meta_hmm)"
        )
    else:
        log_info("  [Pharokka] mode=single | gene caller=PHANOTATE | MMseqs2 + PyHMMER")

    try:
        run_silent(cmd, log_file=rpt_dir / f"{candidate_id}_pharokka.log")
        log_ok("  [Pharokka] OK")
    except Exception as e:
        log_warn(f"  [Pharokka] warning: {e}")


# ── Step 2: Phold run ─────────────────────────────────────────────────────────

def _run_phold(cfg, candidate_id, input_gbk, phold_dir, gbk_out, rpt_dir, force):
    """
    Phold run with -p {candidate_id} → output GBK is {candidate_id}.gbk
    """
    if gbk_out.exists() and gbk_out.stat().st_size > 0 and not force:
        log_info("  [Phold] already exists — skipping  (--force to re-run)")
        return

    if not input_gbk.exists():
        log_warn("  [Phold] skipped — Pharokka GBK not found")
        return

    phold_db   = cfg.databases.phold
    use_gpu    = getattr(cfg.annotate, "phold_gpu",        True)
    batch_size = getattr(cfg.annotate, "phold_batch_size", 1)

    cmd = [
        "phold", "run",
        "-i", str(input_gbk),
        "-o", str(phold_dir),
        "-p", candidate_id,
        "-d", str(phold_db),
        "-t", str(cfg.threads),
        "--hyps",
        "--finetune",
        "-s", _PHOLD_SENSITIVITY,
        "-e", _PHOLD_EVALUE,
        "--card_vfdb_evalue", _PHOLD_CARD_EVALUE,
        "--batch_size", str(batch_size),
        "--force",
        "--foldseek_gpu" if use_gpu else "--cpu",
    ]

    log_info(
        f"  [Phold] {('--foldseek_gpu' if use_gpu else '--cpu')} | "
        f"--hyps --finetune | sensitivity={_PHOLD_SENSITIVITY}"
    )

    try:
        run_silent(cmd, log_file=rpt_dir / f"{candidate_id}_phold.log")
        log_ok("  [Phold] OK")
    except Exception as e:
        log_warn(f"  [Phold] warning: {e}")


# ── Step 3: Phynteny Transformer ──────────────────────────────────────────────

def _run_phynteny(cfg, candidate_id, input_gbk, phynteny_dir, rpt_dir, force):
    """
    Phynteny takes Phold GBK as input (confirmed workflow in benchmarking).
    Default output: phynteny_transformer.gbk + phynteny_per_cds_functions.tsv
    """
    if not input_gbk.exists():
        log_warn("  [Phynteny] skipped — Phold GBK not found")
        return

    existing = list(phynteny_dir.glob("*.gbk"))
    if existing and not force:
        log_info(f"  [Phynteny] already exists — skipping  (--force to re-run)")
        return

    models_path = cfg.databases.phynteny / "models"
    if not models_path.exists():
        log_warn(f"  [Phynteny] models not found: {models_path}")

    confidence = getattr(cfg.annotate, "phynteny_confidence", _PHYNTENY_CONFIDENCE)

    cmd = [
        "phynteny_transformer",
        str(input_gbk),
        "-o",                     str(phynteny_dir),
        "--prefix",               candidate_id,
        "-m",                     str(models_path),
        "--confidence-threshold", str(confidence),
        "-f",
    ]

    log_info(f"  [Phynteny] models={models_path.parent.name}/models | confidence ≥{confidence}")

    try:
        run_silent(cmd, log_file=rpt_dir / f"{candidate_id}_phynteny.log")
        log_ok("  [Phynteny] OK")
    except Exception as e:
        log_warn(f"  [Phynteny] warning: {e}")


def _find_phynteny_gbk(phynteny_dir: Path, candidate_id: str) -> Optional[Path]:
    """
    Locate Phynteny output GBK with fallbacks for different naming conventions.
    Default (no prefix): phynteny_transformer.gbk
    With --prefix {id}: possibly {id}.gbk or {id}_phynteny_transformer.gbk
    """
    candidates = [
        phynteny_dir / "phynteny_transformer.gbk",
        phynteny_dir / f"{candidate_id}.gbk",
        phynteny_dir / f"{candidate_id}_phynteny_transformer.gbk",
    ]
    found = next((p for p in candidates if p.exists()), None)
    if found is None:
        # Try any GBK in the directory
        gbks = list(phynteny_dir.glob("*.gbk"))
        found = gbks[0] if gbks else None
    return found


# ── Step 4: Phold plot ────────────────────────────────────────────────────────

def _run_phold_plot(candidate_id, input_gbk, plots_dir, rpt_dir, force):
    """
    Generate circular genome plot(s) from the final annotated GBK.
    Input is the Phynteny GBK (most annotated) — the plot reflects all 3 tiers.
    --all : one plot per contig (multi-contig genomes generate N plots).
    --dpi 600 : publication-quality.
    -p must match prefix used in phold run.
    """
    if not input_gbk or not input_gbk.exists():
        log_warn("  [Phold plot] skipped — no input GBK found")
        return

    existing = list(plots_dir.glob("*.png")) + list(plots_dir.glob("*.svg"))
    if existing and not force:
        log_info(f"  [Phold plot] already exists ({len(existing)} file(s)) — skipping")
        return

    log_info(f"  [Phold plot] input: {input_gbk.name} (Phynteny final GBK)")

    cmd = [
        "phold", "plot",
        "-i",            str(input_gbk),
        "-o",            str(plots_dir),
        "-p",            candidate_id,
        "--all",
        "--dpi",         _PLOT_DPI,
        "--annotations", _PLOT_ANNOTATIONS,
        "--label_size",  _PLOT_LABEL_SIZE,
        "--interval",    _PLOT_INTERVAL,
        "-f",
    ]

    try:
        run_silent(cmd, log_file=rpt_dir / f"{candidate_id}_phold_plot.log")
        plots = list(plots_dir.glob("*.png")) + list(plots_dir.glob("*.svg"))
        log_ok(f"  [Phold plot] {len(plots)} plot(s) → {plots_dir}")
    except Exception as e:
        log_warn(f"  [Phold plot] warning: {e}")


# ── Parsers ───────────────────────────────────────────────────────────────────

def _parse_pharokka_tsv(tsv: Path) -> dict:
    """Parse Pharokka *_cds_functions.tsv for tier 1 annotation counts."""
    if not tsv.exists():
        return {}
    total = annotated = 0
    try:
        with open(tsv) as f:
            header = f.readline().strip().split("\t")
            col    = {h: i for i, h in enumerate(header)}
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
    return {"total_cds": total, "annotated": annotated,
            "hypothetical": hypo, "pct_annotated": pct}


def _parse_phold_per_cds(phold_dir: Path, candidate_id: str) -> dict:
    """
    Count hypothetical proteins upgraded by Phold.
    Phold generates {candidate_id}_per_cds_functions.tsv (similar to Pharokka).
    Falls back to glob search for the TSV.
    """
    tsv_candidates = [
        phold_dir / f"{candidate_id}_per_cds_functions.tsv",
        phold_dir / f"{candidate_id}_phold_per_cds_functions.tsv",
        phold_dir / "phold_per_cds_functions.tsv",
    ]
    tsv = next((p for p in tsv_candidates if p.exists()), None)
    if tsv is None:
        tsv_files = list(phold_dir.glob("*per_cds*.tsv")) + list(phold_dir.glob("*functions.tsv"))
        tsv = tsv_files[0] if tsv_files else None
    if tsv is None:
        return {}

    upgraded = 0
    try:
        with open(tsv) as f:
            header = f.readline().strip().split("\t")
            col    = {h: i for i, h in enumerate(header)}
            for line in f:
                if not line.strip():
                    continue
                row   = line.strip().split("\t")
                phrog = row[col.get("phrog", 0)] if row else ""
                cat   = row[col.get("category", 1)] if len(row) > 1 else ""
                if (phrog and phrog not in ("No_PHROG", "", "NA")) or \
                   (cat and cat.lower() not in ("unknown function", "hypothetical protein", "")):
                    upgraded += 1
    except Exception:
        return {}
    return {"upgraded": upgraded, "tsv": str(tsv)}


def _parse_phynteny_per_cds(phynteny_dir: Path, candidate_id: str) -> dict:
    """
    Count proteins upgraded by Phynteny via phynteny_per_cds_functions.tsv.
    This is the Phynteny-equivalent of pharokka_cds_functions.tsv.
    """
    tsv_candidates = [
        phynteny_dir / f"{candidate_id}_phynteny_per_cds_functions.tsv",
        phynteny_dir / "phynteny_per_cds_functions.tsv",
        phynteny_dir / f"{candidate_id}_per_cds_functions.tsv",
    ]
    tsv = next((p for p in tsv_candidates if p.exists()), None)
    if tsv is None:
        tsv_files = list(phynteny_dir.glob("*per_cds*.tsv"))
        tsv = tsv_files[0] if tsv_files else None
    if tsv is None:
        return {}

    upgraded = 0
    confidence = _PHYNTENY_CONFIDENCE
    try:
        with open(tsv) as f:
            header = f.readline().strip().split("\t")
            col    = {h.lower().strip(): i for i, h in enumerate(header)}
            conf_col = next(
                (col[k] for k in ("confidence", "score", "phynteny_score") if k in col),
                None
            )
            cat_col = next(
                (col[k] for k in ("category", "function", "phrog_category") if k in col),
                None
            )
            for line in f:
                if not line.strip():
                    continue
                row = line.strip().split("\t")
                cat_ok  = False
                conf_ok = conf_col is None  # no conf col → accept any annotation
                if cat_col is not None and cat_col < len(row):
                    cat_ok = row[cat_col].strip().lower() not in (
                        "", "unknown function", "hypothetical protein", "na", "none"
                    )
                if conf_col is not None and conf_col < len(row):
                    try:
                        conf_ok = float(row[conf_col]) >= confidence
                    except (ValueError, TypeError):
                        conf_ok = False
                if cat_ok and conf_ok:
                    upgraded += 1
    except Exception:
        return {}
    return {"upgraded": upgraded}


# ── Rich display helpers ──────────────────────────────────────────────────────

def _print_completion_panel(
    candidate_id, final_gbk, pharokka_dir, plots_dir, phynteny_dir,
    pharokka_stats, phold_stats, phynteny_stats,
):
    total   = pharokka_stats.get("total_cds",    0)
    p1      = pharokka_stats.get("annotated",    0)
    p2      = phold_stats.get("upgraded",        0)
    p3      = phynteny_stats.get("upgraded",     0)
    final_n = p1 + p2 + p3
    final_h = total - final_n
    final_p = f"{final_n / total * 100:.1f}%" if total else "?"
    plots   = list(plots_dir.glob("*.png")) + list(plots_dir.glob("*.svg"))

    text = Text()
    text.append("Annotation improvement:\n", style="bold")
    text.append(
        f"  Tier 1 Pharokka : {p1:>4}/{total} ({pharokka_stats.get('pct_annotated', '?')})\n",
        style="cyan",
    )
    text.append(
        f"  Tier 2 Phold    : +{p2:>3}  → {p1+p2:>4}/{total} "
        f"({(p1+p2)/total*100:.1f}% if total else '?')\n",
        style="cyan",
    ) if total else None
    text.append(
        f"  Tier 3 Phynteny : +{p3:>3}  → {final_n:>4}/{total} ({final_p})\n",
        style="bold green",
    )
    text.append(
        f"  Still unknown   : {final_h:>4}/{total}\n\n",
        style="dim white",
    )
    text.append("GBK (final) : ", style="dim white")
    text.append(str(final_gbk) + "\n", style="white")
    text.append("Plots       : ", style="dim white")
    text.append(f"{plots_dir}  ({len(plots)} file(s))\n", style="white")
    text.append("Phynteny    : ", style="dim white")
    text.append(str(phynteny_dir), style="white")

    console.print(Panel(
        text,
        title=f"[bold cyan]Annotation complete — {candidate_id}[/bold cyan]",
        border_style="cyan",
        padding=(0, 2),
        width=90,
    ))


# ── TSV summary ───────────────────────────────────────────────────────────────

_TSV_HEADERS = [
    "candidate_id", "sample_id", "genome_bp",
    "total_cds",
    "tier1_pharokka", "tier1_pct",
    "tier2_phold_delta",
    "tier3_phynteny_delta",
    "final_annotated", "final_pct", "final_hypothetical",
]


def _save_tsv(
    candidate_id, sample_id, genome,
    pharokka_stats, phold_stats, phynteny_stats, path,
):
    rows: dict[str, dict] = {}
    if path.exists():
        with open(path) as f:
            old_hdrs = f.readline().rstrip("\n").split("\t")
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if cols and cols[0]:
                    rows[cols[0]] = dict(zip(old_hdrs, cols))

    s     = fasta_stats(genome)
    total = pharokka_stats.get("total_cds",  0)
    p1    = pharokka_stats.get("annotated",  0)
    p2    = phold_stats.get("upgraded",      0)
    p3    = phynteny_stats.get("upgraded",   0)
    final = p1 + p2 + p3
    hypo  = total - final

    rows[candidate_id] = {
        "candidate_id":        candidate_id,
        "sample_id":           sample_id,
        "genome_bp":           str(s.get("total_bp", 0)),
        "total_cds":           str(total),
        "tier1_pharokka":      str(p1),
        "tier1_pct":           pharokka_stats.get("pct_annotated", ""),
        "tier2_phold_delta":   str(p2),
        "tier3_phynteny_delta":str(p3),
        "final_annotated":     str(final),
        "final_pct":           f"{final/total*100:.1f}%" if total else "",
        "final_hypothetical":  str(hypo),
    }

    with open(path, "w") as f:
        f.write("\t".join(_TSV_HEADERS) + "\n")
        for row in rows.values():
            f.write("\t".join(row.get(h, "") for h in _TSV_HEADERS) + "\n")

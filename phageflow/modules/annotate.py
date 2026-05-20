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

Tier 3 — Phynteny Transformer (Grigson et al. 2025, bioRxiv):
    Input : Phold GBK
    Adds /phynteny_category, /phynteny_score, /phynteny_confidence qualifiers
    per CDS. Delta is counted from these qualifiers (confidence ≥ threshold).
    Does NOT modify /product or /function — those reflect the Phold state.

Parser strategy (all tiers):
    All annotation counts are derived from the GBK output files, not from
    tool-specific TSVs.  This is more reliable because:
      - Pharokka _cds_functions.tsv is a per-PHROG/category summary (not per-CDS)
      - Phold per_cds TSV naming varies across versions
      - Phynteny writes its predictions as GBK qualifiers, not in a plain TSV
    Parsing the GBK directly gives per-CDS ground truth regardless of version.
"""

from __future__ import annotations
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
from phageflow.utils.logger import log_step, log_info, log_ok, log_warn, log_error, console
from phageflow.utils.tools import require_tools, run_silent, human_size, mkdirs, fasta_stats

STEP  = "06_annotation"
TOOLS = ["pharokka.py", "phold", "phynteny_transformer"]

# ── Phold parameters ──────────────────────────────────────────────────────────
_PHOLD_SENSITIVITY = "9.5"
_PHOLD_EVALUE      = "1e-3"
_PHOLD_CARD_EVALUE = "1e-10"

# ── Phold plot parameters ─────────────────────────────────────────────────────
_PLOT_DPI         = "600"
_PLOT_ANNOTATIONS = "1"
_PLOT_LABEL_SIZE  = "8"
_PLOT_INTERVAL    = "5000"

# ── Phynteny parameters ───────────────────────────────────────────────────────
_PHYNTENY_CONFIDENCE = 0.7

# ── Annotation category sets ──────────────────────────────────────────────────
_UNKNOWN_CATS = frozenset({
    "unknown function", "hypothetical protein", "", "na", "none",
    "unknown", "uncharacterized protein",
})


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

    gbk_pharokka = pharokka_dir / f"{candidate_id}.gbk"
    gbk_phold    = phold_dir    / f"{candidate_id}.gbk"

    pharokka_stats: dict = {}
    phold_stats:    dict = {}
    phynteny_stats: dict = {}

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
        pharokka_stats = _parse_pharokka_gbk(gbk_pharokka)
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
        phold_stats = _parse_phold_gbk_delta(gbk_phold, gbk_pharokka)
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
        phynteny_stats = _parse_phynteny_gbk_delta(gbk_phynteny, phynteny_c)
        if phynteny_stats:
            phold_cum = (pharokka_stats.get("annotated", 0) +
                         phold_stats.get("upgraded", 0))
            cum = phold_cum + phynteny_stats["upgraded"]
            tot = pharokka_stats.get("total_cds", 0)
            pct = f"{cum / tot * 100:.1f}%" if tot else "?"
            log_ok(
                f"  [Tier 3 Phynteny] Δ+{phynteny_stats['upgraded']} predicted → "
                f"cumulative {cum}/{tot} ({pct})"
            )
        progress.advance(task)

        # 4/4 — Phold plot
        progress.update(task, description="[4/4] Phold plot — circular genome map")
        plot_input = gbk_phynteny if gbk_phynteny else gbk_phold
        _run_phold_plot(candidate_id, plot_input, plots_dir, rpt_dir, force)
        progress.advance(task)

    final_gbk = gbk_phynteny or (gbk_phold if gbk_phold.exists() else gbk_pharokka)

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
        "--dnaapler",
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
    """Locate Phynteny output GBK — checks known naming patterns before glob."""
    candidates = [
        phynteny_dir / "phynteny_transformer.gbk",
        phynteny_dir / "phynteny.gbk",
        phynteny_dir / f"{candidate_id}.gbk",
        phynteny_dir / f"{candidate_id}_phynteny_transformer.gbk",
    ]
    found = next((p for p in candidates if p.exists()), None)
    if found is None:
        gbks = list(phynteny_dir.glob("*.gbk"))
        found = gbks[0] if gbks else None
    return found


# ── Step 4: Phold plot ────────────────────────────────────────────────────────

def _run_phold_plot(candidate_id, input_gbk, plots_dir, rpt_dir, force):
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


# ── GBK parsers ───────────────────────────────────────────────────────────────

def _gbk_cds_qualifiers(gbk_path: Path) -> list[dict]:
    """
    Return a list of qualifier dicts, one per CDS feature in the GBK.

    Extracted keys: locus_tag, product, function, phrog,
                    phynteny_category, phynteny_confidence, phynteny_score.

    Rationale: parsing the GBK is more reliable than reading tool-specific
    TSVs because:
      - Pharokka _cds_functions.tsv is a per-category summary (not per-CDS)
      - Phold per_cds TSV naming conventions vary by version
      - Phynteny writes predictions as /phynteny_* GBK qualifiers only
    """
    if not gbk_path or not gbk_path.exists():
        return []

    try:
        text = gbk_path.read_text(errors="replace")
    except Exception:
        return []

    results = []
    for block in text.split("     CDS ")[1:]:
        def _q(name: str) -> str:
            m = re.search(rf'/{name}="((?:[^"\\]|\\.|\n\s+)*)"', block)
            if not m:
                return ""
            # Collapse wrapped lines: GBK wraps long values across indented lines
            return re.sub(r'\n\s+', ' ', m.group(1)).strip()

        results.append({
            "locus_tag":           _q("locus_tag"),
            "product":             _q("product").lower(),
            "function":            _q("function").lower(),
            "phrog":               _q("phrog"),
            "phynteny_category":   _q("phynteny_category").lower(),
            "phynteny_confidence": _q("phynteny_confidence"),
            "phynteny_score":      _q("phynteny_score"),
        })
    return results


def _is_annotated(cds: dict) -> bool:
    """Return True if the CDS has a real (non-hypothetical) annotation."""
    return (
        cds["product"]  not in _UNKNOWN_CATS or
        cds["function"] not in _UNKNOWN_CATS
    )


def _parse_pharokka_gbk(gbk_path: Path) -> dict:
    """
    Parse the Pharokka GBK to get per-CDS annotation counts.

    Replaces _parse_pharokka_tsv: the _cds_functions.tsv produced by
    Pharokka v1.9+ is a per-PHROG/category summary (one row per unique
    PHROG category), not a per-CDS table. Parsing the GBK directly gives
    the correct total_cds, annotated, and hypothetical counts.
    """
    cds_list = _gbk_cds_qualifiers(gbk_path)
    if not cds_list:
        return {}

    total     = len(cds_list)
    annotated = sum(1 for c in cds_list if _is_annotated(c))
    hypo      = total - annotated
    pct       = f"{annotated / total * 100:.1f}%" if total else "0.0%"

    return {
        "total_cds":    total,
        "annotated":    annotated,
        "hypothetical": hypo,
        "pct_annotated": pct,
    }


def _parse_phold_gbk_delta(phold_gbk: Path, pharokka_gbk: Path) -> dict:
    """
    Count proteins newly annotated by Phold (delta vs Pharokka).

    Delta = CDS where /product changed FROM "hypothetical protein"
    TO a real function in the Phold GBK.  Phold with --hyps only processes
    hypothetical proteins, so any product change from hypothetical → annotated
    is a Phold upgrade.

    Reclassification of already-annotated CDS (e.g. "DNA" → "DNA, RNA and
    nucleotide metabolism") is tracked separately for transparency but is NOT
    counted in 'upgraded' because the protein was already functionally annotated.
    """
    pk_list = _gbk_cds_qualifiers(pharokka_gbk)
    ph_list = _gbk_cds_qualifiers(phold_gbk)

    if not pk_list or not ph_list:
        return {}

    pk_map = {c["locus_tag"]: c for c in pk_list if c["locus_tag"]}
    ph_map = {c["locus_tag"]: c for c in ph_list if c["locus_tag"]}

    upgraded       = 0   # hypothetical → annotated
    reclassified   = 0   # annotated → differently annotated (category change)

    for locus, ph_cds in ph_map.items():
        if locus not in pk_map:
            continue
        pk_cds = pk_map[locus]
        was_hypo  = not _is_annotated(pk_cds)
        now_ann   = _is_annotated(ph_cds)

        if was_hypo and now_ann:
            upgraded += 1
        elif (not was_hypo) and (pk_cds["function"] != ph_cds["function"]):
            reclassified += 1

    if reclassified:
        log_info(
            f"  [Phold] {reclassified} CDS reclassified (function category updated, "
            f"not counted in Δ since product was already annotated)"
        )

    total_phold = sum(1 for c in ph_list if _is_annotated(c))
    return {
        "upgraded":      upgraded,
        "reclassified":  reclassified,
        "total_annotated": total_phold,
    }


def _parse_phynteny_gbk_delta(gbk_path: Optional[Path], confidence: float) -> dict:
    """
    Count Phynteny predictions from /phynteny_category qualifiers.

    Phynteny Transformer does NOT modify /product or /function — it adds
    three new qualifiers per CDS:
        /phynteny_category   — predicted PHROG category
        /phynteny_confidence — prediction confidence (0–1)
        /phynteny_score      — raw model score

    Delta = CDS that were still hypothetical after Phold AND received a
    /phynteny_category with confidence ≥ threshold.
    """
    if not gbk_path:
        return {}

    cds_list = _gbk_cds_qualifiers(gbk_path)
    if not cds_list:
        return {}

    upgraded = 0
    for cds in cds_list:
        cat  = cds["phynteny_category"]
        conf_s = cds["phynteny_confidence"]
        if not cat or cat in _UNKNOWN_CATS:
            continue
        try:
            conf = float(conf_s)
        except (ValueError, TypeError):
            conf = 0.0
        # Count only CDS that are still hypothetical in the /product qualifier
        # (Phold didn't already annotate them) and got a confident Phynteny call
        if cds["product"] in _UNKNOWN_CATS and conf >= confidence:
            upgraded += 1

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
    final_h = max(0, total - final_n)
    plots   = list(plots_dir.glob("*.png")) + list(plots_dir.glob("*.svg"))

    text = Text()
    text.append("Annotation improvement:\n", style="bold")
    text.append(
        f"  Tier 1 Pharokka : {p1:>4}/{total} ({pharokka_stats.get('pct_annotated', '?')})\n",
        style="cyan",
    )

    if total:
        pct2 = f"{(p1 + p2) / total * 100:.1f}%"
        text.append(
            f"  Tier 2 Phold    : +{p2:>3}  → {p1 + p2:>4}/{total} ({pct2})\n",
            style="cyan",
        )
        pct3 = f"{final_n / total * 100:.1f}%"
        text.append(
            f"  Tier 3 Phynteny : +{p3:>3}  → {final_n:>4}/{total} ({pct3})\n",
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
    hypo  = max(0, total - final)

    rows[candidate_id] = {
        "candidate_id":         candidate_id,
        "sample_id":            sample_id,
        "genome_bp":            str(s.get("total_bp", 0)),
        "total_cds":            str(total),
        "tier1_pharokka":       str(p1),
        "tier1_pct":            pharokka_stats.get("pct_annotated", ""),
        "tier2_phold_delta":    str(p2),
        "tier3_phynteny_delta": str(p3),
        "final_annotated":      str(final),
        "final_pct":            f"{final/total*100:.1f}%" if total else "",
        "final_hypothetical":   str(hypo),
    }

    with open(path, "w") as f:
        f.write("\t".join(_TSV_HEADERS) + "\n")
        for row in rows.values():
            f.write("\t".join(row.get(h, "") for h in _TSV_HEADERS) + "\n")

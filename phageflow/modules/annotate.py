"""PhageFlow Module 06 — Structural and functional annotation.

Three-tier pipeline: Pharokka → Phold (run + plot) → Phynteny Transformer

Output structure — one directory per candidate genome:
    results/06_annotation/{candidate_id}/
        pharokka/    gene calling + PHROG annotation
        phold/       structure-based annotation (GBK, TSV)
        plots/       circular genome plots (PNG, 600 dpi)
        phynteny/    ESM2 function predictions

Tier 1 — Pharokka v1.9+ (Bouras et al. 2023, Bioinformatics 39:btac776):
    Gene calling + PHROG sequence-based functional annotation.

    Single-contig mode (default):
        PHANOTATE: optimal phage gene caller; out-performs Prodigal for
        dense, overlapping ORFs (McNair et al. 2019, Bioinformatics 35:4537).
        MMseqs2 + PyHMMER against PHROGs: maximum annotation depth.
        --dnaapler (flag only in v1.4+): reorients genome to start at
        terminase large subunit (Casjens & Thuman-Commike 2011, Virology 418:4).

    Multi-contig mode (-m --meta_hmm):
        prodigal-gv: handles diverse genetic codes in fragmented assemblies.
        --meta_hmm: overrides default --mmseqs2_only in meta mode, enabling
        both MMseqs2 and PyHMMER for equivalent annotation depth to single mode.

Tier 2 — Phold v1.2+ (Bouras et al. 2024, Bioinformatics):
    Structure-based annotation via ProstT5 + Foldseek, followed by circular
    genome plot generation.

    --hyps           : processes only hypothetical proteins (Bouras et al. 2024)
    --finetune       : ProstT5 fine-tuned on phage proteins (best phage coverage)
    -s 9.5           : maximum Foldseek sensitivity (Steinegger et al. 2017,
                       Nat Methods 14:1026)
    --card_vfdb_evalue 1e-10 : strict threshold for AMR/virulence hits
    --foldseek_gpu / --cpu   : configurable via annotate.phold_gpu in config

    phold plot --all --dpi 600 : publication-quality circular genome plot.
    --all plots every contig independently (essential for multi-contig genomes).

Tier 3 — Phynteny Transformer:
    ESM2 + transformer function prediction. Third independent annotation layer.
    -m databases/phynteny_db/models : local pre-trained models.
    --confidence-threshold 0.7 : conservative threshold for assignments.
"""

from __future__ import annotations
import shutil
from pathlib import Path

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
_PHOLD_SENSITIVITY = "9.5"    # maximum Foldseek sensitivity (Steinegger 2017)
_PHOLD_EVALUE      = "1e-3"   # Foldseek e-value (Phold default)
_PHOLD_CARD_EVALUE = "1e-10"  # strict threshold for AMR/VF hits

# ── Phold plot parameters ─────────────────────────────────────────────────────
_PLOT_DPI         = "600"     # publication-quality resolution
_PLOT_ANNOTATIONS = "1"       # label all predicted functions (0–1)
_PLOT_LABEL_SIZE  = "8"       # annotation label font size
_PLOT_INTERVAL    = "5000"    # axis tick interval (bp)

# ── Phynteny parameters ───────────────────────────────────────────────────────
_PHYNTENY_CONFIDENCE = 0.7    # conservative threshold for function assignment


# ── Public entry point ────────────────────────────────────────────────────────

def run(
    cfg:       Config,
    sample_id: str,
    genome:    Path,
    force:     bool = False,
) -> Path:
    """
    Annotate a single phage candidate genome.

    Each candidate gets its own output directory named after the genome stem
    (e.g. Herelleviridae_candidate_001/).

    Parameters
    ----------
    genome : FASTA from quality module (annotation_ready/ directory)

    Returns
    -------
    final_gbk : Path to the final annotated GenBank (Phold output preferred)
    """
    require_tools(*TOOLS)

    genome       = Path(genome)
    candidate_id = genome.stem           # e.g. Herelleviridae_candidate_001

    out_dir      = cfg.results(STEP)
    rpt_dir      = cfg.reports(STEP)
    cdir         = out_dir / candidate_id  # one directory per candidate genome
    mkdirs(cdir, rpt_dir)

    log_step(f"Module 06 — annotation [{candidate_id}] · {cfg.threads} threads")

    if not genome.exists():
        log_error(f"  Genome not found: {genome}")
        raise FileNotFoundError(genome)

    s     = fasta_stats(genome)
    n_ctg = s["n"]
    mode  = "meta" if n_ctg > 1 else "single"

    log_info(
        f"  Genome  : {human_size(genome)} | {n_ctg} contig(s) | "
        f"largest={s['largest_bp']:,} bp | mode={mode}"
    )

    use_gpu    = getattr(cfg.annotate, "phold_gpu",           True)
    batch_size = getattr(cfg.annotate, "phold_batch_size",    1)
    phynteny_c = getattr(cfg.annotate, "phynteny_confidence", _PHYNTENY_CONFIDENCE)
    gpu_tag    = "--foldseek_gpu" if use_gpu else "--cpu"

    log_info(f"  [1/4] Pharokka  : PHANOTATE/prodigal-gv + PHROGs MMseqs2 + PyHMMER")
    log_info(f"  [2/4] Phold run : ProstT5 + Foldseek (--hyps --finetune {gpu_tag})")
    log_info(f"  [3/4] Phold plot: circular genome (--all --dpi {_PLOT_DPI})")
    log_info(f"  [4/4] Phynteny  : ESM2 + transformer (confidence ≥{phynteny_c})")

    pharokka_dir = cdir / "pharokka"
    phold_dir    = cdir / "phold"
    plots_dir    = cdir / "plots"
    phynteny_dir = cdir / "phynteny"
    mkdirs(pharokka_dir, phold_dir, plots_dir, phynteny_dir)

    gbk_pharokka = pharokka_dir / f"{candidate_id}.gbk"
    gbk_phold    = phold_dir    / f"{candidate_id}_phold.gbk"

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
                f"  [Pharokka] {pharokka_stats['total_cds']} CDS | "
                f"annotated={pharokka_stats['annotated']} "
                f"({pharokka_stats['pct_annotated']}) | "
                f"hypothetical={pharokka_stats['hypothetical']}"
            )
        progress.advance(task)

        # 2/4 — Phold run
        progress.update(task, description="[2/4] Phold run — structure-based upgrade")
        _run_phold(cfg, candidate_id, gbk_pharokka, phold_dir, gbk_phold, rpt_dir, force)
        phold_stats = _parse_phold_tsv(phold_dir / f"{candidate_id}_phold.tsv")
        if phold_stats:
            log_ok(
                f"  [Phold] {phold_stats['upgraded']} hypothetical protein(s) upgraded"
            )
        progress.advance(task)

        # 3/4 — Phold plot
        progress.update(task, description="[3/4] Phold plot — circular genome map")
        _run_phold_plot(candidate_id, gbk_phold, plots_dir, rpt_dir, force)
        progress.advance(task)

        # 4/4 — Phynteny Transformer
        progress.update(task, description="[4/4] Phynteny  — ESM2 function prediction")
        input_gbk = gbk_phold if gbk_phold.exists() and gbk_phold.stat().st_size > 0 \
                    else gbk_pharokka
        _run_phynteny(cfg, candidate_id, input_gbk, phynteny_dir, rpt_dir, force)
        phynteny_stats = _parse_phynteny_output(phynteny_dir)
        if phynteny_stats:
            log_ok(
                f"  [Phynteny] {phynteny_stats['upgraded']} additional protein(s) assigned "
                f"(confidence ≥{phynteny_c})"
            )
        progress.advance(task)

    # ── Final output ──────────────────────────────────────────────────────────
    final_gbk = gbk_phold if gbk_phold.exists() and gbk_phold.stat().st_size > 0 \
                else gbk_pharokka

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
    log_ok(f"  Phynteny     : {phynteny_dir}")
    log_info("Next: phageflow safety --sample-id <id> --genome <genome.fasta>")

    return final_gbk


# ── Step 1: Pharokka ──────────────────────────────────────────────────────────

def _run_pharokka(
    cfg:          Config,
    candidate_id: str,
    genome:       Path,
    pharokka_dir: Path,
    gbk_out:      Path,
    n_ctg:        int,
    rpt_dir:      Path,
    force:        bool,
) -> None:
    """
    Run Pharokka gene calling and PHROG annotation.

    Single-contig : PHANOTATE + MMseqs2 + PyHMMER + --dnaapler
    Multi-contig  : prodigal-gv + --meta_hmm (enables PyHMMER in meta mode)

    --dnaapler is a boolean flag in v1.4+ — no argument required.
    """
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
        log_info(
            "  [Pharokka] mode=single | gene caller=PHANOTATE | MMseqs2 + PyHMMER"
        )

    try:
        run_silent(cmd, log_file=rpt_dir / f"{candidate_id}_pharokka.log")
        log_ok("  [Pharokka] OK")
    except Exception as e:
        log_warn(f"  [Pharokka] warning: {e}")


# ── Step 2: Phold run ─────────────────────────────────────────────────────────

def _run_phold(
    cfg:          Config,
    candidate_id: str,
    input_gbk:    Path,
    phold_dir:    Path,
    gbk_out:      Path,
    rpt_dir:      Path,
    force:        bool,
) -> None:
    """
    Run Phold structure-based annotation on the Pharokka GenBank.

    --hyps       : only hypothetical proteins (Bouras et al. 2024 Supp)
    --finetune   : phage-finetuned ProstT5 model
    -s 9.5       : maximum Foldseek sensitivity (Steinegger et al. 2017)
    --card_vfdb_evalue 1e-10 : strict AMR/VF threshold
    """
    if gbk_out.exists() and gbk_out.stat().st_size > 0 and not force:
        log_info("  [Phold] already exists — skipping  (--force to re-run)")
        return

    if not input_gbk.exists():
        log_warn("  [Phold] skipped — Pharokka GBK not found")
        return

    phold_db   = cfg.databases.phold
    if not phold_db.exists():
        log_warn(f"  [Phold] database not found: {phold_db}")

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
    ]

    if use_gpu:
        cmd.append("--foldseek_gpu")
        log_info(
            f"  [Phold] --foldseek_gpu | --finetune | --hyps | "
            f"sensitivity={_PHOLD_SENSITIVITY} | batch_size={batch_size}"
        )
    else:
        cmd.append("--cpu")
        log_info(
            f"  [Phold] --cpu | --finetune | --hyps | "
            f"sensitivity={_PHOLD_SENSITIVITY}"
        )

    try:
        run_silent(cmd, log_file=rpt_dir / f"{candidate_id}_phold.log")
        log_ok("  [Phold] OK")
    except Exception as e:
        log_warn(f"  [Phold] warning: {e}")


# ── Step 3: Phold plot ────────────────────────────────────────────────────────

def _run_phold_plot(
    candidate_id: str,
    gbk_phold:    Path,
    plots_dir:    Path,
    rpt_dir:      Path,
    force:        bool,
) -> None:
    """
    Generate circular genome plot(s) from Phold GenBank output.

    --all  : plot every contig independently — required for multi-contig
             genomes (Herelleviridae, Ackermannviridae fragmented assemblies).
    --dpi 600 : publication-quality resolution (standard for journal figures).
    --annotations 1 : label all predicted functions.

    The prefix (-p) must match the prefix used in phold run.
    """
    if not gbk_phold.exists():
        log_warn("  [Phold plot] skipped — phold GBK not found")
        return

    existing = list(plots_dir.glob("*.png")) + list(plots_dir.glob("*.svg"))
    if existing and not force:
        log_info(f"  [Phold plot] already exists ({len(existing)} file(s)) — skipping")
        return

    cmd = [
        "phold", "plot",
        "-i",            str(gbk_phold),
        "-o",            str(plots_dir),
        "-p",            candidate_id,   # must match phold run prefix
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


# ── Step 4: Phynteny Transformer ──────────────────────────────────────────────

def _run_phynteny(
    cfg:          Config,
    candidate_id: str,
    input_gbk:    Path,
    phynteny_dir: Path,
    rpt_dir:      Path,
    force:        bool,
) -> None:
    """
    Run Phynteny Transformer ESM2-based function prediction.

    -m databases/phynteny_db/models : path to local pre-trained models.
    --confidence-threshold 0.7 : conservative cut-off for function assignments.
    """
    if not input_gbk.exists():
        log_warn("  [Phynteny] skipped — no input GBK found")
        return

    existing = list(phynteny_dir.glob("*.tsv"))
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

    log_info(
        f"  [Phynteny] models={models_path} | confidence ≥{confidence}"
    )

    try:
        run_silent(cmd, log_file=rpt_dir / f"{candidate_id}_phynteny.log")
        log_ok("  [Phynteny] OK")
    except Exception as e:
        log_warn(f"  [Phynteny] warning: {e}")


# ── Parsers ───────────────────────────────────────────────────────────────────

def _parse_pharokka_tsv(tsv: Path) -> dict:
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
                if cat and cat.lower() not in (
                    "unknown function", "hypothetical protein", ""
                ):
                    annotated += 1
    except Exception:
        return {}
    hypo = total - annotated
    pct  = f"{annotated / total * 100:.1f}%" if total else "0.0%"
    return {
        "total_cds":     total,
        "annotated":     annotated,
        "hypothetical":  hypo,
        "pct_annotated": pct,
    }


def _parse_phold_tsv(tsv: Path) -> dict:
    if not tsv.exists():
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
                if phrog and phrog not in ("No_PHROG", "", "NA"):
                    upgraded += 1
    except Exception:
        return {}
    return {"upgraded": upgraded}


def _parse_phynteny_output(phynteny_dir: Path) -> dict:
    if not phynteny_dir.exists():
        return {}
    tsv_files = sorted(phynteny_dir.glob("*.tsv"))
    if not tsv_files:
        return {}

    upgraded = total = 0
    for tsv in tsv_files:
        try:
            with open(tsv) as f:
                header = f.readline().strip().split("\t")
                col    = {h.lower().strip(): i for i, h in enumerate(header)}
                conf_col = next(
                    (col[k] for k in ("confidence", "score", "probability") if k in col),
                    None,
                )
                func_col = next(
                    (col[k] for k in ("function", "prediction", "category",
                                      "phrog", "product", "annotation") if k in col),
                    None,
                )
                for line in f:
                    if not line.strip():
                        continue
                    row = line.strip().split("\t")
                    total += 1
                    func_ok = conf_ok = False
                    if func_col is not None and func_col < len(row):
                        func_ok = row[func_col].strip().lower() not in (
                            "", "unknown function", "hypothetical protein", "na", "none"
                        )
                    if conf_col is not None and conf_col < len(row):
                        try:
                            conf_ok = float(row[conf_col]) >= _PHYNTENY_CONFIDENCE
                        except (ValueError, TypeError):
                            conf_ok = False
                    else:
                        conf_ok = True  # no conf column → accept any function
                    if func_ok and conf_ok:
                        upgraded += 1
        except Exception:
            continue
        if total > 0:
            break
    return {"upgraded": upgraded, "total": total}


# ── Rich display helpers ──────────────────────────────────────────────────────

def _print_completion_panel(
    candidate_id:   str,
    final_gbk:      Path,
    pharokka_dir:   Path,
    plots_dir:      Path,
    phynteny_dir:   Path,
    pharokka_stats: dict,
    phold_stats:    dict,
    phynteny_stats: dict,
) -> None:
    total_cds   = pharokka_stats.get("total_cds",    "?")
    pct_annot   = pharokka_stats.get("pct_annotated", "?")
    phold_up    = phold_stats.get("upgraded",          0)
    phynteny_up = phynteny_stats.get("upgraded",       0)
    plots       = list(plots_dir.glob("*.png")) + list(plots_dir.glob("*.svg"))

    text = Text()
    text.append("✓ ", style="bold green")
    text.append(
        f"Pharokka: {total_cds} CDS ({pct_annot})  |  "
        f"Phold ↑{phold_up}  |  Phynteny ↑{phynteny_up}  |  "
        f"Plots: {len(plots)}\n\n",
        style="cyan",
    )
    text.append("GBK (final) : ", style="dim white")
    text.append(str(final_gbk) + "\n", style="white")
    text.append("GFF         : ", style="dim white")
    text.append(str(pharokka_dir / (candidate_id + ".gff")) + "\n", style="white")
    text.append("Plots       : ", style="dim white")
    text.append(str(plots_dir) + "\n", style="white")
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
    "total_cds", "annotated", "pct_annotated", "hypothetical",
    "phold_upgraded", "phynteny_upgraded",
]


def _save_tsv(
    candidate_id:   str,
    sample_id:      str,
    genome:         Path,
    pharokka_stats: dict,
    phold_stats:    dict,
    phynteny_stats: dict,
    path:           Path,
) -> None:
    rows: dict[str, dict] = {}
    if path.exists():
        with open(path) as f:
            old_hdrs = f.readline().rstrip("\n").split("\t")
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if cols and cols[0]:
                    rows[cols[0]] = dict(zip(old_hdrs, cols))

    s = fasta_stats(genome)
    rows[candidate_id] = {
        "candidate_id":     candidate_id,
        "sample_id":        sample_id,
        "genome_bp":        str(s.get("total_bp", 0)),
        "total_cds":        str(pharokka_stats.get("total_cds",    "")),
        "annotated":        str(pharokka_stats.get("annotated",    "")),
        "pct_annotated":    pharokka_stats.get("pct_annotated",    ""),
        "hypothetical":     str(pharokka_stats.get("hypothetical", "")),
        "phold_upgraded":   str(phold_stats.get("upgraded",        "")),
        "phynteny_upgraded":str(phynteny_stats.get("upgraded",     "")),
    }

    with open(path, "w") as f:
        f.write("\t".join(_TSV_HEADERS) + "\n")
        for row in rows.values():
            f.write("\t".join(row.get(h, "") for h in _TSV_HEADERS) + "\n")

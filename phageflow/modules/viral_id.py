"""PhageFlow Module 04 — Viral identification with geNomad.

Tool:
    geNomad v1.8+ (Camargo et al. 2023, Nature Biotechnology):
        End-to-end classification of contigs as virus / plasmid / chromosome.
        Detects integrated proviruses and assigns viral taxonomy.

        --splits 8         : splits input for parallel processing;
                             each split processed independently then merged.
        --min-score 0.7    : virus score threshold (default); configurable via config.yaml.
                             Camargo et al. 2023: 0.7 yields ~97% precision on diverse viruses.
        --min-hallmarks 0  : minimum viral hallmark genes required for classification;
                             configurable via config.yaml (0 = score-only mode).
        --cleanup          : removes large intermediate files after run,
                             retaining only summary and FASTA outputs.

    Taxonomy:
        geNomad assigns ICTV-aligned taxonomy: Realm > Kingdom > Phylum >
        Class > Order > Family > Genus.  The most specific resolved level
        is reported as best_taxon (unclassified if no assignment).

Input  : results/03_assembly/combined/{sample}_contigs_nr.fasta
Output : results/04_viral_id/{sample}_virus.fna  ← CheckV input
         reports/04_viral_id/genomad_summary.tsv
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
from phageflow.utils.logger import log_step, log_info, log_ok, log_warn, console
from phageflow.utils.tools import require_tools, run_silent, mkdirs, fasta_stats, human_size

STEP  = "04_viral_id"
TOOLS = ["genomad"]

# ── Thresholds ────────────────────────────────────────────────────────────────
_VIRAL_FRACTION_WARN = 0.05   # < 5% viral contigs → unusual for purified phage
_SCORE_GOOD          = 0.90   # virus score: confident viral classification
_SCORE_WARN          = 0.70   # virus score: acceptable but borderline


# ── Public entry point ────────────────────────────────────────────────────────

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

    virus_fna = out_dir / f"{sample_id}_virus.fna"

    log_step(f"Module 04 — viral-id [{sample_id}]")

    if not contigs.exists():
        raise FileNotFoundError(f"[{sample_id}] contigs not found: {contigs}")

    s = fasta_stats(contigs)
    _print_input_table(sample_id, contigs, s)
    log_info(
        f"  geNomad end-to-end  (Camargo et al. 2023)  |  "
        f"min-score={cfg.genomad.min_score}  min-virus-hallmarks={cfg.genomad.min_hallmarks}"
    )
    log_info("  Classifying : virus / plasmid / chromosome + provirus")

    # ── geNomad working paths (fixed naming convention) ───────────────────────
    sdir   = out_dir / sample_id
    prefix = contigs.stem
    v_fna  = sdir / f"{prefix}_summary" / f"{prefix}_virus.fna"
    v_tsv  = sdir / f"{prefix}_summary" / f"{prefix}_virus_summary.tsv"

    # ── Pipeline steps ────────────────────────────────────────────────────────
    with Progress(
        SpinnerColumn(spinner_name="dots2"),
        TextColumn("[bold cyan]{task.description:<48}"),
        BarColumn(bar_width=20, complete_style="cyan", finished_style="green"),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        console=console,
        transient=True,
    ) as progress:
        task = progress.add_task("Initializing...", total=2)

        # 1/2 — geNomad end-to-end
        progress.update(task, description="[1/2] geNomad  — end-to-end classification")
        _run_genomad(cfg, sample_id, contigs, sdir, db, rpt_dir, force)
        progress.advance(task)

        # 2/2 — parse results + copy output
        progress.update(task, description="[2/2] geNomad  — parsing & copying results")
        row = _parse_genomad_output(sample_id, s["n"], sdir, prefix, v_fna, v_tsv)
        if v_fna.exists() and v_fna.stat().st_size > 0:
            shutil.copy(v_fna, virus_fna)
            log_ok("  \\[geNomad] results parsed")
        else:
            log_warn(
                "  \\[geNomad] virus FASTA not found — "
                "check log in reports/04_viral_id/"
            )
        progress.advance(task)

    # ── Metrics, display, save ─────────────────────────────────────────────────
    _validate_output(sample_id, virus_fna, row)
    _save_tsv(row, rpt_dir / "genomad_summary.tsv")
    _print_summary_table(sample_id, row)
    _check_warnings(row)
    _print_completion_panel(sample_id, virus_fna, rpt_dir, row)

    log_step(f"Module 04 completed ✓  [{sample_id}]")
    log_ok(f"  Virus FASTA : {virus_fna}")
    log_info(
        f"  Next: phageflow quality --sample-id {sample_id} "
        f"--virus-fna {virus_fna}"
    )
    return virus_fna


# ── Step implementations ──────────────────────────────────────────────────────

def _run_genomad(
    cfg:       Config,
    sample_id: str,
    contigs:   Path,
    sdir:      Path,
    db:        Path,
    rpt_dir:   Path,
    force:     bool,
) -> None:
    """Run geNomad end-to-end classification."""
    prefix = contigs.stem
    v_fna  = sdir / f"{prefix}_summary" / f"{prefix}_virus.fna"

    if v_fna.exists() and v_fna.stat().st_size > 0 and not force:
        log_info("  \\[geNomad] already processed — skipping  (--force to re-run)")
        return

    mkdirs(sdir)
    run_silent([
        "genomad", "end-to-end",
        "--cleanup",
        "--splits",        "8",
        "--min-score",     str(cfg.genomad.min_score),
        "--min-virus-hallmarks", str(cfg.genomad.min_hallmarks),
        "--threads",       str(cfg.threads),
        str(contigs), str(sdir), str(db),
    ], log_file=rpt_dir / f"{sample_id}_genomad.log", check=False)

    # geNomad exits with code 2 in some valid runs; verify by output file
    if v_fna.exists() and v_fna.stat().st_size > 0:
        log_ok("  \\[geNomad] classification complete")
    else:
        log_warn(
            "  \\[geNomad] no virus FASTA produced — "
            f"check {rpt_dir / f'{sample_id}_genomad.log'}"
        )


# ── Validation ────────────────────────────────────────────────────────────────

def _validate_output(sample_id: str, virus_fna: Path, row: dict) -> None:
    """Warn if no viral contigs were identified, confirm OK otherwise."""
    if not virus_fna.exists() or virus_fna.stat().st_size == 0:
        log_warn(f"  validate · virus FASTA missing or empty: {virus_fna}")
        return

    n = row.get("viral_ctg", 0)
    if n == 0:
        log_warn(
            "  validate · 0 viral contigs identified — "
            "lower genomad.min_score in config.yaml or check assembly output."
        )
    else:
        log_ok(f"  validate · {n} viral contig(s) identified — OK")


# ── Parsers ───────────────────────────────────────────────────────────────────

def _parse_genomad_output(
    sample_id: str,
    n_input:   int,
    sdir:      Path,
    prefix:    str,
    v_fna:     Path,
    v_tsv:     Path,
) -> dict:
    """
    Parse geNomad output into a summary dict.

    Receives n_input from the already-computed fasta_stats to avoid
    redundant I/O on the (potentially large) contigs FASTA.
    """
    # Viral FASTA stats (single read, results reused below)
    v_stats  = fasta_stats(v_fna) if v_fna.exists() else {"n": 0, "total_bp": 0}
    n_viral  = v_stats["n"]
    viral_bp = v_stats["total_bp"]

    # Plasmid / provirus counts
    p_tsv  = sdir / f"{prefix}_summary"         / f"{prefix}_plasmid_summary.tsv"
    pr_tsv = sdir / f"{prefix}_find_proviruses"  / f"{prefix}_provirus.tsv"
    n_plas = _count_tsv_rows(p_tsv)
    n_prov = _count_tsv_rows(pr_tsv)

    # Best virus score + most specific taxonomy from virus_summary.tsv
    best_taxon = "unclassified"
    best_score = 0.0

    if v_tsv.exists():
        try:
            with open(v_tsv) as f:
                header = f.readline().strip().split("\t")
                col    = {h: i for i, h in enumerate(header)}
                for line in f:
                    row = line.strip().split("\t")
                    if not row:
                        continue
                    try:
                        score = float(row[col.get("virus_score", 6)])
                    except (ValueError, IndexError):
                        continue
                    if score > best_score:
                        best_score = score
                        # Extract the most specific non-empty taxonomy level
                        tax_idx = col.get("taxonomy", 10)
                        if len(row) > tax_idx:
                            levels = [
                                lvl.strip()
                                for lvl in row[tax_idx].split(";")
                                if lvl.strip() and lvl.strip().lower() != "unclassified"
                            ]
                            best_taxon = levels[-1] if levels else "unclassified"
        except Exception:
            pass

    return {
        "sample":      sample_id,
        "input_ctg":   n_input,
        "viral_ctg":   n_viral,
        "viral_bp":    viral_bp,
        "plasmid":     n_plas,
        "provirus":    n_prov,
        "best_taxon":  best_taxon,
        "best_score":  round(best_score, 4),
    }


def _count_tsv_rows(path: Path) -> int:
    """Count data rows in a TSV (excludes blank, comment, and header lines)."""
    if not path.exists():
        return 0
    try:
        with open(path) as f:
            return sum(
                1 for line in f
                if line.strip()
                and not line.startswith("#")
                and not line.startswith("seq_name")
            )
    except Exception:
        return 0


# ── Rich display helpers ──────────────────────────────────────────────────────

def _color_score(score: float) -> str:
    """Return Rich-colored string for a geNomad virus score."""
    s = f"{score:.3f}"
    if score >= _SCORE_GOOD: return f"[bold green]{s}[/bold green]"
    if score >= _SCORE_WARN: return f"[bold yellow]{s}[/bold yellow]"
    return f"[bold red]{s}[/bold red]"


def _color_rate(val_str: str, good: float, warn: float) -> str:
    """Return Rich-colored string based on numeric thresholds."""
    try:
        v = float(str(val_str).rstrip("%"))
        if v >= good: return f"[bold green]{val_str}[/bold green]"
        if v >= warn: return f"[bold yellow]{val_str}[/bold yellow]"
        return f"[bold red]{val_str}[/bold red]"
    except (ValueError, AttributeError):
        return str(val_str)


def _fmt(val, unit: str = "") -> str:
    """Format metric value; dim 'N/A' or missing entries."""
    if val in ("N/A", "", None, 0, "0"):
        return f"[dim]{val}{unit}[/dim]" if val in (0, "0") else "[dim]N/A[/dim]"
    return f"{val}{unit}"


def _print_input_table(sample_id: str, contigs: Path, stats: dict) -> None:
    log_info(f"  Contigs : {contigs}  ({human_size(contigs)})")
    log_info(
        f"  Input   : {stats['n']} sequence(s) | "
        f"total={stats['total_bp']:,}bp | largest={stats['largest_bp']:,}bp"
    )


def _print_summary_table(sample_id: str, row: dict) -> None:
    """Print geNomad identification metrics as compact log lines."""
    n_in  = row.get("input_ctg",  0)
    n_vir = row.get("viral_ctg",  0)
    n_bp  = row.get("viral_bp",   0)
    n_pls = row.get("plasmid",    0)
    n_prv = row.get("provirus",   0)
    taxon = row.get("best_taxon", "unclassified")
    score = row.get("best_score", 0.0)

    pct_vir = f"{n_vir / n_in * 100:.1f}%" if n_in else "N/A"
    colored_pct = _color_rate(pct_vir, 50.0, 5.0) if n_in else "[dim]N/A[/dim]"

    log_ok(
        f"  Contigs  : {n_in} input → {n_vir} viral "
        f"({colored_pct})  |  {n_bp:,}bp total viral sequence"
    )
    log_ok(
        f"  Other    : plasmid={_fmt(n_pls)}  provirus={_fmt(n_prv)}"
    )
    log_ok(
        f"  Best hit : {taxon}  |  score={_color_score(score)}"
    )


def _check_warnings(row: dict) -> None:
    """Emit contextual warnings based on viral identification metrics."""
    n_in  = row.get("input_ctg", 0)
    n_vir = row.get("viral_ctg", 0)

    if n_vir == 0:
        log_warn(
            "  No viral contigs identified — "
            "try lowering genomad.min_score in config.yaml (current default: 0.7)."
        )
        return

    if n_in and (n_vir / n_in) < _VIRAL_FRACTION_WARN:
        log_warn(
            f"  Only {n_vir}/{n_in} contigs classified as viral "
            f"({n_vir / n_in * 100:.1f}%) — "
            "verify host removal step or check for contamination."
        )

    if row.get("provirus", 0) > 0:
        log_warn(
            f"  {row['provirus']} integrated provirus(es) detected — "
            "review lifecycle prediction carefully (Module 08)."
        )

    if 0 < row.get("best_score", 0.0) < _SCORE_WARN:
        log_warn(
            f"  Best virus score {row['best_score']:.3f} < {_SCORE_WARN} — "
            "borderline classification; consider manual inspection of contigs."
        )


def _print_completion_panel(
    sample_id: str, virus_fna: Path, rpt_dir: Path, row: dict
) -> None:
    n_in  = row.get("input_ctg",   "?")
    n_vir = row.get("viral_ctg",   "?")
    taxon = row.get("best_taxon",  "unclassified")
    score = row.get("best_score",  0.0)
    n_pls = row.get("plasmid",     0)
    n_prv = row.get("provirus",    0)

    text = Text()
    text.append("✓ ", style="bold green")
    text.append(f"{n_in} contigs  →  ", style="dim white")
    text.append(f"{n_vir}", style="bold green")
    text.append(f" viral  (plasmid={n_pls}  provirus={n_prv})\n", style="cyan")
    text.append(f"  best hit : {taxon}  (score={score:.3f})\n\n", style="dim white")
    text.append("Virus FASTA : ", style="dim white")
    text.append(str(virus_fna) + "\n", style="white")
    text.append("Summary     : ", style="dim white")
    text.append(str(rpt_dir / "genomad_summary.tsv"), style="white")

    console.print(Panel(
        text,
        title=f"[bold cyan]Viral-id complete — {sample_id}[/bold cyan]",
        border_style="cyan",
        padding=(0, 2),
        width=90,
    ))


# ── TSV summary ───────────────────────────────────────────────────────────────

_TSV_HEADERS = [
    "sample", "input_ctg", "viral_ctg", "viral_bp",
    "plasmid", "provirus", "best_taxon", "best_score",
]


def _save_tsv(row: dict, path: Path) -> None:
    """
    Write / update the geNomad summary TSV.

    Existing rows are preserved and migrated to the current schema
    (missing columns filled with empty strings), so re-running with
    an older TSV does not lose previous samples.
    """
    rows: dict[str, dict] = {}

    if path.exists():
        with open(path) as f:
            old_hdrs = f.readline().rstrip("\n").split("\t")
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if cols and cols[0]:
                    rows[cols[0]] = dict(zip(old_hdrs, cols))

    rows[row["sample"]] = {h: str(row.get(h, "")) for h in _TSV_HEADERS}

    with open(path, "w") as f:
        f.write("\t".join(_TSV_HEADERS) + "\n")
        for r in rows.values():
            f.write("\t".join(r.get(h, "") for h in _TSV_HEADERS) + "\n")

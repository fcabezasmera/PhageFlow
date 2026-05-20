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
    Use: conda run -n genomad phageflow viral-id ...
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
    _print_input_table(contigs, s)
    log_info(f"  geNomad min-score : {cfg.genomad.min_score}  (Camargo et al. 2023)")
    log_info("  Classifying       : virus / plasmid / chromosome + provirus detection")

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
        progress.update(task, description="[2/2] geNomad  — parsing results")
        row = _parse_genomad_output(sample_id, contigs, sdir, prefix, v_fna, v_tsv)
        if v_fna.exists() and v_fna.stat().st_size > 0:
            shutil.copy(v_fna, virus_fna)
            log_ok("  \\[geNomad] results parsed")
        else:
            log_warn("  \\[geNomad] virus FASTA not found — check log in reports/04_viral_id/")
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
        f"--virus-fna {virus_fna}  (activate phageflow env)"
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
        log_info("  \\[geNomad] already processed — skipping")
        return

    mkdirs(sdir)
    try:
        run_silent([
            "genomad", "end-to-end",
            "--cleanup",
            "--splits",    "8",
            "--min-score", str(cfg.genomad.min_score),
            "--threads",   str(cfg.threads),
            str(contigs), str(sdir), str(db),
        ], log_file=rpt_dir / f"{sample_id}_genomad.log")
        log_ok("  \\[geNomad] classification complete")
    except Exception as e:
        log_warn(f"  \\[geNomad] warning: {e}")


# ── Validation ────────────────────────────────────────────────────────────────

def _validate_output(sample_id: str, virus_fna: Path, row: dict) -> None:
    """Warn if no viral contigs were identified."""
    if not virus_fna.exists() or virus_fna.stat().st_size == 0:
        log_warn(f"  validate · virus FASTA missing or empty: {virus_fna}")
        return

    n = row.get("viral_ctg", 0)
    if n == 0:
        log_warn(
            f"  validate · 0 viral contigs identified — "
            "check geNomad log or lower --min-score in config.yaml"
        )
    else:
        log_ok(f"  validate · {n} viral contig(s) identified — OK")


# ── Parsers ───────────────────────────────────────────────────────────────────

def _parse_genomad_output(
    sample_id: str,
    fa:        Path,
    sdir:      Path,
    prefix:    str,
    v_fna:     Path,
    v_tsv:     Path,
) -> dict:
    """Parse geNomad output into a summary dict."""
    n_input  = fasta_stats(fa)["n"]
    n_viral  = fasta_stats(v_fna)["n"]        if v_fna.exists() else 0
    viral_bp = fasta_stats(v_fna)["total_bp"] if v_fna.exists() else 0

    p_tsv  = sdir / f"{prefix}_summary"          / f"{prefix}_plasmid_summary.tsv"
    pr_tsv = sdir / f"{prefix}_find_proviruses"  / f"{prefix}_provirus.tsv"

    n_plas = _count_tsv_rows(p_tsv)
    n_prov = _count_tsv_rows(pr_tsv)

    best_fam  = "unclassified"
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
                        if score > best_score:
                            best_score = score
                            tax   = row[col.get("taxonomy", 10)] if len(row) > col.get("taxonomy", 10) else ""
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
        "best_score":  round(best_score, 4),
    }


def _count_tsv_rows(path: Path) -> int:
    """Count data rows in a TSV (excludes blank, comment, and header lines)."""
    if not path.exists():
        return 0
    try:
        return sum(
            1 for line in open(path)
            if line.strip()
            and not line.startswith("#")
            and not line.startswith("seq_name")
        )
    except Exception:
        return 0


# ── Rich display helpers ──────────────────────────────────────────────────────

def _fmt(val, unit: str = "") -> str:
    """Format metric value; dim 'N/A' or missing entries."""
    if val in ("N/A", "", None):
        return "[dim]N/A[/dim]"
    return f"{val}{unit}"


def _print_input_table(contigs: Path, stats: dict) -> None:
    log_info(f"  Contigs : {contigs}  ({human_size(contigs)})")
    log_info(
        f"  Input   : {stats['n']} sequences | "
        f"total={stats['total_bp']}bp | largest={stats['largest_bp']}bp"
    )


def _print_summary_table(sample_id: str, row: dict) -> None:
    """Print geNomad identification metrics as compact log lines."""
    n_in  = row.get("input_ctg",  0)
    n_vir = row.get("viral_ctg",  0)
    n_bp  = row.get("viral_bp",   0)
    n_pls = row.get("plasmid",    0)
    n_prv = row.get("provirus",   0)
    fam   = row.get("best_family", "unclassified")
    score = row.get("best_score",  0.0)

    pct_vir = f"{n_vir / n_in * 100:.1f}%" if n_in else "N/A"

    log_ok(
        f"  Viral    : {_fmt(n_vir)} contig(s)  ({_fmt(n_bp, 'bp')})  "
        f"· {pct_vir} of input"
    )
    log_ok(
        f"  Plasmid  : {_fmt(n_pls)}  |  "
        f"Provirus : {_fmt(n_prv)}  |  "
        f"Best family : {_fmt(fam)}  (score={score:.3f})"
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

    if n_in and n_vir / n_in < _VIRAL_FRACTION_WARN:
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


def _print_completion_panel(
    sample_id: str, virus_fna: Path, rpt_dir: Path, row: dict
) -> None:
    n_in  = row.get("input_ctg",   "?")
    n_vir = row.get("viral_ctg",   "?")
    fam   = row.get("best_family", "unclassified")

    text = Text()
    text.append("✓ ", style="bold green")
    text.append(f"{n_in} contigs  →  ", style="dim white")
    text.append(f"{n_vir}", style="bold green")
    text.append(f" viral  |  best family: {fam}\n\n", style="cyan")
    text.append("Virus FASTA : ", style="dim white")
    text.append(str(virus_fna) + "\n", style="white")
    text.append("Summary     : ", style="dim white")
    text.append(str(rpt_dir / "genomad_summary.tsv"), style="white")

    console.print(Panel(
        text,
        title=f"[bold cyan]Viral-id complete — {sample_id}[/bold cyan]",
        border_style="cyan",
        padding=(0, 2),
        width=72,
    ))


# ── TSV summary ───────────────────────────────────────────────────────────────

def _save_tsv(row: dict, path: Path) -> None:
    """Write / update the geNomad summary TSV."""
    headers = ["sample", "input_ctg", "viral_ctg", "viral_bp",
               "plasmid", "provirus", "best_family"]
    rows: dict[str, list] = {}

    if path.exists():
        with open(path) as f:
            old_hdrs = f.readline().rstrip("\n").split("\t")
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if cols and cols[0]:
                    rows[cols[0]] = dict(zip(old_hdrs, cols))

    rows[row["sample"]] = [str(row.get(h, "")) for h in headers]

    with open(path, "w") as f:
        f.write("\t".join(headers) + "\n")
        for r in rows.values():
            f.write("\t".join(r if isinstance(r, list) else [r.get(h, "") for h in headers]) + "\n")

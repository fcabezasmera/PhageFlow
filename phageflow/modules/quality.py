"""PhageFlow Module 05 — Genome quality assessment and selection.

Tool:
    CheckV v1.0+ (Nayfach et al. 2021, Nature Biotechnology):
        Estimates completeness, contamination, and topology (DTR/ITR)
        for viral genomes using a database of complete viral genomes.

        end_to_end   : runs all CheckV steps in one command.
        --remove_tmp : removes large intermediate files after run.

Quality tiers (MIUVIG standard):
    Complete / High-quality  ≥ min_completeness  → annotation_ready/
    Medium-quality           < min_completeness  → drafts/
    Low-quality / Not-determined                 → discarded

Input  : results/04_viral_id/{sample}_virus.fna
Output : results/05_quality/annotation_ready/{sample}_HQ.fasta  ← annotation input
         results/05_quality/drafts/{sample}_draft.fasta
         reports/05_quality/checkv_summary.tsv
"""

from __future__ import annotations
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

STEP  = "05_quality"
TOOLS = ["checkv"]

# ── Thresholds ────────────────────────────────────────────────────────────────
_COMPLETENESS_GOOD   = 90.0   # completeness: green
_COMPLETENESS_WARN   = 70.0   # completeness: yellow; below → worth a warning


# ── Public entry point ────────────────────────────────────────────────────────

def run(
    cfg:       Config,
    sample_id: str,
    virus_fna: Path,
    force:     bool = False,
) -> Path:
    """
    Assess genome quality with CheckV and select HQ genomes.

    Parameters
    ----------
    sample_id : sample identifier
    virus_fna : viral contigs FASTA (output of viral-id step)

    Returns
    -------
    hq_fasta : Path to annotation-ready HQ genome FASTA
               (empty/missing if no HQ genomes found)
    """
    require_tools(*TOOLS)

    virus_fna = Path(virus_fna)
    out_dir   = cfg.results(STEP)
    rpt_dir   = cfg.reports(STEP)
    final_dir = out_dir / "annotation_ready"
    draft_dir = out_dir / "drafts"
    mkdirs(out_dir, rpt_dir, final_dir, draft_dir)

    hq_fasta = final_dir / f"{sample_id}_HQ.fasta"

    log_step(f"Module 05 — quality [{sample_id}]")

    if not virus_fna.exists():
        log_error(f"  Virus FASTA not found: {virus_fna}")
        raise FileNotFoundError(virus_fna)

    s = fasta_stats(virus_fna)
    _print_input_table(sample_id, virus_fna, s)
    log_info(
        f"  CheckV end-to-end  (Nayfach et al. 2021)  |  "
        f"HQ threshold: completeness ≥ {cfg.checkv.min_completeness}%"
    )

    if not cfg.databases.checkv.exists():
        log_warn(f"  CheckV database not found: {cfg.databases.checkv}")

    # ── Pipeline steps ────────────────────────────────────────────────────────
    sdir      = out_dir / sample_id
    qfile     = sdir / "quality_summary.tsv"
    row       = {}
    hq_seqs   = []
    draft_seqs = []

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

        # 1/2 — CheckV end-to-end
        progress.update(task, description="[1/2] CheckV   — end-to-end assessment")
        _run_checkv(cfg, sample_id, virus_fna, sdir, rpt_dir, force)
        progress.advance(task)

        # 2/2 — Parse and select genomes
        progress.update(task, description="[2/2] CheckV   — selecting HQ genomes")
        if qfile.exists():
            seqs = _load_fasta(virus_fna)
            row, hq_seqs, draft_seqs = _parse_checkv(
                sample_id, qfile, seqs, cfg.checkv.min_completeness
            )
            _write_fasta(hq_seqs,    hq_fasta)
            _write_fasta(draft_seqs, draft_dir / f"{sample_id}_draft.fasta")
            log_ok("  \\[CheckV] genome selection complete")
        else:
            log_warn(
                "  \\[CheckV] no output found — "
                f"check reports/05_quality/{sample_id}_checkv.log"
            )
        progress.advance(task)

    # ── Metrics, display, save ─────────────────────────────────────────────────
    _validate_output(sample_id, hq_fasta, hq_seqs, draft_seqs)
    if row:
        _save_tsv(row, rpt_dir / "checkv_summary.tsv")
        _print_summary_table(sample_id, row)
        _check_warnings(row, cfg.checkv.min_completeness)
    _print_completion_panel(sample_id, hq_fasta, draft_dir, rpt_dir, row, hq_seqs, draft_seqs)

    log_step(f"Module 05 completed ✓  [{sample_id}]")
    log_info("Next: phageflow annotate --sample-id <id> --genome <HQ.fasta>")

    return hq_fasta


# ── Step implementations ──────────────────────────────────────────────────────

def _run_checkv(
    cfg:       Config,
    sample_id: str,
    virus_fna: Path,
    sdir:      Path,
    rpt_dir:   Path,
    force:     bool,
) -> None:
    """Run CheckV end-to-end assessment."""
    qfile = sdir / "quality_summary.tsv"

    if qfile.exists() and qfile.stat().st_size > 0 and not force:
        log_info("  \\[CheckV] already processed — skipping  (--force to re-run)")
        return

    mkdirs(sdir)
    try:
        run_silent([
            "checkv", "end_to_end",
            str(virus_fna), str(sdir),
            "-d", str(cfg.databases.checkv),
            "-t", str(cfg.threads),
            "--remove_tmp",
        ], log_file=rpt_dir / f"{sample_id}_checkv.log")
        log_ok("  \\[CheckV] assessment complete")
    except Exception as e:
        log_warn(f"  \\[CheckV] warning: {e}")


# ── Validation ────────────────────────────────────────────────────────────────

def _validate_output(
    sample_id:  str,
    hq_fasta:   Path,
    hq_seqs:    list,
    draft_seqs: list,
) -> None:
    """Warn if no HQ genomes were found; confirm OK otherwise."""
    if hq_seqs:
        log_ok(f"  validate · {len(hq_seqs)} HQ genome(s) annotation-ready — OK")
    elif draft_seqs:
        log_warn(
            f"  validate · 0 HQ genomes — {len(draft_seqs)} draft(s) in drafts/. "
            "Lower checkv.min_completeness in config.yaml to include them."
        )
    else:
        log_warn(
            f"  validate · 0 genomes passed quality filters for {sample_id}. "
            "Check geNomad viral classification (Module 04)."
        )


# ── Helpers ───────────────────────────────────────────────────────────────────

def _load_fasta(path: Path) -> dict:
    """Load FASTA into {seq_id: sequence} dict."""
    seqs: dict[str, str] = {}
    hdr = seq = ""
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if hdr:
                    seqs[hdr] = seq
                hdr = line[1:].split()[0]
                seq = ""
            else:
                seq += line
    if hdr:
        seqs[hdr] = seq
    return seqs


def _parse_checkv(
    sample_id:        str,
    qfile:            Path,
    seqs:             dict,
    min_completeness: float,
) -> tuple:
    """
    Parse CheckV quality_summary.tsv and split sequences into HQ / draft tiers.

    Returns (summary_dict, hq_seqs, draft_seqs).
    Stores best_completeness as float for reliable comparison downstream.
    """
    complete = hq = mq = lq = 0
    best_len          = 0
    best_completeness = 0.0
    best_quality_tier = ""
    hq_seqs:    list = []
    draft_seqs: list = []

    with open(qfile) as f:
        header = f.readline().strip().split("\t")
        col    = {h: i for i, h in enumerate(header)}
        for line in f:
            row = line.strip().split("\t")
            if not row or len(row) < 3:
                continue

            cid  = row[col.get("contig_id",     0)]
            clen = row[col.get("contig_length",  1)]
            qual = row[col.get("checkv_quality", 7)]
            comp = row[col.get("completeness",   9)]
            seq  = seqs.get(cid, "")

            if qual == "Complete":         complete += 1
            elif qual == "High-quality":   hq       += 1
            elif qual == "Medium-quality": mq       += 1
            else:                          lq       += 1

            # Track best genome by length
            try:
                if int(clen) > best_len:
                    best_len          = int(clen)
                    best_completeness = float(comp) if comp not in ("", "NA", "N/A") else 0.0
                    best_quality_tier = qual
            except (ValueError, TypeError):
                pass

            if seq:
                if qual in ("Complete", "High-quality"):
                    hq_seqs.append((cid, seq, clen, comp, qual))
                elif qual == "Medium-quality":
                    draft_seqs.append((cid, seq, clen, comp, qual))

    summary = {
        "sample":             sample_id,
        "complete":           complete,
        "hq":                 hq,
        "mq":                 mq,
        "lq":                 lq,
        "best_length_bp":     best_len,
        "best_completeness":  best_completeness,
        "best_quality_tier":  best_quality_tier,
    }
    return summary, hq_seqs, draft_seqs


def _write_fasta(seqs_list: list, out_path: Path) -> None:
    if not seqs_list:
        return
    with open(out_path, "w") as f:
        for cid, seq, *_ in seqs_list:
            f.write(f">{cid}\n{seq}\n")


# ── Rich display helpers ──────────────────────────────────────────────────────

def _color_completeness(val: float) -> str:
    """Return Rich-colored string for a completeness percentage."""
    s = f"{val:.1f}%"
    if val >= _COMPLETENESS_GOOD: return f"[bold green]{s}[/bold green]"
    if val >= _COMPLETENESS_WARN: return f"[bold yellow]{s}[/bold yellow]"
    return f"[bold red]{s}[/bold red]"


def _fmt(val, unit: str = "") -> str:
    """Format metric value; dim missing/zero entries."""
    if val in ("N/A", "", None):
        return "[dim]N/A[/dim]"
    return f"{val}{unit}"


def _print_input_table(sample_id: str, virus_fna: Path, stats: dict) -> None:
    log_info(f"  Virus FASTA : {virus_fna}  ({human_size(virus_fna)})")
    log_info(
        f"  Input       : {stats['n']} contig(s) | "
        f"total={stats['total_bp']:,}bp | largest={stats['largest_bp']:,}bp"
    )


def _print_summary_table(sample_id: str, row: dict) -> None:
    """Print CheckV quality tier counts as compact log lines."""
    complete = row.get("complete", 0)
    hq       = row.get("hq",       0)
    mq       = row.get("mq",       0)
    lq       = row.get("lq",       0)
    blen     = row.get("best_length_bp",    0)
    bcomp    = row.get("best_completeness", 0.0)
    btier    = row.get("best_quality_tier", "")

    log_ok(
        f"  Tiers    : Complete={complete}  High-quality={hq}  "
        f"Medium-quality={mq}  Low/ND={lq}"
    )
    log_ok(
        f"  Best     : {blen:,}bp  |  completeness={_color_completeness(bcomp)}  "
        f"|  {_fmt(btier)}"
    )


def _check_warnings(row: dict, min_completeness: float) -> None:
    """Emit contextual warnings based on quality assessment results."""
    total_hq = row.get("complete", 0) + row.get("hq", 0)

    if total_hq == 0 and row.get("mq", 0) > 0:
        log_warn(
            f"  No Complete/HQ genomes — {row['mq']} Medium-quality draft(s) saved. "
            f"Lower checkv.min_completeness (currently {min_completeness}%) to annotate them."
        )
    elif total_hq == 0:
        log_warn(
            "  No genomes passed quality filters — "
            "check geNomad output or assembly contiguity."
        )

    bcomp = row.get("best_completeness", 0.0)
    if 0 < bcomp < _COMPLETENESS_WARN:
        log_warn(
            f"  Best completeness {bcomp:.1f}% < {_COMPLETENESS_WARN}% — "
            "fragmented genome; annotation and lifecycle prediction may be affected."
        )


def _print_completion_panel(
    sample_id:  str,
    hq_fasta:   Path,
    draft_dir:  Path,
    rpt_dir:    Path,
    row:        dict,
    hq_seqs:    list,
    draft_seqs: list,
) -> None:
    complete = row.get("complete", 0)
    hq       = row.get("hq",       0)
    mq       = row.get("mq",       0)
    lq       = row.get("lq",       0)
    total_hq = complete + hq

    text = Text()
    text.append("✓ ", style="bold green")
    text.append(
        f"Complete={complete}  HQ={hq}  MQ={mq}  LQ/ND={lq}  →  ",
        style="dim white",
    )
    text.append(f"{total_hq}", style="bold green")
    text.append(" genome(s) annotation-ready\n\n", style="cyan")

    if hq_seqs:
        text.append("HQ FASTA  : ", style="dim white")
        text.append(str(hq_fasta) + "\n", style="white")
    if draft_seqs:
        text.append("Drafts    : ", style="dim white")
        text.append(str(draft_dir / f"{sample_id}_draft.fasta") + "\n", style="white")
    text.append("Summary   : ", style="dim white")
    text.append(str(rpt_dir / "checkv_summary.tsv"), style="white")

    console.print(Panel(
        text,
        title=f"[bold cyan]Quality complete — {sample_id}[/bold cyan]",
        border_style="cyan",
        padding=(0, 2),
        width=90,
    ))


# ── TSV summary ───────────────────────────────────────────────────────────────

_TSV_HEADERS = [
    "sample", "complete", "hq", "mq", "lq",
    "best_length_bp", "best_completeness", "best_quality_tier",
]


def _save_tsv(row: dict, path: Path) -> None:
    """
    Write / update the CheckV summary TSV.

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

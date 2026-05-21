"""PhageFlow Module 04 — Viral identification with geNomad.

Tool:
    geNomad v1.8+ (Camargo et al. 2023, Nature Biotechnology):
        End-to-end classification of contigs as virus / plasmid / chromosome.
        Detects integrated proviruses and assigns viral taxonomy.

        --splits N         : Adaptive: min(8, n_contigs). Camargo et al. 2023.
        --min-score        : configurable (default 0.7, ~97% precision).
        --min-hallmarks    : configurable (default 0, score-only mode).
        --cleanup          : removes annotate/ intermediates; keeps summary/
                             including aggregated_classification.tsv (v1.8+).

Borderline large-contig rescue (purified_phage mode):
    geNomad scores are calibrated on database sequences. Novel phage lineages
    with no close database relatives consistently receive scores of 0.4-0.6
    despite being genuine phage. For purified preparations, any contig
    >= rescue_min_length_bp with score in [rescue_min_score, min_score) is
    almost certainly phage and is added to the viral FASTA with a WARNING.

    Parameters (config.yaml):
        genomad.rescue_min_score:     0.4  (floor — below this = discard)
        genomad.rescue_min_length_bp: 10000 (bp)

    Mechanism: reads {prefix}_aggregated_classification.tsv from the
    summary/ directory (kept by --cleanup in geNomad v1.8+). If the file
    is absent (older geNomad or manual deletion), rescue is skipped with
    an informative message — not an error.

Output path convention:
    geNomad names outputs using the STEM of the input FASTA (contigs.stem).
    The resolved taxonomy TSV path is written to genomad_summary.tsv so
    quality.py can locate it regardless of input filename.

Input  : results/03_assembly/combined/{sample}_contigs_nr.fasta
Output : results/04_viral_id/{sample}_virus.fna
         reports/04_viral_id/genomad_summary.tsv
"""

from __future__ import annotations
import shutil
import subprocess
from pathlib import Path

from rich.panel import Panel
from rich.text import Text
from rich.progress import (
    Progress, SpinnerColumn, TextColumn,
    BarColumn, MofNCompleteColumn, TimeElapsedColumn,
)

from phageflow.utils.config import Config
from phageflow.utils.logger import log_step, log_info, log_ok, log_warn, console
from phageflow.utils.tools import (
    require_tools, run_silent, mkdirs, fasta_stats, human_size,
)

STEP  = "04_viral_id"
TOOLS = ["genomad"]

_VIRAL_FRACTION_WARN = 0.05
_SCORE_GOOD          = 0.90
_SCORE_WARN          = 0.70
_SPLITS_MAX          = 8
_COL_SCORE           = "virus_score"
_COL_TAXONOMY        = "taxonomy"


# -- Public entry point --------------------------------------------------------

def run(
    cfg:       Config,
    sample_id: str,
    contigs:   Path,
    force:     bool = False,
) -> Path:
    """
    Identify viral contigs with geNomad + borderline large-contig rescue.

    Returns
    -------
    virus_fna : Path to viral contigs FASTA (input for CheckV / quality step)
    """
    require_tools(*TOOLS)

    contigs = Path(contigs)
    out_dir = cfg.results(STEP)
    rpt_dir = cfg.reports(STEP)
    mkdirs(out_dir, rpt_dir)

    if not cfg.databases.genomad.exists():
        log_warn(f"  geNomad database not found: {cfg.databases.genomad}")

    virus_fna = out_dir / f"{sample_id}_virus.fna"

    log_step(f"Module 04 -- viral-id [{sample_id}]")

    if not contigs.exists():
        raise FileNotFoundError(f"[{sample_id}] contigs not found: {contigs}")
    if contigs.stat().st_size == 0:
        raise ValueError(f"[{sample_id}] contigs FASTA is empty: {contigs}")

    s = fasta_stats(contigs)
    if s["n"] == 0:
        raise ValueError(f"[{sample_id}] no sequences in contigs FASTA: {contigs}")

    _check_duplicate_headers(contigs, sample_id)

    n_splits = min(_SPLITS_MAX, max(1, s["n"]))
    prefix   = contigs.stem

    _print_input_table(sample_id, contigs, s)
    log_info(
        f"  geNomad end-to-end  (Camargo et al. 2023)  |  "
        f"min-score={cfg.genomad.min_score}  "
        f"min-virus-hallmarks={cfg.genomad.min_hallmarks}  "
        f"splits={n_splits}/{_SPLITS_MAX}"
    )
    log_info(
        f"  Borderline rescue: score>={cfg.genomad.rescue_min_score} "
        f"AND length>={cfg.genomad.rescue_min_length_bp:,}bp "
        f"(purified_phage mode)"
    )
    log_info(f"  Output prefix: '{prefix}'")
    log_info("  Classifying : virus / plasmid / chromosome + provirus")

    sdir  = out_dir / sample_id
    v_fna = sdir / f"{prefix}_summary" / f"{prefix}_virus.fna"
    v_tsv = sdir / f"{prefix}_summary" / f"{prefix}_virus_summary.tsv"

    with Progress(
        SpinnerColumn(spinner_name="dots2"),
        TextColumn("[bold cyan]{task.description:<52}"),
        BarColumn(bar_width=20, complete_style="cyan", finished_style="green"),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        console=console,
        transient=True,
    ) as progress:
        task = progress.add_task("Initializing...", total=3)

        # 1/3 -- geNomad
        progress.update(task, description="[1/3] geNomad  -- end-to-end classification")
        _run_genomad(cfg, sample_id, contigs, sdir, rpt_dir, prefix, n_splits, force)
        progress.advance(task)

        # 2/3 -- Borderline large-contig rescue
        progress.update(task, description="[2/3] geNomad  -- borderline rescue")
        n_rescued = 0
        if v_fna.exists() and v_fna.stat().st_size > 0:
            n_rescued = _rescue_borderline_contigs(
                cfg, contigs, sdir, prefix, v_fna, rpt_dir,
            )
        progress.advance(task)

        # 3/3 -- Parse + copy
        progress.update(task, description="[3/3] geNomad  -- parsing & copying results")
        row = _parse_genomad_output(
            cfg, sample_id, s["n"], sdir, prefix, v_fna, v_tsv, n_rescued,
        )
        if v_fna.exists() and v_fna.stat().st_size > 0:
            try:
                shutil.copy(v_fna, virus_fna)
            except OSError as e:
                log_warn(f"  [geNomad] copy failed: {e} -- using {v_fna}")
                virus_fna = v_fna
            log_ok("  [geNomad] results parsed")
        else:
            _diagnose_genomad_failure(cfg, sdir, prefix, rpt_dir, sample_id, s["n"])
        progress.advance(task)

    _validate_output(sample_id, virus_fna, row)
    _save_tsv(row, rpt_dir / "genomad_summary.tsv")
    _print_summary_table(sample_id, row)
    _check_warnings(row)
    _print_completion_panel(sample_id, virus_fna, rpt_dir, row)

    log_step(f"Module 04 completed v  [{sample_id}]")
    log_ok(f"  Virus FASTA : {virus_fna}")
    log_info(
        f"  Next: phageflow quality --sample-id {sample_id} "
        f"--virus-fna {virus_fna}"
    )
    return virus_fna


# -- Input validation ----------------------------------------------------------

def _check_duplicate_headers(contigs: Path, sample_id: str) -> None:
    seen: set[str] = set()
    duplicates: list[str] = []
    with open(contigs) as f:
        for line in f:
            if line.startswith(">"):
                hdr = line[1:].split()[0].strip()
                if hdr in seen:
                    duplicates.append(hdr)
                else:
                    seen.add(hdr)
    if duplicates:
        log_warn(
            f"  [validate] {len(duplicates)} duplicate FASTA header(s) in "
            f"{contigs.name} -- geNomad may fail. "
            f"First: '{duplicates[0]}'. Re-run Module 03 with --force."
        )


# -- Step 1: geNomad -----------------------------------------------------------

def _run_genomad(
    cfg:       Config,
    sample_id: str,
    contigs:   Path,
    sdir:      Path,
    rpt_dir:   Path,
    prefix:    str,
    n_splits:  int,
    force:     bool,
) -> None:
    v_fna = sdir / f"{prefix}_summary" / f"{prefix}_virus.fna"
    if v_fna.exists() and v_fna.stat().st_size > 0 and not force:
        log_info("  [geNomad] already processed -- skipping  (--force to re-run)")
        return
    mkdirs(sdir)
    run_silent([
        "genomad", "end-to-end",
        "--cleanup",
        "--splits",              str(n_splits),
        "--min-score",           str(cfg.genomad.min_score),
        "--min-virus-hallmarks", str(cfg.genomad.min_hallmarks),
        "--threads",             str(cfg.threads),
        str(contigs), str(sdir), str(cfg.databases.genomad),
    ], log_file=rpt_dir / f"{sample_id}_genomad.log", check=False)

    if v_fna.exists() and v_fna.stat().st_size > 0:
        log_ok("  [geNomad] classification complete")


# -- Step 2: Borderline large-contig rescue ------------------------------------

def _find_aggregated_tsv(sdir: Path, prefix: str) -> Path | None:
    """
    Locate aggregated_classification.tsv (all sequences with scores).

    geNomad v1.8+ places this in {prefix}_summary/ which is kept by --cleanup.
    Older versions may place it directly in sdir.
    """
    candidates = [
        sdir / f"{prefix}_summary" / f"{prefix}_aggregated_classification.tsv",
        sdir / f"{prefix}_aggregated_classification.tsv",
    ]
    found = next((p for p in candidates if p.exists()), None)
    if found is None:
        hits = list(sdir.rglob("*aggregated_classification.tsv"))
        found = hits[0] if hits else None
    return found


def _load_fasta_local(path: Path) -> dict[str, str]:
    """Load FASTA sequences keyed by first-word header ID."""
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


def _rescue_borderline_contigs(
    cfg:       Config,
    contigs:   Path,
    sdir:      Path,
    prefix:    str,
    virus_fna: Path,
    rpt_dir:   Path,
) -> int:
    """
    Append large borderline-scored contigs to the viral FASTA.

    Contigs with score in [rescue_min_score, min_score) AND length >=
    rescue_min_length_bp are added with explicit WARNING. In purified phage
    preparations these are almost certainly novel phage lineages with no
    close database relatives, which systematically receive lower geNomad
    scores (Camargo et al. 2023).

    Source: {prefix}_aggregated_classification.tsv (all contigs with scores,
    kept by --cleanup in geNomad v1.8+). If absent, rescue is skipped.

    Returns number of rescued contigs.
    """
    rescue_min_score = cfg.genomad.rescue_min_score
    rescue_min_len   = cfg.genomad.rescue_min_length_bp
    min_score        = cfg.genomad.min_score

    agg_tsv = _find_aggregated_tsv(sdir, prefix)
    if not agg_tsv:
        log_info(
            "  [borderline rescue] aggregated_classification.tsv not found -- "
            "rescue skipped. geNomad v1.8+ keeps this in "
            f"{prefix}_summary/. Re-run with --force if needed."
        )
        return 0

    # IDs already included in virus FASTA (passed min_score threshold)
    existing_viral: set[str] = set(_load_fasta_local(virus_fna).keys())

    # All sequences from the input contigs FASTA (for retrieving sequences)
    all_seqs = _load_fasta_local(contigs)

    candidates: list[tuple[str, int, float]] = []
    try:
        with open(agg_tsv) as f:
            header = f.readline().strip().split("\t")
            col    = {h: i for i, h in enumerate(header)}

            id_idx  = col.get("seq_name", 0)
            len_idx = col.get("length", 1)
            sc_idx  = col.get(_COL_SCORE, -1)

            if sc_idx == -1:
                log_warn(
                    f"  [borderline rescue] column '{_COL_SCORE}' not found in "
                    f"{agg_tsv.name} -- rescue skipped."
                )
                return 0

            for line in f:
                parts = line.strip().split("\t")
                if not parts or len(parts) <= max(id_idx, len_idx, sc_idx):
                    continue
                seq_id = parts[id_idx].split("|")[0].strip()
                try:
                    length = int(parts[len_idx])
                    score  = float(parts[sc_idx])
                except (ValueError, IndexError):
                    continue

                if (
                    seq_id not in existing_viral
                    and rescue_min_score <= score < min_score
                    and length >= rescue_min_len
                ):
                    candidates.append((seq_id, length, score))

    except Exception as e:
        log_warn(f"  [borderline rescue] could not parse {agg_tsv.name}: {e}")
        return 0

    if not candidates:
        log_info(
            f"  [borderline rescue] 0 candidates found "
            f"(score [{rescue_min_score}, {min_score}), "
            f"length >={rescue_min_len:,}bp)"
        )
        return 0

    n_added = 0
    with open(virus_fna, "a") as f_out:
        for seq_id, length, score in candidates:
            seq = all_seqs.get(seq_id, "")
            if not seq:
                log_warn(
                    f"  [borderline rescue] {seq_id}: sequence not found in "
                    f"{contigs.name} -- skipping"
                )
                continue
            log_warn(
                f"  [borderline rescue] {seq_id}: {length:,}bp  "
                f"score={score:.3f} (< min_score {min_score}, "
                f">= rescue_min_score {rescue_min_score}) "
                f"-- added to viral FASTA (purified_phage, novel lineage)"
            )
            f_out.write(f">{seq_id}\n{seq}\n")
            n_added += 1

    if n_added:
        log_ok(
            f"  [borderline rescue] {n_added} contig(s) added -- "
            "verify in downstream CheckV + geNomad taxonomy output"
        )
    return n_added


# -- Failure diagnosis ---------------------------------------------------------

def _diagnose_genomad_failure(
    cfg:       Config,
    sdir:      Path,
    prefix:    str,
    rpt_dir:   Path,
    sample_id: str,
    n_contigs: int,
) -> None:
    summary_dir = sdir / f"{prefix}_summary"
    log_path    = rpt_dir / f"{sample_id}_genomad.log"
    if summary_dir.exists() and any(summary_dir.iterdir()):
        log_warn(
            f"  [geNomad] 0 viral contigs from {n_contigs} input contigs "
            f"(min-score={cfg.genomad.min_score}). "
            f"Options: lower genomad.min_score in config.yaml or set "
            f"min_virus_hallmarks: 0."
        )
    else:
        log_warn(
            f"  [geNomad] execution failed -- no summary at {summary_dir}. "
            f"Check: {log_path}"
        )


# -- Validation ----------------------------------------------------------------

def _validate_output(sample_id: str, virus_fna: Path, row: dict) -> None:
    if not virus_fna.exists() or virus_fna.stat().st_size == 0:
        log_warn(f"  validate . virus FASTA missing or empty: {virus_fna}")
        return
    n = row.get("viral_ctg", 0)
    if n == 0:
        log_warn(
            "  validate . 0 viral contigs -- lower genomad.min_score "
            "or check assembly output."
        )
    else:
        log_ok(f"  validate . {n} viral contig(s) identified -- OK")


# -- Parsers -------------------------------------------------------------------

def _parse_genomad_output(
    cfg:        Config,
    sample_id:  str,
    n_input:    int,
    sdir:       Path,
    prefix:     str,
    v_fna:      Path,
    v_tsv:      Path,
    n_rescued:  int,
) -> dict:
    v_stats  = fasta_stats(v_fna) if v_fna.exists() else {"n": 0, "total_bp": 0}
    n_viral  = v_stats["n"]
    viral_bp = v_stats["total_bp"]

    p_tsv  = sdir / f"{prefix}_summary" / f"{prefix}_plasmid_summary.tsv"
    n_plas = _count_tsv_rows(p_tsv)

    pr_tsv = sdir / f"{prefix}_find_proviruses" / f"{prefix}_provirus.tsv"
    if pr_tsv.exists():
        n_prov = _count_tsv_rows(pr_tsv)
    else:
        n_prov = 0
        if (sdir / f"{prefix}_summary").exists():
            log_info(
                "  [geNomad] find_proviruses/ not found (--cleanup or "
                "no proviruses) -- provirus count set to 0."
            )

    best_taxon, best_score = _parse_best_hit(v_tsv)

    return {
        "sample":              sample_id,
        "input_ctg":           n_input,
        "viral_ctg":           n_viral,
        "viral_bp":            viral_bp,
        "plasmid":             n_plas,
        "provirus":            n_prov,
        "best_taxon":          best_taxon,
        "best_score":          round(best_score, 4),
        "rescued_borderline":  n_rescued,
        "genomad_tax_tsv":     str(v_tsv),
    }


def _parse_best_hit(v_tsv: Path) -> tuple[str, float]:
    if not v_tsv.exists():
        return "unclassified", 0.0
    best_taxon = "unclassified"
    best_score = 0.0
    try:
        with open(v_tsv) as f:
            raw_header = f.readline().strip()
            if not raw_header:
                return best_taxon, best_score
            header = raw_header.split("\t")
            col    = {h: i for i, h in enumerate(header)}
            for required in (_COL_SCORE, _COL_TAXONOMY):
                if required not in col:
                    log_warn(
                        f"  [geNomad] column '{required}' missing in "
                        f"{v_tsv.name} -- version mismatch?"
                    )
                    return best_taxon, best_score
            score_idx = col[_COL_SCORE]
            tax_idx   = col[_COL_TAXONOMY]
            for line in f:
                row = line.strip().split("\t")
                if not row or len(row) <= max(score_idx, tax_idx):
                    continue
                try:
                    score = float(row[score_idx])
                except (ValueError, IndexError):
                    continue
                if score > best_score:
                    best_score = score
                    levels = [
                        lvl.strip()
                        for lvl in row[tax_idx].split(";")
                        if lvl.strip()
                        and lvl.strip().lower() != "unclassified"
                    ]
                    best_taxon = levels[-1] if levels else "unclassified"
    except Exception as e:
        log_warn(f"  [geNomad] could not parse {v_tsv.name}: {e}")
    return best_taxon, best_score


def _count_tsv_rows(path: Path) -> int:
    if not path.exists():
        return 0
    try:
        with open(path) as f:
            lines = [l for l in f if l.strip() and not l.startswith("#")]
        return max(0, len(lines) - 1)
    except Exception:
        return 0


# -- Rich display helpers -------------------------------------------------------

def _color_score(score: float) -> str:
    s = f"{score:.3f}"
    if score >= _SCORE_GOOD: return f"[bold green]{s}[/bold green]"
    if score >= _SCORE_WARN: return f"[bold yellow]{s}[/bold yellow]"
    return f"[bold red]{s}[/bold red]"


def _color_rate(val_str: str, good: float, warn: float) -> str:
    try:
        v = float(str(val_str).rstrip("%"))
        if v >= good: return f"[bold green]{val_str}[/bold green]"
        if v >= warn: return f"[bold yellow]{val_str}[/bold yellow]"
        return f"[bold red]{val_str}[/bold red]"
    except (ValueError, AttributeError):
        return str(val_str)


def _fmt(val, unit: str = "") -> str:
    if val in ("N/A", "", None):
        return "[dim]N/A[/dim]"
    return f"{val}{unit}"


def _print_input_table(sample_id: str, contigs: Path, stats: dict) -> None:
    log_info(f"  Contigs : {contigs}  ({human_size(contigs)})")
    log_info(
        f"  Input   : {stats['n']} sequence(s) | "
        f"total={stats['total_bp']:,}bp | largest={stats['largest_bp']:,}bp"
    )


def _print_summary_table(sample_id: str, row: dict) -> None:
    n_in      = row.get("input_ctg",         0)
    n_vir     = row.get("viral_ctg",         0)
    n_bp      = row.get("viral_bp",          0)
    n_pls     = row.get("plasmid",           0)
    n_prv     = row.get("provirus",          0)
    n_resc    = row.get("rescued_borderline",0)
    taxon     = row.get("best_taxon",        "unclassified")
    score     = row.get("best_score",        0.0)

    pct_vir     = f"{n_vir / n_in * 100:.1f}%" if n_in else "N/A"
    colored_pct = _color_rate(pct_vir, 50.0, 5.0) if n_in else "[dim]N/A[/dim]"

    log_ok(
        f"  Contigs  : {n_in} input -> {n_vir} viral "
        f"({colored_pct})  |  {n_bp:,}bp total viral"
    )
    if n_resc:
        log_ok(f"  Rescued  : {n_resc} borderline large contig(s) added")
    log_ok(f"  Other    : plasmid={_fmt(n_pls)}  provirus={_fmt(n_prv)}")
    log_ok(f"  Best hit : {taxon}  |  score={_color_score(score)}")


def _check_warnings(row: dict) -> None:
    n_in  = row.get("input_ctg", 0)
    n_vir = row.get("viral_ctg", 0)
    if n_vir == 0:
        log_warn(
            "  No viral contigs -- lower genomad.min_score in config.yaml."
        )
        return
    if n_in and (n_vir / n_in) < _VIRAL_FRACTION_WARN:
        log_warn(
            f"  Only {n_vir}/{n_in} contigs viral ({n_vir / n_in * 100:.1f}%) "
            "-- verify host removal or check for contamination."
        )
    if row.get("provirus", 0) > 0:
        log_warn(
            f"  {row['provirus']} provirus(es) detected -- "
            "review lifecycle prediction (Module 08)."
        )
    if 0 < row.get("best_score", 0.0) < _SCORE_WARN:
        log_warn(
            f"  Best score {row['best_score']:.3f} < {_SCORE_WARN} -- "
            "borderline; consider manual inspection."
        )


def _print_completion_panel(
    sample_id: str, virus_fna: Path, rpt_dir: Path, row: dict
) -> None:
    n_in  = row.get("input_ctg",  "?")
    n_vir = row.get("viral_ctg",  "?")
    taxon = row.get("best_taxon", "unclassified")
    score = row.get("best_score", 0.0)
    n_pls = row.get("plasmid",    0)
    n_prv = row.get("provirus",   0)
    n_rsc = row.get("rescued_borderline", 0)

    text = Text()
    text.append("v ", style="bold green")
    text.append(f"{n_in} contigs  ->  ", style="dim white")
    text.append(f"{n_vir}", style="bold green")
    text.append(
        f" viral  (plasmid={n_pls}  provirus={n_prv}"
        + (f"  rescued={n_rsc}" if n_rsc else "") + ")\n",
        style="cyan",
    )
    text.append(f"  best hit : {taxon}  (score={score:.3f})\n\n", style="dim white")
    text.append("Virus FASTA : ", style="dim white")
    text.append(str(virus_fna) + "\n", style="white")
    text.append("Summary     : ", style="dim white")
    text.append(str(rpt_dir / "genomad_summary.tsv"), style="white")

    console.print(Panel(
        text,
        title=f"[bold cyan]Viral-id complete -- {sample_id}[/bold cyan]",
        border_style="cyan",
        padding=(0, 2),
        width=90,
    ))


# -- TSV summary ---------------------------------------------------------------

_TSV_HEADERS = [
    "sample", "input_ctg", "viral_ctg", "viral_bp",
    "plasmid", "provirus", "best_taxon", "best_score",
    "rescued_borderline", "genomad_tax_tsv",
]


def _save_tsv(row: dict, path: Path) -> None:
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

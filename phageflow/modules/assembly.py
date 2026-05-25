"""Module"""
from __future__ import annotations
import shutil
import subprocess
from pathlib import Path
from typing import Optional

from rich.panel import Panel
from rich.progress import (
    Progress, SpinnerColumn, TextColumn,
    BarColumn, MofNCompleteColumn, TimeElapsedColumn,
)

from phageflow.utils.config import Config
from phageflow.utils.logger import (
    log_step, log_info, log_ok, log_warn, log_error, console,
)
from phageflow.utils.tools import require_tools, run_silent, human_size, mkdirs

STEP  = "03_assembly"
TOOLS = ["spades.py", "megahit", "cd-hit-est"]

_N50_WARN        = 5_000
_MAX_CONTIG_WARN = 500_000

def run(
    cfg:       Config,
    sample_id: str,
    r1:        Path,
    r2:        Path,
    s1:        Optional[Path] = None,
    force:     bool           = False,
) -> Path:
    """Assemble phage reads for a single sample.

    Parameters
    ----------
    r1, r2 : paired FASTQ from host-removal step
    s1     : singleton FASTQ from bwa-mem2 (optional, improves DTR/ITR coverage)

    Returns
    -------
    Path to NR combined contigs FASTA
    """
    r1 = Path(r1); r2 = Path(r2)
    out_dir     = cfg.results(STEP)
    rpt_dir     = cfg.reports(STEP) / sample_id
    spades_dir  = out_dir / f"{sample_id}_spades"
    megahit_dir = out_dir / f"{sample_id}_megahit"
    contigs_out = out_dir / f"{sample_id}_contigs.fasta"
    log_f       = rpt_dir / f"{sample_id}_assembly.log"
    mkdirs(out_dir, rpt_dir)

    log_step(f"Module 03 — assembly  [{sample_id}]")
    log_info(f"  R1 : {r1}  ({human_size(r1)})")
    log_info(f"  R2 : {r2}  ({human_size(r2)})")

    use_s1 = bool(s1 and Path(s1).exists() and Path(s1).stat().st_size > 0)
    if use_s1:
        log_info(f"  S1 : {s1}  ({human_size(Path(s1))})  → SPAdes --s1")
    else:
        log_info("  S1 : not provided or empty — DTR/ITR boundary reads absent")

    if contigs_out.exists() and contigs_out.stat().st_size > 0 and not force:
        log_info("  Already assembled — skipping  (--force to re-run)")
        return contigs_out

    require_tools(*TOOLS)
    active_warnings: list[str] = []

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

        progress.update(task, description="[1/4] SPAdes  — isolate mode + singletons")
        spades_ok = _run_spades(
            r1, r2, s1 if use_s1 else None,
            spades_dir, log_f, cfg.threads, cfg.memory_gb, force,
        )
        if not spades_ok:
            msg = "SPAdes failed — assembly will use MEGAHIT output only."
            log_warn(f"  {msg}"); active_warnings.append(msg)
        progress.advance(task)

        progress.update(task, description="[2/4] MEGAHIT — no-mercy high-coverage mode")
        megahit_ok = _run_megahit(
            r1, r2, megahit_dir, sample_id,
            cfg.assembly.kmers, log_f, cfg.threads, force,
        )
        if not megahit_ok:
            if not spades_ok:
                log_error("  Both SPAdes and MEGAHIT failed — no contigs produced.")
                raise RuntimeError(
                    f"Assembly failed for {sample_id}: both assemblers returned non-zero."
                )
            msg = "MEGAHIT failed — using SPAdes contigs only."
            log_warn(f"  {msg}"); active_warnings.append(msg)
        progress.advance(task)

        progress.update(task, description="[3/4] cd-hit-est — combine + NR reduction")
        stats = _combine_and_dereplicate(
            spades_dir if spades_ok else None,
            megahit_dir / f"{sample_id}.contigs.fa" if megahit_ok else None,
            contigs_out, cfg.assembly.min_length,
            cfg.threads, log_f, active_warnings,
        )
        progress.advance(task)

        progress.update(task, description="[4/4] Computing assembly stats")
        _check_assembly_warnings(stats, active_warnings)
        _save_tsv({**stats, "sample": sample_id,
                   "kmers": cfg.assembly.kmers,
                   "singletons_used": str(use_s1),
                   "spades_ok": str(spades_ok),
                   "megahit_ok": str(megahit_ok)},
                  rpt_dir / "assembly_summary.tsv")
        _print_completion_panel(
            sample_id, contigs_out, rpt_dir, stats,
            spades_ok, megahit_ok, use_s1, active_warnings,
        )
        progress.advance(task)

    log_step(f"Module 03 completed ✓  [{sample_id}]")
    return contigs_out

# ── SPAdes ────────────────────────────────────────────────────────────────────

def _run_spades(
    r1: Path, r2: Path, s1: Optional[Path],
    out_dir: Path, log_f: Path,
    threads: int, memory_gb: int, force: bool,
) -> bool:
    contigs_fa = out_dir / "contigs.fasta"
    if contigs_fa.exists() and contigs_fa.stat().st_size > 0 and not force:
        log_info("  [SPAdes] already assembled — skipping")
        return True

    if out_dir.exists() and force:
        shutil.rmtree(out_dir, ignore_errors=True)

    out_dir.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "spades.py",
        "--pe1-1", str(r1),
        "--pe1-2", str(r2),
        "-k", "21,33,55,77,99,127",
        "--isolate",
        "--phred-offset", "33",
        "-t", str(threads),
        "-m", str(memory_gb),
        "-o", str(out_dir),
    ]
    if s1 and Path(s1).exists():
        cmd += ["--s1", str(s1)]

    try:
        run_silent(cmd, log_file=log_f)
    except Exception as e:
        log_warn(f"  [SPAdes] non-zero exit: {e}")
        return False

    if not contigs_fa.exists() or contigs_fa.stat().st_size == 0:
        log_warn("  [SPAdes] contigs.fasta missing or empty after run")
        return False

    s = _fasta_stats(contigs_fa)
    log_ok(f"  [SPAdes] {s['n']:,} contigs  N50={s['n50_bp']:,} bp  "
           f"largest={s['largest_bp']:,} bp")
    _cleanup_spades(out_dir)
    return True

def _cleanup_spades(out_dir: Path) -> None:
    """Remove SPAdes intermediates; keep contigs.fasta, assembly_graph.gfa, spades.log."""
    keep = {"contigs.fasta", "assembly_graph.gfa", "assembly_graph.fastg",
            "spades.log", "params.txt"}
    for child in list(out_dir.iterdir()):
        if child.name not in keep:
            if child.is_dir():
                shutil.rmtree(child, ignore_errors=True)
            else:
                child.unlink(missing_ok=True)
    log_info("  [SPAdes] intermediates removed")

# ── MEGAHIT ───────────────────────────────────────────────────────────────────

def _megahit_supports_klist() -> bool:
    """Check if installed MEGAHIT supports --k-list (requires ≥ 1.2.9)."""
    try:
        r = subprocess.run(
            ["megahit", "--help"],
            capture_output=True, text=True, timeout=10,
        )
        return "--k-list" in (r.stdout + r.stderr)
    except Exception:
        return False

def _run_megahit(
    r1: Path, r2: Path,
    out_dir: Path, prefix: str,
    kmers: str, log_f: Path,
    threads: int, force: bool,
) -> bool:
    final_contigs = out_dir / f"{prefix}.contigs.fa"
    if final_contigs.exists() and final_contigs.stat().st_size > 0 and not force:
        log_info("  [MEGAHIT] already assembled — skipping")
        return True

    if out_dir.exists():
        shutil.rmtree(out_dir, ignore_errors=True)
    out_dir.parent.mkdir(parents=True, exist_ok=True)

    kmer_vals = [int(k) for k in kmers.split(",") if k.strip().isdigit()]
    if _megahit_supports_klist():
        kmer_flags = ["--k-list", kmers]
    else:
        k_min  = str(min(kmer_vals))
        k_max  = str(max(kmer_vals))
        diffs  = [b - a for a, b in zip(kmer_vals, kmer_vals[1:])]
        k_step = str(min(diffs)) if diffs else "12"
        kmer_flags = ["--k-min", k_min, "--k-max", k_max, "--k-step", k_step]
        log_info(f"  [MEGAHIT] legacy k-mer flags: {' '.join(kmer_flags)}")

    cmd = [
        "megahit",
        "-1", str(r1), "-2", str(r2),
        *kmer_flags,
        "--min-count", "2",
        "--no-mercy",
        "--min-contig-len", "200",
        "-t", str(threads),
        "-o", str(out_dir),
        "--out-prefix", prefix,
    ]
    try:
        run_silent(cmd, log_file=log_f)
    except Exception as e:
        log_warn(f"  [MEGAHIT] non-zero exit: {e}")
        return False

    if not final_contigs.exists() or final_contigs.stat().st_size == 0:
        log_warn("  [MEGAHIT] contigs file missing or empty after run")
        return False

    s = _fasta_stats(final_contigs)
    log_ok(f"  [MEGAHIT] {s['n']:,} contigs  N50={s['n50_bp']:,} bp  "
           f"largest={s['largest_bp']:,} bp")
    return True

# ── Combine + filter + NR ─────────────────────────────────────────────────────

def _combine_and_dereplicate(
    spades_dir:      Optional[Path],
    megahit_fa:      Optional[Path],
    contigs_out:     Path,
    min_length:      int,
    threads:         int,
    log_f:           Path,
    active_warnings: list[str],
) -> dict:
    """Combine SPAdes + MEGAHIT contigs, filter by length, run cd-hit-est NR.

    cd-hit-est: -c 1.00 -aS 0.85 -n 8 (Fu et al. 2012)
    """
    inputs: list[tuple[str, Path]] = []
    if spades_dir:
        sp_fa = spades_dir / "contigs.fasta"
        if sp_fa.exists() and sp_fa.stat().st_size > 0:
            inputs.append(("SPAdes", sp_fa))
    if megahit_fa and megahit_fa.exists() and megahit_fa.stat().st_size > 0:
        inputs.append(("MEGAHIT", megahit_fa))

    if not inputs:
        raise RuntimeError("No contig FASTA files found from either assembler.")

    combined_fa  = contigs_out.parent / f"{contigs_out.stem}_combined_raw.fasta"
    cdhit_prefix = contigs_out.parent / f"{contigs_out.stem}_cdhit"

    n_raw: dict[str, int] = {}
    with open(combined_fa, "w") as out:
        for source, fa in inputs:
            n = _write_tagged_fasta(fa, out, tag=source, min_len=min_length)
            n_raw[source] = n
            log_info(f"  [{source}] {n:,} contigs ≥{min_length} bp added to pool")

    total_raw = sum(n_raw.values())
    if total_raw == 0:
        msg = f"No contigs ≥{min_length} bp from any assembler."
        log_warn(f"  {msg}"); active_warnings.append(msg)
        combined_fa.unlink(missing_ok=True)
        contigs_out.write_text("")
        return {"n_spades": "0", "n_megahit": "0", "n_combined": "0",
                "n_final": "0", "total_bp": "0", "largest_bp": "0", "n50_bp": "0"}

    run_silent(
        ["cd-hit-est",
         "-i", str(combined_fa), "-o", str(cdhit_prefix),
         "-c", "1.00", "-aS", "0.85", "-n", "8",
         "-d", "0", "-T", str(threads), "-M", "0", "-p", "1"],
        log_file=log_f,
    )

    if cdhit_prefix.exists():
        cdhit_prefix.rename(contigs_out)
    else:
        log_warn("  [cd-hit-est] output not found — using unfiltered combined FASTA")
        combined_fa.rename(contigs_out)

    for p in [combined_fa,
              Path(str(cdhit_prefix) + ".clstr"),
              Path(str(contigs_out)  + ".clstr")]:
        p.unlink(missing_ok=True)

    stats = _fasta_stats(contigs_out)
    log_ok(f"  [NR] {total_raw:,} → {stats['n']:,} contigs after deduplication")

    return {
        "n_spades":   str(n_raw.get("SPAdes",  0)),
        "n_megahit":  str(n_raw.get("MEGAHIT", 0)),
        "n_combined": str(total_raw),
        "n_final":    str(stats["n"]),
        "total_bp":   str(stats["total_bp"]),
        "largest_bp": str(stats["largest_bp"]),
        "n50_bp":     str(stats["n50_bp"]),
    }

def _write_tagged_fasta(in_fa: Path, out_fh, tag: str, min_len: int) -> int:
    """Write contigs ≥ min_len to out_fh, prefixing headers with [tag]."""
    n = 0
    header = ""
    seq_parts: list[str] = []

    def _flush():
        nonlocal n
        if not header:
            return
        seq = "".join(seq_parts)
        if len(seq) >= min_len:
            out_fh.write(f">{tag}|{header[1:].strip()}\n{seq}\n")
            n += 1

    with open(in_fa) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                _flush()
                header = line
                seq_parts = []
            else:
                seq_parts.append(line)
        _flush()
    return n

# ── Sequence utilities ────────────────────────────────────────────────────────

def _fasta_stats(path: Path) -> dict:
    lengths: list[int] = []
    current = 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if current > 0:
                    lengths.append(current)
                current = 0
            else:
                current += len(line.strip())
        if current > 0:
            lengths.append(current)
    if not lengths:
        return {"n": 0, "total_bp": 0, "largest_bp": 0, "n50_bp": 0}
    return {
        "n":          len(lengths),
        "total_bp":   sum(lengths),
        "largest_bp": max(lengths),
        "n50_bp":     _n50(lengths),
    }

def _n50(lengths: list[int]) -> int:
    sorted_l = sorted(lengths, reverse=True)
    half     = sum(sorted_l) / 2
    cumsum   = 0
    for l in sorted_l:
        cumsum += l
        if cumsum >= half:
            return l
    return 0

# ── Warnings ──────────────────────────────────────────────────────────────────

def _check_assembly_warnings(stats: dict, active_warnings: list[str]) -> None:
    try:
        n50 = int(stats.get("n50_bp", 0))
        if n50 == 0:
            msg = "N50 = 0 — no contigs produced."
            log_warn(f"  {msg}"); active_warnings.append(msg)
        elif n50 < _N50_WARN:
            msg = (f"N50 = {n50:,} bp < {_N50_WARN:,} bp — "
                   "fragmented assembly. Check host removal and coverage depth.")
            log_warn(f"  {msg}"); active_warnings.append(msg)
    except (ValueError, TypeError):
        pass
    try:
        largest = int(stats.get("largest_bp", 0))
        if largest > _MAX_CONTIG_WARN:
            msg = (f"Largest contig = {largest:,} bp > {_MAX_CONTIG_WARN:,} bp — "
                   "unusually large. Possible residual host contamination.")
            log_warn(f"  {msg}"); active_warnings.append(msg)
    except (ValueError, TypeError):
        pass

# ── Completion panel ──────────────────────────────────────────────────────────

def _print_completion_panel(
    sample_id:  str,
    contigs_out: Path,
    rpt_dir:    Path,
    stats:      dict,
    spades_ok:  bool,
    megahit_ok: bool,
    singletons: bool,
    active_warnings: list[str],
) -> None:
    def _f(val) -> str:
        v = str(val)
        return "[dim]N/A[/dim]" if v in ("0", "", "None", "N/A") else f"{int(v):,}"

    def _green(val) -> str:
        try:
            v = int(val)
            return (f"[bold green]{v:,}[/bold green]" if v > 0
                    else f"[bold red]{v:,}[/bold red]")
        except: return str(val)

    spades_str  = "[green]✓[/green]" if spades_ok  else "[red]✗ failed[/red]"
    megahit_str = "[green]✓[/green]" if megahit_ok else "[red]✗ failed[/red]"
    s1_str      = "[green]✓[/green]" if singletons else "[dim]none[/dim]"

    lines = [
        "[bold]Assembly statistics[/bold]",
        f"  [cyan]Final contigs  :[/cyan] {_green(stats.get('n_final', 0))} "
        f"  [dim](raw={_f(stats.get('n_combined'))}  "
        f"SPAdes={_f(stats.get('n_spades'))}  "
        f"MEGAHIT={_f(stats.get('n_megahit'))})[/dim]",
        f"  [cyan]N50            :[/cyan] {_f(stats.get('n50_bp'))} bp",
        f"  [cyan]Largest contig :[/cyan] {_f(stats.get('largest_bp'))} bp",
        f"  [cyan]Total assembly :[/cyan] {_f(stats.get('total_bp'))} bp",
        "",
        "[bold]Assemblers[/bold]",
        f"  SPAdes : {spades_str}  MEGAHIT : {megahit_str}  Singletons : {s1_str}",
        "",
        "[bold]Output[/bold]",
        f"  NR contigs : {contigs_out}",
        f"  Summary    : {rpt_dir / 'assembly_summary.tsv'}",
    ]
    if active_warnings:
        lines += ["", "[bold yellow]Warnings[/bold yellow]"]
        lines += [f"  [yellow]⚠ {w}[/yellow]" for w in active_warnings]

    console.print(Panel(
        "\n".join(lines),
        title=f"[bold cyan]Assembly complete — {sample_id}[/bold cyan]",
        border_style="cyan", padding=(0, 2), width=120,
    ))

# ── TSV summary ───────────────────────────────────────────────────────────────

_TSV_HEADERS = [
    "sample", "spades_ok", "megahit_ok", "singletons_used",
    "n_spades", "n_megahit", "n_combined", "n_final",
    "total_bp", "largest_bp", "n50_bp", "kmers",
]

def _save_tsv(row: dict, path: Path) -> None:
    rows: dict[str, dict] = {}
    if path.exists():
        with open(path) as f:
            old_h = f.readline().rstrip("\n").split("\t")
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if cols and cols[0]: rows[cols[0]] = dict(zip(old_h, cols))
    rows[row["sample"]] = {h: str(row.get(h, "")) for h in _TSV_HEADERS}
    with open(path, "w") as f:
        f.write("\t".join(_TSV_HEADERS) + "\n")
        for r in rows.values():
            f.write("\t".join(r.get(h, "") for h in _TSV_HEADERS) + "\n")

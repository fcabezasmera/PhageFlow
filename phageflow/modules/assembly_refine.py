"""Module"""

from __future__ import annotations
import csv
import shutil
import subprocess
from collections import defaultdict
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
from phageflow.utils.tools import require_tools, run_silent, mkdirs

STEP = "05b_refine"

# CheckV quality tier ordering (higher = better)
_TIER_ORDER = {
    "Complete": 4, "High-quality": 3, "Medium-quality": 2,
    "Low-quality": 1, "Not-determined": 0, "": 0,
}

# ── Public entry point ────────────────────────────────────────────────────────

def run(
    cfg:       Config,
    sample_id: str,
    r1:        Path,
    r2:        Path,
    s1:        Optional[Path] = None,
    force:     bool = False,
) -> Path:
    """Iterative assembly refinement for a single sample.

    Parameters
    ----------
    r1, r2 : post-host-removal reads (from 02_host_removal step)
    s1     : singletons from bwa-mem2 host-removal (optional)

    Returns
    -------
    Path to annotation_ready/phages/ (may contain refined candidates)
    """
    r1 = Path(r1); r2 = Path(r2)
    quality_dir = cfg.results("05_quality") / sample_id
    ann_phages  = quality_dir / "annotation_ready" / "phages"
    pre_refine  = ann_phages / "pre_refine"
    rpt_dir     = cfg.reports(STEP) / sample_id
    log_f       = rpt_dir / f"{sample_id}_refine.log"
    mkdirs(rpt_dir)

    log_step(f"Module 05b — assembly refinement  [{sample_id}]")

    candidates = sorted(ann_phages.glob("*.fasta"))
    if not candidates:
        log_warn(
            f"  No annotation_ready candidates found in {ann_phages}. "
            "Run quality.py first."
        )
        return ann_phages

    log_info(f"  Candidates : {len(candidates)} in {ann_phages}")
    log_info(f"  R1 : {r1}")
    log_info(f"  R2 : {r2}")

    if (pre_refine / candidates[0].name).exists() and not force:
        log_info("  Already refined — skipping  (--force to re-run)")
        return ann_phages

    require_tools("bwa-mem2", "samtools", "spades.py", "checkv")

    active_warnings: list[str] = []

    with Progress(
        SpinnerColumn(spinner_name="dots2"),
        TextColumn("[bold cyan]{task.description:<62}"),
        BarColumn(bar_width=20, complete_style="cyan", finished_style="green"),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        console=console,
        transient=True,
    ) as progress:
        task = progress.add_task("Initializing...", total=5)

        # ── Step 1: Build combined reference ──────────────────────────────────
        progress.update(task, description="[1/5] Building combined reference from candidates")
        combined = rpt_dir / "combined_candidates.fasta"
        _build_combined_reference(candidates, combined)
        progress.advance(task)

        # ── Step 2: Map reads, extract unmapped + semi-mapped ─────────────────
        progress.update(task, description="[2/5] Mapping reads → extracting unmapped + bridge reads")
        unmapped_r1 = rpt_dir / "unmapped_R1.fastq.gz"
        unmapped_r2 = rpt_dir / "unmapped_R2.fastq.gz"
        semi_mapped = rpt_dir / "semi_mapped.fastq.gz"
        n_unmapped, n_semi = _extract_unmapped_reads(
            combined, r1, r2, unmapped_r1, unmapped_r2, semi_mapped,
            rpt_dir, log_f, cfg.threads,
        )
        log_info(
            f"  [map] unmapped pairs: {n_unmapped:,}  "
            f"semi-mapped (bridge): {n_semi:,}"
        )
        if n_unmapped + n_semi < 100:
            log_warn(
                "  Very few unmapped reads — candidates likely cover the full genome. "
                "Refinement may not improve assembly."
            )
        progress.advance(task)

        # ── Step 3: SPAdes with trusted-contigs ───────────────────────────────
        progress.update(task, description="[3/5] SPAdes --trusted-contigs refinement")
        spades_dir = rpt_dir / "spades_refine"
        _run_spades_refine(
            combined, unmapped_r1, unmapped_r2, semi_mapped, s1,
            spades_dir, log_f, cfg, force,
        )
        progress.advance(task)

        # ── Step 4: CheckV on refined assembly ────────────────────────────────
        progress.update(task, description="[4/5] CheckV — evaluating refined assembly")
        checkv_dir = rpt_dir / "checkv_refine"
        refined_contigs = spades_dir / "contigs.fasta"
        if not refined_contigs.exists() or refined_contigs.stat().st_size == 0:
            log_warn("  SPAdes produced no contigs — refinement failed")
            progress.advance(task); progress.advance(task)
            return ann_phages

        _run_checkv(refined_contigs, checkv_dir, cfg, log_f, force)
        progress.advance(task)

        # ── Step 5: Compare and replace improved candidates ───────────────────
        progress.update(task, description="[5/5] Comparing quality — replacing improved candidates")
        n_replaced, n_new = _apply_refinements(
            candidates, ann_phages, pre_refine,
            refined_contigs, checkv_dir,
            cfg, active_warnings,
        )
        log_ok(
            f"  [refine] {n_replaced} candidate(s) improved  "
            f"{n_new} new candidate(s) added"
        )
        progress.advance(task)

    _print_completion_panel(sample_id, ann_phages, n_replaced, n_new, active_warnings)
    log_step(f"Module 05b completed ✓  [{sample_id}]")
    return ann_phages

# ── Reference building ────────────────────────────────────────────────────────

def _build_combined_reference(candidates: list[Path], combined: Path) -> None:
    """Concatenate all annotation_ready candidates into a single reference FASTA."""
    with open(combined, "w") as out:
        for fa in candidates:
            cand_id = fa.stem
            seq_parts: list[str] = []
            with open(fa) as f:
                for line in f:
                    if line.startswith(">"):
                        if seq_parts:
                            out.write("".join(seq_parts) + "\n")
                        out.write(f">{cand_id}\n")
                        seq_parts = []
                    else:
                        seq_parts.append(line.rstrip())
            if seq_parts:
                out.write("".join(seq_parts) + "\n")
    log_info(f"  [ref] Combined {len(candidates)} candidates → {combined.name}")

# ── Read extraction ───────────────────────────────────────────────────────────

def _extract_unmapped_reads(
    combined:    Path,
    r1: Path, r2: Path,
    unmapped_r1: Path, unmapped_r2: Path,
    semi_mapped: Path,
    tmp_dir:     Path,
    log_f:       Path,
    threads:     int,
) -> tuple[int, int]:
    """Map reads to combined reference; extract unmapped + semi-mapped.

    Unmapped pairs  (-f 12):  both mates unmapped → NOT in any candidate
    Semi-mapped     (-F 2 -f 1 -F 12): one mate mapped, one not
                              → bridge reads spanning candidate ends

    These two sets together provide SPAdes with material to extend beyond
    the existing candidate boundaries and bridge gaps between fragments.
    """
    sort_prefix = tmp_dir / "sort_tmp"
    bam         = tmp_dir / "mapped.bam"
    st_view     = min(threads, 8)
    st_sort     = min(threads, 8)

    run_silent(["bwa-mem2", "index", str(combined)], log_file=log_f)

    # Align + sort
    run_silent(
        f"bwa-mem2 mem -t {threads} -M {combined} {r1} {r2} "
        f"| samtools sort -@ {st_sort} -n -u -T {sort_prefix} - "
        f"-o {bam}",
        log_file=log_f, shell=True,
    )

    # Extract unmapped pairs (both mates unmapped: -f 12)
    run_silent(
        f"samtools fastq -@ {min(threads,4)} -N -f 12 {bam} "
        f"-1 {unmapped_r1} -2 {unmapped_r2} -0 /dev/null",
        log_file=log_f, shell=True,
    )

    # Extract semi-mapped (paired, NOT both properly mapped, NOT both unmapped)
    # These are bridge reads: one mate maps to a candidate, the other doesn't
    run_silent(
        f"samtools fastq -@ {min(threads,4)} -N -f 1 -F 14 {bam} "
        f"--singleton {semi_mapped} -0 /dev/null -1 /dev/null -2 /dev/null",
        log_file=log_f, shell=True, check=False,
    )

    n_unmapped = _count_fastq(unmapped_r1)
    n_semi     = _count_fastq(semi_mapped)
    return n_unmapped, n_semi

def _count_fastq(path: Path) -> int:
    if not path.exists() or path.stat().st_size == 0:
        return 0
    try:
        import shutil as _sh
        decomp = "pigz -dc" if _sh.which("pigz") else "gzip -dc"
        r = subprocess.run(
            ["bash", "-c", f"{decomp} {path} | wc -l"],
            capture_output=True, text=True, timeout=60,
        )
        return int(r.stdout.strip()) // 4
    except Exception:
        return 0

# ── SPAdes refinement ─────────────────────────────────────────────────────────

def _run_spades_refine(
    combined:    Path,
    unmapped_r1: Path,
    unmapped_r2: Path,
    semi_mapped: Path,
    s1_original: Optional[Path],
    out_dir:     Path,
    log_f:       Path,
    cfg:         Config,
    force:       bool,
) -> None:
    """Run SPAdes with trusted-contigs + unmapped reads.

    --only-assembler: skip BayesHammer error correction (appropriate at
      high coverage and with trusted-contigs as scaffolding anchors).
      Note: --isolate is NOT used here (incompatible with --trusted-contigs).
    --trusted-contigs: SPAdes treats these as reliable sequence anchors.
      The assembler will try to extend through the gaps using the reads.
    Bankevich et al. 2012, J Comput Biol 19:455.
    """
    contigs_out = out_dir / "contigs.fasta"
    if contigs_out.exists() and not force:
        log_info("  [SPAdes-refine] already done — skipping")
        return

    if out_dir.exists():
        shutil.rmtree(out_dir, ignore_errors=True)
    out_dir.parent.mkdir(parents=True, exist_ok=True)

    if not unmapped_r1.exists() or _count_fastq(unmapped_r1) == 0:
        log_warn("  [SPAdes-refine] no unmapped reads — skipping SPAdes refinement")
        return

    cmd = [
        "spades.py",
        "--only-assembler",
        "--trusted-contigs", str(combined),
        "--pe1-1", str(unmapped_r1),
        "--pe1-2", str(unmapped_r2),
        "-k", ",".join(k for k in cfg.assembly.kmers.split(",") if int(k) <= 127),
        "-t", str(cfg.threads),
        "-m", str(cfg.memory_gb),
        "-o", str(out_dir),
    ]

    # Add singletons: semi-mapped bridge reads + original host-removal singletons
    singleton_parts: list[str] = []
    if semi_mapped.exists() and _count_fastq(semi_mapped) > 0:
        singleton_parts.append(str(semi_mapped))
    if s1_original and Path(s1_original).exists():
        singleton_parts.append(str(s1_original))

    if singleton_parts:
        # Merge singletons into one file for SPAdes --s1
        merged_s1 = out_dir.parent / "merged_s1.fastq.gz"
        import shutil as _sh
        pigz = "pigz" if _sh.which("pigz") else "gzip"
        merge_cmd = (
            " ".join(f"zcat {p}" for p in singleton_parts)
            + f" | {pigz} > {merged_s1}"
        )
        subprocess.run(["bash", "-c", merge_cmd], check=False, timeout=300)
        if merged_s1.exists() and merged_s1.stat().st_size > 0:
            cmd += ["--s1", str(merged_s1)]

    try:
        run_silent(cmd, log_file=log_f)
        log_ok("  [SPAdes-refine] assembly complete")
        _cleanup_spades(out_dir)
    except Exception as e:
        log_warn(f"  [SPAdes-refine] failed: {e}")

def _cleanup_spades(out_dir: Path) -> None:
    keep = {"contigs.fasta", "assembly_graph.gfa", "assembly_graph.fastg", "spades.log"}
    for child in list(out_dir.iterdir()):
        if child.name not in keep:
            if child.is_dir():
                shutil.rmtree(child, ignore_errors=True)
            else:
                child.unlink(missing_ok=True)

# ── CheckV evaluation ─────────────────────────────────────────────────────────

def _run_checkv(
    contigs: Path, checkv_dir: Path, cfg: Config, log_f: Path, force: bool,
) -> None:
    summary = checkv_dir / "quality_summary.tsv"
    if summary.exists() and not force:
        return
    db = cfg.databases.checkv
    run_silent(
        ["checkv", "end_to_end", str(contigs), str(checkv_dir),
         "-d", str(db), "-t", str(cfg.threads), "--remove_tmp"],
        log_file=log_f,
    )
    log_ok("  [CheckV-refine] quality assessment complete")

def _load_checkv_quality(checkv_dir: Path) -> dict:
    """Load CheckV quality_summary.tsv → {contig_id: {quality, completeness, length}}"""
    summary = checkv_dir / "quality_summary.tsv"
    if not summary.exists():
        return {}
    result = {}
    with open(summary, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            cid = row.get("contig_id", "").strip()
            result[cid] = {
                "quality":     row.get("checkv_quality", ""),
                "completeness": float(row.get("completeness") or 0) if (row.get("completeness") or "0") not in ("NA", "na", "") else 0.0,
                "length":      int(row.get("contig_length") or 0),
                "termini":     row.get("termini_type", ""),
            }
    return result

# ── Comparison and replacement ────────────────────────────────────────────────

def _apply_refinements(
    original_candidates: list[Path],
    ann_phages:          Path,
    pre_refine:          Path,
    refined_contigs:     Path,
    checkv_dir:          Path,
    cfg:                 Config,
    active_warnings:     list[str],
) -> tuple[int, int]:
    """Compare refined assembly with originals; replace when improved.

    A refined contig is "better" when:
      1. CheckV quality tier is strictly higher than the best existing candidate
      2. OR: refined contig is Complete/HQ and best existing is MQ/lower
      3. OR: refined contig length is >10% longer and quality tier is equal/better

    Originals are backed up to pre_refine/ before any replacement.

    Returns (n_replaced, n_new).
    """
    refined_quality = _load_checkv_quality(checkv_dir)
    if not refined_quality:
        log_warn("  [refine] CheckV produced no results for refined assembly")
        return 0, 0

    # Load refined FASTA sequences
    refined_seqs = _load_fasta(refined_contigs)

    # Current best tier per original candidate
    orig_quality = {}
    for fa in original_candidates:
        orig_quality[fa.stem] = _infer_original_quality(fa)

    best_orig_tier = max(
        (_TIER_ORDER.get(q.get("quality", ""), 0) for q in orig_quality.values()),
        default=0,
    )

    # Backup originals
    pre_refine.mkdir(parents=True, exist_ok=True)
    for fa in original_candidates:
        backup = pre_refine / fa.name
        if not backup.exists():
            shutil.copy2(fa, backup)

    n_replaced = 0
    n_new = 0

    # Find refined contigs that are strictly better than all originals
    improvements: list[tuple[str, dict]] = []
    for contig_id, q in refined_quality.items():
        refined_tier = _TIER_ORDER.get(q["quality"], 0)
        if refined_tier > best_orig_tier:
            improvements.append((contig_id, q))
        elif refined_tier == best_orig_tier and q["length"] > _best_orig_length(orig_quality) * 1.1:
            improvements.append((contig_id, q))

    if not improvements:
        log_info("  [refine] No improvements found — keeping original candidates")
        return 0, 0

    log_ok(f"  [refine] {len(improvements)} refined contig(s) are improvements")

    # Sort improvements by quality tier then completeness
    improvements.sort(key=lambda x: (-_TIER_ORDER.get(x[1]["quality"], 0),
                                     -x[1]["completeness"]))

    # Replace best original with best refined, add others as new candidates
    for i, (contig_id, q) in enumerate(improvements):
        seq = refined_seqs.get(contig_id, "")
        if not seq:
            continue
        if i == 0 and original_candidates:
            # Replace the top original candidate
            target = ann_phages / original_candidates[0].name
            _write_candidate_fasta(target, original_candidates[0].stem, seq, q)
            log_ok(
                f"  [refine] {original_candidates[0].stem} replaced by "
                f"{contig_id} ({q['quality']}, {q['completeness']:.1f}%)"
            )
            n_replaced += 1
        else:
            # Add as new candidate with _refined suffix
            stem    = f"{original_candidates[0].stem}_refined_{i:02d}" if original_candidates else f"Phage_refined_{i:02d}"
            target  = ann_phages / f"{stem}.fasta"
            _write_candidate_fasta(target, stem, seq, q)
            log_info(f"  [refine] New candidate added: {stem} ({q['quality']})")
            n_new += 1

    return n_replaced, n_new

def _infer_original_quality(fa: Path) -> dict:
    """Read quality from FASTA header comment written by quality.py."""
    with open(fa) as f:
        header = f.readline()
    quality = ""
    completeness = 0.0
    length = 0
    for part in header.split("["):
        part = part.rstrip("]").strip()
        if part.startswith("checkv="):
            quality = part.split("=", 1)[1]
        elif part.startswith("completeness="):
            try:
                completeness = float(part.split("=", 1)[1].rstrip("%"))
            except ValueError:
                pass
        elif part.startswith("length="):
            try:
                length = int(part.split("=", 1)[1])
            except ValueError:
                pass
    return {"quality": quality, "completeness": completeness, "length": length}

def _best_orig_length(orig_quality: dict) -> int:
    lengths = [q.get("length", 0) for q in orig_quality.values()]
    return max(lengths) if lengths else 0

def _load_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    current_id = ""
    parts: list[str] = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if current_id:
                    seqs[current_id] = "".join(parts)
                current_id = line[1:].split()[0].strip()
                parts = []
            else:
                parts.append(line.rstrip())
        if current_id:
            seqs[current_id] = "".join(parts)
    return seqs

def _write_candidate_fasta(path: Path, cand_id: str, seq: str, q: dict) -> None:
    with open(path, "w") as f:
        f.write(
            f">{cand_id}  "
            f"[length={q['length']}]  "
            f"[completeness={q['completeness']:.1f}%]  "
            f"[topology={q['termini']}]  "
            f"[checkv={q['quality']}]  "
            f"[source=assembly_refine]\n"
        )
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")

# ── Completion panel ──────────────────────────────────────────────────────────

def _print_completion_panel(
    sample_id: str,
    ann_phages: Path,
    n_replaced: int,
    n_new: int,
    active_warnings: list[str],
) -> None:
    candidates = list(ann_phages.glob("*.fasta"))
    lines = [
        "[bold]Refinement results[/bold]",
        f"  [cyan]Candidates replaced :[/cyan] {n_replaced}",
        f"  [cyan]New candidates added:[/cyan] {n_new}",
        f"  [cyan]Total in ann_ready  :[/cyan] {len(candidates)}",
        "",
        "[bold]Annotation-ready candidates[/bold]",
    ]
    for fa in sorted(candidates):
        lines.append(f"  [green]✓[/green] {fa.stem}")
    lines += [
        "",
        f"  annotation_ready/ : {ann_phages}",
    ]
    if active_warnings:
        lines += ["", "[bold yellow]Warnings[/bold yellow]"]
        lines += [f"  [yellow]⚠ {w}[/yellow]" for w in active_warnings]

    console.print(Panel(
        "\n".join(lines),
        title=f"[bold cyan]Assembly refinement — {sample_id}[/bold cyan]",
        border_style="cyan", padding=(0, 2), width=120,
    ))

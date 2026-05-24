"""PhageFlow Module 04 — Viral identification (geNomad).

Goal: identify viral contigs from the assembly, classify their topology,
determine genetic code, and output a filtered viral FASTA for quality.py.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Tool: geNomad end-to-end
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

geNomad combines k-mer marker gene search (MMseqs2) and a neural-network
classifier (Camargo et al. 2023, Nat Biotechnol 41:1783).

Key flags used:

  --enable-score-calibration --composition virome
      Adjusts scores for virome-dominated samples (purified phage).
      Reduces false negatives in novel lineages without close DB relatives.

  --lenient-taxonomy
      Resolves taxonomy below family rank (subfamily, genus, species).
      Required for biologically meaningful co-binning in quality.py.
      Without this, distinct genera within a family co-bin together.

  --disable-find-proviruses
      Skips provirus detection. For purified phage preparations, bacterial
      chromosomes should be absent. If proviruses appear, it indicates
      host-removal failure rather than genuine lysogeny.

  --sensitivity 4.2
      MMseqs2 marker search sensitivity (default). Increase to 5.7–7.5 for
      more distant relatives at proportional time cost.

  --splits 0
      Automatic MMseqs2 memory management (geNomad default).
      NOT a parallelism flag. Increase only if MMseqs2 OOMs.

  --cleanup
      Remove large intermediate files (annotate/ subdir) after classification.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Filtering
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Main tier   : virus_score ≥ min_score (default 0.7, ~97% precision)
  Rescue tier : min_score > score ≥ rescue_min_score (0.4)
                AND length ≥ rescue_min_length_bp (10 kb)
                Novel phages with no close DB relatives score 0.4–0.6
                despite being genuine phages (Camargo et al. 2023).

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Metadata propagation (critical for downstream modules)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  topology
    Source: _virus_summary.tsv → column 'topology'
    Values: DTR · ITR · No terminal repeats · Provirus · NA
    Used by:
      quality.py  → propagates to rename_map.tsv for downstream decisions
      annotate.py → selects circular vs linear plot; gates dnaapler

  genetic_code
    Source: _virus_summary.tsv → column 'genetic_code'
            Fallback: _virus_genes.tsv → mode per contig
    Values: 11 (standard) · 15 (crAss-like: TGA=Trp)
    Used by:
      annotate.py → --genetic_code flag in Pharokka; ensures correct ORF
                    prediction for CrAss-like phages (Crassvirales)
                    Koonin et al. 2020, mBio 11:e00278-20

  taxonomy (family + finest level)
    Source: _virus_summary.tsv → column 'taxonomy'
    Used by:
      quality.py → co-binning at family level for candidate naming;
                   finest level for biologically coherent grouping

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Output files
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  results/04_viral_id/
      {sample_id}_virus.fna          ← viral contigs → quality.py --virus-fna
      {sample_id}_metadata.tsv       ← topology + genetic_code + score per contig

  reports/04_viral_id/{sample_id}/
      genomad/                       ← geNomad output (TSVs after --cleanup)
          {stem}/
              {stem}_virus_summary.tsv
              {stem}_virus_genes.tsv
              {stem}_taxonomy.tsv
      viral_id_summary.tsv
      {sample_id}_viral_id.log
"""

from __future__ import annotations
import csv
import shutil
from collections import Counter
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

STEP  = "04_viral_id"
TOOLS = ["genomad"]


# ── Public entry point ────────────────────────────────────────────────────────

def run(
    cfg:       Config,
    sample_id: str,
    contigs:   Path,
    force:     bool = False,
) -> tuple[Path, Path]:
    """Classify viral contigs with geNomad for a single sample.

    Parameters
    ----------
    contigs : NR contigs FASTA from assembly step

    Returns
    -------
    (virus_fna, metadata_tsv)
        virus_fna    : path to filtered viral FASTA
        metadata_tsv : path to per-contig metadata (topology, genetic_code, score)
    """
    contigs = Path(contigs)
    out_dir  = cfg.results(STEP)
    rpt_dir  = cfg.reports(STEP) / sample_id
    genomad_dir = rpt_dir / "genomad"   # geNomad writes here
    log_f    = rpt_dir / f"{sample_id}_viral_id.log"
    mkdirs(out_dir, rpt_dir, genomad_dir)

    virus_fna    = out_dir / f"{sample_id}_virus.fna"
    metadata_tsv = out_dir / f"{sample_id}_metadata.tsv"

    log_step(f"Module 04 — viral identification  [{sample_id}]")
    log_info(f"  Contigs : {contigs}  ({_count_fasta(contigs):,} sequences)")

    if virus_fna.exists() and virus_fna.stat().st_size > 0 and not force:
        log_info("  Already classified — skipping  (--force to re-run)")
        return virus_fna, metadata_tsv

    db = cfg.databases.genomad
    if not db or not Path(db).exists():
        log_error(
            "  geNomad database not found. "
            "Set databases.genomad in config.yaml or run: "
            "genomad download-database databases/"
        )
        raise FileNotFoundError(f"geNomad database not found: {db}")

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
        task = progress.add_task("Initializing...", total=2)

        # ── Step 1: geNomad classification ────────────────────────────────────
        progress.update(task, description="[1/2] geNomad  — end-to-end classification")
        _run_genomad(contigs, genomad_dir, Path(db), cfg, log_f, force)
        progress.advance(task)

        # ── Step 2: Parse + filter + write outputs ────────────────────────────
        progress.update(task, description="[2/2] Parsing results — topology · genetic code · taxonomy")
        stats = _process_results(
            contigs, genomad_dir, virus_fna, metadata_tsv,
            cfg, sample_id, active_warnings,
        )
        progress.advance(task)

    _save_tsv({**stats, "sample": sample_id}, rpt_dir / "viral_id_summary.tsv")
    _print_completion_panel(sample_id, virus_fna, metadata_tsv,
                            rpt_dir, stats, active_warnings)
    log_step(f"Module 04 completed ✓  [{sample_id}]")
    return virus_fna, metadata_tsv


# ── geNomad runner ────────────────────────────────────────────────────────────

def _run_genomad(
    contigs:     Path,
    genomad_dir: Path,
    db:          Path,
    cfg:         Config,
    log_f:       Path,
    force:       bool,
) -> None:
    """Run geNomad end-to-end on contigs.

    geNomad creates a subdirectory named after the input file stem inside
    genomad_dir. The virus_summary.tsv and virus_genes.tsv produced there
    are the primary outputs.
    """
    stem = contigs.stem
    expected_summary = genomad_dir / stem / f"{stem}_virus_summary.tsv"

    if expected_summary.exists() and not force:
        log_info("  [geNomad] output already exists — skipping")
        return

    cmd = [
        "genomad", "end-to-end",
        str(contigs),
        str(genomad_dir),
        str(db),
        "--threads", str(cfg.threads),
        "--quiet",
        "--cleanup",                   # remove large annotate/ intermediate dir
        "--sensitivity", str(cfg.genomad.sensitivity),
        "--min-score", str(cfg.genomad.min_score),
        "--min-virus-hallmarks", str(cfg.genomad.min_hallmarks),
    ]

    if cfg.genomad.enable_score_calibration:
        cmd += ["--enable-score-calibration",
                "--composition", cfg.genomad.composition]
    if cfg.genomad.lenient_taxonomy:
        cmd.append("--lenient-taxonomy")
    if cfg.genomad.disable_find_proviruses:
        cmd.append("--disable-find-proviruses")
    if cfg.genomad.splits > 0:
        cmd += ["--splits", str(cfg.genomad.splits)]

    log_info(
        f"  [geNomad] calibration={cfg.genomad.enable_score_calibration} "
        f"composition={cfg.genomad.composition} "
        f"lenient-tax={cfg.genomad.lenient_taxonomy} "
        f"sensitivity={cfg.genomad.sensitivity}"
    )
    run_silent(cmd, log_file=log_f)
    log_ok("  [geNomad] classification complete")


# ── Output processing ─────────────────────────────────────────────────────────

def _process_results(
    contigs:      Path,
    genomad_dir:  Path,
    virus_fna:    Path,
    metadata_tsv: Path,
    cfg:          Config,
    sample_id:    str,
    active_warnings: list[str],
) -> dict:
    """Parse geNomad outputs, filter, extract FASTAs, write metadata TSV."""

    # Locate geNomad output files
    virus_summary_f, virus_genes_f = _find_genomad_outputs(genomad_dir)

    # Parse virus_summary.tsv — scores, topology, taxonomy, genetic_code
    viral_hits = _parse_virus_summary(
        virus_summary_f,
        cfg.genomad.min_score,
        cfg.genomad.rescue_min_score,
        cfg.genomad.rescue_min_length_bp,
    )

    if not viral_hits:
        msg = "No viral contigs identified by geNomad."
        log_warn(f"  {msg}")
        active_warnings.append(msg)
        virus_fna.write_text("")
        metadata_tsv.write_text("")
        return {"n_total": "0", "n_main": "0", "n_rescued": "0",
                "n_viral": "0", "top_family": "N/A", "has_crasslike": "False"}

    # Fallback: get genetic_code from _virus_genes.tsv for any contigs
    # where virus_summary didn't report it (geNomad <1.8 may omit the column)
    codes_from_genes = _parse_genetic_codes_from_genes(
        virus_genes_f, set(viral_hits.keys())
    )
    for seq, info in viral_hits.items():
        if info.get("genetic_code", 0) == 0 and seq in codes_from_genes:
            info["genetic_code"] = codes_from_genes[seq]

    # Count tiers
    n_main    = sum(1 for v in viral_hits.values() if v["status"] == "main")
    n_rescued = sum(1 for v in viral_hits.values() if v["status"] == "rescued")
    n_total   = _count_fasta(contigs)

    log_ok(
        f"  [filter] {len(viral_hits):,} viral contigs "
        f"({n_main} main ≥{cfg.genomad.min_score}  "
        f"{n_rescued} rescued)"
    )

    if n_rescued:
        log_info(
            f"  [rescue] {n_rescued} borderline contigs retained "
            f"(score {cfg.genomad.rescue_min_score}–{cfg.genomad.min_score}, "
            f"length ≥{cfg.genomad.rescue_min_length_bp:,} bp)"
        )

    # Extract viral FASTAs from contigs
    n_written = _extract_viral_fastas(contigs, set(viral_hits.keys()), virus_fna)
    log_ok(f"  [extract] {n_written:,} sequences written to {virus_fna.name}")

    # Taxonomy summary
    families  = _collect_families(viral_hits)
    top_family = families[0][0] if families else "Unknown"
    has_crasslike = any(
        int(v.get("genetic_code", 11)) == 15 for v in viral_hits.values()
    )
    if has_crasslike:
        log_info(
            "  [taxonomy] CrAss-like phage detected (genetic_code=15, TGA=Trp) "
            "— Pharokka will use --genetic_code 15"
        )

    # Propagate topology summary to log
    topologies = Counter(v["topology"] for v in viral_hits.values())
    topo_str = "  ".join(f"{t}={n}" for t, n in topologies.most_common())
    log_info(f"  [topology] {topo_str}")

    # Write metadata TSV (critical link to quality.py and annotate.py)
    _write_metadata(viral_hits, metadata_tsv)
    log_ok(f"  [metadata] written → {metadata_tsv.name}")

    naming_levels   = _collect_naming_levels(viral_hits)
    top_naming      = naming_levels[0][0] if naming_levels else "Phage"
    n_known_family  = sum(1 for v in viral_hits.values() if v["family"] != "Unknown")

    return {
        "n_input_contigs": str(n_total),
        "n_viral":         str(len(viral_hits)),
        "n_main":          str(n_main),
        "n_rescued":       str(n_rescued),
        "top_family":      top_family,
        "top_naming_level": top_naming,
        "n_known_family":  str(n_known_family),
        "has_crasslike":   str(has_crasslike),
        "topologies":      topo_str,
    }


# ── geNomad output discovery ──────────────────────────────────────────────────

def _find_genomad_outputs(genomad_dir: Path) -> tuple[Path, Path]:
    """Locate virus_summary and virus_genes TSVs inside geNomad output dir.

    geNomad creates a subdirectory named after the input file stem.
    We use glob to avoid hardcoding the stem.
    """
    summaries = sorted(genomad_dir.rglob("*_virus_summary.tsv"))
    if not summaries:
        raise FileNotFoundError(
            f"geNomad virus_summary.tsv not found in {genomad_dir}. "
            "Check geNomad log for errors."
        )
    virus_summary_f = summaries[0]
    genomad_out     = virus_summary_f.parent
    stem            = virus_summary_f.stem.replace("_virus_summary", "")
    virus_genes_f   = genomad_out / f"{stem}_virus_genes.tsv"
    return virus_summary_f, virus_genes_f


# ── Parsing ───────────────────────────────────────────────────────────────────

def _parse_virus_summary(
    path:                Path,
    min_score:           float,
    rescue_min_score:    float,
    rescue_min_length_bp: int,
) -> dict:
    """Parse geNomad virus_summary.tsv.

    Returns {seq_name: {score, topology, genetic_code, n_hallmarks,
                        taxonomy, family, finest, status, length}}

    Filtering logic
    ---------------
    main     : score ≥ min_score (0.7 by default, ~97% precision)
    rescued  : rescue_min_score ≤ score < min_score AND length ≥ rescue_min_length_bp
               Novel lineages without close DB relatives score 0.4–0.6.
               Rescue prevents losing genuine phages not in geNomad DB.
               Reference: Camargo et al. 2023, Nat Biotechnol 41:1783.
    discarded: everything else
    """
    if not path.exists():
        raise FileNotFoundError(f"virus_summary.tsv not found: {path}")

    hits: dict = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            seq_name = row.get("seq_name", "").strip()
            if not seq_name:
                continue

            try:
                score  = float(row.get("virus_score", 0) or 0)
                length = int(row.get("length", 0) or 0)
            except (ValueError, TypeError):
                continue

            # Determine tier
            if score >= min_score:
                status = "main"
            elif score >= rescue_min_score and length >= rescue_min_length_bp:
                status = "rescued"
            else:
                continue   # below threshold — discard

            topology     = (row.get("topology")     or "NA").strip()
            n_hallmarks  = int(row.get("n_hallmarks", 0) or 0)
            taxonomy_str = (row.get("taxonomy")     or "").strip()
            tax          = _parse_taxonomy(taxonomy_str)

            # genetic_code: may be in virus_summary (geNomad ≥1.8) or absent
            try:
                gc = int(row.get("genetic_code", 11) or 11)
            except (ValueError, TypeError):
                gc = 11

            hits[seq_name] = {
                "score":        score,
                "length":       length,
                "topology":     topology,
                "n_hallmarks":  n_hallmarks,
                "genetic_code": gc,
                "taxonomy":     taxonomy_str,
                "family":       tax["family"],
                "finest":       tax["finest"],
                "naming_level": tax["naming_level"],
                "status":       status,
            }
    return hits


def _parse_genetic_codes_from_genes(path: Path, viral_seqs: set) -> dict:
    """Parse _virus_genes.tsv for genetic_code per contig (mode across genes).

    Fallback for geNomad versions that don't report genetic_code in
    virus_summary.tsv. All genes of a given contig use the same genetic
    code in geNomad, so the mode = the true value.
    """
    if not path.exists():
        return {}
    codes: dict[str, list] = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            seq = row.get("seq_name", "").strip()
            # virus_genes.tsv has format "seq_name|start|stop" for the gene ID;
            # the contig name is the first part
            contig = seq.split("|")[0] if "|" in seq else seq
            if contig not in viral_seqs:
                continue
            try:
                gc = int(row.get("genetic_code", 11) or 11)
            except (ValueError, TypeError):
                gc = 11
            codes.setdefault(contig, []).append(gc)
    return {
        contig: Counter(vals).most_common(1)[0][0]
        for contig, vals in codes.items()
    }


def _parse_taxonomy(taxonomy_str: str) -> dict:
    """Extract family and finest level from geNomad semicolon-separated lineage.

    geNomad taxonomy format (with --lenient-taxonomy):
        Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;
        Caudovirales;Drexlerviridae;Lambdavirus

    ICTV rank suffixes (used to avoid mis-classifying higher ranks as family):
        Realm   : -viria        e.g. Duplodnaviria
        Kingdom : -virae        e.g. Heunggongvirae
        Phylum  : -viricota     e.g. Uroviricota
        Class   : -viricetes    e.g. Caudoviricetes
        Order   : -virales      e.g. Caudovirales
        Family  : -viridae      e.g. Drexlerviridae  ← only this level
        Genus   : -virus        e.g. Lambdavirus

    Family: ONLY levels ending in 'viridae' (not -viria, -viricota, etc.)
            → 'Unknown' when classification does not reach family rank

    Finest: last non-empty level (genus/species with --lenient-taxonomy)
            → used as naming_level fallback in quality.py when family is Unknown

    naming_level: family if known, else finest non-Viruses level, else 'Phage'
                  → candidate names: Drexlerviridae_candidate_001
                                     Caudovirales_candidate_001   (order-level)
                                     Phage_candidate_001          (unclassified)
    """
    if not taxonomy_str:
        return {"family": "Unknown", "finest": "Unknown", "naming_level": "Phage"}

    levels = [l.strip() for l in taxonomy_str.split(";") if l.strip()]

    # Family: strictly -viridae suffix (ICTV standard)
    family = "Unknown"
    for level in reversed(levels):
        if level.endswith("viridae"):
            family = level
            break

    finest = levels[-1] if levels else "Unknown"

    # naming_level: best available rank for candidate naming
    # Priority: family > finest (if not a root/realm) > 'Phage'
    if family != "Unknown":
        naming_level = family
    elif finest not in ("Unknown", "Viruses", ""):
        naming_level = finest
    else:
        naming_level = "Phage"

    return {"family": family, "finest": finest, "naming_level": naming_level}


# ── FASTA extraction ──────────────────────────────────────────────────────────

def _extract_viral_fastas(
    contigs: Path, viral_seqs: set, output: Path
) -> int:
    """Extract viral contigs from assembly FASTA by sequence ID.

    geNomad's seq_name matches the first whitespace-delimited token of the
    FASTA header (i.e. the sequence ID, not the full description line).
    """
    n = 0
    write = False
    with open(contigs) as f_in, open(output, "w") as f_out:
        for line in f_in:
            if line.startswith(">"):
                header_id = line[1:].split()[0].strip()
                write = header_id in viral_seqs
                if write:
                    n += 1
            if write:
                f_out.write(line)
    return n


# ── Metadata TSV ──────────────────────────────────────────────────────────────

_META_HEADERS = [
    "seq_name", "length", "virus_score", "n_hallmarks",
    "topology", "genetic_code", "family", "finest", "naming_level", "taxonomy", "status",
]


def _write_metadata(hits: dict, path: Path) -> None:
    """Write per-contig metadata TSV.

    This file is the critical data link between viral_id → quality → annotate:
      topology     → rename_map.tsv topology column
      genetic_code → Pharokka --genetic_code
      family       → candidate naming ({Family}_candidate_N)
      finest       → co-binning granularity (--lenient-taxonomy)
    """
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=_META_HEADERS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for seq, info in sorted(hits.items(), key=lambda x: -x[1]["score"]):
            writer.writerow({
                "seq_name":    seq,
                "length":      info["length"],
                "virus_score": f"{info['score']:.4f}",
                "n_hallmarks": info["n_hallmarks"],
                "topology":    info["topology"],
                "genetic_code": info["genetic_code"],
                "family":      info["family"],
                "finest":      info["finest"],
                "naming_level": info.get("naming_level", info["family"] if info["family"] != "Unknown" else "Phage"),
                "taxonomy":    info["taxonomy"],
                "status":      info["status"],
            })


# ── Helpers ───────────────────────────────────────────────────────────────────

def _count_fasta(path: Path) -> int:
    """Count sequences in a FASTA file."""
    if not path.exists() or path.stat().st_size == 0:
        return 0
    n = 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                n += 1
    return n


def _collect_families(hits: dict) -> list[tuple[str, int]]:
    """Return sorted list of (family, count) by frequency."""
    return Counter(v["family"] for v in hits.values()).most_common()


def _collect_naming_levels(hits: dict) -> list[tuple[str, int]]:
    """Return sorted list of (naming_level, count) by frequency.

    naming_level is the best available rank for candidate naming:
    family if known, finest otherwise, 'Phage' if unclassified.
    """
    return Counter(v.get("naming_level", "Phage") for v in hits.values()).most_common()


# ── Completion panel ──────────────────────────────────────────────────────────

def _print_completion_panel(
    sample_id:    str,
    virus_fna:    Path,
    metadata_tsv: Path,
    rpt_dir:      Path,
    stats:        dict,
    active_warnings: list[str],
) -> None:
    def _f(val) -> str:
        return "[dim]N/A[/dim]" if str(val) in ("", "None", "N/A") else str(val)

    n_input  = _f(stats.get("n_input_contigs", "N/A"))
    n_viral  = stats.get("n_viral",   "0")
    n_main   = stats.get("n_main",    "0")
    n_resc   = stats.get("n_rescued", "0")
    family   = _f(stats.get("top_family", "N/A"))
    topos    = _f(stats.get("topologies", "N/A"))
    crasslike = stats.get("has_crasslike", "False") == "True"

    try:
        pct = f"{int(n_viral) / int(n_input) * 100:.1f}%"
    except (ValueError, ZeroDivisionError):
        pct = "N/A"

    # Collect naming_level info for display
    top_naming = _f(stats.get("top_naming_level", family))
    family_known = stats.get("n_known_family", "0") != "0"

    lines = [
        "[bold]Viral contigs[/bold]",
        f"  [cyan]Input contigs  :[/cyan] {n_input}",
        f"  [cyan]Viral retained :[/cyan] [bold green]{n_viral}[/bold green] "
        f"({pct})  [dim]main={n_main}  rescued={n_resc}[/dim]",
        f"  [cyan]Top family     :[/cyan] {family}"
        + ("" if family_known else "  [dim](family unresolved — using finest level for naming)[/dim]"),
        f"  [cyan]Naming level   :[/cyan] {top_naming}",
        f"  [cyan]Topology       :[/cyan] {topos}",
    ]
    if crasslike:
        lines.append(
            "  [yellow]⚡ CrAss-like detected[/yellow]  "
            "[dim]genetic_code=15 → Pharokka --genetic_code 15[/dim]"
        )
    lines += [
        "",
        "[bold]Output files[/bold]",
        f"  Viral FASTA : {virus_fna}",
        f"  Metadata    : {metadata_tsv}",
        f"  Summary     : {rpt_dir / 'viral_id_summary.tsv'}",
    ]
    if active_warnings:
        lines += ["", "[bold yellow]Warnings[/bold yellow]"]
        lines += [f"  [yellow]⚠ {w}[/yellow]" for w in active_warnings]

    console.print(Panel(
        "\n".join(lines),
        title=f"[bold cyan]Viral ID complete — {sample_id}[/bold cyan]",
        border_style="cyan", padding=(0, 2), width=120,
    ))


# ── TSV summary ───────────────────────────────────────────────────────────────

_TSV_HEADERS = [
    "sample", "n_input_contigs", "n_viral", "n_main", "n_rescued",
    "top_family", "top_naming_level", "n_known_family", "has_crasslike", "topologies",
]


def _save_tsv(row: dict, path: Path) -> None:
    rows: dict[str, dict] = {}
    if path.exists():
        with open(path) as f:
            old_h = f.readline().rstrip("\n").split("\t")
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if cols and cols[0]:
                    rows[cols[0]] = dict(zip(old_h, cols))
    rows[row["sample"]] = {h: str(row.get(h, "")) for h in _TSV_HEADERS}
    with open(path, "w") as f:
        f.write("\t".join(_TSV_HEADERS) + "\n")
        for r in rows.values():
            f.write("\t".join(r.get(h, "") for h in _TSV_HEADERS) + "\n")

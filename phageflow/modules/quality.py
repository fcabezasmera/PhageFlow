"""PhageFlow Module 05 — Genome quality assessment and selection (CheckV).

Goal: evaluate viral contigs from viral_id, assign quality tiers based on
completeness and genome integrity, dereplicate near-identical genomes, and
produce annotation-ready FASTA files for annotate.py.

Tool: CheckV end-to-end (Nayfach et al. 2021, Nat Biotechnol 39:578)
  CheckV estimates genome completeness via two complementary methods:
  AAI-based    : amino acid identity against curated viral genome DB.
                 Most precise when a close reference exists in the DB.
  HMM-based    : gene-content lower bound using viral protein HMMs.
                 Used for novel lineages without close DB references.
                 Tends to UNDERESTIMATE completeness for divergent viruses;
                 a contig reporting 40-49% completeness via HMM may be
                 genuinely complete but divergent from all DB entries.

Quality tiers (priority order)
  complete     : CheckV = Complete (DTR or ITR detected in sequence)
                 → annotation_ready/phages/
  high-quality : completeness ≥ 90%
                 → annotation_ready/phages/
  medium-quality: completeness ≥ min_completeness (default 50%)
                 MIUViG standard for medium-quality draft genomes.
                 Roux et al. 2019 (Nat Biotechnol 37:505).
                 → annotation_ready/phages/
  large-nd     : Not-determined, length ≥ large_nd_rescue_bp (20 kb),
                 ≥ 1 viral gene. Captures large divergent phages with no
                 CheckV reference. Jumbo phages (>200 kb) frequently fall
                 in this tier (Al-Shayeb et al. 2020, Nature 578:425).
                 → annotation_ready/phages/
  lq-draft     : length ≥ length_rescue (10 kb), ≥ min_viral_genes (1).
                 Candidate for co-binning by taxonomy.
                 → drafts/ or bin-rescue
  bin-rescue   : LQ drafts of the same naming_level with combined
                 length ≥ min_bin_rescue_bp (30 kb). Captures fragmented
                 genomes where no single contig meets quality thresholds.
                 → annotation_ready/phages/ (multi-contig FASTA)
  discarded    : length < min_contig_bp or no viral genes.

Topology resolution (CheckV > geNomad, confirmed > predicted)
  CheckV termini_type takes precedence over geNomad topology prediction:
  CheckV confirms DTR/ITR by detecting ≥20 bp terminal repeats in the
  actual assembled sequence. geNomad predicts topology from sequence
  composition, which is less reliable for novel lineages.
  Priority: CheckV DTR → DTR (circular, most tailed dsDNA phages)
            CheckV ITR → ITR (linear with ITR, e.g. T7-like phages)
            CheckV NA  → geNomad topology or "No terminal repeats"
  Downstream: DTR → circular genome map + dnaapler reorientation
              ITR → linear genome map + dnaapler terL detection
              other → linear genome map

Dereplication strategy (ANI-based, strain-level)
  Applied to single-contig annotation_ready candidates only.
  Large genomes (≥20 kb): blastn all-vs-all (ANI ≥ 95%, coverage ≥ 80%).
    blastn detects near-identical non-overlapping fragments that mash
    misses due to disjoint k-mer content (see quality.py BLAST discussion).
  Small genomes (<20 kb): minimap2 asm5 all-vs-all (identity ≥ 95%,
    max coverage ≥ 80%). mash is unreliable at this size range due to
    sparse k-mer sketch sampling.
  Fallbacks: mash → cd-hit-est if blastn/minimap2 unavailable.
  Threshold: 95% ANI = ICTV species boundary (Turner et al. 2021,
    Arch Virol 166:2633). Clusters at 95% retain strain-level diversity.
  Within each cluster: longest contig is kept as representative.

Special case warnings (never discards)
  Concatemers : kmer_freq > 1.5 OR genome_copies > 1.5.
                High-coverage lytic replication generates DNA concatemers.
                The sequence is valid; the contig may represent N linked
                genome copies. Common in purified phage preparations.
  Host genes  : host_genes > 0 indicates possible residual contamination
                that survived host removal. Inspect contamination.tsv.

Coverage validation (read recruitment)
  Post-dereplication, post-host-removal reads are mapped back to
  annotation_ready candidates using bwa-mem2.
  Metrics: mean_depth, breadth (≥1x positions), CV (depth uniformity).
  breadth < 95%: uncovered region — possible assembly gap or chimera.
  CV > 1.5     : uneven depth — possible concatemer or chimeric assembly.
  Roux et al. 2019 (eLife 8:e42923) — read recruitment as genome integrity
  validation; Nayfach et al. 2021 (Nat Biotechnol 39:578).

Output layout
  results/05_quality/{sample_id}/
    annotation_ready/phages/   {naming_level}_candidate_{N:03d}.fasta
    annotation_ready/proviruses/
    drafts/                    {naming_level}_draft_{N:03d}.fasta
    checkv/                    raw CheckV output
  reports/05_quality/{sample_id}/
    rename_map.tsv             contig→candidate + all metadata (critical link)
    quality_summary.tsv        module-level statistics
    coverage/                  read recruitment BAM + depth TSV
"""
from __future__ import annotations
import csv
import math
import shutil
import subprocess
import tempfile
from collections import Counter, defaultdict
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

STEP  = "05_quality"
TOOLS = ["checkv"]

# Quality tier constants
_TIER_COMPLETE  = "complete"
_TIER_HQ        = "high-quality"
_TIER_MQ        = "medium-quality"
_TIER_LARGE_ND  = "large-nd"
_TIER_LQ_DRAFT  = "lq-draft"
_TIER_BIN       = "bin-rescue"
_TIER_DEREP     = "dereplicated"
_TIER_DISCARD   = "discarded"

# Tiers that go to annotation_ready (before dereplication / co-binning)
_ANNOTATION_TIERS = {_TIER_COMPLETE, _TIER_HQ, _TIER_MQ, _TIER_LARGE_ND}

# CheckV quality column values
_CKV_COMPLETE = "Complete"
_CKV_HQ       = "High-quality"
_CKV_MQ       = "Medium-quality"
_CKV_LQ       = "Low-quality"
_CKV_ND       = "Not-determined"

# Mash ANI dereplication threshold (>98% ANI → same strain)
# ICTV species boundary is 95% ANI. Turner et al. 2021, Arch Virol 166:2633.
_MASH_DIST_THRESH = 0.02

# BLAST-based dereplication thresholds (used for large genomes, more sensitive than mash)
# 95% ANI = ICTV species boundary (Turner et al. 2021, Arch Virol 166:2633)
# 80% coverage = minimum overlap to consider same genome
_BLAST_ANI_THRESH = 95.0
_BLAST_COV_THRESH = 0.80

# Genome size threshold for dereplication method selection.
# < threshold : minimap2 asm5 (mash unreliable on small genomes — sparse sketch)
# ≥ threshold : mash (fast pairwise on full genome k-mer content)
_SMALL_GENOME_THRESH = 20_000

# ── Public entry point ────────────────────────────────────────────────────────

def run(
    cfg:          Config,
    sample_id:    str,
    virus_fna:    Path,
    force:        bool = False,
    metadata_tsv: Optional[Path] = None,
) -> tuple[Path, Path]:
    """Assess quality, dereplicate, and select annotation-ready viral genomes.

    Returns (ann_phages, drafts_dir) paths.
    """
    virus_fna = Path(virus_fna)
    out_dir   = cfg.results(STEP) / sample_id
    rpt_dir   = cfg.reports(STEP) / sample_id
    checkv_dir = out_dir / "checkv"
    ann_phages = out_dir / "annotation_ready" / "phages"
    ann_prov   = out_dir / "annotation_ready" / "proviruses"
    drafts_dir = out_dir / "drafts"
    log_f      = rpt_dir / f"{sample_id}_quality.log"
    mkdirs(out_dir, rpt_dir, checkv_dir, ann_phages, ann_prov, drafts_dir)

    log_step(f"Module 05 — quality  [{sample_id}]")
    log_info(f"  Viral FASTA : {virus_fna}")

    if not virus_fna.exists() or virus_fna.stat().st_size == 0:
        log_warn("  Viral FASTA is empty — no candidates to evaluate.")
        return ann_phages, drafts_dir

    # Auto-discover metadata from viral_id step if not provided
    if metadata_tsv is None:
        metadata_tsv = cfg.results("04_viral_id") / f"{sample_id}_metadata.tsv"
    if not Path(metadata_tsv).exists():
        log_warn(
            f"  viral_id metadata not found at {metadata_tsv}. "
            "Topology and genetic_code will use CheckV/defaults."
        )
        metadata_tsv = None

    # Check if annotation_ready already has outputs (for idempotency)
    existing = list(ann_phages.glob("*.fasta"))
    if existing and not force:
        log_info(
            f"  {len(existing)} candidate(s) already in annotation_ready — "
            "skipping  (--force to re-run)"
        )
        return ann_phages, drafts_dir

    # Clear output dirs if re-running
    if force:
        for d in [ann_phages, ann_prov, drafts_dir]:
            shutil.rmtree(d, ignore_errors=True)
            d.mkdir(parents=True, exist_ok=True)

    db = cfg.databases.checkv
    if not db or not Path(db).exists():
        log_error(
            "  CheckV database not found. "
            "Set databases.checkv in config.yaml or run: "
            "checkv download_database databases/checkv_db"
        )
        raise FileNotFoundError(f"CheckV database not found: {db}")

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
        task = progress.add_task("Initializing...", total=6)

        # ── Step 1: CheckV ────────────────────────────────────────────────────
        progress.update(task, description="[1/6] CheckV  — completeness + quality assessment")
        _run_checkv(virus_fna, checkv_dir, Path(db), cfg.threads, log_f, force)
        progress.advance(task)

        # ── Step 2: Load results + assign tiers ───────────────────────────────
        progress.update(task, description="[2/6] Parsing quality results + assigning tiers")
        checkv_rows = _load_checkv_results(checkv_dir)
        metadata    = _load_metadata(metadata_tsv) if metadata_tsv else {}
        seqs        = _load_fasta(virus_fna)

        contigs     = _assign_tiers(checkv_rows, metadata, cfg)
        _warn_special_cases(contigs, active_warnings, cfg)
        _warn_naming_level_redundancy(contigs, active_warnings)
        progress.advance(task)

        # ── Step 3: Dereplication ─────────────────────────────────────────────
        progress.update(task, description="[3/6] Dereplicating annotation_ready candidates")
        contigs = _dereplicate(contigs, seqs, cfg.threads, active_warnings)
        progress.advance(task)

        # ── Step 4: Co-bin LQ drafts by taxon ─────────────────────────────────
        progress.update(task, description="[4/6] Co-binning draft contigs by taxonomy")
        contigs = _cobin_drafts(contigs, cfg.checkv.min_bin_rescue_bp)
        progress.advance(task)

        # ── Step 5: Write FASTAs + rename_map ─────────────────────────────────
        progress.update(task, description="[5/6] Writing candidate FASTAs + rename_map")
        stats = _write_outputs(
            contigs, seqs, ann_phages, ann_prov, drafts_dir,
            rpt_dir, sample_id,
        )
        progress.advance(task)

        # ── Step 6: Read recruitment + coverage validation ────────────────────
        progress.update(task, description="[6/6] Read recruitment — coverage validation")
        coverage = _read_recruitment_coverage(
            list(ann_phages.glob("*.fasta")),
            sample_id, cfg, rpt_dir, log_f, active_warnings,
        )
        # Annotate rename_map with coverage metrics if available
        if coverage:
            _annotate_rename_map_coverage(rpt_dir / "rename_map.tsv", coverage)
        progress.advance(task)

    _save_tsv({**stats, "sample": sample_id}, rpt_dir / "quality_summary.tsv")
    _print_completion_panel(
        sample_id, ann_phages, drafts_dir, rpt_dir, stats, active_warnings,
    )
    log_step(f"Module 05 completed ✓  [{sample_id}]")
    return ann_phages, drafts_dir

# ── CheckV ────────────────────────────────────────────────────────────────────

def _run_checkv(
    virus_fna: Path, checkv_dir: Path, db: Path,
    threads: int, log_f: Path, force: bool,
) -> None:
    """Run CheckV end_to_end.

    CheckV estimates genome completeness via:
    1. AAI alignment against a curated viral genome database
    2. HMM-based gene content (lower bound for novel viruses)
    Nayfach et al. 2021, Nat Biotechnol 39:578.
    """
    summary = checkv_dir / "quality_summary.tsv"
    if summary.exists() and not force:
        log_info("  [CheckV] output already exists — skipping")
        return
    run_silent(
        [
            "checkv", "end_to_end",
            str(virus_fna),
            str(checkv_dir),
            "-d", str(db),
            "-t", str(threads),
            "--remove_tmp",
        ],
        log_file=log_f,
    )
    log_ok("  [CheckV] quality assessment complete")

def _load_checkv_results(checkv_dir: Path) -> list[dict]:
    """Load CheckV quality_summary.tsv."""
    summary = checkv_dir / "quality_summary.tsv"
    if not summary.exists():
        raise FileNotFoundError(
            f"CheckV quality_summary.tsv not found in {checkv_dir}."
        )
    rows = []
    with open(summary, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(dict(row))
    return rows

def _load_metadata(path: Path) -> dict:
    """Load viral_id metadata TSV as {seq_name: {topology, genetic_code, ...}}."""
    meta = {}
    if not path or not Path(path).exists():
        return meta
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            meta[row["seq_name"]] = dict(row)
    return meta

def _load_fasta(path: Path) -> dict[str, str]:
    """Load FASTA into {seq_id: sequence} dict."""
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

# ── Tier assignment ───────────────────────────────────────────────────────────

def _assign_tiers(
    checkv_rows: list[dict],
    metadata:    dict,
    cfg:         Config,
) -> dict:
    """Assign quality tier to each contig.

    Returns {seq_name: info_dict} with tier, topology, genetic_code, etc.
    """
    result: dict = {}

    for row in checkv_rows:
        seq_name = row.get("contig_id", "").strip()
        if not seq_name:
            continue

        # Numeric fields with safe defaults
        def _int(key, default=0):
            v = row.get(key, "") or ""
            try: return int(float(v))
            except (ValueError, TypeError): return default

        def _float(key, default=0.0):
            v = row.get(key, "") or ""
            try: return float(v)
            except (ValueError, TypeError): return default

        length           = _int("contig_length")
        viral_genes      = _int("viral_genes")
        host_genes       = _int("host_genes")
        completeness     = _float("completeness")
        kmer_freq        = _float("kmer_freq")
        genome_copies    = _float("genome_copies")
        checkv_quality   = (row.get("checkv_quality") or "").strip()
        termini_type     = (row.get("termini_type")   or "").strip()
        comp_method      = (row.get("completeness_method") or "").strip()
        warnings_str     = (row.get("warnings")       or "").strip()

        # Topology: CheckV (confirmed in sequence) > geNomad (predicted)
        meta = metadata.get(seq_name, {})
        if termini_type == "DTR":
            topology = "DTR"
        elif termini_type == "ITR":
            topology = "ITR"
        else:
            # geNomad topology as fallback context
            genomad_topo = meta.get("topology", "") or ""
            topology = genomad_topo if genomad_topo else "No terminal repeats"

        genetic_code = int(meta.get("genetic_code", 11) or 11)
        naming_level = (meta.get("naming_level") or "Phage").strip() or "Phage"
        virus_score  = _float_safe(meta.get("virus_score", "0"))

        # ── Tier assignment (priority order) ──────────────────────────────────
        if checkv_quality == _CKV_COMPLETE:
            tier = _TIER_COMPLETE
        elif checkv_quality == _CKV_HQ or completeness >= 90.0:
            tier = _TIER_HQ
        elif checkv_quality == _CKV_MQ or (
            completeness >= cfg.checkv.min_completeness and completeness > 0
        ):
            tier = _TIER_MQ
        elif (
            checkv_quality in (_CKV_ND, _CKV_LQ, "")
            and length >= cfg.checkv.large_nd_rescue_bp
            and viral_genes >= 1
        ):
            # Large ND rescue: captures jumbo phages without CheckV references
            tier = _TIER_LARGE_ND
        elif (
            length >= cfg.checkv.min_contig_bp
            and length >= cfg.checkv.length_rescue
            and viral_genes >= cfg.checkv.min_viral_genes
        ):
            tier = _TIER_LQ_DRAFT
        else:
            tier = _TIER_DISCARD

        result[seq_name] = {
            "tier":               tier,
            "topology":           topology,
            "genetic_code":       genetic_code,
            "naming_level":       naming_level,
            "virus_score":        virus_score,
            "checkv_quality":     checkv_quality,
            "completeness":       completeness,
            "completeness_method": comp_method,
            "contig_length":      length,
            "viral_genes":        viral_genes,
            "host_genes":         host_genes,
            "kmer_freq":          kmer_freq,
            "genome_copies":      genome_copies,
            "termini_type":       termini_type,
            "checkv_warnings":    warnings_str,
            # candidate_id assigned later
            "candidate_id":       "",
        }

    return result

def _float_safe(val) -> float:
    try: return float(val)
    except (ValueError, TypeError): return 0.0

# ── Special case warnings ─────────────────────────────────────────────────────

def _warn_special_cases(
    contigs: dict, active_warnings: list[str], cfg: Config,
) -> None:
    """WARN on concatemers and host gene contamination. Never discard."""
    n_concat = 0
    n_host   = 0
    for seq, info in contigs.items():
        if info["tier"] == _TIER_DISCARD:
            continue
        kmer_freq     = info["kmer_freq"]
        genome_copies = info["genome_copies"]
        host_genes    = info["host_genes"]
        if (kmer_freq > cfg.checkv.max_kmer_freq or
                genome_copies > cfg.checkv.max_genome_copies):
            n_concat += 1
        if host_genes > 0:
            n_host += 1

    if n_concat:
        msg = (
            f"{n_concat} contig(s) have kmer_freq>{cfg.checkv.max_kmer_freq} "
            f"or genome_copies>{cfg.checkv.max_genome_copies} — probable "
            "concatemer(s). Sequence is valid; contig may represent multiple "
            "linked genome copies. NOT discarded."
        )
        log_warn(f"  {msg}")
        active_warnings.append(msg)
    if n_host:
        msg = (
            f"{n_host} contig(s) have host_genes > 0 — possible residual "
            "host contamination. Inspect contamination.tsv. NOT discarded."
        )
        log_warn(f"  {msg}")
        active_warnings.append(msg)

# ── Dereplication ─────────────────────────────────────────────────────────────

def _dereplicate(
    contigs: dict,
    seqs:    dict[str, str],
    threads: int,
    active_warnings: list[str],
) -> dict:
    """Dereplicate annotation_ready single-contig candidates at 98% ANI.

    Operates on single-contig annotation_ready contigs only.
    Multi-contig bins and drafts are not dereplicated.

    Method: mash distance (sketch size 10,000, ANI > 98% = dist < 0.02).
    Keeps the longest contig in each cluster.
    Fallback: cd-hit-est at 98% if mash is not available.

    Turner et al. 2021, Arch Virol 166:2633 — 95% ANI = ICTV species boundary.
    """
    singles = [
        seq for seq, info in contigs.items()
        if info["tier"] in _ANNOTATION_TIERS
    ]
    if len(singles) < 2:
        return contigs   # nothing to dereplicate

    # Safety guard: dereplication is designed for a handful of HQ candidates
    # from a purified prep, not for hundreds of short assembly fragments.
    # If there are many candidates, likely a fragmented assembly — skip
    # dereplication to avoid O(N²) mash cost and confusing results.
    # The assembly/host-removal step should be revisited in this case.
    if len(singles) > 200:
        msg = (
            f"{len(singles)} annotation_ready candidates — too many for "
            "pairwise dereplication (expected <20 for purified phage prep). "
            "Likely a fragmented assembly. Dereplication skipped; "
            "review host-removal mode and assembly quality."
        )
        log_warn(f"  [derep] {msg}")
        active_warnings.append(msg)
        return contigs

    # Split by genome size: mash is unreliable on small genomes (<20 kb)
    # because the k-mer sketch undersamples the sequence space.
    # minimap2 asm5 performs actual sequence alignment and is accurate at any size.
    large = [s for s in singles if contigs[s]["contig_length"] >= _SMALL_GENOME_THRESH]
    small = [s for s in singles if contigs[s]["contig_length"] <  _SMALL_GENOME_THRESH]

    clusters: list[list] = []
    has_mash    = bool(shutil.which("mash"))
    has_minimap = bool(shutil.which("minimap2"))

    if len(large) >= 2:
        has_blast = bool(shutil.which("blastn") and shutil.which("makeblastdb"))
        if has_blast:
            clusters += _cluster_blastn(large, seqs, threads)
        elif has_mash:
            log_warn("  [derep] blastn not found — falling back to mash for large genomes.")
            active_warnings.append("blastn not found; used mash for large-genome dereplication (less sensitive)")
            clusters += _cluster_mash(large, seqs, threads)
        else:
            log_warn("  [derep] blastn/mash not found — using cd-hit-est for large genomes.")
            active_warnings.append("blastn/mash not found; used cd-hit-est for large-genome dereplication")
            clusters += _cluster_cdhit(large, seqs, threads)

    if len(small) >= 2:
        if has_minimap:
            log_info(
                f"  [derep] {len(small)} small genome(s) (<{_SMALL_GENOME_THRESH//1000} kb) "
                "— using minimap2 asm5 (mash unreliable at this size)"
            )
            clusters += _cluster_minimap2(small, seqs, threads)
        else:
            log_warn("  [derep] minimap2 not found — using cd-hit-est for small genomes.")
            active_warnings.append("minimap2 not found; used cd-hit-est for small-genome dereplication")
            clusters += _cluster_cdhit(small, seqs, threads)

    n_before = len(singles)
    n_removed = 0
    for cluster in clusters:
        if len(cluster) < 2:
            continue
        keep = max(cluster, key=lambda s: contigs[s]["contig_length"])
        for seq in cluster:
            if seq != keep:
                contigs[seq]["tier"] = _TIER_DEREP
                n_removed += 1

    if n_removed:
        log_info(
            f"  [derep] {n_before} → {n_before - n_removed} candidates "
            f"({n_removed} removed at ANI>{1 - _MASH_DIST_THRESH:.0%})"
        )
    return contigs

class _UnionFind:
    def __init__(self, ids):
        self.parent = {i: i for i in ids}

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, x, y):
        px, py = self.find(x), self.find(y)
        if px != py:
            self.parent[px] = py

    def groups(self) -> list[list]:
        g: dict[str, list] = defaultdict(list)
        for x in self.parent:
            g[self.find(x)].append(x)
        return list(g.values())

def _cluster_mash(singles: list, seqs: dict, threads: int) -> list[list]:
    """Cluster sequences using mash at MASH_DIST_THRESH."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        # Write individual FASTAs for mash -i (per-sequence sketch)
        fasta_combined = tmp / "candidates.fasta"
        with open(fasta_combined, "w") as f:
            for seq in singles:
                if seq in seqs:
                    f.write(f">{seq}\n{seqs[seq]}\n")

        # mash sketch -o PREFIX always appends .msh → final file = PREFIX.msh
        # Pass prefix WITHOUT .msh so the output lands at sketch.msh, not sketch.msh.msh
        sketch_prefix = tmp / "sketch"
        sketch_file   = tmp / "sketch.msh"
        subprocess.run(
            ["mash", "sketch", "-s", "10000", "-i",
             "-o", str(sketch_prefix), str(fasta_combined)],
            capture_output=True, check=True, timeout=120,
        )
        if not sketch_file.exists():
            raise FileNotFoundError(
                f"mash sketch did not produce {sketch_file}. "
                "Check that mash is installed and the input FASTA is valid."
            )
        dist_out = tmp / "distances.tsv"
        subprocess.run(
            ["mash", "dist", "-p", str(min(threads, 8)),
             str(sketch_file), str(sketch_file)],
            stdout=open(dist_out, "w"), capture_output=False,
            check=True, timeout=300,
        )

        uf = _UnionFind(singles)
        with open(dist_out) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 3:
                    continue
                ref, qry = parts[0], parts[1]
                try:
                    dist = float(parts[2])
                except ValueError:
                    continue
                # Extract seq ID from mash output (may include path)
                ref_id = Path(ref).stem
                qry_id = Path(qry).stem
                if ref_id in singles and qry_id in singles and dist < _MASH_DIST_THRESH:
                    uf.union(ref_id, qry_id)

    return uf.groups()

def _cluster_blastn(singles: list, seqs: dict, threads: int) -> list[list]:
    """Cluster sequences using blastn all-vs-all alignment.

    More sensitive than mash for detecting near-identical sequences,
    including non-overlapping fragments of the same genome whose k-mer
    content is disjoint (mash would report high distance for these).

    Two genomes are considered the same strain if:
      identity (pident) >= 95%  AND  coverage >= 80% of shorter sequence

    The coverage criterion handles assembly artifacts at termini where
    one contig is slightly longer than another but represents the same genome.

    Turner et al. 2021, Arch Virol 166:2633 — 95% ANI = ICTV species boundary.
    """
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        fasta = tmp / "candidates.fasta"
        with open(fasta, "w") as f:
            for seq in singles:
                if seq in seqs:
                    f.write(f">{seq}\n{seqs[seq]}\n")

        db = tmp / "blastdb"
        try:
            subprocess.run(
                ["makeblastdb", "-in", str(fasta), "-dbtype", "nucl",
                 "-out", str(db), "-parse_seqids"],
                capture_output=True, check=True, timeout=60,
            )
        except Exception:
            log_warn("  [derep] makeblastdb failed — falling back to cd-hit-est")
            return _cluster_cdhit(singles, seqs, threads)

        blast_out = tmp / "blast.tsv"
        try:
            subprocess.run(
                ["blastn", "-query", str(fasta), "-db", str(db),
                 "-outfmt", "6 qseqid sseqid pident length qlen slen evalue bitscore",
                 "-evalue", "1e-10",
                 "-num_threads", str(min(threads, 8)),
                 "-out", str(blast_out)],
                capture_output=True, check=True, timeout=600,
            )
        except Exception:
            log_warn("  [derep] blastn failed — falling back to cd-hit-est")
            return _cluster_cdhit(singles, seqs, threads)

        uf = _UnionFind(singles)
        if blast_out.exists():
            with open(blast_out) as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) < 6:
                        continue
                    try:
                        qseqid = parts[0]
                        sseqid = parts[1]
                        pident = float(parts[2])
                        aln_len = int(parts[3])
                        qlen   = int(parts[4])
                        slen   = int(parts[5])
                    except (ValueError, IndexError):
                        continue

                    if qseqid == sseqid:
                        continue
                    if qseqid not in singles or sseqid not in singles:
                        continue

                    min_len  = min(qlen, slen)
                    coverage = aln_len / min_len if min_len > 0 else 0.0

                    if pident >= _BLAST_ANI_THRESH and coverage >= _BLAST_COV_THRESH:
                        uf.union(qseqid, sseqid)

        return uf.groups()


def _cluster_minimap2(singles: list, seqs: dict, threads: int) -> list[list]:
    """Dereplicate small genomes (<20 kb) using minimap2 all-vs-all alignment.

    mash is unreliable for genomes below ~15–20 kb because the 10,000-hash
    sketch over-samples the entire k-mer space, producing inaccurate distance
    estimates. minimap2 asm5 performs actual sequence alignment and is accurate
    at any genome size, including circular ssDNA phages (Inoviridae, Microviridae)
    and small tailed phages.

    Two genomes are considered the same if:
      identity ≥ 95%  AND  max(coverage_query, coverage_target) ≥ 80%

    The coverage criterion handles cases where one genome is longer due to
    assembly artifacts at the termini — as long as 80% of one genome aligns
    to the other at high identity, they are the same phage.
    """
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        fasta = tmp / "small.fasta"
        with open(fasta, "w") as f:
            for seq in singles:
                if seq in seqs:
                    f.write(f">{seq}\n{seqs[seq]}\n")

        paf = tmp / "alignments.paf"
        # asm5: accurate alignment for assemblies with <5% divergence.
        # -c: output CIGAR string (needed for accurate match count).
        result = subprocess.run(
            ["minimap2", "-c", "-x", "asm5",
             "-t", str(min(threads, 8)),
             str(fasta), str(fasta)],
            stdout=open(paf, "w"),
            stderr=subprocess.DEVNULL,
            timeout=180,
        )
        if result.returncode != 0 or not paf.exists():
            log_warn("  [derep] minimap2 failed — falling back to cd-hit-est for small genomes")
            return _cluster_cdhit(singles, seqs, threads)

        uf = _UnionFind(singles)
        with open(paf) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 11:
                    continue
                try:
                    qname  = parts[0]; qlen = int(parts[1])
                    qstart = int(parts[2]); qend = int(parts[3])
                    tname  = parts[5]; tlen = int(parts[6])
                    tstart = int(parts[7]); tend = int(parts[8])
                    matches   = int(parts[9])
                    block_len = int(parts[10])
                except (ValueError, IndexError):
                    continue

                if qname == tname or block_len == 0:
                    continue

                identity = matches / block_len
                cov_q    = (qend - qstart) / qlen if qlen else 0
                cov_t    = (tend - tstart) / tlen if tlen else 0

                if identity >= 0.95 and max(cov_q, cov_t) >= 0.80:
                    if qname in singles and tname in singles:
                        uf.union(qname, tname)

        return uf.groups()

def _cluster_cdhit(singles: list, seqs: dict, threads: int) -> list[list]:
    """Fallback: cd-hit-est at 98% for dereplication."""
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        inp = tmp / "input.fasta"
        out = tmp / "nr"
        with open(inp, "w") as f:
            for seq in singles:
                if seq in seqs:
                    f.write(f">{seq}\n{seqs[seq]}\n")
        subprocess.run(
            ["cd-hit-est", "-i", str(inp), "-o", str(out),
             "-c", "0.98", "-aS", "0.85", "-n", "8",
             "-T", str(threads), "-M", "0", "-d", "0"],
            capture_output=True, check=False, timeout=300,
        )
        clstr = Path(str(out) + ".clstr")
        if not clstr.exists():
            return [[s] for s in singles]

        clusters: list[list] = []
        current: list[str] = []
        with open(clstr) as f:
            for line in f:
                if line.startswith(">Cluster"):
                    if current:
                        clusters.append(current)
                    current = []
                else:
                    # Extract seq ID: ">seq_id..."
                    parts = line.strip().split(">")
                    if len(parts) > 1:
                        seq_id = parts[1].split("...")[0].strip()
                        if seq_id in singles:
                            current.append(seq_id)
        if current:
            clusters.append(current)
        return clusters or [[s] for s in singles]

# ── Co-binning ────────────────────────────────────────────────────────────────

def _cobin_drafts(contigs: dict, min_bin_bp: int) -> dict:
    """Co-bin LQ draft contigs of the same naming_level.

    Groups lq-draft contigs by naming_level. If combined length ≥ min_bin_bp,
    upgrades all contigs in the group to bin-rescue (annotation_ready).
    Otherwise they remain as lq-draft (→ drafts/).

    Contigs with naming_level = 'Phage' (unclassified) are co-binned together
    only if all drafts are unclassified; this prevents merging genuinely
    distinct phages that share the 'Phage' label due to absent taxonomy.
    """
    drafts = {
        seq: info for seq, info in contigs.items()
        if info["tier"] == _TIER_LQ_DRAFT
    }
    if not drafts:
        return contigs

    # Naming levels that already have annotation_ready singles.
    # Bin-rescue is suppressed for these: the LQ drafts are redundant
    # fragments of the same genome, not a distinct new candidate.
    already_annotated = {
        info["naming_level"]
        for info in contigs.values()
        if info["tier"] in _ANNOTATION_TIERS
    }

    # Group by naming_level
    groups: dict[str, list[str]] = defaultdict(list)
    for seq, info in drafts.items():
        groups[info["naming_level"]].append(seq)

    for naming_level, members in groups.items():
        if naming_level in already_annotated:
            log_info(
                f"  [co-bin] {naming_level}: annotation_ready single(s) exist "
                f"— {len(members)} LQ draft(s) kept as drafts (likely redundant fragments)"
            )
            continue
        combined_len = sum(contigs[m]["contig_length"] for m in members)
        if combined_len >= min_bin_bp:
            for m in members:
                contigs[m]["tier"] = _TIER_BIN
            log_info(
                f"  [co-bin] {naming_level}: {len(members)} contig(s) "
                f"→ bin-rescue ({combined_len:,} bp combined)"
            )

    return contigs

# ── Write outputs ─────────────────────────────────────────────────────────────

def _write_outputs(
    contigs:    dict,
    seqs:       dict[str, str],
    ann_phages: Path,
    ann_prov:   Path,
    drafts_dir: Path,
    rpt_dir:    Path,
    sample_id:  str,
) -> dict:
    """Assign candidate IDs, write FASTAs, write rename_map.tsv."""

    counters: dict[str, int] = defaultdict(int)

    def _next_id(naming_level: str, suffix: str) -> str:
        counters[f"{naming_level}_{suffix}"] += 1
        n = counters[f"{naming_level}_{suffix}"]
        return f"{naming_level}_{suffix}_{n:03d}"

    # ── Assign candidate IDs ──────────────────────────────────────────────────

    # Annotation-ready singles (complete, HQ, MQ, large-ND)
    # Sorted: complete first, then by completeness desc, then by length desc
    def _sort_key(item):
        tier_order = {
            _TIER_COMPLETE: 0, _TIER_HQ: 1, _TIER_MQ: 2, _TIER_LARGE_ND: 3,
        }
        return (tier_order.get(item[1]["tier"], 99),
                -item[1]["completeness"], -item[1]["contig_length"])

    for seq, info in sorted(contigs.items(), key=_sort_key):
        if info["tier"] in _ANNOTATION_TIERS:
            info["candidate_id"] = _next_id(info["naming_level"], "candidate")

    # Bin-rescue (multi-contig): one candidate_id per naming_level group
    bin_groups: dict[str, list[str]] = defaultdict(list)
    for seq, info in contigs.items():
        if info["tier"] == _TIER_BIN:
            bin_groups[info["naming_level"]].append(seq)

    for naming_level, members in bin_groups.items():
        cand_id = _next_id(naming_level, "candidate")
        for m in members:
            contigs[m]["candidate_id"] = cand_id

    # Drafts
    for seq, info in contigs.items():
        if info["tier"] == _TIER_LQ_DRAFT:
            info["candidate_id"] = _next_id(info["naming_level"], "draft")

    # ── Write FASTAs ──────────────────────────────────────────────────────────

    # Single-contig annotation_ready
    for seq, info in contigs.items():
        if info["tier"] not in _ANNOTATION_TIERS:
            continue
        cand_id = info["candidate_id"]
        fasta   = ann_phages / f"{cand_id}.fasta"
        seq_str = seqs.get(seq, "")
        if not seq_str:
            continue
        with open(fasta, "w") as f:
            # Rename header to candidate_id for clean downstream tools
            f.write(
                f">{cand_id}  "
                f"[original={seq}]  "
                f"[length={info['contig_length']}]  "
                f"[completeness={info['completeness']:.1f}%]  "
                f"[topology={info['topology']}]  "
                f"[checkv={info['checkv_quality']}]\n"
            )
            # 60-char line wrap
            for i in range(0, len(seq_str), 60):
                f.write(seq_str[i:i+60] + "\n")

    # Bin-rescue (multi-contig): all members in one FASTA, original headers
    written_bins: set[str] = set()
    for naming_level, members in bin_groups.items():
        cand_id = contigs[members[0]]["candidate_id"]
        if cand_id in written_bins:
            continue
        written_bins.add(cand_id)
        fasta = ann_phages / f"{cand_id}.fasta"
        with open(fasta, "w") as f:
            for m in members:
                seq_str = seqs.get(m, "")
                if seq_str:
                    f.write(f">{m}  [candidate={cand_id}]\n")
                    for i in range(0, len(seq_str), 60):
                        f.write(seq_str[i:i+60] + "\n")

    # Drafts (single-contig per draft)
    for seq, info in contigs.items():
        if info["tier"] != _TIER_LQ_DRAFT:
            continue
        cand_id = info["candidate_id"]
        fasta   = drafts_dir / f"{cand_id}.fasta"
        seq_str = seqs.get(seq, "")
        if not seq_str:
            continue
        with open(fasta, "w") as f:
            f.write(f">{cand_id}  [original={seq}]  [length={info['contig_length']}]\n")
            for i in range(0, len(seq_str), 60):
                f.write(seq_str[i:i+60] + "\n")

    # ── Write rename_map.tsv ──────────────────────────────────────────────────
    _write_rename_map(contigs, rpt_dir / "rename_map.tsv")

    # ── Compute stats ─────────────────────────────────────────────────────────
    tier_counts = Counter(info["tier"] for info in contigs.values())
    n_ann  = tier_counts.get(_TIER_COMPLETE, 0) + tier_counts.get(_TIER_HQ, 0) + \
             tier_counts.get(_TIER_MQ, 0) + tier_counts.get(_TIER_LARGE_ND, 0)
    n_bin  = len(written_bins)
    n_draft = tier_counts.get(_TIER_LQ_DRAFT, 0)
    n_derep = tier_counts.get(_TIER_DEREP, 0)
    n_disc  = tier_counts.get(_TIER_DISCARD, 0)
    n_candidates = len(list(ann_phages.glob("*.fasta")))

    log_ok(
        f"  [quality] {n_candidates} candidate(s) → annotation_ready  "
        f"[complete={tier_counts.get(_TIER_COMPLETE,0)}  "
        f"HQ={tier_counts.get(_TIER_HQ,0)}  "
        f"MQ={tier_counts.get(_TIER_MQ,0)}  "
        f"large-ND={tier_counts.get(_TIER_LARGE_ND,0)}  "
        f"bin={n_bin}]"
    )

    return {
        "n_input":       str(len(contigs)),
        "n_candidates":  str(n_candidates),
        "n_complete":    str(tier_counts.get(_TIER_COMPLETE,  0)),
        "n_hq":          str(tier_counts.get(_TIER_HQ,        0)),
        "n_mq":          str(tier_counts.get(_TIER_MQ,        0)),
        "n_large_nd":    str(tier_counts.get(_TIER_LARGE_ND,  0)),
        "n_bin":         str(n_bin),
        "n_draft":       str(n_draft),
        "n_dereplicated": str(n_derep),
        "n_discarded":   str(n_disc),
    }

# ── Rename map ────────────────────────────────────────────────────────────────

_RENAME_HEADERS = [
    "seq_name", "candidate_id", "tier",
    "topology", "genetic_code", "naming_level",
    "checkv_quality", "completeness", "completeness_method",
    "contig_length", "viral_genes", "host_genes",
    "kmer_freq", "genome_copies", "termini_type", "checkv_warnings",
]

def _write_rename_map(contigs: dict, path: Path) -> None:
    """Write rename_map.tsv — the critical metadata link to annotate.py.

    Contains all per-contig metadata needed downstream:
      topology     → annotate.py: circular/linear plot, dnaapler decision
      genetic_code → annotate.py: Pharokka --genetic_code
      naming_level → quality naming, co-binning context
      candidate_id → the output FASTA filename (without .fasta)
      tier         → which output directory the candidate went to
    """
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=_RENAME_HEADERS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for seq, info in sorted(contigs.items()):
            writer.writerow({
                "seq_name":          seq,
                "candidate_id":      info.get("candidate_id", ""),
                "tier":              info["tier"],
                "topology":          info["topology"],
                "genetic_code":      info["genetic_code"],
                "naming_level":      info["naming_level"],
                "checkv_quality":    info["checkv_quality"],
                "completeness":      f"{info['completeness']:.2f}",
                "completeness_method": info.get("completeness_method", ""),
                "contig_length":     info["contig_length"],
                "viral_genes":       info["viral_genes"],
                "host_genes":        info["host_genes"],
                "kmer_freq":         f"{info['kmer_freq']:.3f}",
                "genome_copies":     f"{info['genome_copies']:.3f}",
                "termini_type":      info["termini_type"],
                "checkv_warnings":   info.get("checkv_warnings", ""),
            })

# ── Redundancy warning ───────────────────────────────────────────────────────

def _warn_naming_level_redundancy(contigs: dict, active_warnings: list[str]) -> None:
    """Warn when multiple candidates of the same naming_level look redundant.

    If the combined length of all annotation_ready candidates in a naming_level
    exceeds 1.5× the length of the longest single candidate, the candidates
    are likely overlapping fragments of the same genome rather than distinct
    phages. This is a diagnostic warning — no candidates are discarded.

    Root cause: mash distance between non-overlapping fragments of the same
    large genome appears high (k-mer content is disjoint), so pairwise ANI
    dereplication does not cluster them. The redundancy manifests as multiple
    candidates from the same family representing the same phage genome.
    """
    groups: dict[str, list] = defaultdict(list)
    for seq, info in contigs.items():
        if info["tier"] in _ANNOTATION_TIERS or info["tier"] == _TIER_BIN:
            groups[info["naming_level"]].append(info["contig_length"])

    for naming_level, lengths in groups.items():
        if len(lengths) < 2:
            continue
        total   = sum(lengths)
        largest = max(lengths)
        if total > 1.5 * largest:
            msg = (
                f"{naming_level}: {len(lengths)} candidate(s), combined "
                f"{total:,} bp ({total / largest:.1f}× largest). "
                "Possible redundant fragments of the same genome — "
                "dereplication could not cluster non-overlapping contig pairs. "
                "Review annotation_ready candidates and select the most complete."
            )
            log_warn(f"  [redundancy] {msg}")
            active_warnings.append(msg)

# ── Completion panel ──────────────────────────────────────────────────────────

def _print_completion_panel(
    sample_id:  str,
    ann_phages: Path,
    drafts_dir: Path,
    rpt_dir:    Path,
    stats:      dict,
    active_warnings: list[str],
) -> None:
    def _f(key) -> str:
        v = str(stats.get(key, "0"))
        return v if v else "0"

    candidates = list(ann_phages.glob("*.fasta"))

    lines = [
        "[bold]Quality summary[/bold]",
        f"  [cyan]Input contigs     :[/cyan] {_f('n_input')}",
        f"  [cyan]Candidates        :[/cyan] [bold green]{_f('n_candidates')}[/bold green]"
        f"  [dim](complete={_f('n_complete')}  HQ={_f('n_hq')}  "
        f"MQ={_f('n_mq')}  large-ND={_f('n_large_nd')}  bin={_f('n_bin')})[/dim]",
        f"  [cyan]Drafts            :[/cyan] {_f('n_draft')}",
        f"  [cyan]Dereplicated      :[/cyan] {_f('n_dereplicated')}",
        f"  [cyan]Discarded         :[/cyan] {_f('n_discarded')}",
    ]

    if candidates:
        lines += ["", "[bold]Annotation-ready candidates[/bold]"]
        for fa in sorted(candidates):
            lines.append(f"  [green]✓[/green] {fa.stem}")

    lines += [
        "",
        "[bold]Output directories[/bold]",
        f"  annotation_ready/ : {ann_phages}",
        f"  drafts/           : {drafts_dir}",
        f"  rename_map.tsv    : {rpt_dir / 'rename_map.tsv'}",
    ]

    if active_warnings:
        lines += ["", "[bold yellow]Warnings[/bold yellow]"]
        lines += [f"  [yellow]⚠ {w}[/yellow]" for w in active_warnings]

    console.print(Panel(
        "\n".join(lines),
        title=f"[bold cyan]Quality complete — {sample_id}[/bold cyan]",
        border_style="cyan", padding=(0, 2), width=120,
    ))

# ── Read recruitment + coverage ──────────────────────────────────────────────

def _read_recruitment_coverage(
    candidates:      list,
    sample_id:       str,
    cfg,
    rpt_dir:         Path,
    log_f:           Path,
    active_warnings: list[str],
) -> dict:
    """Map post-host-removal reads to annotation_ready candidates.

    Computes per-candidate:
      mean_depth : mean read depth across all positions
      breadth    : fraction of positions covered at ≥1x
      cv         : coefficient of variation of depth (std/mean)

    Roux et al. 2019 (eLife 8:e42923) — read recruitment is the primary
    validation of phage genome integrity.
    Nayfach et al. 2021 (Nat Biotechnol 39:578) — read mapping for circular
    contig confirmation.

    Reads auto-discovered from results/02_host_removal/{sample_id}_R*.
    Silently skipped if reads not found.
    """
    if not candidates:
        return {}

    r1 = cfg.results("02_host_removal") / f"{sample_id}_R1.fastq.gz"
    r2 = cfg.results("02_host_removal") / f"{sample_id}_R2.fastq.gz"
    if not r1.exists() or not r2.exists():
        log_info(
            f"  [coverage] reads not found at {r1.parent} — "
            "coverage validation skipped"
        )
        return {}

    cov_dir = rpt_dir / "coverage"
    cov_dir.mkdir(parents=True, exist_ok=True)
    combined  = cov_dir / "combined_candidates.fasta"
    bam       = cov_dir / "reads_vs_candidates.bam"
    depth_tsv = cov_dir / "depth.tsv"
    threads   = cfg.threads

    # Combined reference: use candidate_id as header for easy depth parsing
    with open(combined, "w") as out:
        for fa in candidates:
            cand_id = fa.stem
            seq_parts: list[str] = []
            with open(fa) as f:
                for line in f:
                    if line.startswith(">"):
                        if seq_parts:
                            out.write("".join(seq_parts))
                            out.write("\n")
                        out.write(f">{cand_id}\n")
                        seq_parts = []
                    else:
                        seq_parts.append(line.rstrip())
            if seq_parts:
                out.write("".join(seq_parts))
                out.write("\n")

    try:
        run_silent(["bwa-mem2", "index", str(combined)], log_file=log_f)
        run_silent(
            # -M: mark split reads as secondary (consistent with host_removal).
            # Default parameters used here (not permissive): coverage validation
            # measures genuine mapping of reads to their source genome.
            f"bwa-mem2 mem -t {threads} -M {combined} {r1} {r2} "
            f"| samtools view -@ {min(threads,8)} -bF 2304 -u "
            f"| samtools sort -@ {min(threads, 8)} -o {bam}",
            log_file=log_f, shell=True,
        )
        run_silent(["samtools", "index", str(bam)], log_file=log_f)
        subprocess.run(
            # -Q 20: skip alignments with mapQ < 20 to avoid
            # spurious coverage from low-confidence multimappers.
            ["samtools", "depth", "-a", "-Q", "20", str(bam)],
            stdout=open(depth_tsv, "w"),
            stderr=subprocess.DEVNULL,
            check=True, timeout=600,
        )
    except Exception as e:
        log_warn(f"  [coverage] alignment failed: {e} — coverage validation skipped")
        return {}

    depths_by_cand: dict[str, list[int]] = defaultdict(list)
    with open(depth_tsv) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            try:
                depths_by_cand[parts[0]].append(int(parts[2]))
            except ValueError:
                continue

    results: dict[str, dict] = {}
    for cand_id, depths in depths_by_cand.items():
        if not depths:
            continue
        n        = len(depths)
        mean_d   = sum(depths) / n
        breadth  = sum(1 for d in depths if d >= 1) / n
        variance = sum((d - mean_d) ** 2 for d in depths) / n
        cv       = math.sqrt(variance) / mean_d if mean_d > 0 else 0.0

        results[cand_id] = {
            "mean_depth": f"{mean_d:.1f}",
            "breadth":    f"{breadth:.4f}",
            "cv":         f"{cv:.3f}",
        }
        log_info(
            f"  [coverage] {cand_id}: "
            f"depth={mean_d:.0f}x  breadth={breadth * 100:.1f}%  CV={cv:.2f}"
        )
        if breadth < 0.95:
            msg = (
                f"{cand_id}: breadth={breadth * 100:.1f}% (<95%) — "
                "uncovered region(s). Possible assembly gap or chimera."
            )
            log_warn(f"  [coverage] {msg}")
            active_warnings.append(msg)
        elif cv > 1.5:
            msg = (
                f"{cand_id}: coverage CV={cv:.2f} (>1.5) — "
                "uneven depth. Possible concatemer or chimeric assembly."
            )
            log_warn(f"  [coverage] {msg}")
            active_warnings.append(msg)
        else:
            log_ok(f"  [coverage] {cand_id}: uniform coverage — genome integrity OK")

    return results

def _annotate_rename_map_coverage(rename_map: Path, coverage: dict) -> None:
    """Add mean_depth, breadth, cv columns to existing rename_map.tsv."""
    if not rename_map.exists() or not coverage:
        return
    with open(rename_map) as f:
        lines = f.readlines()
    if not lines:
        return

    headers = lines[0].rstrip("\n").split("\t")
    if "mean_depth" not in headers:
        headers += ["mean_depth", "breadth", "cv"]
        lines[0] = "\t".join(headers) + "\n"

    out_lines = [lines[0]]
    for line in lines[1:]:
        parts = line.rstrip("\n").split("\t")
        row   = dict(zip(headers[:len(parts)], parts))
        cand  = row.get("candidate_id", "")
        cov   = coverage.get(cand, {})
        row["mean_depth"] = cov.get("mean_depth", "")
        row["breadth"]    = cov.get("breadth", "")
        row["cv"]         = cov.get("cv", "")
        out_lines.append("\t".join(row.get(h, "") for h in headers) + "\n")

    with open(rename_map, "w") as f:
        f.writelines(out_lines)

# ── TSV summary ───────────────────────────────────────────────────────────────

_TSV_HEADERS = [
    "sample", "n_input", "n_candidates",
    "n_complete", "n_hq", "n_mq", "n_large_nd", "n_bin",
    "n_draft", "n_dereplicated", "n_discarded",
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

"""PhageFlow Module 05 — Genome quality assessment and selection.

Goal: maximize recovery of Complete and High-quality phage genomes.

Tools:
    CheckV v1.0+ (Nayfach et al. 2021, Nature Biotechnology 39:578)
    mash v2.3   (Ondov et al. 2016, Genome Biology 17:132)   -- all contigs when available
    cd-hit-est  (Fu et al. 2012, Bioinformatics 28:3150)     -- fallback only
    seqkit      (Shen et al. 2016, PLoS ONE 11:e0163962)

Output structure (per-sample folders):
    results/05_quality/annotation_ready/phages/{sample_id}/
    results/05_quality/annotation_ready/proviruses/{sample_id}/
    results/05_quality/drafts/{sample_id}/
    All output paths are namespaced by sample_id to prevent cross-sample
    filename collisions (e.g. two samples both producing
    Ackermannviridae_candidate_001).

Dereplication strategy:
    mash is used for ALL contig sizes when available.
    Rationale: bacteriophage genomes are frequently circular. Different
    assemblers (SPAdes, MEGAHIT) linearize the same circular genome at
    different positions, producing circular permutations. cd-hit-est uses
    linear alignment and cannot detect permutations: individual HSPs each
    cover only a fraction of the sequence, failing the -aS threshold even
    when combined coverage is 100%. Inoviridae (ssDNA, ~6 kb) is the
    canonical example, but the same applies to any circular phage genome.
    mash computes MinHash distances from k-mer composition, which is
    invariant to rotation — two permutations of the same circular sequence
    have distance ≈ 0 and are correctly clustered (Ondov et al. 2016).

    Threshold: _MASH_ANI_THRESH = 0.02 (= ANI > 98%; Turner et al. 2021,
    Arch Virol 166:2633)

    Fallback (mash not in PATH): cd-hit-est at 98% ANI for all sizes.
    A _LARGE_CONTIG_BP split (50 kb) is applied only in the cd-hit-est
    fallback because cd-hit-est's alignment band cannot reliably span
    large-genome indels between assemblers (Fu et al. 2012).

Taxonomy path:
    Reads genomad_tax_tsv from viral-id's genomad_summary.tsv.
    Falls back to convention-based path for backward compatibility.

Completeness in rename_map:
    'ND' for Not-determined genomes (not '0.0', which is misleading).
"""

from __future__ import annotations
import shutil
import subprocess
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from rich.panel import Panel
from rich.text import Text
from rich.progress import (
    Progress, SpinnerColumn, TextColumn,
    BarColumn, MofNCompleteColumn, TimeElapsedColumn,
)

from phageflow.utils.config import Config
from phageflow.utils.logger import (
    log_step, log_info, log_ok, log_warn, log_error, console,
)
from phageflow.utils.tools import (
    require_tools, run_silent, human_size, mkdirs, fasta_stats,
)

STEP  = "05_quality"
TOOLS = ["checkv", "cd-hit-est", "seqkit"]

_COMPLETENESS_GOOD = 90.0
_COMPLETENESS_WARN = 50.0
_CONTAM_WARN       = 5.0
_HOST_GENE_WARN    = 3

_CDHIT_ID = "0.98"
_CDHIT_AS = "0.85"
_CDHIT_N  = "8"

# Dereplication thresholds
# _LARGE_CONTIG_BP is only used in the cd-hit-est FALLBACK path.
# When mash is available it is NOT used — mash handles all sizes uniformly.
_LARGE_CONTIG_BP  = 50_000   # cd-hit-est fallback only: >= this -> wider band
_MASH_ANI_THRESH  = 0.02     # distance < 0.02 = ANI > 98% (Turner et al. 2021)
_MASH_SKETCH_SIZE = "10000"  # MinHash sketch size

_TIER_ORDER = [
    "Complete", "High-quality", "Medium-quality", "Low-quality", "Not-determined"
]
_HIGH_CONFIDENCE_TIERS = frozenset({"Complete", "High-quality"})


# -- Data structure ------------------------------------------------------------

@dataclass
class ContigRecord:
    contig_id:         str
    sequence:          str
    length_bp:         int
    completeness:      float
    completeness_str:  str    # 'ND' for Not-determined, '85.3' otherwise
    quality_tier:      str
    is_provirus:       bool
    viral_genes:       int
    is_multicontig:    bool = False
    sub_records:       list = field(default_factory=list)


# -- Public entry point --------------------------------------------------------

def run(
    cfg:       Config,
    sample_id: str,
    virus_fna: Path,
    force:     bool = False,
) -> list[Path]:
    """
    Assess genome quality, rescue, dereplicate, co-bin, and format candidates.

    Output paths are namespaced under {sample_id}/ to prevent cross-sample
    filename collisions.

    Returns
    -------
    hq_fastas : list of annotation-ready FASTA paths
    """
    require_tools(*TOOLS)

    virus_fna = Path(virus_fna)
    out_dir   = cfg.results(STEP)
    rpt_dir   = cfg.reports(STEP)

    # Per-sample output folders -- prevent cross-sample filename collisions
    phage_dir = out_dir / "annotation_ready" / "phages"     / sample_id
    prov_dir  = out_dir / "annotation_ready" / "proviruses" / sample_id
    draft_dir = out_dir / "drafts"                          / sample_id
    tmp_dir   = out_dir / "tmp"                             / sample_id
    mkdirs(out_dir, rpt_dir, phage_dir, prov_dir, draft_dir, tmp_dir)

    log_step(f"Module 05 -- quality [{sample_id}]")

    if not virus_fna.exists():
        log_error(f"  Virus FASTA not found: {virus_fna}")
        raise FileNotFoundError(virus_fna)
    if virus_fna.stat().st_size == 0:
        log_warn(f"  Virus FASTA is empty -- no genomes to assess: {virus_fna}")
        return []

    s = fasta_stats(virus_fna)
    _print_input_table(sample_id, virus_fna, s)
    log_info(
        f"  CheckV  (Nayfach et al. 2021)  |  "
        f"HQ threshold >={cfg.checkv.min_completeness}%  |  "
        f"Complete/HQ always promoted"
    )
    log_info(
        f"  Rescue  : LQ >=1 gene | ND >={cfg.checkv.length_rescue:,}bp -> drafts/ | "
        f"ND >={cfg.checkv.large_nd_rescue_bp:,}bp -> annotation_ready/"
    )

    import shutil as _sh
    _dedup_tool = "mash" if _sh.which("mash") else "cd-hit-est (fallback)"
    log_info(
        f"  Dedup   : {_dedup_tool} {_MASH_ANI_THRESH*100:.0f}% dist threshold "
        f"(all sizes, circular-safe) | Ondov et al. 2016"
    )
    log_info(
        f"  Co-bin  : >=30kb shared-taxonomy drafts -> annotation_ready/ | "
        f"output: annotation_ready/{{phages,proviruses}}/{sample_id}/"
    )

    if not cfg.databases.checkv.exists():
        log_warn(f"  CheckV database not found: {cfg.databases.checkv}")

    genomad_tsv = _find_genomad_tsv(cfg, sample_id)
    tax_map     = _load_genomad_taxonomy(genomad_tsv)
    if not tax_map:
        log_warn(
            "  geNomad taxonomy not loaded -- co-binning disabled. "
            "Run viral-id before quality."
        )

    sdir  = out_dir / sample_id
    qfile = sdir / "quality_summary.tsv"

    summary:       dict               = {}
    hq_records:    list[ContigRecord] = []
    draft_records: list[ContigRecord] = []
    promoted_bins: list[ContigRecord] = []
    hq_fastas:     list[Path]         = []

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

        # 1/4 -- CheckV
        progress.update(task, description="[1/4] CheckV   -- end-to-end assessment")
        _run_checkv(cfg, sample_id, virus_fna, sdir, rpt_dir, force)
        progress.advance(task)

        # 2/4 -- Tier selection + rescue
        progress.update(task, description="[2/4] CheckV   -- tier selection + rescue")
        if qfile.exists():
            seqs      = _load_fasta(virus_fna)
            prov_seqs = _load_provirus_seqs(sdir / "proviruses.fna")
            summary, hq_records, draft_records = _parse_checkv(
                sample_id, qfile, seqs, prov_seqs, cfg.checkv,
            )
            promoted_bins, draft_records = _cobin_draft_rescue(
                draft_records, tax_map, sample_id,
                cfg.checkv.min_bin_rescue_bp,
            )
            summary["rescued_draft_bins"] = len(promoted_bins)
            _write_fasta_records(
                draft_records, draft_dir / f"{sample_id}_draft.fasta",
            )
            log_ok("  [CheckV] genome selection complete")
        else:
            log_warn(
                "  [CheckV] quality_summary.tsv not found -- "
                f"check reports/{STEP}/{sample_id}_checkv.log"
            )
        progress.advance(task)

        # 3/4 -- Dereplication (mash for all sizes; cd-hit-est fallback)
        progress.update(
            task,
            description="[3/4] dedup    -- mash all-size (circular-safe)",
        )
        if hq_records:
            hq_records, n_before, n_after = _dereplicate(
                hq_records, tmp_dir, sample_id, rpt_dir, cfg.threads, force,
            )
            log_ok(
                f"  [dedup] {n_before} -> {n_after} single-contig genomes  "
                f"(-{n_before - n_after} redundant)"
            )
        hq_records.extend(promoted_bins)
        if promoted_bins:
            log_ok(f"  [draft rescue] {len(promoted_bins)} bin(s) promoted")
        progress.advance(task)

        # 4/4 -- MQ co-bin + rename + seqkit
        progress.update(task, description="[4/4] seqkit   -- co-binning, renaming & formatting")
        if hq_records:
            hq_records = _cobin_by_taxon(hq_records, tax_map, sample_id)
            rename_map = _build_rename_map(hq_records, tax_map)
            _save_rename_map(rename_map, rpt_dir / f"{sample_id}_rename_map.tsv")
            hq_fastas = _write_candidates(
                hq_records, rename_map, phage_dir, prov_dir, tmp_dir,
                rpt_dir / f"{sample_id}_seqkit.log",
            )
            n_multi = sum(1 for r in hq_records if r.is_multicontig)
            n_phage = sum(1 for r in hq_records if not r.is_provirus)
            n_prov  = sum(1 for r in hq_records if r.is_provirus)
            log_ok(
                f"  [seqkit] {len(hq_fastas)} file(s) -> "
                f"{phage_dir.parent.name}/{sample_id}/  "
                f"(phage={n_phage}  provirus={n_prov}  multi-contig={n_multi})"
            )
        progress.advance(task)

    shutil.rmtree(tmp_dir, ignore_errors=True)

    _validate_output(sample_id, hq_fastas, draft_records)
    _print_summary_table(sample_id, summary, hq_fastas, draft_records)
    _check_warnings(summary, cfg.checkv.min_completeness)
    if summary:
        _save_tsv(summary, rpt_dir / "checkv_summary.tsv")
    _print_completion_panel(
        sample_id, phage_dir, prov_dir, draft_dir, rpt_dir,
        summary, hq_fastas, draft_records,
    )

    log_step(f"Module 05 completed v  [{sample_id}]")
    log_info("Next: phageflow annotate --sample-id <id> --genome <candidate.fasta>")
    return hq_fastas


# -- Step 1: CheckV ------------------------------------------------------------

def _run_checkv(cfg, sample_id, virus_fna, sdir, rpt_dir, force):
    qfile = sdir / "quality_summary.tsv"
    if qfile.exists() and qfile.stat().st_size > 0 and not force:
        log_info("  [CheckV] already processed -- skipping  (--force to re-run)")
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
        log_ok("  [CheckV] assessment complete")
    except Exception as e:
        log_warn(f"  [CheckV] warning: {e}")


# -- Step 2: Tier selection + rescue -------------------------------------------

def _parse_checkv(sample_id, qfile, seqs, prov_seqs, chkv):
    tiers: dict[str, int] = {t: 0 for t in _TIER_ORDER}
    rescued_lq = rescued_nd = rescued_nd_density = rescued_large_nd = 0
    total_viral_genes = total_host_genes = 0
    best_len = 0
    best_completeness = 0.0
    best_quality_tier = ""
    max_contamination = 0.0
    hq_records:    list[ContigRecord] = []
    draft_records: list[ContigRecord] = []

    try:
        with open(qfile) as f:
            header = f.readline().strip().split("\t")
            col    = {h: i for i, h in enumerate(header)}

            for line in f:
                row = line.strip().split("\t")
                if not row or len(row) < 3:
                    continue

                cid      = _col(row, col, "contig_id",      0)
                clen_s   = _col(row, col, "contig_length",  1)
                qual     = _col(row, col, "checkv_quality", 7)
                comp_s   = _col(row, col, "completeness",   9)
                contm_s  = _col(row, col, "contamination", 10)
                vgn_s    = _col(row, col, "viral_genes",    4)
                hgn_s    = _col(row, col, "host_genes",     5)
                prov_col = _col(row, col, "provirus",        2)

                is_provirus = (
                    prov_col.strip().lower() == "yes"
                    or "|provirus" in cid.lower()
                )

                clen_int = int(clen_s or 0)
                vgn_int  = int(vgn_s  or 0)

                if clen_int < chkv.min_contig_bp:
                    log_info(
                        f"  [filter] {cid} discarded -- "
                        f"{clen_int}bp < min_contig_bp {chkv.min_contig_bp}bp"
                    )
                    continue

                gene_density = vgn_int / (clen_int / 1000) if clen_int > 0 else 0.0
                tiers[qual if qual in tiers else "Not-determined"] += 1

                try:   total_viral_genes += vgn_int
                except (ValueError, TypeError): pass
                try:   total_host_genes  += int(hgn_s)
                except (ValueError, TypeError): pass
                try:   max_contamination = max(max_contamination, float(contm_s))
                except (ValueError, TypeError): pass

                is_nd = qual == "Not-determined" or comp_s in ("", "NA", "N/A")
                try:
                    comp_val = 0.0 if is_nd else float(comp_s)
                except (ValueError, TypeError):
                    comp_val = 0.0
                    is_nd = True
                comp_str = "ND" if is_nd else f"{comp_val:.1f}"

                if clen_int > best_len:
                    best_len          = clen_int
                    best_completeness = comp_val
                    best_quality_tier = qual

                seq = seqs.get(cid, "")
                if is_provirus and cid in prov_seqs:
                    seq = prov_seqs[cid]
                if not seq:
                    continue

                record = ContigRecord(
                    contig_id        = cid,
                    sequence         = seq,
                    length_bp        = clen_int,
                    completeness     = comp_val,
                    completeness_str = comp_str,
                    quality_tier     = qual,
                    is_provirus      = is_provirus,
                    viral_genes      = vgn_int,
                )

                if qual in ("Complete", "High-quality"):
                    hq_records.append(record)
                elif qual == "Medium-quality":
                    if comp_val >= chkv.min_completeness:
                        hq_records.append(record)
                    else:
                        draft_records.append(record)
                elif qual == "Low-quality":
                    if vgn_int >= chkv.min_viral_genes or gene_density >= chkv.min_gene_density:
                        log_warn(
                            f"  [LQ rescue] {cid} -- "
                            f"{vgn_int} viral gene(s), "
                            f"density={gene_density:.2f}/kb -> drafts/"
                        )
                        draft_records.append(record)
                        rescued_lq += 1
                elif qual == "Not-determined":
                    if clen_int >= chkv.large_nd_rescue_bp and vgn_int >= 1:
                        log_warn(
                            f"  [large-ND rescue] {cid} -- "
                            f"{clen_int:,}bp, {vgn_int} gene(s) -> annotation_ready/"
                        )
                        hq_records.append(record)
                        rescued_large_nd += 1
                    elif clen_int >= chkv.length_rescue and vgn_int >= 1:
                        log_warn(
                            f"  [ND rescue] {cid} -- "
                            f"{clen_int:,}bp, {vgn_int} gene(s) -> drafts/"
                        )
                        draft_records.append(record)
                        rescued_nd += 1
                    elif gene_density >= chkv.min_gene_density:
                        log_warn(
                            f"  [density-ND rescue] {cid} -- "
                            f"density={gene_density:.2f}/kb -> drafts/"
                        )
                        draft_records.append(record)
                        rescued_nd_density += 1

    except Exception as e:
        log_warn(f"  [CheckV] could not parse {qfile}: {e}")

    rescued_total = rescued_lq + rescued_nd + rescued_nd_density + rescued_large_nd
    if rescued_total:
        log_ok(
            f"  [rescue] LQ={rescued_lq}  ND-moderate={rescued_nd}  "
            f"ND-density={rescued_nd_density}  ND-large->HQ={rescued_large_nd}"
        )

    return {
        "sample":             sample_id,
        "complete":           tiers.get("Complete",       0),
        "hq":                 tiers.get("High-quality",   0),
        "mq":                 tiers.get("Medium-quality", 0),
        "lq":                 tiers.get("Low-quality",    0),
        "nd":                 tiers.get("Not-determined", 0),
        "rescued_lq":         rescued_lq,
        "rescued_nd":         rescued_nd,
        "rescued_nd_density": rescued_nd_density,
        "rescued_large_nd":   rescued_large_nd,
        "rescued_draft_bins": 0,
        "hq_promoted":        len(hq_records),
        "draft_saved":        len(draft_records),
        "best_length_bp":     best_len,
        "best_completeness":  best_completeness,
        "best_quality_tier":  best_quality_tier,
        "total_viral_genes":  total_viral_genes,
        "total_host_genes":   total_host_genes,
        "max_contamination":  round(max_contamination, 2),
    }, hq_records, draft_records


# -- Draft co-bin rescue -------------------------------------------------------

def _cobin_draft_rescue(draft_records, tax_map, sample_id, min_bin_rescue_bp):
    taxon_buckets: dict[str, list[ContigRecord]] = defaultdict(list)
    no_taxon:      list[ContigRecord]             = []

    for rec in draft_records:
        taxon = tax_map.get(rec.contig_id)
        if not taxon:
            log_info(
                f"  [draft-bin] {rec.contig_id}: no taxonomy -> "
                "individual draft (Camargo et al. 2023)"
            )
            no_taxon.append(rec)
            continue
        taxon_buckets[taxon].append(rec)

    promoted:  list[ContigRecord] = []
    remaining: list[ContigRecord] = list(no_taxon)

    for taxon, records in taxon_buckets.items():
        total_len = sum(r.length_bp for r in records)
        total_vgn = sum(r.viral_genes for r in records)
        if len(records) >= 2 and total_len >= min_bin_rescue_bp:
            log_warn(
                f"  [draft-bin rescue] {taxon}: {len(records)} contigs, "
                f"{total_len:,}bp, {total_vgn} viral genes "
                f"-> annotation_ready/"
            )
            bin_id = "__bin__" + "||".join(r.contig_id for r in records)
            promoted.append(ContigRecord(
                contig_id        = bin_id,
                sequence         = "",
                length_bp        = total_len,
                completeness     = max(r.completeness for r in records),
                completeness_str = "ND-bin",
                quality_tier     = _best_quality_tier(r.quality_tier for r in records),
                is_provirus      = any(r.is_provirus for r in records),
                viral_genes      = total_vgn,
                is_multicontig   = True,
                sub_records      = records,
            ))
        else:
            remaining.extend(records)

    return promoted, remaining


# -- Step 3: Dereplication (mash for all sizes; cd-hit-est fallback) -----------

def _dereplicate(
    records: list[ContigRecord],
    tmp_dir: Path,
    sample_id: str,
    rpt_dir: Path,
    threads: int,
    force: bool,
) -> tuple[list[ContigRecord], int, int]:
    """
    Dereplicate HQ single-contig genomes.

    Strategy (in priority order):
    1. mash (all sizes) — preferred when mash is in PATH.
       mash computes distances from k-mer composition (MinHash sketches),
       which is invariant to sequence orientation and starting position.
       This correctly handles circular phage genomes whose linearizations
       by SPAdes and MEGAHIT begin at different positions (circular
       permutations). cd-hit-est cannot detect these because individual
       HSPs each cover only a fraction of the sequence.
       Reference: Ondov et al. (2016) Genome Biology 17:132.

    2. cd-hit-est (fallback when mash not in PATH) — split by size:
       < _LARGE_CONTIG_BP (50 kb): standard parameters.
       >= _LARGE_CONTIG_BP: wider alignment band (-band 500).
       Both at 98% ANI (Turner et al. 2021, Arch Virol 166:2633).
       WARNING: cd-hit-est will miss circular permutations.

    Both target distance < 0.02 = ANI > 98%, above the ICTV species
    boundary of 95%, to remove assembly duplicates while preserving
    biological strain diversity.
    """
    n_before = len(records)
    if n_before <= 1:
        return records, n_before, n_before

    import shutil as _sh
    if _sh.which("mash"):
        log_info(
            "  [dedup] using mash for all sizes (circular permutation-safe)"
        )
        result = _mash_dereplicate(records, tmp_dir, sample_id, rpt_dir)
    else:
        log_warn(
            "  [dedup] mash not found -- falling back to cd-hit-est. "
            "Circular permutations (e.g. Inoviridae) may NOT be deduplicated. "
            "Install mash: mamba install -n phageflow bioconda::mash"
        )
        small = [r for r in records if r.length_bp <  _LARGE_CONTIG_BP]
        large = [r for r in records if r.length_bp >= _LARGE_CONTIG_BP]
        result = (
            _cdhit_dereplicate(small, tmp_dir, sample_id, rpt_dir, threads)
            if len(small) > 1 else small
        ) + (
            _cdhit_dereplicate_large(large, tmp_dir, sample_id, rpt_dir, threads)
            if len(large) > 1 else large
        )

    return result, n_before, len(result)


def _cdhit_dereplicate(
    records: list[ContigRecord],
    tmp_dir: Path,
    sample_id: str,
    rpt_dir: Path,
    threads: int,
) -> list[ContigRecord]:
    """cd-hit-est dereplication for contigs < 50kb (fallback only)."""
    if len(records) <= 1:
        return records

    tmp_in  = tmp_dir / f"{sample_id}_small_raw.fasta"
    tmp_out = tmp_dir / f"{sample_id}_small_derep.fasta"
    _write_fasta_records(records, tmp_in)

    try:
        run_silent([
            "cd-hit-est",
            "-i",  str(tmp_in), "-o", str(tmp_out),
            "-c",  _CDHIT_ID, "-aS", _CDHIT_AS,
            "-G",  "0", "-n", _CDHIT_N, "-d", "0",
            "-sc", "1", "-T", str(threads), "-M", "4000",
        ], log_file=rpt_dir / f"{sample_id}_cdhit_small.log")
    except Exception as e:
        log_warn(f"  [cd-hit-est] small contig dedup warning: {e}")
        return records
    finally:
        Path(str(tmp_out) + ".clstr").unlink(missing_ok=True)

    rep_ids = set(_load_fasta(tmp_out).keys())
    return [r for r in records if r.contig_id in rep_ids]


def _cdhit_dereplicate_large(
    records: list[ContigRecord],
    tmp_dir: Path,
    sample_id: str,
    rpt_dir: Path,
    threads: int,
) -> list[ContigRecord]:
    """cd-hit-est dereplication for contigs >= 50kb (fallback only, wider band)."""
    if len(records) <= 1:
        return records

    tmp_in  = tmp_dir / f"{sample_id}_large_raw.fasta"
    tmp_out = tmp_dir / f"{sample_id}_large_derep.fasta"
    _write_fasta_records(records, tmp_in)

    try:
        run_silent([
            "cd-hit-est",
            "-i", str(tmp_in), "-o", str(tmp_out),
            "-c", _CDHIT_ID, "-aS", _CDHIT_AS,
            "-G", "0", "-n", _CDHIT_N, "-d", "0",
            "-sc", "1", "-T", str(threads), "-M", "8000",
            "-band", "500",   # wider band for large sequences
        ], log_file=rpt_dir / f"{sample_id}_cdhit_large.log")
    except Exception as e:
        log_warn(f"  [cd-hit-est] large contig dedup warning: {e}")
        return records
    finally:
        Path(str(tmp_out) + ".clstr").unlink(missing_ok=True)

    rep_ids = set(_load_fasta(tmp_out).keys())
    return [r for r in records if r.contig_id in rep_ids]


def _mash_dereplicate(
    records: list[ContigRecord],
    tmp_dir: Path,
    sample_id: str,
    rpt_dir: Path,
) -> list[ContigRecord]:
    """
    mash-based dereplication for ALL contig sizes.

    mash sketch -i sketches each sequence from the multi-FASTA independently.
    mash dist sketch.msh sketch.msh gives all pairwise distances.
    Union-Find clusters sequences with distance < _MASH_ANI_THRESH (0.02 = ANI > 98%);
    the longest representative is kept from each cluster.

    Why mash for all sizes (including small phages like Inoviridae ~6 kb):
        mash distance = 1 - Jaccard(k-mer sets), computed from MinHash sketches.
        k-mer composition is invariant to the starting position of a circular
        sequence linearization. Two SPAdes/MEGAHIT assemblies of the same
        circular phage genome linearized at different positions share all k-mers
        (ignoring assembly artifacts), yielding distance ≈ 0.
        cd-hit-est uses linear alignment and fails for circular permutations:
        each HSP covers only a fraction of the sequence, failing -aS thresholds.

    Reference: Ondov et al. (2016) Genome Biology 17:132.
    """
    if len(records) <= 1:
        return records

    combined  = tmp_dir / f"{sample_id}_all.fasta"
    sketch_pf = tmp_dir / f"{sample_id}_all"     # mash appends .msh
    sketch    = Path(str(sketch_pf) + ".msh")

    _write_fasta_records(records, combined)

    try:
        run_silent(
            ["mash", "sketch", "-i",
             "-s", _MASH_SKETCH_SIZE,
             "-o", str(sketch_pf),
             str(combined)],
            log_file=rpt_dir / f"{sample_id}_mash.log",
        )
    except Exception as e:
        log_warn(
            f"  [mash sketch] warning: {e} -- "
            "falling back to cd-hit-est (circular permutations may not be caught)"
        )
        return records

    try:
        result = subprocess.run(
            ["mash", "dist", str(sketch), str(sketch)],
            capture_output=True, text=True, check=False,
        )
        dist_output = result.stdout
    except Exception as e:
        log_warn(f"  [mash dist] warning: {e} -- skipping dereplication")
        return records

    # Union-Find clustering
    id_to_record = {r.contig_id: r for r in records}
    parent: dict[str, str] = {r.contig_id: r.contig_id for r in records}

    def _find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def _union(x: str, y: str) -> None:
        px, py = _find(x), _find(y)
        if px == py:
            return
        # Keep the longer sequence as cluster root
        if id_to_record[px].length_bp >= id_to_record[py].length_bp:
            parent[py] = px
        else:
            parent[px] = py

    n_pairs = 0
    for line in dist_output.splitlines():
        parts = line.strip().split("\t")
        if len(parts) < 3:
            continue
        # mash output: ref_path query_path distance p-value shared/total
        # seq name is the last whitespace-separated token in the path field
        ref = parts[0].split()[-1]
        qry = parts[1].split()[-1]
        try:
            dist = float(parts[2])
        except (ValueError, IndexError):
            continue
        if ref == qry:
            continue
        if ref in parent and qry in parent and dist < _MASH_ANI_THRESH:
            _union(ref, qry)
            n_pairs += 1

    if n_pairs:
        log_info(
            f"  [mash] {n_pairs} near-identical pair(s) found "
            f"(dist < {_MASH_ANI_THRESH}, ANI > {(1 - _MASH_ANI_THRESH)*100:.0f}%)"
        )

    representatives = {_find(r.contig_id) for r in records}
    return [r for r in records if r.contig_id in representatives]


# -- Step 4a: MQ co-binning ----------------------------------------------------

def _cobin_by_taxon(records, tax_map, sample_id):
    already_multi     = [r for r in records if r.is_multicontig]
    singles           = [r for r in records if not r.is_multicontig]
    high_conf_singles = [r for r in singles if r.quality_tier in _HIGH_CONFIDENCE_TIERS]
    mq_singles        = [r for r in singles if r.quality_tier not in _HIGH_CONFIDENCE_TIERS]

    taxon_buckets: dict[str, list[ContigRecord]] = defaultdict(list)
    no_taxon: list[ContigRecord] = []

    for rec in mq_singles:
        taxon = tax_map.get(rec.contig_id)
        if not taxon:
            no_taxon.append(rec)
            continue
        taxon_buckets[taxon].append(rec)

    result = list(already_multi) + list(high_conf_singles) + list(no_taxon)

    if high_conf_singles:
        log_info(
            f"  [co-bin] {len(high_conf_singles)} Complete/HQ kept as "
            "individual candidates (Turner et al. 2021: 95% ANI species boundary)"
        )

    for taxon, recs in taxon_buckets.items():
        if len(recs) == 1:
            result.append(recs[0])
        else:
            total_len = sum(r.length_bp for r in recs)
            total_vgn = sum(r.viral_genes for r in recs)
            log_info(
                f"  [co-bin] {taxon}: {len(recs)} MQ contigs "
                f"({total_len:,}bp) -> multi-contig candidate"
            )
            bin_id = "__bin__" + "||".join(r.contig_id for r in recs)
            result.append(ContigRecord(
                contig_id        = bin_id,
                sequence         = "",
                length_bp        = total_len,
                completeness     = max(r.completeness for r in recs),
                completeness_str = f"{max(r.completeness for r in recs):.1f}",
                quality_tier     = _best_quality_tier(r.quality_tier for r in recs),
                is_provirus      = any(r.is_provirus for r in recs),
                viral_genes      = total_vgn,
                is_multicontig   = True,
                sub_records      = recs,
            ))

    return result


# -- Step 4b: Rename + write ---------------------------------------------------

def _build_rename_map(records, tax_map):
    counters:   dict[str, int] = defaultdict(int)
    rename_map: list[dict]     = []

    for rec in records:
        if rec.is_multicontig:
            taxon    = (tax_map.get(rec.sub_records[0].contig_id) or "Unknown"
                        if rec.sub_records else "Unknown")
            n_ctg    = len(rec.sub_records)
            orig_ids = ",".join(r.contig_id for r in rec.sub_records)
        else:
            taxon    = tax_map.get(rec.contig_id) or "Unknown"
            n_ctg    = 1
            orig_ids = rec.contig_id

        counters[taxon] += 1
        new_id = f"{taxon}_candidate_{counters[taxon]:03d}"

        log_info(
            f"  [rename] {orig_ids[:60]} -> {new_id}  "
            f"({taxon}, {n_ctg} contig(s), {rec.completeness_str})"
        )
        rename_map.append({
            "new_id":       new_id,
            "original_id":  orig_ids,
            "taxon":        taxon,
            "length_bp":    rec.length_bp,
            "completeness": rec.completeness_str,
            "quality_tier": rec.quality_tier,
            "is_provirus":  rec.is_provirus,
            "n_contigs":    n_ctg,
        })

    return rename_map


def _write_candidates(records, rename_map, phage_dir, prov_dir, tmp_dir, log_file):
    written: list[Path] = []
    for rec, mapping in zip(records, rename_map):
        target_dir = prov_dir if rec.is_provirus else phage_dir
        out_fasta  = target_dir / f"{mapping['new_id']}.fasta"
        tmp_fasta  = tmp_dir    / f"{mapping['new_id']}_raw.tmp"
        new_id     = mapping["new_id"]

        if rec.is_multicontig:
            pairs = [
                (f"{new_id}_ctg{i + 1:03d}", r.sequence)
                for i, r in enumerate(rec.sub_records)
                if r.sequence
            ]
        else:
            pairs = [(new_id, rec.sequence)]

        _write_fasta_pairs(pairs, tmp_fasta)
        _format_seqkit(tmp_fasta, out_fasta, log_file)
        tmp_fasta.unlink(missing_ok=True)
        written.append(out_fasta)

    return written


# -- Taxonomy helpers ----------------------------------------------------------

def _find_genomad_tsv(cfg: Config, sample_id: str) -> Path:
    """
    Locate geNomad virus_summary.tsv.
    Primary: read genomad_tax_tsv from viral-id's genomad_summary.tsv.
    Fallback: convention-based path.
    """
    genomad_rpt = cfg.reports("04_viral_id") / "genomad_summary.tsv"
    if genomad_rpt.exists():
        try:
            with open(genomad_rpt) as f:
                header = f.readline().strip().split("\t")
                col    = {h: i for i, h in enumerate(header)}
                if "genomad_tax_tsv" in col:
                    idx = col["genomad_tax_tsv"]
                    for line in f:
                        parts = line.strip().split("\t")
                        if parts and parts[0] == sample_id and idx < len(parts):
                            tsv_path = Path(parts[idx])
                            if tsv_path.exists():
                                return tsv_path
                            log_info(
                                f"  [taxonomy] reported path not found: "
                                f"{tsv_path} -- trying fallback"
                            )
        except Exception as e:
            log_info(f"  [taxonomy] could not read genomad_summary.tsv: {e}")

    prefix   = f"{sample_id}_contigs_nr"
    fallback = (
        cfg.results("04_viral_id")
        / sample_id
        / f"{prefix}_summary"
        / f"{prefix}_virus_summary.tsv"
    )
    if not fallback.exists():
        log_warn(
            f"  [taxonomy] geNomad TSV not found at {fallback}. "
            "Co-binning disabled. Run viral-id before quality."
        )
    return fallback


def _load_genomad_taxonomy(tsv: Path) -> dict[str, str | None]:
    tax_map: dict[str, str | None] = {}
    if not tsv or not tsv.exists():
        return tax_map
    try:
        with open(tsv) as f:
            raw_header = f.readline().strip()
            if not raw_header:
                return tax_map
            header  = raw_header.split("\t")
            col     = {h: i for i, h in enumerate(header)}
            id_col  = next(
                (c for c in ("seq_name", "contig_id", "sequence_name") if c in col),
                None,
            )
            if id_col is None:
                log_warn(
                    f"  [taxonomy] no ID column in {tsv.name} "
                    f"({header[:6]}). Co-binning disabled."
                )
                return tax_map
            id_idx  = col[id_col]
            tax_idx = col.get("taxonomy", -1)

            for line in f:
                if not line.strip():
                    continue
                parts = line.strip().split("\t")
                if id_idx >= len(parts):
                    continue
                cid   = parts[id_idx].split("|")[0].strip()
                taxon: str | None = None
                if 0 <= tax_idx < len(parts):
                    raw = parts[tax_idx].strip()
                    if raw.upper() not in ("NA", "N/A", ""):
                        levels = [
                            lvl.strip() for lvl in raw.split(";")
                            if lvl.strip()
                            and "unclassified" not in lvl.strip().lower()
                        ]
                        if levels:
                            taxon = (
                                levels[-1]
                                .replace(" ", "_").replace("/", "_")
                                .replace("\\", "_").replace(":", "")
                                .replace(";", "")
                            )
                tax_map[cid] = taxon

    except Exception as e:
        log_warn(f"  [taxonomy] could not parse {tsv}: {e}")

    n_classified = sum(1 for v in tax_map.values() if v)
    log_ok(
        f"  [taxonomy] {len(tax_map)} contigs, "
        f"{n_classified} with family-level classification"
    )
    return tax_map


def _load_provirus_seqs(prov_fna: Path) -> dict[str, str]:
    if not prov_fna.exists():
        return {}
    seqs = _load_fasta(prov_fna)
    if seqs:
        log_info(f"  [CheckV] {len(seqs)} provirus sequence(s) loaded")
    return seqs


def _best_quality_tier(tiers) -> str:
    tier_list = list(tiers)
    for t in _TIER_ORDER:
        if t in tier_list:
            return t
    return tier_list[0] if tier_list else "Not-determined"


# -- FASTA I/O -----------------------------------------------------------------

def _load_fasta(path: Path) -> dict[str, str]:
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


def _write_fasta_records(records: list[ContigRecord], out_path: Path) -> None:
    with open(out_path, "w") as f:
        for rec in records:
            if rec.sequence and not rec.is_multicontig:
                f.write(f">{rec.contig_id}\n{rec.sequence}\n")


def _write_fasta_pairs(pairs: list[tuple[str, str]], out_path: Path) -> None:
    with open(out_path, "w") as f:
        for hdr, seq in pairs:
            if seq:
                f.write(f">{hdr}\n{seq}\n")


def _col(row: list, col: dict, name: str, fallback: int) -> str:
    idx = col.get(name, fallback)
    return row[idx] if idx < len(row) else ""


def _format_seqkit(in_fasta: Path, out_fasta: Path, log_f: Path) -> None:
    log_f.parent.mkdir(parents=True, exist_ok=True)
    try:
        with open(out_fasta, "w") as fout, open(log_f, "a") as ferr:
            result = subprocess.run(
                ["seqkit", "seq", "-w", "60", str(in_fasta)],
                stdout=fout, stderr=ferr,
            )
        if result.returncode != 0:
            log_warn("  [seqkit] non-zero exit -- copying without line-wrap")
            shutil.copy(in_fasta, out_fasta)
    except Exception as e:
        log_warn(f"  [seqkit] warning: {e} -- copying without line-wrap")
        shutil.copy(in_fasta, out_fasta)


def _save_rename_map(rename_map: list[dict], path: Path) -> None:
    headers = [
        "new_id", "original_id", "taxon", "length_bp",
        "completeness", "quality_tier", "n_contigs",
    ]
    with open(path, "w") as f:
        f.write("\t".join(headers) + "\n")
        for row in rename_map:
            f.write("\t".join(str(row.get(h, "")) for h in headers) + "\n")


# -- Validation ----------------------------------------------------------------

def _validate_output(
    sample_id: str,
    hq_fastas: list[Path],
    draft_records: list[ContigRecord],
) -> None:
    if hq_fastas:
        log_ok(f"  validate . {len(hq_fastas)} HQ genome(s) annotation-ready -- OK")
    elif draft_records:
        log_warn(
            f"  validate . 0 HQ genomes -- {len(draft_records)} draft(s) saved. "
            "Lower checkv.min_completeness or check min_bin_rescue_bp."
        )
    else:
        log_warn(
            f"  validate . 0 genomes passed quality filters for {sample_id}. "
            "Check geNomad output (Module 04) and assembly contiguity."
        )


# -- Rich display helpers -------------------------------------------------------

def _color_completeness(val: float) -> str:
    s = f"{val:.1f}%"
    if val >= _COMPLETENESS_GOOD: return f"[bold green]{s}[/bold green]"
    if val >= _COMPLETENESS_WARN: return f"[bold yellow]{s}[/bold yellow]"
    return f"[bold red]{s}[/bold red]"


def _color_contamination(val: float) -> str:
    s = f"{val:.1f}%"
    if val == 0.0:          return f"[bold green]{s}[/bold green]"
    if val <= _CONTAM_WARN: return f"[bold yellow]{s}[/bold yellow]"
    return f"[bold red]{s}[/bold red]"


def _fmt(val, unit: str = "") -> str:
    if val in ("N/A", "", None):
        return "[dim]N/A[/dim]"
    return f"{val}{unit}"


def _print_input_table(sample_id: str, virus_fna: Path, stats: dict) -> None:
    log_info(f"  Virus FASTA : {virus_fna}  ({human_size(virus_fna)})")
    log_info(
        f"  Input       : {stats['n']} contig(s) | "
        f"total={stats['total_bp']:,}bp | largest={stats['largest_bp']:,}bp | "
        f"N50={stats['n50']:,}bp"
    )


def _print_summary_table(sample_id, summary, hq_fastas, draft_records):
    if not summary:
        return
    log_ok(
        f"  Tiers    : Complete={summary.get('complete', 0)}  "
        f"HQ={summary.get('hq', 0)}  MQ={summary.get('mq', 0)}  "
        f"LQ={summary.get('lq', 0)}  ND={summary.get('nd', 0)}"
    )
    rescued = sum(summary.get(k, 0) for k in (
        "rescued_lq", "rescued_nd", "rescued_nd_density",
        "rescued_large_nd", "rescued_draft_bins",
    ))
    if rescued:
        log_ok(
            f"  Rescued  : LQ={summary.get('rescued_lq', 0)}  "
            f"ND-mod={summary.get('rescued_nd', 0)}  "
            f"ND-dens={summary.get('rescued_nd_density', 0)}  "
            f"ND-large->HQ={summary.get('rescued_large_nd', 0)}  "
            f"draft-bins={summary.get('rescued_draft_bins', 0)}"
        )
    log_ok(
        f"  Selected : {len(hq_fastas)} -> annotation_ready/  |  "
        f"{len(draft_records)} -> drafts/"
    )
    log_ok(
        f"  Best     : {summary.get('best_length_bp', 0):,}bp  |  "
        f"completeness={_color_completeness(summary.get('best_completeness', 0.0))}  |  "
        f"{_fmt(summary.get('best_quality_tier', ''))}"
    )
    log_ok(
        f"  Genes    : viral={_fmt(summary.get('total_viral_genes', 0))}  "
        f"host={_fmt(summary.get('total_host_genes', 0))}  |  "
        f"max contamination="
        f"{_color_contamination(summary.get('max_contamination', 0.0))}"
    )


def _check_warnings(summary: dict, min_completeness: float) -> None:
    if not summary:
        return
    if summary.get("hq_promoted", 0) == 0 and summary.get("draft_saved", 0) > 0:
        log_warn(
            f"  No genomes at {min_completeness}% threshold -- "
            f"{summary['draft_saved']} draft(s) saved. "
            "Lower checkv.min_completeness or check min_bin_rescue_bp."
        )
    elif summary.get("hq_promoted", 0) == 0 and not summary.get("rescued_draft_bins", 0):
        log_warn("  No genomes passed quality filters -- check geNomad output.")
    if 0 < summary.get("best_completeness", 0.0) < _COMPLETENESS_WARN:
        log_warn(
            f"  Best completeness {summary['best_completeness']:.1f}% "
            f"< {_COMPLETENESS_WARN}% -- fragmented genome."
        )
    if summary.get("total_host_genes", 0) > _HOST_GENE_WARN:
        log_warn(
            f"  {summary['total_host_genes']} host gene(s) -- "
            "possible contamination. Check Module 02."
        )
    if summary.get("max_contamination", 0.0) > _CONTAM_WARN:
        log_warn(
            f"  Max contamination {summary['max_contamination']:.1f}% "
            f"> {_CONTAM_WARN}% -- inspect quality_summary.tsv."
        )


def _print_completion_panel(
    sample_id, phage_dir, prov_dir, draft_dir, rpt_dir,
    summary, hq_fastas, draft_records,
) -> None:
    text = Text()
    text.append("v ", style="bold green")
    text.append(
        f"Complete={summary.get('complete', 0)}  HQ={summary.get('hq', 0)}  "
        f"MQ={summary.get('mq', 0)}  LQ={summary.get('lq', 0)}  "
        f"ND={summary.get('nd', 0)}  ->  ",
        style="dim white",
    )
    text.append(f"{len(hq_fastas)}", style="bold green")
    text.append(" genome(s) annotation-ready\n\n", style="cyan")
    if hq_fastas:
        text.append("Phages Dir : ", style="dim white")
        text.append(str(phage_dir) + "\n", style="white")
        text.append("Proviruses : ", style="dim white")
        text.append(str(prov_dir) + "\n", style="white")
        text.append("Rename map : ", style="dim white")
        text.append(str(rpt_dir / f"{sample_id}_rename_map.tsv") + "\n", style="white")
    if draft_records:
        text.append("Drafts     : ", style="dim white")
        text.append(str(draft_dir / f"{sample_id}_draft.fasta") + "\n", style="white")
    text.append("Summary    : ", style="dim white")
    text.append(str(rpt_dir / "checkv_summary.tsv"), style="white")

    console.print(Panel(
        text,
        title=f"[bold cyan]Quality complete -- {sample_id}[/bold cyan]",
        border_style="cyan",
        padding=(0, 2),
        width=90,
    ))


# -- TSV summary ---------------------------------------------------------------

_TSV_HEADERS = [
    "sample",
    "complete", "hq", "mq", "lq", "nd",
    "rescued_lq", "rescued_nd", "rescued_nd_density", "rescued_large_nd",
    "rescued_draft_bins", "hq_promoted", "draft_saved",
    "best_length_bp", "best_completeness", "best_quality_tier",
    "total_viral_genes", "total_host_genes", "max_contamination",
]


def _save_tsv(summary: dict, path: Path) -> None:
    rows: dict[str, dict] = {}
    if path.exists():
        with open(path) as f:
            old_hdrs = f.readline().rstrip("\n").split("\t")
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if cols and cols[0]:
                    rows[cols[0]] = dict(zip(old_hdrs, cols))
    rows[summary["sample"]] = {h: str(summary.get(h, "")) for h in _TSV_HEADERS}
    with open(path, "w") as f:
        f.write("\t".join(_TSV_HEADERS) + "\n")
        for row in rows.values():
            f.write("\t".join(row.get(h, "") for h in _TSV_HEADERS) + "\n")

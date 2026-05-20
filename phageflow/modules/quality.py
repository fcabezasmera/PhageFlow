"""PhageFlow Module 05 — Genome quality assessment and selection.

Tools and parameters (literature-based for purified phage preparations):

CheckV v1.0+ (Nayfach et al. 2021, Nature Biotechnology 39:578–585):
    End-to-end viral genome quality assessment.

    Quality tiers:
        Complete        : circular/DTR/ITR-confirmed genome (~100% completeness)
        High-quality    : ≥90% completeness, <5% contamination
        Medium-quality  : 50–90% completeness
        Low-quality     : <50% completeness with ≥1 viral gene
        Not-determined  : no database relatives for completeness estimation

    Rescue paths in purified_phage mode:
    ┌──────────────────────────────────────────────────────┬────────────────────┐
    │ Condition                                            │ Destination        │
    ├──────────────────────────────────────────────────────┼────────────────────┤
    │ Complete / High-quality                              │ annotation_ready/  │
    │ MQ ≥ min_completeness (default 50%)                 │ annotation_ready/  │
    │ MQ < min_completeness                                │ drafts/            │
    │ LQ + ≥ min_viral_genes (1) OR density ≥ 0.5 /kb    │ drafts/            │
    │ ND + ≥ large_nd_rescue_bp (30 kb) + ≥ 1 viral gene │ annotation_ready/ ↑│
    │ ND + ≥ length_rescue (10 kb) + ≥ 1 viral gene       │ drafts/            │
    │ ND + gene density ≥ 0.5 /kb                         │ drafts/            │
    │ < min_contig_bp (1,500 bp)                           │ discarded          │
    └──────────────────────────────────────────────────────┴────────────────────┘

    Rationale for large-ND rescue (≥30 kb):
        CheckV completeness relies on HMM profiles derived from reference
        sequences. Herelleviridae (>100 kb) and Ackermannviridae (>150 kb)
        lacking close database relatives are consistently Not-determined
        regardless of actual completeness. In a purified phage preparation,
        a contig ≥30 kb with ≥1 viral gene is almost certainly phage.
        (Camargo et al. 2023, Nat Biotechnol; Adriaenssens & Brister 2017, Viruses)

    Rationale for gene-density threshold (0.5 genes/kb):
        Phage genomes average 1.5–2 genes/kb; 0.5/kb is a conservative lower
        bound that filters assembly artefacts while retaining large-ND sequences
        with sparse hallmark coverage. (Roux et al. 2019, eLife)

Multi-contig bin assembly — confidence criteria:
    Contigs are grouped into a single candidate genome when ALL of the
    following hold:
      1. Each contig independently classified as viral by geNomad (score ≥0.7),
         giving ~97% per-contig precision (Camargo et al. 2023).
      2. All contigs share the same geNomad family-level taxonomy, indicating
         a common viral lineage.
      3. Combined length ≥ min_bin_rescue_bp (30 kb) for draft rescue, or
         co-binned at any combined length for HQ genomes.
      4. At least 1 viral gene across the bin (viral marker confirmed).
    In a purified phage preparation, co-elution of two phages from the same
    family is uncommon, making shared-taxonomy co-binning reliable.
    Limitation: two co-purified phages of the same family would be merged.
    Downstream CheckV contamination and Pharokka annotation can reveal this.
    (Buttimer et al. 2020, Front Microbiol 11:710;
     Adriaenssens et al. 2020, Viruses 12:955)

Draft co-bin rescue:
    Contigs individually below threshold sharing a geNomad taxonomy are grouped.
    If combined length ≥ min_bin_rescue_bp, the bin is promoted to
    annotation_ready/ for Pharokka --meta processing.

cd-hit-est (Fu et al. 2012, Bioinformatics 28:3150–3152):
    -c 0.98: intra-sample dereplication at 98% ANI (single-contig only).
    ICTV species boundary = 95% ANI (Turner et al. 2021, Arch Virol 166:2633);
    98% removes assembly duplicates while preserving biological strain diversity.

seqkit (Shen et al. 2016, PLoS ONE 11:e0163962):
    seq -w 60: standard FASTA line-wrapping for downstream tool compatibility.
"""

from __future__ import annotations
import shutil
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
import subprocess

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
TOOLS = ["checkv", "cd-hit-est", "seqkit"]

# ── Quality thresholds ────────────────────────────────────────────────────────
_COMPLETENESS_GOOD = 90.0   # HQ lower bound (Nayfach et al. 2021)
_COMPLETENESS_WARN = 50.0   # MQ lower bound (Nayfach et al. 2021)
_CONTAM_WARN       = 5.0    # maximum acceptable contamination % (Nayfach et al. 2021)
_HOST_GENE_WARN    = 3      # host gene count triggering a warning

# cd-hit-est dereplication parameters (Turner et al. 2021, Arch Virol 166:2633)
_CDHIT_ID = "0.98"   # 98% ANI — removes assembly duplicates, preserves strains
_CDHIT_AS = "0.85"   # minimum alignment fraction of shorter sequence
_CDHIT_N  = "8"      # word length (correct for ≥0.90 identity)

# CheckV tier ranking (most to least confident)
_TIER_ORDER = ["Complete", "High-quality", "Medium-quality", "Low-quality", "Not-determined"]


# ── Data structure ────────────────────────────────────────────────────────────

@dataclass
class ContigRecord:
    """Represents a viral contig or assembled multi-contig genomic bin.

    For single-contig records: sequence holds DNA, sub_records is empty.
    For multi-contig bins (is_multicontig=True): sequence is an empty string
    and sub_records holds the constituent ContigRecord objects. Their sequences
    are written as a multi-FASTA for downstream Pharokka --meta processing.

    FASTA headers written per record:
        Single  : >{new_id}                        e.g. >Herelleviridae_candidate_001
        Multi   : >{new_id}_ctg001, _ctg002 ...    one line per constituent contig
    """
    contig_id:      str
    sequence:       str    # DNA; empty string for multi-contig bins
    length_bp:      int
    completeness:   float  # 0.0 when Not-determined
    quality_tier:   str
    is_provirus:    bool
    viral_genes:    int
    is_multicontig: bool = False
    sub_records:    list  = field(default_factory=list)  # List[ContigRecord]


# ── Public entry point ────────────────────────────────────────────────────────

def run(
    cfg:       Config,
    sample_id: str,
    virus_fna: Path,
    force:     bool = False,
) -> list[Path]:
    """
    Assess genome quality with CheckV, rescue, dereplicate, co-bin, and format.

    Pipeline
    --------
    CheckV end-to-end
        → tier selection + rescue paths
        → draft co-bin rescue          (multi-contig bins from drafts/)
        → cd-hit-est 98% ANI           (single-contig HQ only)
        → HQ co-binning by geNomad taxon
        → rename + seqkit -w 60
        → one FASTA per candidate genome

    Returns
    -------
    hq_fastas : list of annotation-ready FASTA paths
    """
    require_tools(*TOOLS)

    virus_fna = Path(virus_fna)
    out_dir   = cfg.results(STEP)
    rpt_dir   = cfg.reports(STEP)
    phage_dir = out_dir / "annotation_ready" / "phages"
    prov_dir  = out_dir / "annotation_ready" / "proviruses"
    draft_dir = out_dir / "drafts"
    tmp_dir   = out_dir / "tmp"
    mkdirs(out_dir, rpt_dir, phage_dir, prov_dir, draft_dir, tmp_dir)

    log_step(f"Module 05 — quality [{sample_id}]")

    if not virus_fna.exists():
        log_error(f"  Virus FASTA not found: {virus_fna}")
        raise FileNotFoundError(virus_fna)

    s = fasta_stats(virus_fna)
    _print_input_table(sample_id, virus_fna, s)
    log_info(
        f"  CheckV end-to-end  (Nayfach et al. 2021)  |  "
        f"HQ threshold ≥{cfg.checkv.min_completeness}%  "
        f"(Complete/HQ always promoted)"
    )
    log_info(
        f"  min_contig_bp={cfg.checkv.min_contig_bp} bp  |  "
        f"LQ rescue: ≥{cfg.checkv.min_viral_genes} viral gene(s) or "
        f"density ≥{cfg.checkv.min_gene_density}/kb  |  "
        f"ND moderate: ≥{cfg.checkv.length_rescue:,} bp → drafts/  |  "
        f"ND large: ≥{cfg.checkv.large_nd_rescue_bp:,} bp → annotation_ready/"
    )
    log_info(
        f"  Draft co-bin rescue: ≥{cfg.checkv.min_bin_rescue_bp:,} bp → annotation_ready/  |  "
        f"cd-hit-est {_CDHIT_ID} ANI  (Turner et al. 2021)"
    )

    if not cfg.databases.checkv.exists():
        log_warn(f"  CheckV database not found: {cfg.databases.checkv}")

    # Pre-load geNomad taxonomy map for co-binning and renaming
    genomad_tsv = _find_genomad_tsv(cfg, sample_id)
    tax_map     = _load_genomad_taxonomy(genomad_tsv)

    sdir  = out_dir / sample_id
    qfile = sdir / "quality_summary.tsv"

    summary:       dict               = {}
    hq_records:    list[ContigRecord] = []
    draft_records: list[ContigRecord] = []
    promoted_bins: list[ContigRecord] = []
    hq_fastas:     list[Path]         = []

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

        # 1/4 — CheckV end-to-end
        progress.update(task, description="[1/4] CheckV   — end-to-end assessment")
        _run_checkv(cfg, sample_id, virus_fna, sdir, rpt_dir, force)
        progress.advance(task)

        # 2/4 — Tier selection + rescue + draft co-bin rescue
        progress.update(task, description="[2/4] CheckV   — tier selection + rescue")
        if qfile.exists():
            seqs      = _load_fasta(virus_fna)
            prov_seqs = _load_provirus_seqs(sdir / "proviruses.fna")
            summary, hq_records, draft_records = _parse_checkv(
                sample_id, qfile, seqs, prov_seqs, cfg.checkv,
            )
            # Group draft contigs by geNomad taxonomy; promote bins ≥ threshold
            promoted_bins, draft_records = _cobin_draft_rescue(
                draft_records, tax_map, sample_id, cfg.checkv.min_bin_rescue_bp,
            )
            summary["rescued_draft_bins"] = len(promoted_bins)
            _write_fasta_records(draft_records, draft_dir / f"{sample_id}_draft.fasta")
            log_ok("  [CheckV] genome selection complete")
        else:
            log_warn(
                "  [CheckV] no output found — "
                f"check reports/{STEP}/{sample_id}_checkv.log"
            )
        progress.advance(task)

        # 3/4 — cd-hit-est dereplication (single-contig HQ only)
        progress.update(task, description="[3/4] cd-hit-est — dereplication 98% ANI")
        if hq_records:
            hq_records, n_before, n_after = _dereplicate(
                hq_records, tmp_dir, sample_id, rpt_dir, cfg.threads, force,
            )
            log_ok(
                f"  [cd-hit-est] {n_before} → {n_after} single-contig genomes  "
                f"(−{n_before - n_after} redundant at {_CDHIT_ID} ANI)"
            )
        # Append promoted bins after dereplication — unique by construction
        hq_records.extend(promoted_bins)
        if promoted_bins:
            log_ok(f"  [draft rescue] {len(promoted_bins)} bin(s) promoted to annotation pool")
        progress.advance(task)

        # 4/4 — Co-binning by taxon + rename + seqkit format
        progress.update(task, description="[4/4] seqkit   — co-binning, renaming & formatting")
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
                f"  [seqkit] {len(hq_fastas)} file(s) written  "
                f"(phage={n_phage}  provirus={n_prov}  multi-contig={n_multi})"
            )
        progress.advance(task)

    shutil.rmtree(tmp_dir, ignore_errors=True)

    # ── Metrics, display, save ─────────────────────────────────────────────────
    _validate_output(sample_id, hq_fastas, draft_records)
    _print_summary_table(sample_id, summary, hq_fastas, draft_records)
    _check_warnings(summary, cfg.checkv.min_completeness)
    if summary:
        _save_tsv(summary, rpt_dir / "checkv_summary.tsv")
    _print_completion_panel(
        sample_id, phage_dir, prov_dir, draft_dir, rpt_dir,
        summary, hq_fastas, draft_records,
    )

    log_step(f"Module 05 completed ✓  [{sample_id}]")
    log_info("Next: phageflow annotate --sample-id <id> --genome <candidate.fasta>")
    return hq_fastas


# ── Step 1: CheckV ────────────────────────────────────────────────────────────

def _run_checkv(cfg, sample_id, virus_fna, sdir, rpt_dir, force):
    qfile = sdir / "quality_summary.tsv"
    if qfile.exists() and qfile.stat().st_size > 0 and not force:
        log_info("  [CheckV] already processed — skipping  (--force to re-run)")
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


# ── Step 2: Tier selection + rescue ──────────────────────────────────────────

def _parse_checkv(sample_id, qfile, seqs, prov_seqs, chkv):
    """
    Parse quality_summary.tsv and classify contigs into HQ and draft bins.
    Tier logic: Nayfach et al. 2021 + purified_phage rescue extensions.
    """
    tiers: dict[str, int] = {t: 0 for t in _TIER_ORDER}
    rescued_lq = rescued_nd = rescued_nd_density = rescued_large_nd = 0
    total_viral_genes = total_host_genes = 0
    best_len          = 0
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
                        f"  [filter] {cid} discarded — "
                        f"{clen_int} bp < min_contig_bp {chkv.min_contig_bp} bp"
                    )
                    continue

                gene_density = vgn_int / (clen_int / 1000) if clen_int > 0 else 0.0

                tiers[qual if qual in tiers else "Not-determined"] += 1
                try:   total_viral_genes += vgn_int
                except (ValueError, TypeError): pass
                try:   total_host_genes  += int(hgn_s)
                except (ValueError, TypeError): pass
                try:
                    max_contamination = max(max_contamination, float(contm_s))
                except (ValueError, TypeError): pass

                try:
                    comp_val = float(comp_s) if comp_s not in ("", "NA", "N/A") else 0.0
                except (ValueError, TypeError):
                    comp_val = 0.0

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
                    contig_id    = cid,
                    sequence     = seq,
                    length_bp    = clen_int,
                    completeness = comp_val,
                    quality_tier = qual,
                    is_provirus  = is_provirus,
                    viral_genes  = vgn_int,
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
                            f"  [LQ rescue] {cid} — "
                            f"{vgn_int} viral gene(s), density={gene_density:.2f}/kb → drafts/"
                        )
                        draft_records.append(record)
                        rescued_lq += 1

                elif qual == "Not-determined":
                    if clen_int >= chkv.large_nd_rescue_bp and vgn_int >= 1:
                        log_warn(
                            f"  [large-ND rescue] {cid} — "
                            f"{clen_int:,} bp, {vgn_int} viral gene(s), "
                            f"density={gene_density:.2f}/kb → annotation_ready/"
                        )
                        hq_records.append(record)
                        rescued_large_nd += 1
                    elif clen_int >= chkv.length_rescue and vgn_int >= 1:
                        log_warn(
                            f"  [ND rescue] {cid} — "
                            f"{clen_int:,} bp, {vgn_int} viral gene(s) → drafts/"
                        )
                        draft_records.append(record)
                        rescued_nd += 1
                    elif gene_density >= chkv.min_gene_density and clen_int >= chkv.min_contig_bp:
                        log_warn(
                            f"  [density-ND rescue] {cid} — "
                            f"density={gene_density:.2f}/kb → drafts/"
                        )
                        draft_records.append(record)
                        rescued_nd_density += 1

    except Exception as e:
        log_warn(f"  [CheckV] could not parse {qfile}: {e}")

    rescued_total = rescued_lq + rescued_nd + rescued_nd_density + rescued_large_nd
    if rescued_total:
        log_ok(
            f"  [rescue] LQ={rescued_lq}  ND-moderate={rescued_nd}  "
            f"ND-density={rescued_nd_density}  ND-large→HQ={rescued_large_nd}"
        )

    summary = {
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
    }
    return summary, hq_records, draft_records


# ── Draft co-bin rescue ───────────────────────────────────────────────────────

def _cobin_draft_rescue(draft_records, tax_map, sample_id, min_bin_rescue_bp):
    """
    Promote draft bins ≥ min_bin_rescue_bp to annotation_ready/.

    Confidence that a bin represents a single phage genome:
      - Each contig independently classified by geNomad ≥0.7 (~97% precision)
      - Shared taxonomy implies a common viral lineage
      - In purified preparations, same-family co-purification is uncommon
      - Downstream CheckV contamination and Pharokka annotation will reveal
        any incorrectly merged genomes

    References: Buttimer et al. 2020, Front Microbiol;
                Adriaenssens et al. 2020, Viruses 12:955;
                Camargo et al. 2023, Nat Biotechnol
    """
    taxon_buckets: dict[str, list[ContigRecord]] = defaultdict(list)
    for rec in draft_records:
        taxon = tax_map.get(rec.contig_id, "Unknown")
        taxon_buckets[taxon].append(rec)

    promoted:  list[ContigRecord] = []
    remaining: list[ContigRecord] = []

    for taxon, records in taxon_buckets.items():
        total_len = sum(r.length_bp for r in records)
        total_vgn = sum(r.viral_genes for r in records)

        if len(records) >= 2 and total_len >= min_bin_rescue_bp:
            log_warn(
                f"  [draft-bin rescue] {taxon}: {len(records)} contigs, "
                f"{total_len:,} bp total, {total_vgn} viral genes "
                f"→ annotation_ready/ (multi-contig candidate)"
            )
            log_info(
                f"  [bin confidence] geNomad ≥0.7 per contig + shared {taxon} "
                f"taxonomy + {total_vgn} viral genes → high confidence phage bin"
            )
            bin_id = "__bin__" + "||".join(r.contig_id for r in records)
            promoted.append(ContigRecord(
                contig_id      = bin_id,
                sequence       = "",
                length_bp      = total_len,
                completeness   = max(r.completeness for r in records),
                quality_tier   = _best_quality_tier(r.quality_tier for r in records),
                is_provirus    = any(r.is_provirus for r in records),
                viral_genes    = total_vgn,
                is_multicontig = True,
                sub_records    = records,
            ))
        else:
            remaining.extend(records)

    return promoted, remaining


# ── Step 3: Dereplication ─────────────────────────────────────────────────────

def _dereplicate(records, tmp_dir, sample_id, rpt_dir, threads, force):
    """Dereplicate single-contig HQ genomes with cd-hit-est at 98% ANI."""
    n_before = len(records)
    if n_before <= 1:
        return records, n_before, n_before

    tmp_in  = tmp_dir / f"{sample_id}_hq_raw.fasta"
    tmp_out = tmp_dir / f"{sample_id}_hq_derep.fasta"
    _write_fasta_records(records, tmp_in)

    try:
        run_silent([
            "cd-hit-est",
            "-i",  str(tmp_in),
            "-o",  str(tmp_out),
            "-c",  _CDHIT_ID,
            "-aS", _CDHIT_AS,
            "-G",  "0",
            "-n",  _CDHIT_N,
            "-d",  "0",
            "-sc", "1",
            "-T",  str(threads),
            "-M",  "4000",
        ], log_file=rpt_dir / f"{sample_id}_cdhit_derep.log")
    except Exception as e:
        log_warn(f"  [cd-hit-est] dereplication warning: {e} — using all genomes")
        return records, n_before, n_before
    finally:
        Path(str(tmp_out) + ".clstr").unlink(missing_ok=True)

    rep_ids = set(_load_fasta(tmp_out).keys())
    derep   = [r for r in records if r.contig_id in rep_ids]
    return derep, n_before, len(derep)


# ── Step 4a: Co-binning by taxon ──────────────────────────────────────────────

def _cobin_by_taxon(records, tax_map, sample_id):
    """
    Group single-contig HQ records sharing a geNomad taxon into multi-contig bins.
    Multi-contig bins already present pass through unchanged.
    """
    already_multi = [r for r in records if r.is_multicontig]
    singles       = [r for r in records if not r.is_multicontig]

    taxon_buckets: dict[str, list[ContigRecord]] = defaultdict(list)
    for rec in singles:
        taxon = tax_map.get(rec.contig_id, "Unknown")
        taxon_buckets[taxon].append(rec)

    result = list(already_multi)

    for taxon, recs in taxon_buckets.items():
        if len(recs) == 1:
            result.append(recs[0])
        else:
            total_len = sum(r.length_bp for r in recs)
            total_vgn = sum(r.viral_genes for r in recs)
            log_info(
                f"  [co-bin] {taxon}: {len(recs)} HQ contigs "
                f"({total_len:,} bp, {total_vgn} viral genes) → multi-contig candidate"
            )
            bin_id = "__bin__" + "||".join(r.contig_id for r in recs)
            result.append(ContigRecord(
                contig_id      = bin_id,
                sequence       = "",
                length_bp      = total_len,
                completeness   = max(r.completeness for r in recs),
                quality_tier   = _best_quality_tier(r.quality_tier for r in recs),
                is_provirus    = any(r.is_provirus for r in recs),
                viral_genes    = total_vgn,
                is_multicontig = True,
                sub_records    = recs,
            ))

    return result


# ── Step 4b: Rename + write ───────────────────────────────────────────────────

def _build_rename_map(records, tax_map):
    """
    Assign new_id to each record.

    Naming convention (all candidates use the same pattern):
        {Family}_candidate_{NNN}

    The rename_map TSV stores n_contigs > 1 for multi-contig bins, which
    downstream modules (annotate) use to select Pharokka mode (--meta vs --single).

    FASTA headers:
        Single  : >{Family}_candidate_{NNN}
        Multi   : >{Family}_candidate_{NNN}_ctg001  (one per contig)
    """
    counters:   dict[str, int] = defaultdict(int)
    rename_map: list[dict]     = []

    for rec in records:
        if rec.is_multicontig:
            taxon = (
                tax_map.get(rec.sub_records[0].contig_id, "Unknown")
                if rec.sub_records else "Unknown"
            )
            n_ctg    = len(rec.sub_records)
            orig_ids = ",".join(r.contig_id for r in rec.sub_records)
        else:
            taxon    = tax_map.get(rec.contig_id, "Unknown")
            n_ctg    = 1
            orig_ids = rec.contig_id

        counters[taxon] += 1
        new_id = f"{taxon}_candidate_{counters[taxon]:03d}"

        log_info(
            f"  [rename] {orig_ids[:60]} → {new_id}  "
            f"({taxon}, {n_ctg} contig(s))"
        )
        rename_map.append({
            "new_id":       new_id,
            "original_id":  orig_ids,
            "taxon":        taxon,
            "length_bp":    rec.length_bp,
            "completeness": f"{rec.completeness:.1f}",
            "quality_tier": rec.quality_tier,
            "is_provirus":  rec.is_provirus,
            "n_contigs":    n_ctg,
        })

    return rename_map


def _write_candidates(records, rename_map, phage_dir, prov_dir, tmp_dir, log_file):
    """
    Write one FASTA per candidate.

    Single-contig → FASTA with one sequence:
        >Herelleviridae_candidate_001

    Multi-contig → FASTA with one sequence per constituent contig:
        >Herelleviridae_candidate_002_ctg001
        >Herelleviridae_candidate_002_ctg002

    Multi-FASTA files are automatically processed by Pharokka in --meta mode
    (detected by counting sequences in the file, not by filename).
    """
    written: list[Path] = []

    for rec, mapping in zip(records, rename_map):
        target_dir = prov_dir if rec.is_provirus else phage_dir
        out_fasta  = target_dir / f"{mapping['new_id']}.fasta"
        tmp_fasta  = tmp_dir    / f"{mapping['new_id']}_raw.tmp"
        new_id     = mapping["new_id"]

        if rec.is_multicontig:
            # One sequence per constituent contig, named {new_id}_ctg001 etc.
            pairs = [
                (f"{new_id}_ctg{i + 1:03d}", r.sequence)
                for i, r in enumerate(rec.sub_records)
                if r.sequence
            ]
        else:
            # Single sequence — clean header, no metadata
            pairs = [(new_id, rec.sequence)]

        _write_fasta_pairs(pairs, tmp_fasta)
        _format_seqkit(tmp_fasta, out_fasta, log_file)
        tmp_fasta.unlink(missing_ok=True)
        written.append(out_fasta)

    return written


# ── Taxonomy helpers ──────────────────────────────────────────────────────────

def _find_genomad_tsv(cfg, sample_id):
    prefix = f"{sample_id}_contigs_nr"
    return (
        cfg.results("04_viral_id")
        / sample_id
        / f"{prefix}_summary"
        / f"{prefix}_virus_summary.tsv"
    )


def _load_genomad_taxonomy(tsv):
    tax_map: dict[str, str] = {}
    if not tsv.exists():
        log_warn(f"  [rename] geNomad taxonomy file not found: {tsv}")
        return tax_map
    try:
        with open(tsv) as f:
            header  = f.readline().strip().split("\t")
            col     = {h: i for i, h in enumerate(header)}
            tax_idx = col.get("taxonomy", -1)
            id_idx  = col.get("seq_name", 0)
            for line in f:
                if not line.strip():
                    continue
                parts     = line.strip().split("\t")
                cid       = parts[id_idx] if id_idx < len(parts) else ""
                cid_clean = cid.split("|")[0]
                taxon     = "Unknown"
                if 0 <= tax_idx < len(parts):
                    raw    = parts[tax_idx].strip()
                    levels = [
                        lvl.strip() for lvl in raw.split(";")
                        if lvl.strip()
                        and "unclassified" not in lvl.strip().lower()
                        and raw.upper() not in ("NA", "N/A", "")
                    ]
                    if levels:
                        taxon = levels[-1]
                taxon = (
                    taxon.replace(" ", "_").replace("/", "_")
                         .replace("\\", "_").replace(":", "").replace(";", "")
                )
                tax_map[cid_clean] = taxon
                if cid != cid_clean:
                    tax_map[cid] = taxon
    except Exception as e:
        log_warn(f"  [rename] could not parse geNomad taxonomy: {e}")
    return tax_map


def _load_provirus_seqs(prov_fna):
    if not prov_fna.exists():
        return {}
    seqs = _load_fasta(prov_fna)
    if seqs:
        log_info(
            f"  [CheckV] {len(seqs)} provirus sequence(s) loaded "
            f"from proviruses.fna (host flanks removed)"
        )
    return seqs


def _best_quality_tier(tiers):
    tier_list = list(tiers)
    for t in _TIER_ORDER:
        if t in tier_list:
            return t
    return tier_list[0] if tier_list else "Not-determined"


# ── FASTA I/O ─────────────────────────────────────────────────────────────────

def _load_fasta(path):
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


def _write_fasta_records(records, out_path):
    if not records:
        return
    with open(out_path, "w") as f:
        for rec in records:
            if rec.sequence and not rec.is_multicontig:
                f.write(f">{rec.contig_id}\n{rec.sequence}\n")


def _write_fasta_pairs(pairs, out_path):
    with open(out_path, "w") as f:
        for hdr, seq in pairs:
            if seq:
                f.write(f">{hdr}\n{seq}\n")


def _col(row, col, name, fallback):
    idx = col.get(name, fallback)
    return row[idx] if idx < len(row) else ""


def _format_seqkit(in_fasta, out_fasta, log_f):
    log_f.parent.mkdir(parents=True, exist_ok=True)
    try:
        with open(out_fasta, "w") as fout, open(log_f, "a") as ferr:
            result = subprocess.run(
                ["seqkit", "seq", "-w", "60", str(in_fasta)],
                stdout=fout, stderr=ferr,
            )
        if result.returncode != 0:
            log_warn("  [seqkit] non-zero exit — copying without line-wrap")
            shutil.copy(in_fasta, out_fasta)
    except Exception as e:
        log_warn(f"  [seqkit] warning: {e} — copying without line-wrap")
        shutil.copy(in_fasta, out_fasta)


def _save_rename_map(rename_map, path):
    headers = [
        "new_id", "original_id", "taxon", "length_bp",
        "completeness", "quality_tier", "n_contigs",
    ]
    with open(path, "w") as f:
        f.write("\t".join(headers) + "\n")
        for row in rename_map:
            f.write("\t".join(str(row.get(h, "")) for h in headers) + "\n")


# ── Validation ────────────────────────────────────────────────────────────────

def _validate_output(sample_id, hq_fastas, draft_records):
    if hq_fastas:
        log_ok(f"  validate · {len(hq_fastas)} HQ genome(s) annotation-ready — OK")
    elif draft_records:
        log_warn(
            f"  validate · 0 HQ genomes — {len(draft_records)} draft(s) saved. "
            "Consider lowering checkv.min_completeness in config.yaml."
        )
    else:
        log_warn(
            f"  validate · 0 genomes passed quality filters for {sample_id}. "
            "Check geNomad output (Module 04) and assembly contiguity."
        )


# ── Rich display helpers ──────────────────────────────────────────────────────

def _color_completeness(val):
    s = f"{val:.1f}%"
    if val >= _COMPLETENESS_GOOD: return f"[bold green]{s}[/bold green]"
    if val >= _COMPLETENESS_WARN: return f"[bold yellow]{s}[/bold yellow]"
    return f"[bold red]{s}[/bold red]"


def _color_contamination(val):
    s = f"{val:.1f}%"
    if val == 0.0:          return f"[bold green]{s}[/bold green]"
    if val <= _CONTAM_WARN: return f"[bold yellow]{s}[/bold yellow]"
    return f"[bold red]{s}[/bold red]"


def _fmt(val, unit=""):
    if val in ("N/A", "", None):
        return "[dim]N/A[/dim]"
    return f"{val}{unit}"


def _print_input_table(sample_id, virus_fna, stats):
    log_info(f"  Virus FASTA : {virus_fna}  ({human_size(virus_fna)})")
    log_info(
        f"  Input       : {stats['n']} contig(s) | "
        f"total={stats['total_bp']:,} bp | largest={stats['largest_bp']:,} bp | "
        f"N50={stats['n50']:,} bp"
    )


def _print_summary_table(sample_id, summary, hq_fastas, draft_records):
    if not summary:
        return
    log_ok(
        f"  Tiers    : Complete={summary.get('complete', 0)}  "
        f"High-quality={summary.get('hq', 0)}  "
        f"Medium-quality={summary.get('mq', 0)}  "
        f"Low-quality={summary.get('lq', 0)}  "
        f"ND={summary.get('nd', 0)}"
    )
    rescued = (
        summary.get("rescued_lq",         0) +
        summary.get("rescued_nd",         0) +
        summary.get("rescued_nd_density", 0) +
        summary.get("rescued_large_nd",   0) +
        summary.get("rescued_draft_bins", 0)
    )
    if rescued:
        log_ok(
            f"  Rescued  : LQ={summary.get('rescued_lq', 0)}  "
            f"ND-mod={summary.get('rescued_nd', 0)}  "
            f"ND-dens={summary.get('rescued_nd_density', 0)}  "
            f"ND-large→HQ={summary.get('rescued_large_nd', 0)}  "
            f"draft-bins→HQ={summary.get('rescued_draft_bins', 0)}"
        )
    log_ok(
        f"  Selected : {len(hq_fastas)} → annotation_ready/  |  "
        f"{len(draft_records)} → drafts/"
    )
    log_ok(
        f"  Best     : {summary.get('best_length_bp', 0):,} bp  |  "
        f"completeness={_color_completeness(summary.get('best_completeness', 0.0))}  |  "
        f"{_fmt(summary.get('best_quality_tier', ''))}"
    )
    log_ok(
        f"  Genes    : viral={_fmt(summary.get('total_viral_genes', 0))}  "
        f"host={_fmt(summary.get('total_host_genes', 0))}  |  "
        f"max contamination={_color_contamination(summary.get('max_contamination', 0.0))}"
    )


def _check_warnings(summary, min_completeness):
    if not summary:
        return
    if summary.get("hq_promoted", 0) == 0 and summary.get("draft_saved", 0) > 0:
        log_warn(
            f"  No genomes met the {min_completeness}% threshold — "
            f"{summary['draft_saved']} draft(s) saved. "
            "Lower checkv.min_completeness or check min_bin_rescue_bp."
        )
    elif summary.get("hq_promoted", 0) == 0 and not summary.get("rescued_draft_bins", 0):
        log_warn(
            "  No genomes passed quality filters — "
            "verify geNomad output or assembly contiguity."
        )
    if 0 < summary.get("best_completeness", 0.0) < _COMPLETENESS_WARN:
        log_warn(
            f"  Best completeness {summary['best_completeness']:.1f}% < {_COMPLETENESS_WARN}% — "
            "fragmented genome; annotation may be affected."
        )
    if summary.get("total_host_genes", 0) > _HOST_GENE_WARN:
        log_warn(
            f"  {summary['total_host_genes']} host gene(s) detected — "
            "possible contamination. Check Module 02 (host-removal)."
        )
    if summary.get("max_contamination", 0.0) > _CONTAM_WARN:
        log_warn(
            f"  Max contamination {summary['max_contamination']:.1f}% > {_CONTAM_WARN}% — "
            "inspect quality_summary.tsv."
        )


def _print_completion_panel(
    sample_id, phage_dir, prov_dir, draft_dir, rpt_dir,
    summary, hq_fastas, draft_records,
):
    text = Text()
    text.append("✓ ", style="bold green")
    text.append(
        f"Complete={summary.get('complete', 0)}  "
        f"HQ={summary.get('hq', 0)}  "
        f"MQ={summary.get('mq', 0)}  "
        f"LQ={summary.get('lq', 0)}  "
        f"ND={summary.get('nd', 0)}  →  ",
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
        title=f"[bold cyan]Quality complete — {sample_id}[/bold cyan]",
        border_style="cyan",
        padding=(0, 2),
        width=90,
    ))


# ── TSV summary ───────────────────────────────────────────────────────────────

_TSV_HEADERS = [
    "sample",
    "complete", "hq", "mq", "lq", "nd",
    "rescued_lq", "rescued_nd", "rescued_nd_density", "rescued_large_nd",
    "rescued_draft_bins",
    "hq_promoted", "draft_saved",
    "best_length_bp", "best_completeness", "best_quality_tier",
    "total_viral_genes", "total_host_genes", "max_contamination",
]


def _save_tsv(summary, path):
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

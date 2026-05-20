"""PhageFlow Module 05 — Genome quality assessment and selection"""

from __future__ import annotations
import shutil
from collections import defaultdict
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

# ── Thresholds ────────────────────────────────────────────────────────────────
_COMPLETENESS_GOOD  = 90.0
_COMPLETENESS_WARN  = 50.0
_CONTAM_WARN        = 5.0
_HOST_GENE_WARN     = 3
_GENE_DENSITY_GOOD  = 1.0   # genes virales/kb: señal fuerte (Roux 2019)

# cd-hit-est: 0.98 ANI intra-muestra (Turner et al. 2021 Arch Virol 166:2633)
_CDHIT_ID = "0.98"
_CDHIT_AS = "0.85"
_CDHIT_N  = "8"

_TIER_ORDER = ["Complete", "High-quality", "Medium-quality", "Low-quality", "Not-determined"]

# Sentinel para identificar multi-contig en cualquier tamaño de tupla
_MULTI_SENTINEL = "__MULTICONTIG__"


# ── Public entry point ────────────────────────────────────────────────────────

def run(
    cfg:       Config,
    sample_id: str,
    virus_fna: Path,
    force:     bool = False,
) -> list[Path]:
    """
    Assess genome quality with CheckV, dereplicate, co-bin, rename, and format.

    Pipeline
    --------
    CheckV → tier selection → draft co-bin rescue → cd-hit-est 98% ANI →
    co-binning HQ por taxón → taxonomy rename → seqkit -w 60 → one FASTA per candidate
    """
    require_tools(*TOOLS)

    virus_fna  = Path(virus_fna)
    out_dir    = cfg.results(STEP)
    rpt_dir    = cfg.reports(STEP)
    phage_dir  = out_dir / "annotation_ready" / "phages"
    prov_dir   = out_dir / "annotation_ready" / "proviruses"
    draft_dir  = out_dir / "drafts"
    tmp_dir    = out_dir / "tmp"
    mkdirs(out_dir, rpt_dir, phage_dir, prov_dir, draft_dir, tmp_dir)

    log_step(f"Module 05 — quality [{sample_id}] · {cfg.threads} threads")

    if not virus_fna.exists():
        log_error(f"  Virus FASTA not found: {virus_fna}")
        raise FileNotFoundError(virus_fna)

    s = fasta_stats(virus_fna)
    _print_input_table(sample_id, virus_fna, s)
    log_info(
        f"  CheckV end-to-end  (Nayfach et al. 2021)  |  "
        f"HQ threshold ≥ {cfg.checkv.min_completeness}%  "
        f"(Complete/HQ always promoted)"
    )
    log_info(
        f"  min_contig_bp={cfg.checkv.min_contig_bp}bp  |  "
        f"LQ rescue: ≥{cfg.checkv.min_viral_genes} gene(s) or "
        f"≥{cfg.checkv.min_gene_density}/kb  |  "
        f"ND rescue: ≥{cfg.checkv.length_rescue:,}bp → drafts/  |  "
        f"large-ND: ≥{cfg.checkv.large_nd_rescue_bp:,}bp → annotation_ready/"
    )
    log_info(
        f"  Draft co-bin rescue: ≥{cfg.checkv.min_bin_rescue_bp:,}bp total → annotation_ready/  |  "
        f"cd-hit-est {_CDHIT_ID} ANI  (Turner et al. 2021)"
    )

    if not cfg.databases.checkv.exists():
        log_warn(f"  CheckV database not found: {cfg.databases.checkv}")

    # Pre-carga del mapa taxonómico (necesario para draft rescue en paso 2)
    genomad_tsv = _find_genomad_tsv(cfg, sample_id)
    tax_map     = _load_genomad_taxonomy(genomad_tsv)

    sdir        = out_dir / sample_id
    qfile       = sdir / "quality_summary.tsv"
    row:         dict = {}
    hq_seqs:     list = []
    draft_seqs:  list = []
    extra_bins:  list = []   # bins promovidos desde draft rescue (multi-contig)
    hq_fastas:   list = []

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

        # 2/4 — Parse tiers + draft co-bin rescue
        progress.update(task,
            description="[2/4] CheckV   — tier selection + draft bin rescue")
        if qfile.exists():
            seqs      = _load_fasta(virus_fna)
            prov_seqs = _load_provirus_seqs(sdir / "proviruses.fna")
            row, hq_seqs, draft_seqs = _parse_checkv(
                sample_id, qfile, seqs, prov_seqs,
                cfg.checkv.min_completeness,
                cfg.checkv.min_viral_genes,
                cfg.checkv.length_rescue,
                cfg.checkv.min_contig_bp,
                cfg.checkv.large_nd_rescue_bp,
                cfg.checkv.min_gene_density,
            )
            # Draft co-bin rescue: agrupa drafts por taxón geNomad y promueve
            # bins ≥ min_bin_rescue_bp a annotation_ready/ (Buttimer 2020;
            # Adriaenssens 2020)
            extra_bins, draft_seqs = _cobin_draft_rescue(
                draft_seqs, tax_map, sample_id,
                cfg.checkv.min_bin_rescue_bp,
            )
            row["rescued_draft_bins"] = len(extra_bins)
            _write_fasta(draft_seqs, draft_dir / f"{sample_id}_draft.fasta")
            log_ok("  \\[CheckV] genome selection complete")
        else:
            log_warn(
                "  \\[CheckV] no output found — "
                f"check reports/05_quality/{sample_id}_checkv.log"
            )
        progress.advance(task)

        # 3/4 — Dereplication: solo single-contig hq_seqs a 0.98 ANI
        progress.update(task, description="[3/4] cd-hit-est — dereplication 98% ANI")
        if hq_seqs:
            hq_seqs, n_before, n_after = _dereplicate(
                hq_seqs, tmp_dir, sample_id, rpt_dir, cfg.threads, force
            )
            log_ok(
                f"  \\[cd-hit-est] {n_before} → {n_after} genomes  "
                f"(−{n_before - n_after} redundant at {_CDHIT_ID} ANI)"
            )
        # Añadir bins del draft rescue DESPUÉS de la dereplicación:
        # los bins multi-contig ya son únicos por construcción (1 taxón = 1 bin).
        hq_seqs.extend(extra_bins)
        if extra_bins:
            log_ok(f"  \\[draft rescue] {len(extra_bins)} bin(s) added to annotation pool")
        progress.advance(task)

        # 4/4 — Co-binning HQ + taxonomy rename + seqkit format
        progress.update(task,
            description="[4/4] seqkit   — co-binning, renaming & formatting")
        if hq_seqs:
            # Co-bin contigs single del mismo taxón (pass-through para multi ya existentes)
            hq_seqs = _cobin_by_taxon(hq_seqs, tax_map, sample_id)

            renamed, rename_map = _rename_candidates(hq_seqs, tax_map)
            _save_rename_map(rename_map, rpt_dir / f"{sample_id}_rename_map.tsv")

            hq_fastas = _write_candidates(
                renamed, rename_map, phage_dir, prov_dir, tmp_dir,
                rpt_dir / f"{sample_id}_seqkit.log",
            )
            n_multi  = sum(1 for r in rename_map if r.get("n_contigs", 1) > 1)
            n_phage  = sum(1 for r in rename_map if not r["is_provirus"])
            n_prov   = sum(1 for r in rename_map if r["is_provirus"])
            log_ok(
                f"  \\[seqkit] {len(hq_fastas)} file(s) written  "
                f"(phage={n_phage}  provirus={n_prov}  multi-contig={n_multi})"
            )
        progress.advance(task)

    shutil.rmtree(tmp_dir, ignore_errors=True)

    # ── Metrics, display, save ─────────────────────────────────────────────────
    _validate_output(sample_id, hq_fastas, draft_seqs)
    _print_summary_table(sample_id, row, hq_fastas, draft_seqs)
    _check_warnings(row, cfg.checkv.min_completeness)
    if row:
        _save_tsv(row, rpt_dir / "checkv_summary.tsv")

    _print_completion_panel(
        sample_id, phage_dir, prov_dir, draft_dir, rpt_dir,
        row, hq_fastas, draft_seqs,
    )

    log_step(f"Module 05 completed ✓  [{sample_id}]")
    log_info("Next: phageflow annotate --sample-id <id> --genome <candidate.fasta>")

    return hq_fastas


# ── Step implementations ──────────────────────────────────────────────────────

def _run_checkv(cfg, sample_id, virus_fna, sdir, rpt_dir, force):
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


def _dereplicate(seqs_list, tmp_dir, sample_id, rpt_dir, threads, force):
    """Dereplica genomas HQ single-contig con cd-hit-est al 0.98 ANI."""
    n_before = len(seqs_list)
    if n_before <= 1:
        return seqs_list, n_before, n_before

    tmp_in  = tmp_dir / f"{sample_id}_hq_raw.fasta"
    tmp_out = tmp_dir / f"{sample_id}_hq_derep.fasta"

    _write_fasta(seqs_list, tmp_in)

    try:
        run_silent([
            "cd-hit-est",
            "-i", str(tmp_in), "-o", str(tmp_out),
            "-c", _CDHIT_ID, "-aS", _CDHIT_AS,
            "-G", "0", "-n", _CDHIT_N, "-d", "0",
            "-sc", "1", "-T", str(threads), "-M", "4000",
        ], log_file=rpt_dir / f"{sample_id}_cdhit_derep.log")
    except Exception as e:
        log_warn(f"  \\[cd-hit-est] dereplication warning: {e} — usando todos los genomas")
        return seqs_list, n_before, n_before
    finally:
        Path(str(tmp_out) + ".clstr").unlink(missing_ok=True)

    rep_ids = set(_load_fasta(tmp_out).keys())
    derep   = [item for item in seqs_list if item[0] in rep_ids]
    return derep, n_before, len(derep)


# ── CheckV parser (v3.2) ──────────────────────────────────────────────────────

def _parse_checkv(
    sample_id:          str,
    qfile:              Path,
    seqs:               dict,
    prov_seqs:          dict,
    min_completeness:   float,
    min_viral_genes:    int,
    length_rescue:      int,
    min_contig_bp:      int,
    large_nd_rescue_bp: int,
    min_gene_density:   float,
) -> tuple:
    """
    Parsea quality_summary.tsv y separa secuencias en HQ / draft.

    Formato de tupla extendido (v3.2):
        single-contig → (cid, seq, clen, comp, qual, is_provirus, vgn_int)  [7]

    Lógica de tiers (actualizada):
    ─────────────────────────────────────────────────────────────────────────
    Complete / High-quality               → annotation_ready/  (siempre)
    Medium-quality ≥ min_completeness     → annotation_ready/
    Medium-quality < min_completeness     → drafts/
    Low-quality + (≥ min_viral_genes OR ≥ min_gene_density genes/kb)
                                          → drafts/  (rescate LQ)
    Not-determined + ≥ large_nd_rescue_bp + ≥ 1 gen viral
                                          → annotation_ready/  (rescate ND grande)
    Not-determined + ≥ length_rescue + ≥ 1 gen viral
                                          → drafts/  (rescate ND moderado)
    Not-determined + gene_density ≥ min_gene_density
                                          → drafts/  (rescate ND por densidad)
    < min_contig_bp                       → descartado
    ─────────────────────────────────────────────────────────────────────────
    """
    complete = hq_count = mq_count = lq_count = nd_count = 0
    rescued_lq = rescued_nd = rescued_large_nd = rescued_nd_density = 0
    total_viral_genes = total_host_genes = 0
    best_len = 0
    best_completeness = 0.0
    best_quality_tier = ""
    max_contamination = 0.0
    hq_seqs:    list = []
    draft_seqs: list = []

    try:
        with open(qfile) as f:
            header = f.readline().strip().split("\t")
            col    = {h: i for i, h in enumerate(header)}

            for line in f:
                row = line.strip().split("\t")
                if not row or len(row) < 3:
                    continue

                cid      = _col(row, col, "contig_id",     0)
                clen     = _col(row, col, "contig_length",  1)
                qual     = _col(row, col, "checkv_quality", 7)
                comp     = _col(row, col, "completeness",   9)
                contm    = _col(row, col, "contamination", 10)
                vgn      = _col(row, col, "viral_genes",    4)
                hgn      = _col(row, col, "host_genes",     5)
                prov_col = _col(row, col, "provirus",       2)

                is_provirus = (
                    prov_col.strip().lower() == "yes"
                    or "|provirus" in cid.lower()
                )

                clen_int = int(clen or 0)
                vgn_int  = int(vgn  or 0)

                # ── Filtro de longitud mínima ─────────────────────────────────
                if clen_int < min_contig_bp:
                    log_info(
                        f"  \\[filter] {cid} descartado — "
                        f"{clen_int}bp < min_contig_bp {min_contig_bp}bp"
                    )
                    continue

                # ── Densidad génica viral ──────────────────────────────────────
                gene_density = vgn_int / (clen_int / 1000) if clen_int > 0 else 0.0

                # ── Contadores globales ───────────────────────────────────────
                if qual == "Complete":          complete  += 1
                elif qual == "High-quality":    hq_count  += 1
                elif qual == "Medium-quality":  mq_count  += 1
                elif qual == "Low-quality":     lq_count  += 1
                else:                           nd_count  += 1

                try: total_viral_genes += vgn_int
                except (ValueError, TypeError): pass
                try: total_host_genes  += int(hgn)
                except (ValueError, TypeError): pass
                try:
                    c = float(contm)
                    if c > max_contamination:
                        max_contamination = c
                except (ValueError, TypeError): pass

                try:
                    c = float(comp) if comp not in ("", "NA", "N/A") else 0.0
                    if clen_int > best_len:
                        best_len          = clen_int
                        best_completeness = c
                        best_quality_tier = qual
                except (ValueError, TypeError):
                    pass

                # Obtener secuencia (trimmed para provirus)
                seq = seqs.get(cid, "")
                if is_provirus and cid in prov_seqs:
                    seq = prov_seqs[cid]
                if not seq:
                    continue

                try:
                    comp_val = float(comp) if comp not in ("", "NA", "N/A") else 0.0
                except (ValueError, TypeError):
                    comp_val = 0.0

                # ── Lógica de tiers (v3.2) ────────────────────────────────────
                tup = (cid, seq, clen, comp, qual, is_provirus, vgn_int)  # 7-elem

                if qual in ("Complete", "High-quality"):
                    hq_seqs.append(tup)

                elif qual == "Medium-quality":
                    if comp_val >= min_completeness:
                        hq_seqs.append(tup)
                    else:
                        draft_seqs.append(tup)

                elif qual == "Low-quality":
                    # Rescate LQ: count absoluto O densidad génica
                    # (Nayfach 2021; Roux 2019 eLife)
                    if vgn_int >= min_viral_genes or gene_density >= min_gene_density:
                        log_warn(
                            f"  \\[LQ rescue] {cid} — "
                            f"{vgn_int} viral gene(s), "
                            f"density={gene_density:.2f}/kb → drafts/"
                        )
                        draft_seqs.append(tup)
                        rescued_lq += 1

                elif qual == "Not-determined":
                    # Rescate ND grande → annotation_ready/
                    # CheckV subestima completeness para Herelleviridae/Ackermannviridae
                    # sin referente en DB. Contigs ND ≥ large_nd_rescue_bp en
                    # preparaciones purificadas son prácticamente siempre fago.
                    # (Nayfach 2021 Supp; Camargo 2023; Adriaenssens 2020 Viruses)
                    if clen_int >= large_nd_rescue_bp and vgn_int >= 1:
                        log_warn(
                            f"  \\[large-ND rescue] {cid} — "
                            f"{clen_int:,}bp ≥ {large_nd_rescue_bp:,}bp, "
                            f"{vgn_int} viral gene(s), "
                            f"density={gene_density:.2f}/kb → annotation_ready/"
                        )
                        hq_seqs.append(tup)
                        rescued_large_nd += 1

                    # Rescate ND moderado → drafts/
                    elif clen_int >= length_rescue and vgn_int >= 1:
                        log_warn(
                            f"  \\[ND rescue] {cid} — "
                            f"{clen_int:,}bp ≥ {length_rescue:,}bp, "
                            f"{vgn_int} viral gene(s) → drafts/"
                        )
                        draft_seqs.append(tup)
                        rescued_nd += 1

                    # Rescate ND por densidad génica → drafts/
                    elif gene_density >= min_gene_density and clen_int >= min_contig_bp:
                        log_warn(
                            f"  \\[density-ND rescue] {cid} — "
                            f"density={gene_density:.2f}/kb ≥ {min_gene_density}/kb, "
                            f"{clen_int:,}bp → drafts/"
                        )
                        draft_seqs.append(tup)
                        rescued_nd_density += 1

    except Exception as e:
        log_warn(f"  \\[CheckV] could not parse {qfile}: {e}")

    rescued_total = rescued_lq + rescued_nd + rescued_nd_density + rescued_large_nd
    if rescued_total:
        log_ok(
            f"  \\[rescue] LQ={rescued_lq}  ND-moderate={rescued_nd}  "
            f"ND-density={rescued_nd_density}  ND-large→HQ={rescued_large_nd}"
        )

    summary = {
        "sample":             sample_id,
        "complete":           complete,
        "hq":                 hq_count,
        "mq":                 mq_count,
        "lq":                 lq_count,
        "nd":                 nd_count,
        "rescued_lq":         rescued_lq,
        "rescued_nd":         rescued_nd,
        "rescued_nd_density": rescued_nd_density,
        "rescued_large_nd":   rescued_large_nd,
        "rescued_draft_bins": 0,   # actualizado en run() tras _cobin_draft_rescue
        "hq_promoted":        len(hq_seqs),
        "draft_saved":        len(draft_seqs),
        "best_length_bp":     best_len,
        "best_completeness":  best_completeness,
        "best_quality_tier":  best_quality_tier,
        "total_viral_genes":  total_viral_genes,
        "total_host_genes":   total_host_genes,
        "max_contamination":  round(max_contamination, 2),
    }
    return summary, hq_seqs, draft_seqs


# ── Draft co-bin rescue (v3.2, NUEVA FUNCIÓN) ────────────────────────────────

def _cobin_draft_rescue(
    draft_seqs:        list,
    tax_map:           dict,
    sample_id:         str,
    min_bin_rescue_bp: int,
) -> tuple[list, list]:
    """
    Agrupa draft sequences por taxón geNomad y promueve bins a annotation_ready/.

    En una preparación purificada, múltiples contigs MQ/LQ/ND del mismo taxón
    casi siempre representan un único genoma fragmentado por el ensamble. Si la
    suma de sus longitudes supera min_bin_rescue_bp, el bin se promueve como
    entrada multi-contig para Pharokka --meta (modo correcto para multi-FASTA).

    Singletons (taxón único o bin demasiado corto) quedan en drafts/.

    Referencias
    -----------
    - Buttimer et al. 2020, Front Microbiol (Herelleviridae fragmented assemblies)
    - Adriaenssens et al. 2020, Viruses 12:955 (Ackermannviridae co-binning)
    - Camargo et al. 2023, Nat Biotechnol (geNomad taxonomy for binning)
    """
    taxon_buckets: dict[str, list] = defaultdict(list)
    for item in draft_seqs:
        cid   = item[0]
        taxon = tax_map.get(cid, "Unknown")
        taxon_buckets[taxon].append(item)

    promoted:  list = []
    remaining: list = []

    for taxon, items in taxon_buckets.items():
        total_len = sum(int(i[2] or 0) for i in items)
        total_vgn = sum(i[6] if len(i) > 6 else 0 for i in items)

        if len(items) >= 2 and total_len >= min_bin_rescue_bp:
            max_comp     = max(
                float(i[3]) if i[3] not in ("", "NA", "N/A") else 0.0
                for i in items
            )
            best_qual    = _best_quality_tier_from_list([i[4] for i in items])
            any_provirus = any(i[5] for i in items)
            joined_cids  = "||".join(i[0] for i in items)

            log_warn(
                f"  \\[draft-bin rescue] {taxon}: {len(items)} contigs, "
                f"{total_len:,}bp total, {total_vgn} viral genes "
                f"→ annotation_ready/ (multi-contig)"
            )
            # 8-element multi-contig tuple (consistent with _cobin_by_taxon)
            promoted.append((
                f"__MULTICONTIG__{joined_cids}",
                _MULTI_SENTINEL,
                str(total_len),
                str(round(max_comp, 1)),
                best_qual,
                any_provirus,
                str(total_vgn),
                items,
            ))
        else:
            remaining.extend(items)

    return promoted, remaining


# ── Co-binning HQ por taxón (v3.2: pass-through para multi existentes) ────────

def _cobin_by_taxon(
    seqs_list: list,
    tax_map:   dict,
    sample_id: str,
) -> list:
    """
    Agrupa contigs single-contig del mismo taxón en bins multi-contig.

    Bins multi-contig ya existentes (del draft rescue o de llamadas previas)
    se pasan sin modificar para evitar re-binning incorrecto.

    Referencias
    -----------
    - Adriaenssens et al. 2020, Viruses 12:955
    - Buttimer et al. 2020, Front Microbiol
    """
    # Separar ya-multi de single-contig
    already_multi = [item for item in seqs_list if item[1] == _MULTI_SENTINEL]
    single        = [item for item in seqs_list if item[1] != _MULTI_SENTINEL]

    taxon_buckets: dict[str, list] = defaultdict(list)
    for item in single:
        cid   = item[0]
        taxon = tax_map.get(cid, "Unknown")
        taxon_buckets[taxon].append(item)

    result = list(already_multi)

    for taxon, items in taxon_buckets.items():
        if len(items) == 1:
            result.append(items[0])
        else:
            total_len    = sum(int(i[2] or 0) for i in items)
            total_vgn    = sum(i[6] if len(i) > 6 else 0 for i in items)
            max_comp     = max(
                float(i[3]) if i[3] not in ("", "NA", "N/A") else 0.0
                for i in items
            )
            best_qual    = _best_quality_tier_from_list([i[4] for i in items])
            any_provirus = any(i[5] for i in items)
            joined_cids  = "||".join(i[0] for i in items)

            log_info(
                f"  \\[co-bin] {taxon}: {len(items)} contigs "
                f"({total_len:,}bp, {total_vgn} viral genes) → multi-contig bin"
            )
            result.append((
                f"__MULTICONTIG__{joined_cids}",
                _MULTI_SENTINEL,
                str(total_len),
                str(round(max_comp, 1)),
                best_qual,
                any_provirus,
                str(total_vgn),  # 7th: total viral genes
                items,           # 8th: orig_items
            ))

    return result


# ── Rename + write (v3.2: maneja tuplas de 7 y 8 elementos) ──────────────────

def _rename_candidates(seqs_list, tax_map):
    """
    Renombra a {Taxon}_candidate_{NNN} o {Taxon}_multicontig_{NNN}.

    Detección de multi-contig: item[1] == "__MULTICONTIG__" (8-elem).
    Single-contig: item[1] es una secuencia de ADN (7-elem).
    """
    counters:     dict[str, int] = defaultdict(int)
    renamed_seqs: list           = []
    rename_map:   list           = []

    for item in seqs_list:
        is_multi = item[1] == _MULTI_SENTINEL

        if is_multi:
            _, _, clen, comp, qual, is_provirus, total_vgn, orig_items = item
            taxon = tax_map.get(orig_items[0][0], "Unknown")
            counters[taxon] += 1
            new_id    = f"{taxon}_multicontig_{counters[taxon]:03d}"
            orig_cids = [i[0] for i in orig_items]
            desc      = f"{new_id} [original_ids={','.join(orig_cids)}]"

            renamed_seqs.append((desc, _MULTI_SENTINEL, clen, comp,
                                  qual, is_provirus, total_vgn, orig_items))
            rename_map.append({
                "new_id":       new_id,
                "original_id":  ",".join(orig_cids),
                "taxon":        taxon,
                "length_bp":    clen,
                "completeness": comp,
                "quality_tier": qual,
                "is_provirus":  is_provirus,
                "n_contigs":    len(orig_items),
            })
            log_info(
                f"  \\[rename] {len(orig_items)} contigs → {new_id} "
                f"({taxon}, multi-contig, {total_vgn} viral genes)"
            )
        else:
            cid, seq, clen, comp, qual, is_provirus, vgn = item
            taxon = tax_map.get(cid, "Unknown")
            counters[taxon] += 1
            new_id = f"{taxon}_candidate_{counters[taxon]:03d}"
            desc   = f"{new_id} [original_id={cid}]"

            renamed_seqs.append((desc, seq, clen, comp, qual, is_provirus, vgn))
            rename_map.append({
                "new_id":       new_id,
                "original_id":  cid,
                "taxon":        taxon,
                "length_bp":    clen,
                "completeness": comp,
                "quality_tier": qual,
                "is_provirus":  is_provirus,
                "n_contigs":    1,
            })
            log_info(f"  \\[rename] {cid} → {new_id}  ({taxon})")

    return renamed_seqs, rename_map


def _write_candidates(renamed_seqs, rename_map, phage_dir, prov_dir, tmp_dir, log_file):
    """Escribe un FASTA por candidato. Multi-contig → multi-FASTA."""
    written_paths = []
    for item, mapping in zip(renamed_seqs, rename_map):
        is_multi   = item[1] == _MULTI_SENTINEL
        target_dir = prov_dir if mapping["is_provirus"] else phage_dir
        out_fasta  = target_dir / f"{mapping['new_id']}.fasta"
        tmp_fasta  = tmp_dir    / f"{mapping['new_id']}_raw.tmp"

        if is_multi:
            desc, _, _, _, _, is_provirus, total_vgn, orig_items = item
            seqs_to_write = [(i[0], i[1]) for i in orig_items if i[1]]
            _write_fasta(seqs_to_write, tmp_fasta)
        else:
            desc, seq, *_ = item
            _write_fasta([(desc, seq)], tmp_fasta)

        _format_seqkit(tmp_fasta, out_fasta, log_file)
        tmp_fasta.unlink(missing_ok=True)
        written_paths.append(out_fasta)

    return written_paths


# ── Nuevas funciones auxiliares ───────────────────────────────────────────────

def _find_genomad_tsv(cfg: Config, sample_id: str) -> Path:
    prefix = f"{sample_id}_contigs_nr"
    return (
        cfg.results("04_viral_id")
        / sample_id
        / f"{prefix}_summary"
        / f"{prefix}_virus_summary.tsv"
    )


def _load_genomad_taxonomy(tsv: Path) -> dict[str, str]:
    """Build {contig_id: best_taxon_label} from geNomad virus_summary.tsv."""
    tax_map: dict[str, str] = {}
    if not tsv.exists():
        log_warn(f"  \\[rename] geNomad taxonomy file not found: {tsv}")
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
                if tax_idx >= 0 and tax_idx < len(parts):
                    raw = parts[tax_idx].strip()
                    if raw and raw.upper() not in ("NA", "N/A", ""):
                        levels = [
                            lvl.strip() for lvl in raw.split(";")
                            if lvl.strip()
                            and "unclassified" not in lvl.strip().lower()
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
        log_warn(f"  \\[rename] could not parse geNomad taxonomy: {e}")
    return tax_map


def _load_provirus_seqs(prov_fna: Path) -> dict:
    """
    Carga secuencias trimmed de proviruses.fna (CheckV).
    Usa secuencias sin flancos de hospedero. (Nayfach 2021 Supp Methods)
    """
    if not prov_fna.exists():
        return {}
    seqs = _load_fasta(prov_fna)
    if seqs:
        log_info(
            f"  \\[CheckV] {len(seqs)} provirus sequence(s) trimmed "
            f"(proviruses.fna) — flancos de hospedero removidos"
        )
    return seqs


def _best_quality_tier_from_list(tiers: list) -> str:
    for t in _TIER_ORDER:
        if t in tiers:
            return t
    return tiers[0] if tiers else "Not-determined"


# ── Formato seqkit ────────────────────────────────────────────────────────────

def _format_seqkit(in_fasta: Path, out_fasta: Path, log_f: Path) -> None:
    log_f.parent.mkdir(parents=True, exist_ok=True)
    try:
        with open(out_fasta, "w") as fout, open(log_f, "a") as ferr:
            result = subprocess.run(
                ["seqkit", "seq", "-w", "60", str(in_fasta)],
                stdout=fout, stderr=ferr,
            )
        if result.returncode != 0:
            log_warn(f"  \\[seqkit] non-zero exit — copiando sin formato")
            shutil.copy(in_fasta, out_fasta)
    except Exception as e:
        log_warn(f"  \\[seqkit] formatting warning: {e} — copiando sin formato")
        shutil.copy(in_fasta, out_fasta)


def _save_rename_map(rename_map: list, path: Path) -> None:
    headers = ["new_id", "original_id", "taxon", "length_bp",
               "completeness", "quality_tier", "n_contigs"]
    with open(path, "w") as f:
        f.write("\t".join(headers) + "\n")
        for row in rename_map:
            f.write("\t".join(str(row.get(h, "")) for h in headers) + "\n")


# ── Helpers FASTA ─────────────────────────────────────────────────────────────

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


def _col(row: list, col: dict, name: str, fallback: int) -> str:
    idx = col.get(name, fallback)
    return row[idx] if idx < len(row) else ""


def _write_fasta(seqs_list: list, out_path: Path) -> None:
    """Escribe FASTA. Compatible con tuplas de cualquier tamaño (usa pos 0 y 1)."""
    if not seqs_list:
        return
    with open(out_path, "w") as f:
        for item in seqs_list:
            cid, seq = item[0], item[1]
            if seq and seq != _MULTI_SENTINEL:
                f.write(f">{cid}\n{seq}\n")


# ── Validation ────────────────────────────────────────────────────────────────

def _validate_output(sample_id, hq_fastas, draft_seqs):
    if hq_fastas:
        log_ok(f"  validate · {len(hq_fastas)} HQ genome(s) annotation-ready — OK")
    elif draft_seqs:
        log_warn(
            f"  validate · 0 HQ genomes — {len(draft_seqs)} draft(s) en drafts/. "
            "Considera bajar checkv.min_completeness en config.yaml."
        )
    else:
        log_warn(
            f"  validate · 0 genomes pasaron filtros para {sample_id}. "
            "Revisa el output de geNomad (Módulo 04)."
        )


# ── Rich display helpers ──────────────────────────────────────────────────────

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


def _print_input_table(sample_id, virus_fna, stats):
    log_info(f"  Virus FASTA : {virus_fna}  ({human_size(virus_fna)})")
    log_info(
        f"  Input       : {stats['n']} contig(s) | "
        f"total={stats['total_bp']:,}bp | largest={stats['largest_bp']:,}bp | "
        f"N50={stats['n50']:,}bp"
    )


def _print_summary_table(sample_id, row, hq_fastas, draft_seqs):
    if not row:
        return
    log_ok(
        f"  Tiers    : Complete={row.get('complete', 0)}  "
        f"High-quality={row.get('hq', 0)}  "
        f"Medium-quality={row.get('mq', 0)}  "
        f"Low-quality={row.get('lq', 0)}  "
        f"ND={row.get('nd', 0)}"
    )
    rescued = (row.get('rescued_lq', 0) + row.get('rescued_nd', 0) +
               row.get('rescued_nd_density', 0) + row.get('rescued_large_nd', 0) +
               row.get('rescued_draft_bins', 0))
    if rescued:
        log_ok(
            f"  Rescued  : LQ={row.get('rescued_lq', 0)}  "
            f"ND-mod={row.get('rescued_nd', 0)}  "
            f"ND-dens={row.get('rescued_nd_density', 0)}  "
            f"ND-large→HQ={row.get('rescued_large_nd', 0)}  "
            f"draft-bins→HQ={row.get('rescued_draft_bins', 0)}"
        )
    log_ok(
        f"  Selected : {len(hq_fastas)} → annotation_ready/  |  "
        f"{len(draft_seqs)} → drafts/"
    )
    log_ok(
        f"  Best     : {row.get('best_length_bp', 0):,}bp  |  "
        f"completeness={_color_completeness(row.get('best_completeness', 0.0))}  |  "
        f"{_fmt(row.get('best_quality_tier', ''))}"
    )
    log_ok(
        f"  Genes    : viral={_fmt(row.get('total_viral_genes', 0))}  "
        f"host={_fmt(row.get('total_host_genes', 0))}  |  "
        f"max contamination={_color_contamination(row.get('max_contamination', 0.0))}"
    )


def _check_warnings(row, min_completeness):
    if not row:
        return
    if row.get("hq_promoted", 0) == 0 and row.get("draft_saved", 0) > 0:
        log_warn(
            f"  No genomes met the {min_completeness}% threshold — "
            f"{row['draft_saved']} draft(s) saved. "
            "Baja checkv.min_completeness o revisa min_bin_rescue_bp."
        )
    elif row.get("hq_promoted", 0) == 0 and not row.get("rescued_draft_bins", 0):
        log_warn(
            "  No genomes passed quality filters — "
            "verifica el output de geNomad o la contiguidad del ensamble."
        )
    if 0 < row.get("best_completeness", 0.0) < _COMPLETENESS_WARN:
        log_warn(
            f"  Best completeness {row['best_completeness']:.1f}% < {_COMPLETENESS_WARN}% — "
            "genoma fragmentado; anotación y predicción de ciclo de vida pueden verse afectados."
        )
    if row.get("total_host_genes", 0) > _HOST_GENE_WARN:
        log_warn(
            f"  {row['total_host_genes']} host gene(s) — posible contaminación. "
            "Revisa el Módulo 02 (host-removal)."
        )
    if row.get("max_contamination", 0.0) > _CONTAM_WARN:
        log_warn(
            f"  Max contamination {row['max_contamination']:.1f}% > {_CONTAM_WARN}% — "
            "inspecciona quality_summary.tsv."
        )


def _print_completion_panel(
    sample_id, phage_dir, prov_dir, draft_dir, rpt_dir,
    row, hq_fastas, draft_seqs,
):
    complete = row.get("complete", 0)
    hq_count = row.get("hq",      0)
    mq_count = row.get("mq",      0)
    lq_count = row.get("lq",      0)
    nd_count = row.get("nd",      0)

    text = Text()
    text.append("✓ ", style="bold green")
    text.append(
        f"Complete={complete}  HQ={hq_count}  MQ={mq_count}  "
        f"LQ={lq_count}  ND={nd_count}  →  ",
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
    if draft_seqs:
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

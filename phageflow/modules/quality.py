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
_COMPLETENESS_GOOD  = 90.0   # completeness: green
_COMPLETENESS_WARN  = 50.0   # completeness: yellow; below → red
_CONTAM_WARN        = 5.0    # contamination %: warn if above
_HOST_GENE_WARN     = 3      # host genes: warn if more than this

# cd-hit-est dereplication — Nayfach et al. 2021 (CheckV / IMG/VR)
_CDHIT_ID  = "0.95"   # 95% ANI — ICTV bacteriophage species boundary
_CDHIT_AS  = "0.85"   # 85% alignment fraction of shorter genome
_CDHIT_N   = "8"      # word length for 0.90–1.00 identity range


# ── Public entry point ────────────────────────────────────────────────────────

def run(
    cfg:       Config,
    sample_id: str,
    virus_fna: Path,
    force:     bool = False,
) -> list[Path]:
    """
    Assess genome quality with CheckV, dereplicate, rename, and format.

    Pipeline
    --------
    CheckV → tier selection → cd-hit-est 95% ANI → taxonomy rename
    → seqkit -w 60 → one FASTA per candidate in phages/ or proviruses/
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
        f"  Dereplication: cd-hit-est  {_CDHIT_ID} ANI / {_CDHIT_AS} AF  "
        f"(Nayfach et al. 2021 · ICTV species boundary)"
    )

    if not cfg.databases.checkv.exists():
        log_warn(f"  CheckV database not found: {cfg.databases.checkv}")

    # ── Working paths ─────────────────────────────────────────────────────────
    sdir       = out_dir / sample_id
    qfile      = sdir / "quality_summary.tsv"
    row:        dict = {}
    hq_seqs:    list = []   # tuples: (cid, seq, clen, comp, qual, is_provirus)
    draft_seqs: list = []
    hq_fastas:  list = []   # final output Paths

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
        task = progress.add_task("Initializing...", total=4)

        # 1/4 — CheckV end-to-end
        progress.update(task, description="[1/4] CheckV   — end-to-end assessment")
        _run_checkv(cfg, sample_id, virus_fna, sdir, rpt_dir, force)
        progress.advance(task)

        # 2/4 — Parse, select, write drafts
        progress.update(task, description="[2/4] CheckV   — selecting HQ genomes")
        if qfile.exists():
            seqs = _load_fasta(virus_fna)
            row, hq_seqs, draft_seqs = _parse_checkv(
                sample_id, qfile, seqs, cfg.checkv.min_completeness
            )
            _write_fasta(draft_seqs, draft_dir / f"{sample_id}_draft.fasta")
            log_ok("  \\[CheckV] genome selection complete")
        else:
            log_warn(
                "  \\[CheckV] no output found — "
                f"check reports/05_quality/{sample_id}_checkv.log"
            )
        progress.advance(task)

        # 3/4 — Dereplication (cd-hit-est 95% ANI)
        progress.update(task, description="[3/4] cd-hit-est — dereplication 95% ANI")
        if hq_seqs:
            hq_seqs, n_before, n_after = _dereplicate(
                hq_seqs, tmp_dir, sample_id, rpt_dir, cfg.threads, force
            )
            log_ok(
                f"  \\[cd-hit-est] {n_before} → {n_after} genomes  "
                f"(−{n_before - n_after} redundant at {_CDHIT_ID} ANI)"
            )
        progress.advance(task)

        # 4/4 — Taxonomy rename + per-file seqkit format
        progress.update(task, description="[4/4] seqkit   — renaming & formatting")
        if hq_seqs:
            genomad_tsv = _find_genomad_tsv(cfg, sample_id)
            tax_map     = _load_genomad_taxonomy(genomad_tsv)
            renamed, rename_map = _rename_candidates(hq_seqs, tax_map)
            _save_rename_map(rename_map, rpt_dir / f"{sample_id}_rename_map.tsv")
            
            hq_fastas = _write_candidates(
                renamed, rename_map, phage_dir, prov_dir, tmp_dir,
                rpt_dir / f"{sample_id}_seqkit.log",
            )
            log_ok(
                f"  \\[seqkit] {len(hq_fastas)} file(s) written  "
                f"({sum(1 for r in rename_map if not r['is_provirus'])} phage  |  "
                f"{sum(1 for r in rename_map if r['is_provirus'])} provirus)"
            )
        progress.advance(task)

    # ── Cleanup ───────────────────────────────────────────────────────────────
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


def _dereplicate(
    seqs_list:  list,
    tmp_dir:    Path,
    sample_id:  str,
    rpt_dir:    Path,
    threads:    int,
    force:      bool,
) -> tuple[list, int, int]:
    """Dereplicate HQ genomes with cd-hit-est at 95% ANI."""
    n_before = len(seqs_list)
    if n_before <= 1:
        return seqs_list, n_before, n_before

    tmp_in  = tmp_dir / f"{sample_id}_hq_raw.fasta"
    tmp_out = tmp_dir / f"{sample_id}_hq_derep.fasta"

    _write_fasta(seqs_list, tmp_in)

    try:
        run_silent([
            "cd-hit-est",
            "-i", str(tmp_in),
            "-o", str(tmp_out),
            "-c", _CDHIT_ID,
            "-aS", _CDHIT_AS,
            "-G",  "0",
            "-n",  _CDHIT_N,
            "-d",  "0",
            "-sc", "1",
            "-T",  str(threads),
            "-M",  "4000",
        ], log_file=rpt_dir / f"{sample_id}_cdhit_derep.log")
    except Exception as e:
        log_warn(f"  \\[cd-hit-est] dereplication warning: {e} — using all genomes")
        return seqs_list, n_before, n_before
    finally:
        Path(str(tmp_out) + ".clstr").unlink(missing_ok=True)

    rep_ids = {h for h, *_ in _load_fasta(tmp_out).items()}
    derep   = [item for item in seqs_list if item[0] in rep_ids]

    return derep, n_before, len(derep)


def _find_genomad_tsv(cfg: Config, sample_id: str) -> Path:
    """Locate the per-contig geNomad virus_summary.tsv for taxonomy lookup."""
    prefix  = f"{sample_id}_contigs_nr"
    genomad = (
        cfg.results("04_viral_id")
        / sample_id
        / f"{prefix}_summary"
        / f"{prefix}_virus_summary.tsv"
    )
    return genomad


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
                parts = line.strip().split("\t")
                cid   = parts[id_idx] if id_idx < len(parts) else ""
                cid_clean = cid.split("|")[0]

                taxon = "Unknown"
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
                    taxon.replace(" ", "_")
                         .replace("/", "_")
                         .replace("\\", "_")
                         .replace(":", "")
                         .replace(";", "")
                )
                tax_map[cid_clean] = taxon
                if cid != cid_clean:
                    tax_map[cid] = taxon

    except Exception as e:
        log_warn(f"  \\[rename] could not parse geNomad taxonomy: {e}")

    return tax_map


def _rename_candidates(
    seqs_list: list,
    tax_map:   dict[str, str],
) -> tuple[list, list[dict]]:
    """Rename sequences to {Taxon}_candidate_{NNN} format."""
    counters: dict[str, int] = defaultdict(int)
    renamed_seqs: list       = []
    rename_map:   list       = []

    for cid, seq, clen, comp, qual, is_provirus in seqs_list:
        taxon = tax_map.get(cid, "Unknown")
        counters[taxon] += 1
        new_id = f"{taxon}_candidate_{counters[taxon]:03d}"
        desc   = f"{new_id} [original_id={cid}]"

        renamed_seqs.append((desc, seq, clen, comp, qual, is_provirus))
        rename_map.append({
            "new_id":       new_id,
            "original_id":  cid,
            "taxon":        taxon,
            "length_bp":    clen,
            "completeness": comp,
            "quality_tier": qual,
            "is_provirus":  is_provirus,
        })
        log_info(f"  \\[rename] {cid}  →  {new_id}  ({taxon})")

    return renamed_seqs, rename_map


def _write_candidates(
    renamed_seqs: list,
    rename_map: list[dict],
    phage_dir: Path,
    prov_dir: Path,
    tmp_dir: Path,
    log_file: Path,
) -> list[Path]:
    """Writes individual candidate sequences using seqkit logic."""
    written_paths = []
    for (desc, seq, _, _, _, is_provirus), mapping in zip(renamed_seqs, rename_map):
        target_dir = prov_dir if is_provirus else phage_dir
        out_fasta = target_dir / f"{mapping['new_id']}.fasta"
        
        # Intermediate raw file before processing with seqkit
        tmp_fasta = tmp_dir / f"{mapping['new_id']}_raw.tmp"
        _write_fasta([(desc, seq)], tmp_fasta)
        
        _format_seqkit(tmp_fasta, out_fasta, log_file)
        tmp_fasta.unlink(missing_ok=True)
        written_paths.append(out_fasta)
        
    return written_paths


def _format_seqkit(in_fasta: Path, out_fasta: Path, log_f: Path) -> None:
    """Wrap FASTA sequences at 60 bp per line with seqkit seq -w 60."""
    log_f.parent.mkdir(parents=True, exist_ok=True)
    try:
        with open(out_fasta, "w") as fout, open(log_f, "a") as ferr:
            result = subprocess.run(
                ["seqkit", "seq", "-w", "60", str(in_fasta)],
                stdout=fout, stderr=ferr,
            )
        if result.returncode != 0:
            log_warn(f"  \\[seqkit] non-zero exit ({result.returncode}) — copying unformatted")
            shutil.copy(in_fasta, out_fasta)
    except Exception as e:
        log_warn(f"  \\[seqkit] formatting warning: {e} — copying unformatted")
        shutil.copy(in_fasta, out_fasta)


# ── Validation ────────────────────────────────────────────────────────────────

def _validate_output(
    sample_id:  str,
    hq_fastas:  list,
    draft_seqs: list,
) -> None:
    """Confirm HQ genomes found, or warn with actionable suggestion."""
    if hq_fastas:
        log_ok(f"  validate · {len(hq_fastas)} HQ genome(s) annotation-ready — OK")
    elif draft_seqs:
        log_warn(
            f"  validate · 0 HQ genomes — {len(draft_seqs)} draft(s) saved to drafts/. "
            "Lower checkv.min_completeness in config.yaml to include them."
        )
    else:
        log_warn(
            f"  validate · 0 genomes passed quality filters for {sample_id}. "
            "Check geNomad viral classification (Module 04)."
        )


# ── Helpers ───────────────────────────────────────────────────────────────────

def _load_fasta(path: Path) -> dict[str, str]:
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
    """Parse CheckV quality_summary.tsv and split sequences into HQ / draft tiers."""
    complete = hq_count = mq_count = lq_count = 0
    total_viral_genes = total_host_genes = 0
    best_len          = 0
    best_completeness = 0.0
    best_quality_tier = ""
    max_contamination = 0.0
    hq_seqs:    list  = []
    draft_seqs: list  = []

    try:
        with open(qfile) as f:
            header = f.readline().strip().split("\t")
            col    = {h: i for i, h in enumerate(header)}

            for line in f:
                row = line.strip().split("\t")
                if not row or len(row) < 3:
                    continue

                cid        = _col(row, col, "contig_id",      0)
                clen       = _col(row, col, "contig_length",   1)
                qual       = _col(row, col, "checkv_quality",   7)
                comp       = _col(row, col, "completeness",    9)
                contm      = _col(row, col, "contamination",  10)
                vgn        = _col(row, col, "viral_genes",     4)
                hgn        = _col(row, col, "host_genes",      5)
                prov_col   = _col(row, col, "provirus",        2)

                is_provirus = (
                    prov_col.strip().lower() == "yes"
                    or "|provirus" in cid.lower()
                )

                seq = seqs.get(cid, "")

                if qual == "Complete":         complete  += 1
                elif qual == "High-quality":   hq_count  += 1
                elif qual == "Medium-quality": mq_count  += 1
                else:                          lq_count  += 1

                try: total_viral_genes += int(vgn)
                except (ValueError, TypeError): pass
                try: total_host_genes  += int(hgn)
                except (ValueError, TypeError): pass
                try:
                    c = float(contm)
                    if c > max_contamination:
                        max_contamination = c
                except (ValueError, TypeError): pass

                try:
                    l = int(clen)
                    c = float(comp) if comp not in ("", "NA", "N/A") else 0.0
                    if l > best_len:
                        best_len          = l
                        best_completeness = c
                        best_quality_tier = qual
                except (ValueError, TypeError): pass

                if not seq:
                    continue

                try:
                    comp_val = float(comp) if comp not in ("", "NA", "N/A") else 0.0
                except (ValueError, TypeError):
                    comp_val = 0.0

                if qual in ("Complete", "High-quality"):
                    hq_seqs.append((cid, seq, clen, comp, qual, is_provirus))
                elif qual == "Medium-quality":
                    if comp_val >= min_completeness:
                        hq_seqs.append((cid, seq, clen, comp, qual, is_provirus))
                    else:
                        draft_seqs.append((cid, seq, clen, comp, qual, is_provirus))

    except Exception as e:
        log_warn(f"  \\[CheckV] could not parse {qfile}: {e}")

    summary = {
        "sample":            sample_id,
        "complete":          complete,
        "hq":                hq_count,
        "mq":                mq_count,
        "lq":                lq_count,
        "hq_promoted":       len(hq_seqs),
        "draft_saved":       len(draft_seqs),
        "best_length_bp":    best_len,
        "best_completeness": best_completeness,
        "best_quality_tier": best_quality_tier,
        "total_viral_genes": total_viral_genes,
        "total_host_genes":  total_host_genes,
        "max_contamination": round(max_contamination, 2),
    }
    return summary, hq_seqs, draft_seqs


def _col(row: list, col: dict, name: str, fallback: int) -> str:
    idx = col.get(name, fallback)
    return row[idx] if idx < len(row) else ""


def _write_fasta(seqs_list: list, out_path: Path) -> None:
    if not seqs_list:
        return
    with open(out_path, "w") as f:
        for cid, seq, *_ in seqs_list:
            f.write(f">{cid}\n{seq}\n")


def _save_rename_map(rename_map: list, path: Path) -> None:
    headers = ["new_id", "original_id", "taxon", "length_bp",
               "completeness", "quality_tier"]
    with open(path, "w") as f:
        f.write("\t".join(headers) + "\n")
        for row in rename_map:
            f.write("\t".join(str(row.get(h, "")) for h in headers) + "\n")


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


def _print_input_table(sample_id: str, virus_fna: Path, stats: dict) -> None:
    log_info(f"  Virus FASTA : {virus_fna}  ({human_size(virus_fna)})")
    log_info(
        f"  Input       : {stats['n']} contig(s) | "
        f"total={stats['total_bp']:,}bp | largest={stats['largest_bp']:,}bp | "
        f"N50={stats['n50']:,}bp"
    )


def _print_summary_table(
    sample_id: str, row: dict, hq_fastas: list, draft_seqs: list
) -> None:
    if not row:
        return
    log_ok(
        f"  Tiers    : Complete={row.get('complete', 0)}  "
        f"High-quality={row.get('hq', 0)}  "
        f"Medium-quality={row.get('mq', 0)}  "
        f"Low/ND={row.get('lq', 0)}"
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


def _check_warnings(row: dict, min_completeness: float) -> None:
    if not row:
        return

    if row.get("hq_promoted", 0) == 0 and row.get("draft_saved", 0) > 0:
        log_warn(
            f"  No genomes met the {min_completeness}% threshold — "
            f"{row['draft_saved']} draft(s) saved. "
            "Lower checkv.min_completeness in config.yaml to annotate them."
        )
    elif row.get("hq_promoted", 0) == 0:
        log_warn(
            "  No genomes passed quality filters — "
            "check geNomad output or assembly contiguity."
        )

    if 0 < row.get("best_completeness", 0.0) < _COMPLETENESS_WARN:
        log_warn(
            f"  Best completeness {row['best_completeness']:.1f}% < {_COMPLETENESS_WARN}% — "
            "fragmented genome; annotation and lifecycle prediction may be affected."
        )

    if row.get("total_host_genes", 0) > _HOST_GENE_WARN:
        log_warn(
            f"  {row['total_host_genes']} host gene(s) detected — possible contamination. "
            "Review Module 02 (host-removal)."
        )

    if row.get("max_contamination", 0.0) > _CONTAM_WARN:
        log_warn(
            f"  Max contamination {row['max_contamination']:.1f}% > {_CONTAM_WARN}% — "
            "genome may contain host sequence. Inspect CheckV quality_summary.tsv."
        )


def _print_completion_panel(
    sample_id:  str,
    phage_dir:  Path,
    prov_dir:   Path,
    draft_dir:  Path,
    rpt_dir:    Path,
    row:        dict,
    hq_fastas:  list,
    draft_seqs: list,
) -> None:
    complete = row.get("complete", 0)
    hq_count = row.get("hq",      0)
    mq_count = row.get("mq",      0)
    lq_count = row.get("lq",      0)

    text = Text()
    text.append("✓ ", style="bold green")
    text.append(
        f"Complete={complete}  HQ={hq_count}  MQ={mq_count}  LQ/ND={lq_count}  →  ",
        style="dim white",
    )
    text.append(f"{len(hq_fastas)}", style="bold green")
    text.append(" genome(s) renamed & annotation-ready\n\n", style="cyan")

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
    "complete", "hq", "mq", "lq",
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

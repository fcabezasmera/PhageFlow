"""PhageFlow Module 05 — Genome quality assessment and selection.

Tool:
    CheckV (Nayfach et al. 2021, Nature Biotechnology):
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

from phageflow.utils.config import Config
from phageflow.utils.logger import log_step, log_info, log_ok, log_warn, log_error
from phageflow.utils.tools import require_tools, run_silent, mkdirs

STEP  = "05_quality"
TOOLS = ["checkv"]


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

    db = cfg.databases.checkv
    if not db.exists():
        log_warn(f"  CheckV database not found: {db}")

    log_step(f"Module 05 — CheckV quality assessment [{sample_id}] · {cfg.threads} threads")

    if not virus_fna.exists():
        log_error(f"  Virus FASTA not found: {virus_fna}")
        raise FileNotFoundError(virus_fna)

    n_input = sum(1 for l in open(virus_fna) if l.startswith(">"))
    log_info(f"  Input   : {n_input} viral contig(s)")
    log_info(f"  HQ threshold : completeness ≥ {cfg.checkv.min_completeness}%")
    log_info("  Nayfach et al. 2021, Nature Biotechnology")

    # ── CheckV end-to-end ─────────────────────────────────────────────────────
    sdir  = out_dir / sample_id
    qfile = sdir / "quality_summary.tsv"

    if qfile.exists() and qfile.stat().st_size > 0 and not force:
        log_info("  Already processed — skipping (use --force to re-run)")
    else:
        mkdirs(sdir)
        log_info("  [CheckV] running end-to-end assessment...")
        try:
            run_silent([
                "checkv", "end_to_end",
                str(virus_fna), str(sdir),
                "-d", str(db),
                "-t", str(cfg.threads),
                "--remove_tmp",
            ], log_file=rpt_dir / f"{sample_id}_checkv.log")
            log_ok("  [CheckV] OK")
        except Exception as e:
            log_warn(f"  [CheckV] warning: {e}")

    # ── Parse and select genomes ──────────────────────────────────────────────
    if not qfile.exists():
        log_warn(f"  No CheckV output found — check log: "
                 f"reports/05_quality/{sample_id}_checkv.log")
        return final_dir / f"{sample_id}_HQ.fasta"

    seqs        = _load_fasta(virus_fna)
    row, hq_seqs, draft_seqs = _parse_checkv(sample_id, qfile, seqs,
                                              cfg.checkv.min_completeness)

    hq_fasta    = final_dir / f"{sample_id}_HQ.fasta"
    draft_fasta = draft_dir / f"{sample_id}_draft.fasta"

    _write_fasta(hq_seqs,    hq_fasta)
    _write_fasta(draft_seqs, draft_fasta)

    log_ok(
        f"  [{sample_id}] "
        f"Complete={row['complete']} HQ={row['hq']} "
        f"MQ={row['mq']} LQ={row['lq']} | "
        f"best={row['best_length']} {row['best_quality']}"
    )

    if hq_seqs:
        log_ok(f"  Annotation-ready : {hq_fasta}  ({len(hq_seqs)} genome(s))")
    else:
        log_warn(f"  No HQ genomes for {sample_id} — only drafts or nothing")
        if draft_seqs:
            log_warn(f"  Draft genomes    : {draft_fasta}  ({len(draft_seqs)} genome(s))")

    _save_tsv(row, rpt_dir / "checkv_summary.tsv")

    log_step(f"Module 05 completed ✓  [{sample_id}]")
    log_info("Next: phageflow annotate --sample-id <id> --genome <HQ.fasta>")

    return hq_fasta


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_fasta(path: Path) -> dict:
    """Load FASTA into {seq_id: sequence} dict."""
    seqs = {}
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
    sample_id: str,
    qfile:     Path,
    seqs:      dict,
    min_completeness: float,
) -> tuple:
    complete = hq = mq = lq = 0
    best_len = 0
    best_quality = ""
    hq_seqs    = []
    draft_seqs = []

    with open(qfile) as f:
        header = f.readline().strip().split("\t")
        col    = {h: i for i, h in enumerate(header)}
        for line in f:
            row = line.strip().split("\t")
            if not row or len(row) < 3:
                continue
            cid  = row[col.get("contig_id",       0)]
            clen = row[col.get("contig_length",    1)]
            qual = row[col.get("checkv_quality",   7)]
            comp = row[col.get("completeness",     9)]
            seq  = seqs.get(cid, "")

            if qual == "Complete":           complete += 1
            elif qual == "High-quality":     hq       += 1
            elif qual == "Medium-quality":   mq       += 1
            else:                            lq       += 1

            try:
                if int(clen) > best_len:
                    best_len     = int(clen)
                    best_quality = f"{comp}% {qual}"
            except ValueError:
                pass

            if seq:
                if qual in ("Complete", "High-quality"):
                    hq_seqs.append((cid, seq, clen, comp, qual))
                elif qual == "Medium-quality":
                    draft_seqs.append((cid, seq, clen, comp, qual))

    summary = {
        "sample":       sample_id,
        "complete":     complete,
        "hq":           hq,
        "mq":           mq,
        "lq":           lq,
        "best_length":  f"{best_len}bp",
        "best_quality": best_quality,
    }
    return summary, hq_seqs, draft_seqs


def _write_fasta(seqs_list: list, out_path: Path) -> None:
    if not seqs_list:
        return
    with open(out_path, "w") as f:
        for cid, seq, *_ in seqs_list:
            f.write(f">{cid}\n{seq}\n")


def _save_tsv(row: dict, path: Path) -> None:
    headers = ["sample", "complete", "hq", "mq", "lq",
               "best_length", "best_quality"]
    rows = {}
    if path.exists():
        with open(path) as f:
            f.readline()
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if cols:
                    rows[cols[0]] = cols
    rows[row["sample"]] = [str(row.get(h, "")) for h in headers]
    with open(path, "w") as f:
        f.write("\t".join(headers) + "\n")
        for r in rows.values():
            f.write("\t".join(r) + "\n")

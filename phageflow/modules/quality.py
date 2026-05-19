"""PhageFlow Module 05 — Genome quality assessment and selection.

Tools:
    CheckV : completeness, contamination, topology (DTR/ITR)

Quality tiers:
    Complete / High-quality  → annotation_ready/  (full annotation)
    Medium-quality           → drafts/             (reported in supplement)
    Low-quality              → discarded from main analysis
"""

from __future__ import annotations
from pathlib import Path
from typing import List

from phageflow.utils.config import Config, Sample
from phageflow.utils.logger import log_step, log_info, log_ok, log_warn
from phageflow.utils.tools import require_tools, run_silent, mkdirs

STEP  = "05_quality"
TOOLS = ["checkv"]


def run(cfg: Config, samples: List[Sample], force: bool = False) -> None:
    require_tools(*TOOLS)

    out_dir   = cfg.results(STEP)
    rpt_dir   = cfg.reports(STEP)
    final_dir = out_dir / "annotation_ready"
    draft_dir = out_dir / "drafts"
    mkdirs(out_dir, rpt_dir, final_dir, draft_dir)

    db = cfg.databases.checkv
    if not db.exists():
        log_warn(f"CheckV database not found: {db}")

    log_step(f"Module 05 — CheckV quality assessment ({len(samples)} samples · {cfg.threads} threads)")

    summary_rows  = []
    selected_rows = []

    for sample in samples:
        fa = cfg.results("04_viral_id") / f"{sample.sample_id}_virus.fna"
        if not fa.exists():
            log_warn(f"[{sample.sample_id}] virus.fna not found")
            continue

        n_input = sum(1 for l in open(fa) if l.startswith(">"))
        log_step(f"[{sample.sample_id}]  {sample.cohort}  |  {n_input} viral contigs")

        sdir   = out_dir / sample.sample_id
        qfile  = sdir / "quality_summary.tsv"

        if qfile.exists() and not force:
            log_info("  Already processed — skipping")
        else:
            mkdirs(sdir)
            try:
                run_silent([
                    "checkv", "end_to_end",
                    str(fa), str(sdir),
                    "-d", str(db),
                    "-t", str(cfg.threads),
                    "--remove_tmp",
                ], log_file=rpt_dir / f"{sample.sample_id}_checkv.log")
            except Exception as e:
                log_warn(f"  CheckV warning: {e}")

        if qfile.exists():
            row, sel = _parse_checkv(sample, fa, qfile, final_dir, draft_dir)
            summary_rows.append(row)
            selected_rows.extend(sel)
            log_ok(
                f"  [{sample.sample_id}] "
                f"Complete={row['complete']} HQ={row['hq']} "
                f"MQ={row['mq']} LQ={row['lq']} | "
                f"best={row['best_length']} {row['best_quality']}"
            )
        else:
            log_warn(f"  [{sample.sample_id}] no CheckV output")

    log_step("CheckV summary")
    _print_table(summary_rows)

    log_step("Selected genomes")
    _print_table(selected_rows)

    _save_tsv(summary_rows, rpt_dir / "checkv_summary.tsv")
    _save_tsv(selected_rows, rpt_dir / "genome_selection.tsv")

    n_hq    = len([f for f in final_dir.glob("*_HQ.fasta")])
    n_draft = len([f for f in draft_dir.glob("*_draft.fasta")])
    log_ok(f"Annotation-ready : {final_dir}/  ({n_hq} genomes)")
    log_ok(f"Drafts           : {draft_dir}/  ({n_draft} genomes)")
    log_step("Module 05 completed ✓")
    log_info("Next: phageflow annotate")


def _parse_checkv(
    sample: Sample, fa: Path, qfile: Path,
    final_dir: Path, draft_dir: Path
) -> tuple:
    seqs = {}
    with open(fa) as f:
        hdr = seq = ""
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if hdr: seqs[hdr.lstrip(">").split()[0]] = seq
                hdr, seq = line, ""
            else: seq += line
        if hdr: seqs[hdr.lstrip(">").split()[0]] = seq

    complete = hq = mq = lq = 0
    best_len = best_q = best_comp = "0"
    hq_seqs = []; draft_seqs = []

    with open(qfile) as f:
        header = f.readline().strip().split("\t")
        col = {h: i for i, h in enumerate(header)}
        for line in f:
            row = line.strip().split("\t")
            cid  = row[col.get("contig_id", 0)]
            clen = row[col.get("contig_length", 1)]
            comp = row[col.get("completeness", 9)]
            qual = row[col.get("checkv_quality", 7)]
            seq  = seqs.get(cid, "")

            if qual == "Complete":       complete += 1
            elif qual == "High-quality": hq       += 1
            elif qual == "Medium-quality": mq     += 1
            else:                          lq     += 1

            try:
                if int(clen) > int(best_len):
                    best_len = clen; best_q = qual; best_comp = comp
            except ValueError:
                pass

            if seq:
                if qual in ("Complete", "High-quality"):
                    hq_seqs.append((cid, seq, clen, comp, qual))
                elif qual == "Medium-quality":
                    draft_seqs.append((cid, seq, clen, comp, qual))

    def _write(seqs_list, out_path):
        if not seqs_list: return
        with open(out_path, "w") as o:
            for cid, seq, *_ in seqs_list:
                o.write(f">{cid}\n{seq}\n")

    _write(hq_seqs,    final_dir / f"{sample.sample_id}_HQ.fasta")
    _write(draft_seqs, draft_dir / f"{sample.sample_id}_draft.fasta")

    summary_row = {
        "sample": sample.sample_id, "cohort": sample.cohort,
        "complete": complete, "hq": hq, "mq": mq, "lq": lq,
        "best_length": f"{best_len}bp",
        "best_quality": f"{best_comp}% {best_q}",
    }
    selected = [
        {"sample": sample.sample_id, "tier": "HQ",
         "contig_id": cid, "length_bp": clen, "completeness": comp, "quality": qual}
        for cid, _, clen, comp, qual in hq_seqs
    ] + [
        {"sample": sample.sample_id, "tier": "draft",
         "contig_id": cid, "length_bp": clen, "completeness": comp, "quality": qual}
        for cid, _, clen, comp, qual in draft_seqs
    ]
    return summary_row, selected


def _print_table(rows):
    if not rows: return
    headers = list(rows[0].keys())
    widths  = {h: len(h) for h in headers}
    for row in rows:
        for h in headers:
            widths[h] = max(widths[h], len(str(row.get(h,""))))
    print("  "+"  ".join(h.ljust(widths[h]) for h in headers))
    print("  "+"-"*(sum(widths.values())+2*len(headers)))
    for row in rows:
        print("  "+"  ".join(str(row.get(h,"")).ljust(widths[h]) for h in headers))


def _save_tsv(rows, path):
    if not rows: return
    headers = list(rows[0].keys())
    with open(path,"w") as f:
        f.write("\t".join(headers)+"\n")
        for row in rows:
            f.write("\t".join(str(row.get(h,"")) for h in headers)+"\n")

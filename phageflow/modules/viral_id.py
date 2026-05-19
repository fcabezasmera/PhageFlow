"""PhageFlow Module 04 — Viral identification with geNomad.

Input  : results/03_assembly/combined/{sample}_contigs_nr.fasta
Output : results/04_viral_id/{sample}/  (virus.fna, taxonomy, plasmid, provirus)
"""

from __future__ import annotations
from pathlib import Path
from typing import List

from phageflow.utils.config import Config, Sample
from phageflow.utils.logger import log_step, log_info, log_ok, log_warn
from phageflow.utils.tools import require_tools, run_silent, mkdirs, fasta_stats

STEP  = "04_viral_id"
TOOLS = ["genomad"]


def run(cfg: Config, samples: List[Sample], force: bool = False) -> None:
    require_tools(*TOOLS)

    out_dir = cfg.results(STEP)
    rpt_dir = cfg.reports(STEP)
    mkdirs(out_dir, rpt_dir)

    db = cfg.databases.genomad
    if not db.exists():
        log_warn(f"geNomad database not found: {db}")

    log_step(f"Module 04 — geNomad viral identification ({len(samples)} samples)")

    summary_rows = []
    for sample in samples:
        fa = cfg.results("03_assembly") / "combined" / f"{sample.sample_id}_contigs_nr.fasta"
        if not fa.exists():
            log_warn(f"[{sample.sample_id}] NR contigs not found: {fa}")
            continue

        n_input = fasta_stats(fa)["n"]
        log_step(f"[{sample.sample_id}]  {sample.cohort}  |  {n_input} NR contigs")

        sdir   = out_dir / sample.sample_id
        prefix = fa.stem
        v_fna  = sdir / f"{prefix}_summary" / f"{prefix}_virus.fna"
        v_tsv  = sdir / f"{prefix}_summary" / f"{prefix}_virus_summary.tsv"

        if v_fna.exists() and not force:
            log_info("  Already processed — skipping")
        else:
            mkdirs(sdir)
            try:
                run_silent([
                    "genomad", "end-to-end",
                    "--cleanup", "--splits", "8",
                    "--min-score", str(cfg.genomad.min_score),
                    "--threads", str(cfg.threads),
                    str(fa), str(sdir), str(db),
                ], log_file=rpt_dir / f"{sample.sample_id}_genomad.log")
            except Exception as e:
                log_warn(f"  geNomad warning: {e}")

        if v_tsv.exists():
            row = _parse_genomad_output(sample, fa, sdir, prefix, v_fna, v_tsv, out_dir)
            summary_rows.append(row)
            log_ok(
                f"  [{sample.sample_id}] viral={row['viral_ctg']} ({row['viral_bp']}bp) "
                f"| plasmid={row['plasmid']} | provirus={row['provirus']} "
                f"| family={row['best_family']}"
            )
        else:
            log_warn(f"  [{sample.sample_id}] no geNomad output")

    log_step("geNomad summary")
    _print_table(summary_rows)

    tsv = rpt_dir / "genomad_summary.tsv"
    _save_tsv(summary_rows, tsv)
    log_ok(f"Virus FASTAs : {out_dir}/*_virus.fna  ← CheckV input")
    log_ok(f"Summary TSV  : {tsv}")
    log_step("Module 04 completed ✓")
    log_info("Next: phageflow quality  (activate phageflow env)")


def _parse_genomad_output(
    sample: Sample, fa: Path, sdir: Path, prefix: str,
    v_fna: Path, v_tsv: Path, out_dir: Path
) -> dict:
    n_input  = fasta_stats(fa)["n"]
    n_viral  = sum(1 for l in open(v_fna) if l.startswith(">")) if v_fna.exists() else 0
    viral_bp = sum(len(l.rstrip()) for l in open(v_fna) if not l.startswith(">")) if v_fna.exists() else 0

    p_tsv   = sdir / f"{prefix}_summary" / f"{prefix}_plasmid_summary.tsv"
    pr_tsv  = sdir / f"{prefix}_find_proviruses" / f"{prefix}_provirus.tsv"
    n_plas  = sum(1 for l in open(p_tsv) if not l.startswith("#")) - 1 if p_tsv.exists() else 0
    n_prov  = sum(1 for l in open(pr_tsv) if not l.startswith("#")) - 1 if pr_tsv.exists() else 0
    n_plas  = max(0, n_plas)
    n_prov  = max(0, n_prov)

    best_fam = "unclassified"
    if v_tsv.exists():
        best_score = 0.0
        with open(v_tsv) as f:
            header = f.readline().strip().split("\t")
            col = {h: i for i, h in enumerate(header)}
            for line in f:
                row = line.strip().split("\t")
                try:
                    score = float(row[col.get("virus_score", 6)])
                    if score > best_score:
                        best_score = score
                        tax = row[col.get("taxonomy", 10)]
                        parts = tax.split(";")
                        for p in reversed(parts):
                            if p.strip():
                                best_fam = p.strip()
                                break
                except (ValueError, IndexError):
                    pass

    if v_fna.exists():
        import shutil
        shutil.copy(v_fna, out_dir / f"{sample.sample_id}_virus.fna")

    return {
        "sample":      sample.sample_id,
        "cohort":      sample.cohort,
        "input_ctg":   n_input,
        "viral_ctg":   n_viral,
        "viral_bp":    viral_bp,
        "plasmid":     n_plas,
        "provirus":    n_prov,
        "best_family": best_fam,
    }


def _print_table(rows):
    if not rows: return
    headers = ["sample","cohort","input_ctg","viral_ctg","viral_bp","plasmid","provirus","best_family"]
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

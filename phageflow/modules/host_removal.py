"""PhageFlow Module 02 — Host read removal.

Three modes (priority order):

  1. --host-file  path
     Local FASTA file, folder of FASTAs, or text file with FASTA paths.
     Uses bwa-mem2 — most precise when host genome is available.
     Reference: Li (2013) arXiv:1303.3997

  2. --accessions / --accessions-file
     Download NCBI Assembly genomes, then align with bwa-mem2.

  3. No arguments → Kraken2 filter (default when no host provided)
     Classify reads; keep unclassified + Viruses (taxid 10239).

     Parameters (literature-based):
       --confidence 0.2
           Reduces bacterial false-positives while keeping novel phage
           as unclassified. CS=0.2-0.4 recommended for metagenomics
           (aBIOTECH 2024; Loomis et al. 2021; Li et al. 2022).
       --minimum-hit-groups 2
           Requires 2+ distinct k-mer groups; reduces spurious hits
           from low-complexity reads (Wood et al. 2019, Genome Biology).

     Read classes retained:
       - Unclassified (taxid 0): novel phage absent from the database.
         For purified phage preps, new phages appear almost entirely here.
       - Viruses (taxid 10239) + all descendants: known phage sequences.

     Discarded:
       - Bacteria (taxid 2), Archaea (taxid 2157), Eukaryota (taxid 2759).

     Extraction uses seqtk subseq — efficient O(n) read ID filtering
     (Gregory et al. 2020, Cell; Krishnamurthy & Wang 2017, Virology).

     Requires: databases.kraken2 set in config.yaml.
"""

from __future__ import annotations
import subprocess
import zipfile
from pathlib import Path
from typing import List, Optional

from phageflow.utils.config import Config
from phageflow.utils.logger import log_step, log_info, log_ok, log_warn, log_error
from phageflow.utils.tools import require_tools, run_silent, human_size, mkdirs

STEP       = "02_host_removal"
FASTA_EXTS = {".fasta", ".fa", ".fna", ".fas"}
TAXID_VIRUSES = 10239


def run(
    cfg:             Config,
    sample_id:       str,
    r1:              Path,
    r2:              Path,
    force:           bool                = False,
    host_file:       Optional[Path]      = None,
    accessions:      Optional[List[str]] = None,
    accessions_file: Optional[Path]      = None,
    kraken_db:       Optional[Path]      = None,
) -> tuple:
    """Remove host reads for a single sample. Returns (r1_out, r2_out)."""

    r1 = Path(r1); r2 = Path(r2)
    out_dir  = cfg.results(STEP)
    rpt_dir  = cfg.reports(STEP)
    host_dir = cfg.workdir / "data" / "host_genomes"
    tmp_dir  = out_dir / "tmp"
    mkdirs(out_dir, rpt_dir, host_dir, tmp_dir)

    log_step(f"Module 02 — host removal [{sample_id}] · {cfg.threads} threads")
    log_info(f"  R1 : {human_size(r1)}")

    r1_out = out_dir / f"{sample_id}_R1.fastq.gz"
    r2_out = out_dir / f"{sample_id}_R2.fastq.gz"

    if r1_out.exists() and r1_out.stat().st_size > 0 and not force:
        log_info("  Already processed — skipping (use --force to re-run)")
        return r1_out, r2_out

    if host_file or accessions or accessions_file:
        require_tools("bwa-mem2", "samtools")
        combined = host_dir / "combined_hosts.fasta"
        fastas   = _resolve_fastas(host_file, accessions, accessions_file,
                                   host_dir, rpt_dir)
        _build_index(fastas, combined, rpt_dir, force)
        stats = _bwa_align(sample_id, r1, r2, r1_out, r2_out,
                           combined, tmp_dir, rpt_dir, cfg.threads)
        log_ok(f"  host={stats['pct_host']}  phage={stats['pct_phage']}  "
               f"reads={stats['reads_phage']}")
        try:
            if float(stats['pct_phage'].rstrip('%')) < 70:
                log_warn("  Phage fraction <70% — possible contamination")
        except ValueError:
            pass
        _save_tsv({**stats, "sample": sample_id, "mode": "bwa-mem2"},
                  rpt_dir / "host_removal_summary.tsv")
    else:
        require_tools("kraken2", "seqtk")
        db = kraken_db or cfg.databases.kraken2
        if not db or not Path(db).exists():
            log_error(
                "Kraken2 database not configured. "
                "Set 'databases.kraken2' in config.yaml "
                "or use --host-file / --accessions."
            )
            return r1_out, r2_out
        stats = _kraken2_filter(sample_id, r1, r2, r1_out, r2_out,
                                Path(db), tmp_dir, rpt_dir, cfg.threads)
        _save_tsv({**stats, "sample": sample_id, "mode": "kraken2-filter"},
                  rpt_dir / "host_removal_summary.tsv")

    import shutil
    shutil.rmtree(tmp_dir, ignore_errors=True)

    log_ok(f"  Phage reads : {r1_out}")
    log_step(f"Module 02 completed ✓  [{sample_id}]")
    log_info("Next: phageflow assembly")
    return r1_out, r2_out


# ---------------------------------------------------------------------------
# Kraken2 filter
# ---------------------------------------------------------------------------

def _kraken2_filter(sample_id, r1, r2, r1_out, r2_out,
                    db, tmp_dir, rpt_dir, threads):
    log_info("  Mode: Kraken2 filter")
    log_info("  Keep: unclassified (novel phage) + Viruses (taxid 10239)")
    log_info("  Params: --confidence 0.2 --minimum-hit-groups 2")

    report  = rpt_dir / f"{sample_id}_k2.report"
    k2_out  = tmp_dir / f"{sample_id}_k2.tsv"
    un_r1   = tmp_dir / f"{sample_id}_unclassified_1.fastq"
    un_r2   = tmp_dir / f"{sample_id}_unclassified_2.fastq"

    # Step 1: classify reads
    log_info("  [1/3] Kraken2 classification...")
    run_silent([
        "kraken2",
        "--db",    str(db),
        "--report", str(report),
        "--output", str(k2_out),
        "--unclassified-out", f"{tmp_dir}/{sample_id}_unclassified#.fastq",
        "--paired",
        "--confidence",        "0.2",
        "--minimum-hit-groups", "2",
        "--threads", str(threads),
        str(r1), str(r2),
    ], log_file=rpt_dir / f"{sample_id}_k2.log", check=False)

    # Step 2: get all viral taxids from report (descendants of taxid 10239)
    log_info("  [2/3] Identifying viral taxids from report...")
    viral_taxids = _parse_viral_taxids(report)
    log_info(f"    {len(viral_taxids)} viral taxid(s) in database")

    # Step 3: extract viral read IDs from Kraken2 output
    viral_ids_file = tmp_dir / f"{sample_id}_viral_ids.txt"
    n_viral_reads  = 0
    if k2_out.exists() and viral_taxids:
        with open(k2_out) as f, open(viral_ids_file, "w") as o:
            for line in f:
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                if parts[0].strip() != "C":
                    continue
                try:
                    taxid = int(parts[2].strip())
                except ValueError:
                    continue
                if taxid in viral_taxids:
                    o.write(parts[1].strip() + "\n")
                    n_viral_reads += 1
    log_info(f"    {n_viral_reads} reads classified as Viruses")

    # Step 4: merge unclassified + viral → output
    log_info("  [3/3] Merging reads with seqtk...")
    for un_fq, src_fq, out_fq in [(un_r1, r1, r1_out), (un_r2, r2, r2_out)]:
        has_viral = viral_ids_file.exists() and n_viral_reads > 0
        if un_fq.exists() and has_viral:
            run_silent(
                f"{{ gzip -c {un_fq}; "
                f"seqtk subseq {src_fq} {viral_ids_file} | gzip; }} > {out_fq}",
                log_file=rpt_dir / f"{sample_id}_k2.log",
                shell=True, check=False
            )
        elif un_fq.exists():
            run_silent(f"gzip -c {un_fq} > {out_fq}",
                       log_file=rpt_dir / f"{sample_id}_k2.log",
                       shell=True, check=False)
        elif has_viral:
            run_silent(
                f"seqtk subseq {src_fq} {viral_ids_file} | gzip > {out_fq}",
                log_file=rpt_dir / f"{sample_id}_k2.log",
                shell=True, check=False)

    # Stats from report
    n_total = n_unclass = n_viral = 0
    if report.exists():
        with open(report) as f:
            for line in f:
                p = line.strip().split("\t")
                if len(p) < 6: continue
                try:
                    count = int(p[1])
                    name  = p[5].strip()
                    if name == "unclassified": n_unclass = count
                    if name == "root":         n_total   = count + n_unclass
                    if name == "Viruses":      n_viral   = count
                except ValueError:
                    pass

    n_kept   = n_unclass + n_viral
    pct_kept = f"{n_kept/n_total*100:.1f}%" if n_total else "N/A"
    pct_host = f"{(n_total-n_kept)/n_total*100:.1f}%" if n_total else "N/A"

    log_ok(f"  Unclassified={n_unclass:,} + Viral={n_viral:,} "
           f"→ kept {n_kept:,} ({pct_kept})")
    if n_total > 0 and n_unclass / n_total < 0.3:
        log_warn("  <30% unclassified — phage may be partially in the DB")

    # Cleanup
    for f in [k2_out, viral_ids_file, un_r1, un_r2]:
        if f and Path(f).exists():
            Path(f).unlink()

    return {"reads_in": str(n_total), "reads_phage": str(n_kept),
            "pct_host": pct_host, "pct_phage": pct_kept}


def _parse_viral_taxids(report: Path) -> set:
    """Extract all taxids descending from Viruses (10239) in Kraken2 report."""
    taxids = set()
    if not report.exists():
        return taxids
    in_viral    = False
    viral_depth = None
    with open(report) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6: continue
            try:
                taxid    = int(parts[4].strip())
                name_raw = parts[5]
                depth    = len(name_raw) - len(name_raw.lstrip())
                name     = name_raw.strip()
            except (ValueError, IndexError):
                continue
            if name == "Viruses" or taxid == TAXID_VIRUSES:
                in_viral = True; viral_depth = depth
                taxids.add(taxid); continue
            if in_viral:
                if depth > viral_depth:
                    taxids.add(taxid)
                else:
                    in_viral = False; viral_depth = None
    return taxids


# ---------------------------------------------------------------------------
# bwa-mem2 helpers
# ---------------------------------------------------------------------------

def _resolve_fastas(host_file, accessions, accessions_file, host_dir, rpt_dir):
    if host_file:
        p = Path(host_file)
        if p.is_dir():
            fastas = [f for f in sorted(p.iterdir()) if f.suffix.lower() in FASTA_EXTS]
            log_info(f"  Mode: local folder ({len(fastas)} FASTAs)"); return fastas
        elif p.suffix.lower() in FASTA_EXTS:
            log_info(f"  Mode: local FASTA ({p.name})"); return [p]
        else:
            valid = [Path(l.strip()) for l in open(p)
                     if l.strip() and not l.startswith("#")
                     and Path(l.strip()).exists()
                     and Path(l.strip()).suffix.lower() in FASTA_EXTS]
            log_info(f"  Mode: FASTA list ({len(valid)} files)"); return valid
    if accessions:
        log_info(f"  Mode: NCBI accessions ({len(accessions)})")
        return _download(accessions, host_dir, rpt_dir)
    if accessions_file and Path(accessions_file).exists():
        accs = [l.strip() for l in open(accessions_file)
                if l.strip() and not l.startswith("#")]
        log_info(f"  Mode: accessions file ({len(accs)})")
        return _download(accs, host_dir, rpt_dir)
    return []


def _download(accs, host_dir, rpt_dir):
    require_tools("datasets")
    dl_dir = host_dir / "ncbi_dataset"
    zip_out = host_dir / "ncbi_download.zip"
    acc_file = host_dir / "accessions.txt"
    acc_file.write_text("\n".join(accs) + "\n")
    log_info(f"  Downloading: {', '.join(accs)}")
    try:
        run_silent(["datasets","download","genome","accession",
                    "--inputfile", str(acc_file), "--include","genome",
                    "--filename", str(zip_out), "--no-progressbar"],
                   log_file=rpt_dir/"ncbi_download.log")
    except Exception as e:
        log_warn(f"  Download warning: {e}")
    fastas = []
    if zip_out.exists():
        with zipfile.ZipFile(zip_out) as z: z.extractall(dl_dir)
        zip_out.unlink()
        fastas = sorted(dl_dir.rglob("*.fna"))
        log_ok(f"  Downloaded {len(fastas)} genome(s)")
    return fastas


def _build_index(fastas, combined, rpt_dir, force):
    if Path(str(combined)+".bwt.2bit.64").exists() and not force:
        log_info("  Index already exists — skipping"); return
    if not fastas:
        log_error("No host FASTAs. Check --host-file or --accessions."); return
    log_info("  Building combined reference...")
    with open(combined,"w") as out:
        for fna in fastas:
            n = 0
            with open(fna) as f:
                for line in f:
                    out.write(f">{fna.stem}_{line[1:]}" if line.startswith(">") else line)
                    n += line.startswith(">")
            log_info(f"    + {fna.name} ({n} seqs, {human_size(fna)})")
    log_info(f"  Indexing ({human_size(combined)})...")
    run_silent(["bwa-mem2","index",str(combined)], log_file=rpt_dir/"bwa_index.log")
    run_silent(["samtools","faidx", str(combined)], log_file=rpt_dir/"bwa_index.log")
    log_ok("  Index ready ✓")


def _bwa_align(sample_id, r1, r2, r1_out, r2_out, combined, tmp_dir, rpt_dir, threads):
    bam = tmp_dir / f"{sample_id}.bam"
    log = rpt_dir / f"{sample_id}_bwa.log"
    run_silent(
        f"bwa-mem2 mem -t {threads} -M {combined} {r1} {r2} "
        f"| samtools view -bS -@ 4 "
        f"| samtools sort -n -@ 4 -T {tmp_dir}/{sample_id}_sort -o {bam}",
        log_file=log, shell=True)
    def _c(flags):
        r = subprocess.run(["samtools","view","-c"]+flags+[str(bam)],
                           capture_output=True, text=True)
        return int(r.stdout.strip() or 0)
    total = _c([]); host_r = _c(["-F","4"]); phage = _c(["-f","12","-F","256"])
    pct_host  = f"{host_r/total*100:.2f}%" if total else "0.00%"
    pct_phage = f"{phage/total*100:.2f}%"  if total else "0.00%"
    run_silent(
        f"samtools view -bS -f 12 -F 256 {bam} "
        f"| samtools fastq -@ 4 -1 {r1_out} -2 {r2_out} "
        f"-0 /dev/null -s /dev/null -n",
        log_file=log, shell=True)
    bam.unlink(missing_ok=True)
    return {"reads_in":str(total),"reads_phage":str(phage),
            "pct_host":pct_host,"pct_phage":pct_phage}


def _save_tsv(row, path):
    headers = ["sample","mode","reads_in","reads_phage","pct_host","pct_phage"]
    rows = {}
    if path.exists():
        with open(path) as f:
            f.readline()
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if cols: rows[cols[0]] = cols
    rows[row["sample"]] = [str(row.get(h,"")) for h in headers]
    with open(path,"w") as f:
        f.write("\t".join(headers)+"\n")
        for r in rows.values(): f.write("\t".join(r)+"\n")

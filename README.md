# PhageFlow

Modular bacteriophage genomics pipeline for Illumina paired-end sequencing of purified phage preparations.

**Goal:** Complete/High-quality phage genomes per sample, fully annotated.

---

## Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Pipeline Modules](#pipeline-modules)
- [Directory Structure](#directory-structure)
- [Configuration Reference](#configuration-reference)
- [Command Reference](#command-reference)
- [Scientific Rationale](#scientific-rationale)

---

## Overview

PhageFlow runs one sample at a time through a fixed module sequence. Each module writes to `results/` and `reports/`, reads its inputs from the previous module's output, and can be re-run independently with `--force`.

```
Raw reads (PE150 Illumina)
    │
    ▼ phageflow qc
01  Quality control + trimming          fastp · FastQC · MultiQC
    │
    ▼ phageflow host-removal
02  Host read removal                   bwa-mem2 (reference) · Kraken2 (no-reference)
    │  ↳ singletons retained for DTR/ITR coverage
    │
    ▼ phageflow assembly
03  De novo assembly                    SPAdes --isolate · MEGAHIT --no-mercy · cd-hit-est
    │
    ▼ phageflow viral-id
04  Viral identification                geNomad
    │
    ▼ phageflow quality
05  Genome quality + selection          CheckV · mash dereplication
    │  ↳ annotation_ready/phages/*.fasta → one per HQ candidate
    │
    ▼ phageflow annotate   (per candidate)
06  Structural annotation               Pharokka → Phold → Phynteny
    │
    ▼ phageflow report     (per candidate)
07  HTML report                         per-candidate summary
```

---

## Requirements

### Software

| Tool | Version | Used by |
|------|---------|---------|
| fastp | ≥ 1.0 | qc |
| FastQC | any | qc |
| MultiQC | any | qc |
| bwa-mem2 | any | host-removal |
| samtools | **≥ 1.15** | host-removal |
| seqtk | any | host-removal |
| kraken2 | any | host-removal (no-reference mode) |
| datasets (NCBI) | any | host-removal (accession download) |
| pigz | any | host-removal (optional, faster compression) |
| SPAdes | ≥ 3.13 (`--isolate`) | assembly |
| MEGAHIT | any | assembly |
| cd-hit-est | any | assembly |
| geNomad | ≥ 1.8 | viral-id |
| CheckV | any | quality |
| mash | any | quality (dereplication) |
| seqkit | any | quality |
| Pharokka | any | annotate |
| Phold | any | annotate |
| phynteny_transformer | any | annotate |
| dnaapler | any | annotate |

> **Critical:** samtools ≥ 1.15 is required. `blast-legacy` (pulled in by `abricate`) installs samtools 0.1.19 as a dependency, which shadows the modern version. Fix:
> ```bash
> mamba install -n phageflow "samtools>=1.15"
> samtools --version | head -1   # must show 1.15+
> ```

### Python

Python ≥ 3.11, plus: `click`, `rich`, `pyyaml`, `biopython`

### Databases

| Database | Module | Setup |
|----------|--------|-------|
| geNomad DB | viral-id | `genomad download-database databases/` |
| CheckV DB | quality | `checkv download_database databases/checkv_db` |
| Pharokka DB | annotate | `pharokka.py --install_databases` |
| Phold DB | annotate | `phold install` |
| Phynteny DB | annotate | download from Phynteny releases |
| Kraken2 DB | host-removal | only needed in no-reference mode |

---

## Installation

```bash
git clone https://github.com/your-lab/PhageFlow.git
cd PhageFlow

mamba env create -f environment.yml
conda activate phageflow
pip install -e .

# Fix samtools version (if needed)
mamba install -n phageflow "samtools>=1.15"

# Verify tools
phageflow check-tools
```

---

## Quick Start

### Single sample — with known host

```bash
# 1. Quality control
phageflow qc \
    --sample-id sampleA \
    --r1 raw/sampleA_R1.fastq.gz \
    --r2 raw/sampleA_R2.fastq.gz \
    -o /data/project/

# 2. Host removal (provide host genome accessions)
phageflow host-removal \
    --sample-id sampleA \
    --r1 /data/project/results/01_qc/sampleA_R1.fastq.gz \
    --r2 /data/project/results/01_qc/sampleA_R2.fastq.gz \
    --accessions GCF_000005845.2,GCF_009928615.1 \
    -o /data/project/

# 3. Assembly
phageflow assembly \
    --sample-id sampleA \
    --r1 /data/project/results/02_host_removal/sampleA_R1.fastq.gz \
    --r2 /data/project/results/02_host_removal/sampleA_R2.fastq.gz \
    --s1 /data/project/results/02_host_removal/sampleA_singletons.fastq.gz \
    -o /data/project/

# 4. Viral identification
phageflow viral-id \
    --sample-id sampleA \
    --contigs /data/project/results/03_assembly/sampleA_contigs.fasta \
    -o /data/project/

# 5. Quality + selection
phageflow quality \
    --sample-id sampleA \
    --virus-fna /data/project/results/04_viral_id/sampleA_virus.fna \
    -o /data/project/

# 6. Annotate each HQ candidate
for genome in /data/project/results/05_quality/sampleA/annotation_ready/phages/*.fasta; do
    candidate=$(basename "${genome}" .fasta)
    phageflow annotate \
        --sample-id sampleA \
        --genome "${genome}" \
        -o /data/project/
done
```

### Single sample — no reference host

```bash
# When you don't know the host, use Kraken2 mode (requires databases.kraken2)
phageflow host-removal \
    --sample-id sampleB \
    --r1 raw/sampleB_R1.fastq.gz \
    --r2 raw/sampleB_R2.fastq.gz \
    -o /data/project/
```

> **Note:** Kraken2 mode does not capture DTR/ITR boundary reads (singletons). If the host is known, always prefer `--accessions` or `--host-file`. See [Host Removal](#02-host-removal) below.

### Using accessions file

```bash
# hosts.txt — one GCF/GCA accession per line
cat hosts.txt
GCF_000005845.2   # E. coli K-12 MG1655
GCF_009928615.1   # Klebsiella pneumoniae

phageflow host-removal \
    --sample-id sampleA \
    --r1 ... --r2 ... \
    --accessions-file hosts.txt \
    -o /data/project/
```

---

## Pipeline Modules

### 01 — QC

Tool chain: `fastp` → `FastQC` → `MultiQC`

Key parameters (see [Scientific Rationale](#scientific-rationale)):

| Parameter | Value | Reason |
|-----------|-------|--------|
| `--qualified_quality_phred` | Q20 | 99% base accuracy |
| `--average_qual` | Q25 | stricter per-read threshold |
| `--unqualified_percent_limit` | 10% | more stringent than fastp default (40%) |
| `--length_required` | 75 bp | minimum for PE150 downstream |
| `--correction` | on | overlap-based PE error correction |
| `--trim_poly_x --poly_x_min_len 10` | 10 bp | removes poly-G artefacts (NextSeq/NovaSeq) |

Outputs:
```
results/01_qc/
    {sample}_R1.fastq.gz
    {sample}_R2.fastq.gz

reports/01_qc/{sample}/
    {sample}_fastp.html
    {sample}_fastp.json
    {sample}_R1_fastqc.html / .zip
    {sample}_R2_fastqc.html / .zip
    multiqc/multiqc_qc.html
    qc_summary.tsv
```

---

### 02 — Host Removal

Two modes based on whether a reference genome is available:

**Mode 1 — bwa-mem2 (reference known, preferred)**

Triggered by: `--host-file`, `--accessions`, or `--accessions-file`

```
bwa-mem2 mem | samtools view -bF 2304 | samtools sort -n | samtools fastq
    -f 4 -F 256
    -1 R1_out  -2 R2_out   (both mates unmapped → phage pairs)
    -s singletons           (one mate unmapped → DTR/ITR boundary reads)
```

Singletons are reads where one mate mapped to the host and one did not. These occur at genome termini (DTR/ITR junctions) and are **critical for CheckV Complete classification**. They are passed to SPAdes `--s1` in the assembly step.

**Level A contamination check** (automatic, requires `databases.kraken2`): after bwa-mem2 extraction, 50k reads are subsampled and classified with Kraken2. A warning is emitted if >5% are bacterial.

**Mode 2 — Kraken2 (no reference)**

Triggered by: no flags provided

Retains unclassified reads + Viruses (taxid 10239) at `--confidence 0.5 --minimum-hit-groups 3`.

> **Limitation:** Kraken2 paired mode processes read pairs as a unit. Boundary reads at DTR/ITR termini may be discarded if one mate has bacterial k-mers. No singletons are produced. Enable `kraken2_postfilter: true` in config.yaml for partial singleton recovery via adaptive bwa-mem2 post-filter.

**Host genome download**

Genomes are downloaded to `reports/02_host_removal/{sample}/host_genomes/` (per-sample, isolated). On subsequent runs, if the bwa-mem2 index already exists, download and re-indexing are skipped automatically. Use `--force` to refresh.

The `always_include_accessions` field in config.yaml lets you add lab-specific strains that are always merged with any user-provided accessions:

```yaml
host_removal:
  always_include_accessions:
    - GCF_000005845.2   # E. coli K-12 MG1655 (propagation host)
```

Outputs:
```
results/02_host_removal/
    {sample}_R1.fastq.gz
    {sample}_R2.fastq.gz
    {sample}_singletons.fastq.gz   (bwa-mem2 mode only)

reports/02_host_removal/{sample}/
    host_genomes/
        combined_hosts.fasta       (reference + bwa-mem2 index)
    {sample}_host_removal.log
    {sample}_k2.report             (Kraken2 mode)
    {sample}_levelA_contamination.report  (bwa-mem2 + K2 DB)
    host_removal_summary.tsv
```

---

### 03 — Assembly

Two-assembler strategy with NR reduction:

**SPAdes** (`--isolate --only-assembler -k 21,33,55,77,99,127`)

Primary assembler. `--isolate` disables BayesHammer read error correction, which is appropriate and preferred at high coverage (>200×) because correction can collapse real SNPs. Singletons passed as `--s1` if available.

**MEGAHIT** (`--no-mercy --min-count 2 --k-list 21,33,55,77,99,127`)

Secondary assembler, always runs. `--no-mercy` disables k-mer rescue at very low depth, which is counterproductive at high coverage. Both SPAdes and MEGAHIT scaffolds are combined.

**cd-hit-est** (`-c 1.00 -aS 0.85`)

Deduplicates near-identical contigs produced by both assemblers. `-aS 0.85` ensures that a shorter contig which is 85% covered by a longer one is removed, eliminating functional duplicates not caught by 100% identity alone.

Large SPAdes intermediate directories (`corrected/`, `K21/`, `K33/`, ...) are deleted after assembly to reclaim disk space. `scaffolds.fasta` and `assembly_graph.gfa` are retained.

Outputs:
```
results/03_assembly/
    {sample}_spades/
        scaffolds.fasta
        assembly_graph.gfa
        spades.log
    {sample}_megahit/
        {sample}.contigs.fa
    {sample}_contigs.fasta         (NR combined output → input to viral-id)

reports/03_assembly/{sample}/
    assembly_summary.tsv
    {sample}_assembly.log
```

---

### 04 — Viral Identification

Tool: `geNomad end-to-end`

Key parameters:

| Parameter | Value | Reason |
|-----------|-------|--------|
| `--min-score` | 0.7 | ~97% precision (Camargo et al. 2023) |
| `--enable-score-calibration` | yes | adjusts scores for virome-dominated samples |
| `--composition` | virome | calibration composition prior for purified phage |
| `--lenient-taxonomy` | yes | genus-level taxonomy → better co-binning |
| `--disable-find-proviruses` | yes | not expected in purified phage |
| `--quiet` | yes | suppress per-contig progress |

Borderline contigs (score 0.4–0.7, length ≥ 10 kb) are rescued in quality.py to avoid losing novel phages without close database relatives.

Outputs:
```
results/04_viral_id/
    {sample}_virus.fna             (viral contigs → input to quality)

reports/04_viral_id/{sample}/
    genomad_summary.tsv
    {sample}_assembly.fna_virus_summary.tsv
```

---

### 05 — Quality + Selection

Tool: `CheckV` + mash dereplication

CheckV quality tiers:

| Tier | Criteria | Destination |
|------|---------|-------------|
| Complete | CheckV Complete (DTR/ITR) | `annotation_ready/phages/` |
| High-quality | completeness ≥ 90% | `annotation_ready/phages/` |
| Medium-quality | completeness ≥ 50% | `annotation_ready/phages/` |
| Low-quality rescue | ≥1 viral gene, length ≥ 10 kb | `drafts/` |
| Large ND rescue | length ≥ 30 kb, ≥1 gene | `annotation_ready/phages/` |
| Bin rescue | same taxon, combined ≥ 30 kb | `annotation_ready/phages/` |

Dereplication removes near-identical candidates (ANI > 98% by mash, retaining the longer contig). The 98% threshold is above the ICTV species boundary of 95% ANI (Turner et al. 2021), ensuring distinct strains are kept.

The `termini_type` column from CheckV (`DTR` / `ITR` / `NA`) is propagated to `rename_map.tsv` as `topology`. This drives downstream decisions in annotate.py (circular plot vs linear, dnaapler application).

Outputs:
```
results/05_quality/{sample}/
    annotation_ready/phages/
        {Family}_candidate_001.fasta
        {Family}_candidate_002.fasta   (if multiple HQ)
    annotation_ready/proviruses/       (if any)
    drafts/
    checkv/

reports/05_quality/{sample}/
    checkv_summary.tsv
    rename_map.tsv                     (contig → candidate mapping + topology)
    assembly_summary.tsv
```

---

### 06 — Annotation

Three-tier annotation cascade:

```
Pharokka   → base GBK: CDS, tRNA, tmRNA, CRISPR, pseudogenes
    ↓
Phold      → structure-based /product updates (Foldseek + ProstT5)
    ↓
Phynteny   → /phynteny_category, /phynteny_score per CDS
    ↓
Merge      → CDS from Phynteny GBK + non-CDS from Pharokka GBK
    ↓
{candidate_id}_annotated.gbk   (canonical output)
```

Dnaapler reorientation (`dnaapler phage`) is applied before annotation when CheckV reports `termini_type = DTR` (circular phage), standardising genomes to begin at the large terminase subunit as per convention (Bouras et al. 2024, JOSS).

`genetic_code` is propagated from geNomad's `_virus_genes.tsv` (code 15 for crAss-like phages, code 11 otherwise) to ensure correct translation in Pharokka.

Outputs:
```
results/06_annotation/{sample}/{candidate_id}/
    {candidate_id}_annotated.gbk       (final merged annotation)
    {candidate_id}_pharokka/           (Pharokka raw output)
    {candidate_id}_phold/              (Phold raw output)
    {candidate_id}_phynteny/           (Phynteny raw output)

reports/06_annotation/{sample}/
    annotation_summary.tsv
    {candidate_id}_delta_report.tsv    (per-gene annotation gain across tiers)
```

---

## Directory Structure

```
PhageFlow/
├── phageflow/
│   ├── cli.py
│   ├── modules/
│   │   ├── qc.py
│   │   ├── host_removal.py
│   │   ├── assembly.py
│   │   ├── viral_id.py
│   │   ├── quality.py
│   │   └── annotate.py
│   ├── utils/
│   │   ├── config.py
│   │   ├── logger.py
│   │   └── tools.py
│   └── config/
│       └── default_config.yaml
├── config/
│   ├── config.yaml           ← your project config
│   └── samples.tsv           ← sample manifest
├── databases/                ← geNomad, CheckV, Pharokka, etc.
├── results/                  ← pipeline outputs  (or -o /path/to/project)
└── reports/                  ← logs, TSVs, HTML  (or --reports-dir /path)

Per run (-o /data/project/):
/data/project/
├── results/
│   ├── 01_qc/
│   ├── 02_host_removal/
│   ├── 03_assembly/
│   ├── 04_viral_id/
│   ├── 05_quality/
│   └── 06_annotation/
└── reports/
    ├── 01_qc/{sample}/
    ├── 02_host_removal/{sample}/
    │   └── host_genomes/    ← per-sample, isolated
    ├── 03_assembly/{sample}/
    ├── 04_viral_id/{sample}/
    ├── 05_quality/{sample}/
    └── 06_annotation/{sample}/
```

---

## Configuration Reference

Generate a config template:
```bash
phageflow config         # writes config/config.yaml
phageflow samples        # writes config/samples.tsv
```

### config.yaml

```yaml
project: PhageFlow
memory_gb: 32

# threads: auto-detected as 90% of logical CPUs
# Uncomment to override:
# threads: 16

databases:
  checkv:   databases/checkv_db/checkv-db-v1.5
  pharokka: databases/pharokka_db
  genomad:  databases/genomad_db
  phold:    databases/phold_db
  phynteny: databases/phynteny_db
  kraken2:  databases/k2_db     # only needed for Kraken2 mode

host_removal:
  # Accessions always merged with --host-file / --accessions.
  # Add lab propagation strains here so they are always included.
  always_include_accessions: []
  # Example:
  #   always_include_accessions:
  #     - GCF_000005845.2   # E. coli K-12 MG1655

  # Level A: WARN if >contamination_warn_pct% bacterial in phage reads.
  # Only runs when databases.kraken2 is configured.
  contamination_warn_pct: 5.0

  # Level B: adaptive bwa-mem2 post-filter in Kraken2 mode (opt-in).
  kraken2_postfilter: false
  postfilter_min_pct: 1.0

assembly:
  kmers:      "21,33,55,77,99,127"
  min_length: 200

genomad:
  min_score:                  0.7
  min_virus_hallmarks:        0
  rescue_min_score:           0.4
  rescue_min_length_bp:       10000
  sensitivity:                4.2
  splits:                     0       # 0 = MMseqs2 auto-manages memory
  disable_find_proviruses:    true
  enable_score_calibration:   true
  composition:                virome
  lenient_taxonomy:           true

checkv:
  min_completeness:    50
  min_viral_genes:     1
  length_rescue:       10000
  min_contig_bp:       1500
  large_nd_rescue_bp:  30000
  min_bin_rescue_bp:   30000
  min_gene_density:    0.5

annotate:
  phold_gpu:           false
  phold_batch_size:    1
  phynteny_confidence: 0.7
  genetic_code:        11
```

### samples.tsv

```tsv
sample_id	r1	r2
sampleA	/data/raw/sampleA_R1.fastq.gz	/data/raw/sampleA_R2.fastq.gz
sampleB	/data/raw/sampleB_R1.fastq.gz	/data/raw/sampleB_R2.fastq.gz
```

---

## Command Reference

All commands accept `-o DIR` (output base directory), `--reports-dir DIR`, `-t N` (threads), and `--force`.

```bash
phageflow qc \
    --sample-id ID \
    --r1 R1.fastq.gz --r2 R2.fastq.gz \
    [-o DIR] [-t N] [--force]

phageflow host-removal \
    --sample-id ID \
    --r1 R1.fastq.gz --r2 R2.fastq.gz \
    [--accessions GCF_1,GCF_2]     \  # bwa-mem2: comma-separated accessions
    [--accessions-file hosts.txt]   \  # bwa-mem2: one accession per line
    [--host-file genome.fasta]      \  # bwa-mem2: local FASTA/folder/list
    [--kraken-db /path/to/k2_db]    \  # Kraken2: explicit DB path
    [-o DIR] [-t N] [--force]

phageflow assembly \
    --sample-id ID \
    --r1 R1.fastq.gz --r2 R2.fastq.gz \
    [--s1 singletons.fastq.gz] \
    [-o DIR] [-t N] [--force]

phageflow viral-id \
    --sample-id ID \
    --contigs contigs.fasta \
    [-o DIR] [-t N] [--force]

phageflow quality \
    --sample-id ID \
    --virus-fna virus.fna \
    [-o DIR] [-t N] [--force]

phageflow annotate \
    --sample-id ID \
    --genome candidate.fasta \
    [-o DIR] [-t N] [--force]

phageflow check-tools   # verify all required tools are available
```

---

## Scientific Rationale

### Why `--isolate --only-assembler` in SPAdes?

At >200× coverage, BayesHammer (SPAdes read error correction) collapses real SNPs and introduces assembly artefacts. `--isolate` disables all pre-processing. Prjibelski et al. 2020 (*Curr Protocols* 70:e102); Roux et al. 2019 (*eLife* 8:e42923).

### Why keep singletons?

In bwa-mem2 mode, singleton reads are mates where one read maps to the host and the other does not. These reads originate at DTR/ITR boundaries (the genome termini). Including them in assembly (`--s1`) improves terminal coverage and directly increases the probability of CheckV classifying the genome as Complete. Nayfach et al. 2021 (*Nat Biotechnol* 39:578): 90% of DTR contigs with estimated completeness meet the high-quality standard.

### Why bwa-mem2 over Kraken2 when a reference is available?

Alignment-based methods have higher specificity (fewer phage reads incorrectly discarded) than k-mer classifiers. Critically, bwa-mem2 correctly identifies and retains singletons at DTR/ITR boundaries; Kraken2 paired mode discards these reads via LCA resolution. References: HoCoRT benchmark (Kracherberger et al. 2023); Nayfach et al. 2021.

### Why cd-hit-est `-aS 0.85` and not `-aS 1.00`?

`-aS 1.00` only removes exact-length duplicates. `-aS 0.85` additionally removes contained sequences: a shorter contig that is entirely contained within a longer one (with ≥85% coverage) is a functional duplicate, even if not a perfect-length match. This situation commonly arises when SPAdes and MEGAHIT agree on the same genomic region but produce slightly different-length scaffolds. Fu et al. 2012 (*Bioinformatics* 28:3150).

### Why 98% ANI for dereplication in quality.py?

The ICTV species boundary for bacteriophages is 95% ANI (Turner et al. 2021, *Arch Virol* 166:2633). Candidates above 98% ANI are the same strain assembled twice (or two near-identical preparations of the same phage). Candidates between 95–98% are distinct strains of the same species and are both retained.

### Genetic code 15 for crAss-like phages

CrAss-like phages (Crassvirales) use an alternative genetic code where TGA encodes tryptophan instead of acting as a stop codon. geNomad's `_virus_genes.tsv` reports the genetic code used per contig. PhageFlow propagates this to Pharokka's `--genetic_code` to ensure correct ORF prediction. Koonin et al. 2020 (*mBio* 11:e00278-20).

---

## Key References

| Reference | DOI | Used for |
|-----------|-----|---------|
| Chen et al. 2018 | 10.1186/s13059-018-1568-0 | fastp QC parameters |
| Wick & Holt 2022 | 10.1099/mgen.0.000788 | QC quality thresholds |
| Bankevich et al. 2012 | 10.1089/cmb.2012.0021 | SPAdes k-mers |
| Prjibelski et al. 2020 | 10.1002/cpbi.102 | SPAdes `--isolate` |
| Roux et al. 2019 | 10.7554/eLife.42923 | `--only-assembler`, MIUViG standards |
| Li et al. 2015 | 10.1093/bioinformatics/btv033 | MEGAHIT `--no-mercy` |
| Fu et al. 2012 | 10.1093/bioinformatics/bts565 | cd-hit-est `-aS 0.85` |
| Nayfach et al. 2021 | 10.1038/s41587-020-00774-7 | CheckV, singletons, DTR/ITR |
| Wood et al. 2019 | 10.1186/s13059-019-1891-0 | Kraken2 confidence parameters |
| Vasimuddin et al. 2019 | 10.1109/IPDPS.2019.00041 | bwa-mem2 |
| Li et al. 2009 | 10.1093/bioinformatics/btp352 | samtools flags |
| Camargo et al. 2023 | 10.1038/s41587-023-01953-y | geNomad parameters |
| Turner et al. 2021 | 10.1007/s00705-021-05156-1 | ICTV ANI species boundary |
| Ondov et al. 2016 | 10.1186/s13059-016-0997-x | mash ANI dereplication |
| Bouras et al. 2023 | 10.1093/bioinformatics/btac776 | Pharokka |
| Bouras et al. 2024 | 10.21105/joss.05968 | dnaapler reorientation |
| Grigson et al. 2025 | bioRxiv 2025 | Phynteny transformer |

---

*PhageFlow — purified phage genomics, from reads to annotated genome.*

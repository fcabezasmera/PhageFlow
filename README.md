# PhageFlow

**Modular bacteriophage genomics pipeline for complete genome recovery from Illumina paired-end sequencing of purified phage preparations.**

[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-0.1.0-green.svg)](CHANGELOG.md)

---

## Authors

| Name | Affiliation |
|------|-------------|
| Fausto Cabezas-Mera | Universidad Tecnológica Metropolitana, Santiago, Chile |
| Estefania Tisalema-Guanopatin | Universidad Tecnológica Metropolitana, Santiago, Chile |
| Dayra Valle | Universidad Tecnológica Metropolitana, Santiago, Chile |
| Antonella Nole | Universidad Tecnológica Metropolitana, Santiago, Chile |

---

## Table of Contents

1. [Overview](#overview)
2. [Pipeline Architecture](#pipeline-architecture)
3. [Requirements](#requirements)
4. [Installation](#installation)
5. [Database Setup](#database-setup)
6. [Quick Start](#quick-start)
7. [Module Reference](#module-reference)
   - [Module 01: Quality Control](#module-01-quality-control)
   - [Module 02: Host Read Removal](#module-02-host-read-removal)
   - [Module 03: De Novo Assembly](#module-03-de-novo-assembly)
   - [Module 04: Viral Identification](#module-04-viral-identification)
   - [Module 05: Genome Quality Assessment](#module-05-genome-quality-assessment)
   - [Module 05b: Assembly Refinement](#module-05b-assembly-refinement-optional)
   - [Module 06: Structural Annotation](#module-06-structural-annotation)
8. [Configuration Reference](#configuration-reference)
9. [Output Structure](#output-structure)
10. [CLI Reference](#cli-reference)
11. [Scientific Basis](#scientific-basis)
12. [References](#references)
13. [License](#license)

---

## Overview

PhageFlow is designed for researchers working with **purified phage preparations** sequenced on Illumina paired-end platforms (PE150). It implements a sequential, module-based workflow that processes raw reads through quality control, host depletion, de novo assembly, viral identification, genome quality assessment, and multi-tier structural annotation.

**Design principles:**

- **Modularity**: each module is independently re-runnable with `--force` to re-execute or `--from-module` to resume
- **Reproducibility**: all parameters are controlled by a single `config/config.yaml`; no hard-coded thresholds
- **Transparency**: every module writes structured TSV reports and logs all external commands
- **Scientific rigour**: all methodological decisions are documented with primary literature references

**Supported sample types:**

| Sample type | Recommended mode | Notes |
|-------------|-----------------|-------|
| Purified phage (single) | bwa-mem2 + SPAdes `--isolate` | Primary use case |
| Purified phage (mixed) | bwa-mem2 + SPAdes `--isolate` | Multiple hosts in `always_include_accessions` |
| Low-coverage preparation | MEGAHIT fallback | Automatic when SPAdes produces no contigs |
| Fragmented genome | assembly-refine | Run after quality step |

---

## Pipeline Architecture

```
Raw reads (PE150 Illumina)
    │
    ▼  phageflow qc
01  Quality control + trimming
    fastp (adapter trim · PE correction · poly-X · complexity filter)
    FastQC (per-read QC)
    MultiQC (aggregate report)
    │
    ▼  phageflow host-removal
02  Host read removal
    bwa-mem2 (alignment-based, streaming, no BAM on disk)
    ↳ singletons retained for DTR/ITR boundary coverage
    Optional: Kraken2 (classification-based, with optional bwa-mem2 post-filter)
    Optional: Level A contamination check (subsample Kraken2 diagnostic)
    │
    ▼  phageflow assembly
03  De novo assembly
    SPAdes --isolate (primary; with singletons as --s1)
    MEGAHIT --no-mercy --min-count 2 (secondary)
    cd-hit-est (NR reduction at 100% identity)
    │
    ▼  phageflow viral-id
04  Viral identification
    geNomad end-to-end (k-mer + neural network)
    ↳ topology inference (DTR / ITR / No terminal repeats)
    ↳ genetic code detection (standard=11, CrAss-like=15)
    ↳ taxonomy (family → finest level, --lenient-taxonomy)
    │
    ▼  phageflow quality
05  Genome quality and selection
    CheckV end-to-end (completeness via AAI + HMM)
    mash / minimap2 (dereplication at >98% ANI)
    co-binning (LQ drafts → annotation_ready when combined ≥30 kb)
    bwa-mem2 (read recruitment coverage validation)
    │
    [optional: phageflow assembly-refine]
05b Iterative assembly refinement
    Reads unmapped to candidates → SPAdes --trusted-contigs
    CheckV re-evaluation → replace improved candidates
    │
    ▼  phageflow annotate
06  Structural annotation (three-tier cascade)
    Tier 1  Pharokka   gene calling (PHANOTATE/prodigal-gv)
                       PHROGs MMseqs2 + PyHMMER
    Tier 2  Phold      ProstT5 + Foldseek structure-based upgrade
    Tier 3  Phynteny   Transformer + ESM2 synteny upgrade
    Plot               phold plot (circular) OR pyGenomeViz (linear multi-track)
    │
    ▼  phageflow safety
    Biosafety screening  CARD (AMR) + VFDB (virulence factors)
    │
    ▼  phageflow report
    Final HTML report per candidate genome
```

---

## Requirements

### Hardware

| Resource | Minimum | Recommended |
|----------|---------|-------------|
| CPU | 8 cores | 24+ cores |
| RAM | 16 GB | 32–64 GB |
| GPU | — | NVIDIA GPU for Phold (ProstT5 + Foldseek) |
| Storage | 50 GB | 200 GB (databases ~80 GB) |

### Software

- Linux (Ubuntu 20.04+ recommended)
- conda / mamba ≥ 23.x
- CUDA ≥ 11.8 (optional, for `phold --foldseek_gpu`)

---

## Installation

```bash
git clone https://github.com/fcabezasmera/PhageFlow.git
cd PhageFlow

# Create environment (large — ~80 tools, ~8 GB)
mamba env create -f environment.yml
conda activate phageflow

# Install PhageFlow package (editable mode for development)
pip install -e .

# Verify installation
phageflow check-tools
```

> **samtools version conflict**: if `abricate` or `blast-legacy` are present,
> force samtools ≥ 1.15 (required for `samtools fastq -N`):
> ```bash
> mamba install -n phageflow "samtools>=1.15"
> ```

---

## Database Setup

PhageFlow requires several external databases. All paths are configured in `config/config.yaml` under the `databases:` section.

### CheckV (~2 GB)

```bash
checkv download_database databases/checkv_db
```

Default path: `databases/checkv_db/checkv-db-v1.5`

### geNomad (~7 GB)

```bash
genomad download-database databases/genomad_db
```

Default path: `databases/genomad_db`

### Pharokka (~2 GB)

```bash
install_databases.py -o databases/pharokka_db
```

Default path: `databases/pharokka_db`

### Phold (~18 GB with ESM2 embeddings)

```bash
phold install -d databases/phold_db
```

Default path: `databases/phold_db`

### Phynteny Transformer (~1 GB)

```bash
phynteny_transformer install -d databases/phynteny_db
```

Default path: `databases/phynteny_db`

### Kraken2 (optional, for classification-based host removal)

Download a pre-built database from [https://benlangmead.github.io/aws-indexes/k2](https://benlangmead.github.io/aws-indexes/k2). Set `databases.kraken2` in config.yaml.

---

## Quick Start

### 1. Initialize a project

```bash
phageflow init /path/to/myproject
```

This creates the project structure:
```
myproject/
  config/config.yaml   ← edit this
  raw/                 ← place your FASTQ files here
  results/             ← pipeline outputs (auto-created)
  reports/             ← logs, TSVs, HTML reports (auto-created)
```

### 2. Place FASTQ files and edit config

```bash
cp /path/to/reads/*.fastq.gz /path/to/myproject/raw/

# Edit database paths and thread count
nano /path/to/myproject/config/config.yaml
```

### 3. Run the full pipeline

```bash
# Auto-detect samples from raw/, download host genome from NCBI
phageflow run \
    --project /path/to/myproject \
    --accessions GCF_000005845.2   # propagation host accession

# Multiple hosts
phageflow run \
    --project /path/to/myproject \
    --accessions GCF_000005845.2,GCF_000013425.1
```

### 4. Run individual modules

```bash
phageflow qc \
    --r1 raw/sample_R1.fastq.gz \
    --r2 raw/sample_R2.fastq.gz

phageflow host-removal \
    --r1 results/01_qc/sample_R1.fastq.gz \
    --r2 results/01_qc/sample_R2.fastq.gz \
    --accessions GCF_000005845.2

phageflow assembly \
    --r1 results/02_host_removal/sample_R1.fastq.gz \
    --r2 results/02_host_removal/sample_R2.fastq.gz \
    --s1 results/02_host_removal/sample_singletons.fastq.gz

phageflow viral-id \
    --contigs results/03_assembly/sample_contigs.fasta

phageflow quality \
    --virus-fna results/04_viral_id/sample_virus.fna

phageflow annotate \
    --genome results/05_quality/sample/annotation_ready/phages/Drexlerviridae_candidate_001.fasta
```

### 5. Resume an interrupted run

```bash
phageflow run --project /path/to/myproject --from-module quality
```

### 6. Check pipeline status

```bash
phageflow status --project /path/to/myproject
```

---

## Module Reference

### Module 01: Quality Control

**Command**: `phageflow qc`

**Tools**: fastp · FastQC · MultiQC

Performs adapter trimming, quality filtering, paired-end overlap correction, poly-X tail removal, and low-complexity filtering. Runs FastQC on trimmed reads (R1+R2 in parallel) and aggregates results with MultiQC.

#### fastp parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `--detect_adapter_for_pe` | auto | Automatic adapter detection (Chen et al. 2018) |
| `--correction` | enabled | PE overlap-based base error correction |
| `--overlap_len_require` | 10 bp | Minimum overlap for correction |
| `--overlap_diff_percent_limit` | 10% | Max mismatch rate in overlap |
| `--cut_right` | enabled | Sliding-window 3′ quality trimming (Bolger et al. 2014) |
| `--cut_right_window_size` | 4 bp | Window size for sliding-window trimming |
| `--cut_right_mean_quality` | Q20 | Quality floor for sliding-window trimming |
| `--qualified_quality_phred` | 20 | Per-base quality threshold (99% accuracy) |
| `--unqualified_percent_limit` | 10% | Max fraction of low-quality bases per read |
| `--average_qual` | 20 | Read-level mean quality floor |
| `--n_base_limit` | 5 | Max N calls per read |
| `--length_required` | 75 bp | Minimum read length for PE150 |
| `--low_complexity_filter` | enabled | Removes homopolymer reads (Roux et al. 2019) |
| `--complexity_threshold` | 15% | Minimum sequence complexity |
| `--trim_poly_x` | enabled | Poly-X tail removal (targets NextSeq/NovaSeq poly-G) |
| `--poly_x_min_len` | 10 bp | Minimum poly-X length to trigger trimming |

> **Note on duplication**: fastp uses k-mer sampling, not coordinate-based deduplication. At high phage coverage (100–5000×), apparent duplication rates of 50–70% are normal (Head et al. 2014). PhageFlow does NOT deduplicate reads.

#### Thresholds and warnings

| Metric | Warning level | Notes |
|--------|--------------|-------|
| Pass rate | < 75% | Expected 75–85% with strict HQ filters (Wick & Holt 2022) |
| Q30 rate | < 75% | Assembly contiguity may be affected |
| Total reads after filtering | < 50,000 | Minimum for 50× coverage on 150 kb genome |
| Duplication rate | > 70% | Normal at high phage coverage; informational only |

#### Outputs

```
results/01_qc/
    {sample}_R1.fastq.gz      ← trimmed reads
    {sample}_R2.fastq.gz

reports/01_qc/{sample}/
    {sample}_fastp.json        ← raw metrics (parsed by downstream modules)
    {sample}_fastp.html        ← fastp HTML report
    {sample}_R1_fastqc.html    ← FastQC per-read report
    {sample}_R2_fastqc.html
    multiqc/multiqc_qc.html    ← aggregated QC report
    qc_summary.tsv             ← structured metrics TSV
```

---

### Module 02: Host Read Removal

**Command**: `phageflow host-removal`

**Tools**: bwa-mem2 · samtools ≥ 1.15 (primary) or Kraken2 · seqtk (alternative)

Removes host reads from the phage preparation. The primary mode uses alignment-based removal (bwa-mem2), which is more specific than k-mer classification and correctly identifies and retains DTR/ITR boundary reads as singletons for downstream assembly.

#### bwa-mem2 mode (default, recommended)

```
bwa-mem2 mem | samtools view -bF 2304 | samtools sort -n | samtools fastq
    -f 4 -F 256    → unmapped reads (phage)
    -1/-2          → both unmapped → phage pairs
    -s             → one unmapped  → DTR/ITR singletons → SPAdes --s1
    -0 /dev/null   → mapped (host) → discard
```

The streaming pipeline avoids writing a BAM file to disk. The `-F 2304` flag excludes secondary and supplementary alignments. Singleton recovery is critical for phages with DTR/ITR termini — these reads span the junction between the end and start of the linear/circular map and are discarded by paired-mode tools (Nayfach et al. 2021).

**Host references** can be provided as:
- NCBI accession(s): `--accessions GCF_000005845.2` (auto-downloaded via `datasets`)
- Local FASTA: `--host-file /path/to/host.fasta`
- Folder of FASTAs: `--host-file /path/to/genomes/`
- Text file of paths: `--host-file paths.txt`

Multiple accessions are concatenated into a single combined reference for one-pass alignment.

#### Kraken2 mode (alternative)

Used when no reference genome is available. Retains unclassified reads + all viral reads (taxid 10239).

> **Limitation**: Kraken2 paired mode discards singletons at DTR/ITR boundaries. Enable `kraken2_postfilter: true` in config.yaml to add a bwa-mem2 post-filter step that recovers these reads.

#### Level A contamination check

After bwa-mem2 removal, a subsample of 50,000 phage reads is classified with Kraken2 as a diagnostic. Any bacterial genus > 5% (configurable via `contamination_warn_pct`) triggers a warning. Reads are NOT removed — this is diagnostic only.

#### Outputs

```
results/02_host_removal/
    {sample}_R1.fastq.gz          ← phage paired reads
    {sample}_R2.fastq.gz
    {sample}_singletons.fastq.gz  ← DTR/ITR boundary reads

reports/02_host_removal/{sample}/
    host_genomes/
        combined_hosts.fasta      ← concatenated reference
        combined_hosts.fasta.*    ← bwa-mem2 index
    {sample}_levelA_contamination.report  ← Kraken2 diagnostic (if DB configured)
    host_removal_summary.tsv
    {sample}_host_removal.log
```

---

### Module 03: De Novo Assembly

**Command**: `phageflow assembly`

**Tools**: SPAdes · MEGAHIT · cd-hit-est

Two-assembler strategy with non-redundant (NR) reduction. Both assemblers run independently; their contigs are combined, filtered by minimum length, and deduplicated with cd-hit-est.

#### SPAdes (primary)

```
spades.py --pe1-1 R1 --pe1-2 R2 --s1 singletons
          --isolate -k 21,33,55,77,99,127 --phred-offset 33
```

`--isolate` is the recommended preset for high-coverage isolates. It disables BayesHammer error correction (unnecessary and harmful at >200× coverage, where it collapses real SNPs) and optimises graph construction for near-complete genomes (Prjibelski et al. 2020). The `--s1 singletons` flag passes DTR/ITR boundary reads recovered by host-removal, improving terminal coverage and CheckV Complete classification.

> **Important**: `--only-assembler` is implied by `--isolate` and must **not** be passed explicitly — SPAdes raises a compatibility error if both flags are provided.

#### MEGAHIT (secondary)

```
megahit -1 R1 -2 R2 --k-list 21,33,55,77,99,127,141
        --no-mercy --min-count 2 --min-contig-len 200
```

`--no-mercy` disables low-depth k-mer rescue, reducing chimeric contigs at high coverage. `--min-count 2` discards k-mers seen only once (sequencing errors). MEGAHIT does not support singleton input; singletons are omitted.

> **Note**: `--k-list` requires MEGAHIT ≥ 1.2.9. PhageFlow auto-detects the installed version and falls back to `--k-min/--k-max/--k-step` for older versions.

#### cd-hit-est NR reduction

```
cd-hit-est -c 1.00 -aS 0.85 -n 8 -d 0 -M 0 -p 1
```

Removes exact duplicates and near-identical contigs (≥85% coverage at 100% identity). Reduces redundancy from assembler agreement on the same genomic region.

#### Warnings

| Condition | Warning |
|-----------|---------|
| N50 < 5,000 bp | Fragmented assembly — check host removal and read depth |
| Largest contig > 500,000 bp | Unusually large — possible host contamination or concatenated assembly |
| Both assemblers fail | Pipeline aborts for this sample |

#### Outputs

```
results/03_assembly/
    {sample}_spades/
        contigs.fasta           ← SPAdes contigs (scaffolds.fasta NOT used)
        assembly_graph.gfa      ← assembly graph
        spades.log
    {sample}_megahit/
        {sample}.contigs.fa     ← MEGAHIT contigs
    {sample}_contigs.fasta      ← NR combined output → viral-id input

reports/03_assembly/{sample}/
    assembly_summary.tsv
    {sample}_assembly.log
```

> **Why not scaffolds.fasta?** SPAdes scaffolds contain N-runs at gap positions. These break CheckV DTR/ITR detection (which requires an uninterrupted terminal repeat sequence) and interfere with genome circularity assessment.

---

### Module 04: Viral Identification

**Command**: `phageflow viral-id`

**Tools**: geNomad

Classifies contigs as viral or non-viral using geNomad's combined k-mer marker gene search (MMseqs2) and neural-network classifier. Extracts taxonomy, topology, and genetic code for each viral contig.

#### geNomad flags

| Flag | Value | Rationale |
|------|-------|-----------|
| `--enable-score-calibration` | enabled | Adjusts scores for virome-dominated samples |
| `--composition` | virome | Reduces false negatives for novel lineages |
| `--lenient-taxonomy` | enabled | Resolves below family rank (genus, species) |
| `--disable-find-proviruses` | enabled | Proviruses not expected in purified preps |
| `--sensitivity` | 4.2 | Default MMseqs2 sensitivity |
| `--cleanup` | enabled | Removes large `annotate/` intermediate directory |

#### Filtering tiers

| Tier | Condition | Rationale |
|------|-----------|-----------|
| Main | `virus_score ≥ 0.7` | ~97% precision (Camargo et al. 2023) |
| Rescued | `0.4 ≤ score < 0.7` AND `length ≥ 10 kb` | Novel lineages without DB relatives score 0.4–0.6 |
| Discarded | Everything else | |

#### Metadata propagated downstream

| Field | Source | Used by |
|-------|--------|---------|
| `topology` | geNomad `_virus_summary.tsv` | quality.py → rename_map; annotate.py → plot mode, dnaapler |
| `genetic_code` | geNomad (11=standard, 15=CrAss-like) | annotate.py → Pharokka `--genetic_code` |
| `naming_level` | Taxonomy parse (family > finest > 'Phage') | quality.py → candidate file naming |

#### Outputs

```
results/04_viral_id/
    {sample}_virus.fna          ← viral contigs → quality input
    {sample}_metadata.tsv       ← topology + genetic_code + score per contig

reports/04_viral_id/{sample}/
    genomad/
        {stem}/
            {stem}_virus_summary.tsv
            {stem}_virus_genes.tsv
    viral_id_summary.tsv
    {sample}_viral_id.log
```

---

### Module 05: Genome Quality Assessment

**Command**: `phageflow quality`

**Tools**: CheckV · mash · minimap2 (or cd-hit-est fallback)

Evaluates viral contigs with CheckV, assigns quality tiers, dereplicates near-identical candidates, co-bins low-quality draft contigs by taxonomy, validates genome integrity with read recruitment coverage, and writes annotation-ready FASTAs.

#### CheckV quality tiers

| Tier | Condition | Output directory |
|------|-----------|-----------------|
| `complete` | CheckV = Complete (DTR or ITR detected) | `annotation_ready/phages/` |
| `high-quality` | completeness ≥ 90% | `annotation_ready/phages/` |
| `medium-quality` | completeness ≥ 50% (configurable) | `annotation_ready/phages/` |
| `large-nd` | Not-determined, length ≥ 30 kb, ≥ 1 viral gene | `annotation_ready/phages/` |
| `bin-rescue` | LQ drafts co-binned by taxonomy, combined ≥ 30 kb | `annotation_ready/phages/` |
| `lq-draft` | length ≥ 10 kb, ≥ 1 viral gene | `drafts/` |
| `discarded` | below all thresholds | not written |

#### Topology resolution (priority order)

CheckV `termini_type` takes precedence over geNomad topology because CheckV confirms terminal repeats in the actual assembled sequence (≥20 bp detected repeat), while geNomad predicts topology from sequence composition.

```
CheckV DTR   → topology = "DTR"  (circular, most tailed dsDNA phages)
CheckV ITR   → topology = "ITR"  (linear with ITR, e.g. T7-like Autographiviridae)
CheckV NA    → use geNomad topology as context ("No terminal repeats" or "NA")
```

#### Dereplication

Applied to single-contig annotation_ready candidates only, at >98% ANI (Turner et al. 2021 ICTV species boundary = 95% ANI).

| Genome size | Method |
|-------------|--------|
| ≥ 20 kb | mash distance (sketch size 10,000) |
| < 20 kb | minimap2 asm5 all-vs-all (mash unreliable on small genomes) |
| Fallback | cd-hit-est at 98% if neither available |

#### Candidate naming scheme

Candidates are named after their taxonomic level: `{naming_level}_candidate_{NNN}.fasta`

Examples: `Drexlerviridae_candidate_001.fasta`, `Ackermannviridae_candidate_002.fasta`, `Phage_candidate_001.fasta` (unclassified)

#### Coverage validation

Post-assignment, phage reads are mapped back to annotation_ready candidates (bwa-mem2) to compute mean depth, breadth, and coefficient of variation (CV). Warnings are raised for:
- Breadth < 95% (potential assembly gap or chimera)
- CV > 1.5 (uneven coverage, possible concatemer)

#### Outputs

```
results/05_quality/{sample}/
    annotation_ready/phages/
        {naming_level}_candidate_001.fasta   ← annotation input
        {naming_level}_candidate_002.fasta
    drafts/
        {naming_level}_draft_001.fasta       ← low-quality candidates
    checkv/
        quality_summary.tsv
        contamination.tsv
        completeness.tsv

reports/05_quality/{sample}/
    rename_map.tsv          ← contig→candidate with full metadata
    quality_summary.tsv
    coverage/
        depth.tsv
    {sample}_quality.log
```

---

### Module 05b: Assembly Refinement (Optional)

**Command**: `phageflow assembly-refine`

Run **after** `quality` when annotation_ready contains fragmented assemblies (MQ candidates, multiple candidates from same family, breadth < 95%).

**Strategy**: extracts reads that do not map to existing candidates (unmapped pairs) and reads that bridge a candidate end and a gap (semi-mapped "bridge" reads), then re-assembles with SPAdes using `--trusted-contigs` to anchor the graph at known correct sequences and extend through gaps.

```bash
phageflow assembly-refine \
    --r1 results/02_host_removal/sample_R1.fastq.gz \
    --r2 results/02_host_removal/sample_R2.fastq.gz \
    --s1 results/02_host_removal/sample_singletons.fastq.gz
```

Refined candidates replace originals when CheckV quality tier is strictly higher or genome is ≥10% longer at equal tier. Original candidates are backed up to `annotation_ready/phages/pre_refine/`.

---

### Module 06: Structural Annotation

**Command**: `phageflow annotate`

**Tools**: Pharokka → Phold → Phynteny Transformer

Three-tier annotation cascade with per-CDS delta tracking from tool TSVs.

#### Tier 1: Pharokka

Gene calling with PHANOTATE (single-contig) or prodigal-gv (multi-contig), PHROGs database search via MMseqs2 and PyHMMER. For single-contig genomes, `--dnaapler` reorients the sequence to canonical origin (terminase large subunit). For multi-contig genomes, `--meta --meta_hmm` are used; `--dnaapler` is not added (incompatible with meta mode).

#### Tier 2: Phold

Structure-based functional annotation via ProstT5 (protein language model embeddings) + Foldseek (structural alignment). Upgrades hypothetical proteins from Tier 1 by finding structural homologs with known function. Key flags: `--hyps` (process hypothetical proteins only), `--finetune` (phage-finetuned ProstT5 model), `--foldseek_gpu` (GPU acceleration).

`annotation_confidence` (high / medium / low) is the primary quality indicator per phold documentation.

#### Tier 3: Phynteny Transformer

Synteny-aware functional prediction using Transformer + ESM2. Adds `/phynteny_category`, `/phynteny_score`, `/phynteny_confidence` qualifiers to the GBK. Does NOT modify `/product` — predictions are additive.

> **Important**: The per-CDS TSV output filename contains a confirmed typo in the tool: `phynteny_per_cds_funcions.tsv` (not `functions`). PhageFlow handles this with a fallback search.

#### Plot strategy

| Genome | Method | Rationale |
|--------|--------|-----------|
| Single-contig | phold plot (circular PNG + SVG) | Reads `/product`, Foldseek qualifiers, confidence |
| Multi-contig | pyGenomeViz (linear multi-track) | Single figure for all contigs; phold generates N separate maps |

#### Delta tracking

All three tier TSVs are joined on `locus_tag` to produce a `{candidate_id}_delta_report.tsv` recording which tier annotated each CDS and the full evidence chain.

#### Outputs

```
results/06_annotation/{candidate_id}/
    pharokka/
        {candidate_id}.gbk      ← Tier 1 GBK
        {candidate_id}.gff      ← GFF3 for external tools
        {candidate_id}_cds_final_merged_output.tsv
    phold/
        {candidate_id}.gbk      ← Tier 2 GBK (canonical output)
        {candidate_id}_per_cds_predictions.tsv
    phynteny/
        phynteny_transformer.gbk  ← Tier 3 GBK (/phynteny_* only)
        phynteny_per_cds_funcions.tsv
    plots/
        {candidate_id}.png      ← circular or linear genome map
        {candidate_id}.svg

reports/06_annotation/
    annotation_summary.tsv      ← aggregate stats across all candidates
    {candidate_id}/
        {candidate_id}_delta_report.tsv
        {candidate_id}_pharokka.log
        {candidate_id}_phold.log
        {candidate_id}_phynteny.log
        {candidate_id}_plot.log
```

---

## Configuration Reference

All parameters are set in `config/config.yaml`. Run `phageflow config` to generate a template.

### QC (`qc:`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `average_qual` | 20 | fastp `--average_qual`: read-level mean quality floor |
| `qualified_quality_phred` | 20 | Per-base quality threshold |
| `unqualified_percent_limit` | 10 | Max % low-quality bases per read |
| `length_required` | 75 | Minimum read length (bp) |
| `cut_right_window_size` | 4 | Sliding-window size (bp) |
| `cut_right_mean_quality` | 20 | Sliding-window quality floor |
| `correction` | true | PE overlap correction |
| `overlap_len_require` | 10 | Minimum overlap for PE correction (bp) |
| `low_complexity_filter` | true | Enable low-complexity filter |
| `complexity_threshold` | 15 | Minimum complexity (%) |
| `trim_poly_x` | true | Poly-X tail removal |
| `poly_x_min_len` | 10 | Minimum poly-X length to trim (bp) |
| `n_base_limit` | 5 | Max N calls per read |

### Host Removal (`host_removal:`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `always_include_accessions` | `[]` | NCBI accessions always added to host reference |
| `contamination_warn_pct` | 5.0 | Level A warning threshold (% bacterial reads) |
| `kraken2_postfilter` | false | Enable bwa-mem2 post-filter after Kraken2 |
| `postfilter_min_pct` | 1.0 | Minimum bacterial % to trigger Level B download |

### Assembly (`assembly:`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `kmers` | `"21,33,55,77,99,127,141"` | k-mer sizes for SPAdes and MEGAHIT |
| `min_length` | 200 | Minimum contig length (bp) for NR pool |
| `iterative_refinement` | false | Run assembly-refine automatically after quality |

### Viral Identification (`genomad:`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_score` | 0.7 | Primary viral score threshold (~97% precision) |
| `rescue_min_score` | 0.4 | Lower threshold for rescue tier |
| `rescue_min_length_bp` | 10000 | Minimum length for rescue tier (bp) |
| `sensitivity` | 4.2 | MMseqs2 search sensitivity |
| `enable_score_calibration` | true | Score calibration for virome samples |
| `composition` | virome | Sample composition for calibration |
| `lenient_taxonomy` | true | Resolve below family rank |
| `disable_find_proviruses` | true | Skip provirus detection |

### Quality Assessment (`checkv:`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_completeness` | 50 | MQ threshold (%) |
| `min_viral_genes` | 1 | Minimum viral genes for draft rescue |
| `length_rescue` | 10000 | Minimum length for LQ draft (bp) |
| `min_contig_bp` | 1500 | Global minimum contig length (bp) |
| `large_nd_rescue_bp` | 30000 | Length threshold for large-ND rescue (bp) |
| `min_bin_rescue_bp` | 30000 | Combined length threshold for bin rescue (bp) |
| `max_kmer_freq` | 1.5 | Concatemer warning threshold |
| `max_genome_copies` | 1.5 | Concatemer warning threshold |

### Annotation (`annotate:`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `phold_gpu` | true | Use `--foldseek_gpu` (disable for CPU-only nodes) |
| `phold_finetune` | true | Phage-finetuned ProstT5 model |
| `phold_batch_size` | 1 | ProstT5 batch size (increase on GPU: 4–8) |
| `phold_sensitivity` | 9.5 | Foldseek sensitivity |
| `phynteny_confidence` | 0.8 | Phynteny prediction confidence threshold |
| `genetic_code` | 11 | Default genetic code (overridden per-sample from geNomad) |

---

## Output Structure

```
<project>/
├── results/
│   ├── 01_qc/
│   │   ├── {sample}_R1.fastq.gz
│   │   └── {sample}_R2.fastq.gz
│   ├── 02_host_removal/
│   │   ├── {sample}_R1.fastq.gz
│   │   ├── {sample}_R2.fastq.gz
│   │   └── {sample}_singletons.fastq.gz
│   ├── 03_assembly/
│   │   ├── {sample}_spades/
│   │   ├── {sample}_megahit/
│   │   └── {sample}_contigs.fasta         ← NR combined contigs
│   ├── 04_viral_id/
│   │   ├── {sample}_virus.fna             ← viral contigs
│   │   └── {sample}_metadata.tsv          ← topology + genetic_code
│   ├── 05_quality/
│   │   └── {sample}/
│   │       ├── annotation_ready/
│   │       │   └── phages/
│   │       │       └── {naming_level}_candidate_NNN.fasta   ← annotation input
│   │       ├── drafts/
│   │       └── checkv/
│   └── 06_annotation/
│       └── {candidate_id}/
│           ├── pharokka/                  ← Tier 1 (Pharokka)
│           ├── phold/                     ← Tier 2 (Phold) ← canonical GBK
│           ├── phynteny/                  ← Tier 3 (Phynteny)
│           └── plots/                     ← genome map (PNG + SVG)
└── reports/
    ├── 01_qc/{sample}/
    │   ├── {sample}_fastp.json
    │   ├── {sample}_fastp.html
    │   ├── multiqc/multiqc_qc.html
    │   └── qc_summary.tsv
    ├── 02_host_removal/{sample}/
    ├── 03_assembly/{sample}/
    ├── 04_viral_id/{sample}/
    ├── 05_quality/{sample}/
    │   └── rename_map.tsv                 ← critical metadata link
    ├── 06_annotation/
    │   ├── annotation_summary.tsv         ← aggregate (all candidates)
    │   └── {candidate_id}/
    │       └── {candidate_id}_delta_report.tsv
    └── pipeline_status.tsv                ← module completion tracking
```

### Key files for downstream use

| File | Description | Used for |
|------|-------------|----------|
| `results/05_quality/{sample}/annotation_ready/phages/*.fasta` | Annotation-ready candidate genomes | Annotate, safety screening, NCBI submission |
| `reports/05_quality/{sample}/rename_map.tsv` | Full metadata per contig (topology, completeness, genetic_code) | Interpreting annotation results |
| `results/06_annotation/{candidate}/phold/{candidate}.gbk` | Canonical annotated GenBank (Tier 2) | Genome browser, comparative genomics |
| `reports/06_annotation/annotation_summary.tsv` | Per-candidate annotation statistics | Manuscript preparation |
| `reports/06_annotation/{candidate}/{candidate}_delta_report.tsv` | Per-CDS annotation provenance | Quality filtering, functional analysis |

---

## CLI Reference

```
phageflow [OPTIONS] COMMAND [ARGS]

Commands:
  qc               Quality control and trimming (fastp + FastQC + MultiQC)
  host-removal     Remove host reads (bwa-mem2 or Kraken2)
  assembly         De novo assembly (SPAdes + MEGAHIT + cd-hit-est)
  assembly-refine  Iterative assembly refinement (optional; after quality)
  viral-id         Viral identification (geNomad)
  quality          Genome quality and selection (CheckV)
  annotate         Structural annotation (Pharokka → Phold → Phynteny)
  safety           Biosafety screening (CARD + VFDB)
  report           Generate final HTML report per candidate genome
  run              Execute all modules for all samples in a directory
  check-tools      Verify all required tools are installed
  config           Copy default config.yaml template
  init             Initialise a new project directory
  status           Show pipeline status for all samples in a project

Options (all commands):
  -c, --config     Path to config.yaml [default: config/config.yaml]
  -w, --workdir    Pipeline working directory
  -o, --output-dir Base output directory (results/ and reports/ created inside)
  --reports-dir    Override reports directory independently of -o
  -t, --threads    CPU threads [default: auto 90% of logical CPUs]
  --force          Force re-run even if outputs already exist
  --project        Project directory (from phageflow init); overrides --config/--workdir
```

---

## Scientific Basis

Methodological decisions are documented inline in each module. Key decisions:

**Why SPAdes `--isolate`?**
At >200× coverage, BayesHammer error correction collapses real SNPs and introduces artefacts in the assembly graph. The `--isolate` preset disables correction entirely and is optimised for near-complete genomes (Prjibelski et al. 2020).

**Why retain singletons?**
Reads spanning the junction of a circular genome map to both ends of a linearised assembly. bwa-mem2 in paired mode discards these as unmapped mates; PhageFlow retains them as singletons and passes them to SPAdes via `--s1`. This improves terminal coverage and is critical for CheckV DTR/ITR detection and Complete classification (Nayfach et al. 2021).

**Why CheckV over geNomad for topology?**
geNomad predicts topology from k-mer composition; CheckV detects DTR/ITR in the actual assembled sequence with a repeat-finding algorithm requiring ≥20 bp confirmed repeat. CheckV topology is therefore more reliable for the assembly quality determination. PhageFlow uses CheckV as the authoritative source and geNomad as fallback.

**Why three annotation tiers?**
Pharokka (MMseqs2 + PyHMMER) annotates ~30–50% of CDS on first pass. Phold structure-based annotation upgrades an additional ~10–20% of hypothetical proteins by finding structural homologs regardless of sequence divergence. Phynteny adds synteny context to predict function for remaining unknowns. Each tier is tracked independently via per-CDS delta reports.

---

## References

- Bankevich A et al. (2012) SPAdes: a new genome assembly algorithm. *J Comput Biol* 19:455
- Bouras G et al. (2023) Pharokka: a fast scalable bacteriophage annotation tool. *Bioinformatics* 39:btac776
- Bouras G et al. (2025) Phold: structure-based functional annotation of phage proteins. *bioRxiv* 2025.08.05.668817
- Bolger AM et al. (2014) Trimmomatic: a flexible trimmer for Illumina sequence data. *Bioinformatics* 30:2114
- Camargo AP et al. (2023) Identification of mobile genetic elements with geNomad. *Nat Biotechnol* 41:1783
- Chen S et al. (2018) fastp: an ultra-fast all-in-one FASTQ preprocessor. *Genome Biology* 19:274
- Fu L et al. (2012) CD-HIT: accelerated for clustering the next-generation sequencing data. *Bioinformatics* 28:3150
- Grigson SR et al. (2025) Phynteny: a synteny-based approach to phage annotation. *bioRxiv* 2025.07.28.667340
- Head SR et al. (2014) Library construction for next-generation sequencing: overviews and challenges. *Biotechniques* 56:61
- Li D et al. (2015) MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly. *Bioinformatics* 31:1674
- Nayfach S et al. (2021) CheckV assesses the quality and completeness of metagenome-assembled viral genomes. *Nat Biotechnol* 39:578
- Prjibelski A et al. (2020) Using SPAdes de novo assembler. *Curr Protoc Bioinf* 70:e102
- Roux S et al. (2019) Minimum information about an uncultivated virus genome (MIUViG). *Nat Biotechnol* 37:29
- Turner D et al. (2021) Abolishment of morphology-based taxa and change to binomial species names: 2022 taxonomy update of the ICTV bacterial viruses subcommittee. *Arch Virol* 166:2633
- Wick R & Holt K (2022) Benchmarking of long-read assemblers for prokaryote whole genome sequencing. *Microb Genomics* 8:mgen000788

---

## License

MIT © 2025 Cabezas-Mera et al. — see [LICENSE](LICENSE).

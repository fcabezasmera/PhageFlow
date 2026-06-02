# PhageFlow

**Modular bacteriophage genomics pipeline for complete and high-quality viral genome recovery from Illumina paired-end sequencing.**

[![Python 3.11](https://img.shields.io/badge/python-3.11-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-1.0.0-green.svg)](CHANGELOG.md)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://anaconda.org/bioconda/phageflow)

---

## Authors

| Name | Affiliation |
|------|-------------|
| Fausto Cabezas-Mera | Universidad Tecnológica Metropolitana, Santiago, Chile |
| Estefania Tisalema-Guanopatin | Universidad Tecnológica Metropolitana, Santiago, Chile |
| Dayra Valle | Universidad Internacional SEK, Quito, Ecuador |
| Antonella Nole | Universidad Internacional SEK, Quito, Ecuador |
| Katty Coral Carrillo | Universidad Internacional SEK, Quito, Ecuador |

---

## Table of Contents

1. [Overview](#overview)
2. [Pipeline Architecture](#pipeline-architecture)
3. [Requirements](#requirements)
4. [Installation](#installation)
5. [Database Setup](#database-setup)
6. [Quick Start](#quick-start)
7. [Module Reference](#module-reference)
8. [Configuration Reference](#configuration-reference)
9. [Output Structure](#output-structure)
10. [CLI Reference](#cli-reference)
11. [Scientific Basis](#scientific-basis)
12. [References](#references)
13. [License](#license)

---

## Overview

PhageFlow recovers complete and high-quality (HQ) phage genomes from Illumina paired-end sequencing. It supports **purified phage preparations**, **mixed virome samples**, and **environmental metagenomes**. No prior knowledge of sample composition is required — PhageFlow makes no assumptions about the non-host fraction until viral identification (Module 04).

**Design principles:**

- **Composition-agnostic**: modules 01–03 do not assume the input is viral
- **Modularity**: each module is independently re-runnable (`--force`) or resumable (`--from-module`)
- **Reproducibility**: all parameters controlled by a single `config/config.yaml`
- **Evidence-based**: all methodological decisions documented with primary literature
- **Read-length aware**: auto-detects PE150/PE250/PE300 and adjusts QC and assembly parameters

**Supported sample types:**

| Sample type | Host removal mode | Notes |
|-------------|-------------------|-------|
| Purified phage (single host) | bwa-mem2 + accessions | Recommended; most accurate |
| Purified phage (multiple hosts) | bwa-mem2 + multiple accessions | Concatenated reference |
| Virome / enriched preparation | Kraken2 | No host genome required |
| Environmental metagenome | Kraken2 | unclassified + viral reads retained |
| No host removal needed | pass-through | Reads forwarded directly to assembly |

---

## Pipeline Architecture

```
Raw reads (PE150 / PE250 / PE300 Illumina)
    │
    ▼  phageflow qc
01  Quality control + trimming
    fastp — adapter trim · PE correction · poly-X · complexity filter
            auto-detects read length; adjusts length_required and overlap
    FastQC — per-read QC reports
    MultiQC — aggregated QC report
    │
    ▼  phageflow host-removal
02  Host read removal
    bwa-mem2 (recommended) — permissive alignment (-A 1 -B 2 -O 2,2)
    ↳ singletons retained for DTR/ITR boundary coverage → SPAdes --s1
    ↳ Level A contamination check (Kraken2 diagnostic on 50k subsample)
    Kraken2 (alternative) — classification-based; unclassified + viral retained
    │
    ▼  phageflow assembly
03  De novo assembly
    SPAdes standard mode — paired reads + singletons (--s1)
    MEGAHIT — complementary assembler (different graph algorithm)
    cd-hit-est — NR reduction (100% identity, 85% coverage)
    k-mer range auto-adjusted: PE150 → k≤127; PE250 → k≤241; PE300 → k≤281
    │
    ▼  phageflow coverage
03b Coverage profiling (CoverM)
    mean · trimmed_mean (5-95%) · covered_bases · variance
    Warns: low coverage (<5x), ultra-high (>1000x), high CV (>2.0)
    │
    ▼  phageflow viral-id
04  Viral identification (geNomad)
    Neural network + MMseqs2 marker gene search
    ↳ topology inference (DTR / ITR / No terminal repeats)
    ↳ genetic code detection (standard=11, CrAss-like=15)
    ↳ taxonomy (lenient — resolves below family rank)
    Rescue tier: borderline contigs (score 0.4–0.7) ≥ 3 kb retained
    │
    ▼  phageflow quality
05  Genome quality and selection (CheckV)
    Quality tiers: Complete → HQ (≥90%) → MQ (≥50%) → large-ND (≥20 kb) → LQ-draft
    Dereplication: blastn (large ≥20 kb) + minimap2 (small <20 kb) at 95% ANI
    ↳ Circular rotation detection: cumulative multi-hit blastn coverage
    Co-binning: LQ drafts by taxonomy → bin-rescue if combined ≥30 kb
    bwa-mem2 read recruitment: depth · breadth · CV per candidate
    │
    ▼  phageflow annotate
06  Structural annotation (two-tier cascade)
    Tier 1  Pharokka  — PHANOTATE (single) / prodigal-gv (meta)
                        MMseqs2 + PyHMMER vs PHROGs
                        --coding_table from genome header (CrAss-like → 15)
                        --dnaapler reorientation (single-contig only)
    Tier 2  Phold     — ProstT5 3Di tokens + Foldseek structural search
                        --hyps: upgrades hypothetical proteins only
                        --finetune: phage-finetuned ProstT5 model
    Plot    phold plot (circular, single-contig) / pyGenomeViz (linear, multi-contig)
    │
    ▼  phageflow resistance
07  Resistance + biosafety screening
    AMR       CARD (Pharokka seq + Phold struct, dual-evidence)
    Virulence VFDB (Pharokka seq + Phold struct, dual-evidence)
    ACR       Anti-CRISPR proteins (Phold ACR database)
    Defense   DefenseFinder systems (Phold)
    TA        Toxin-antitoxin (Phold NetFlax)
    biosafety_flag: YES if any AMR or virulence hits
```

---

## Requirements

### Hardware

| Resource | Minimum | Recommended |
|----------|---------|-------------|
| CPU | 8 cores | 24+ cores |
| RAM | 32 GB | 64 GB |
| GPU | — | NVIDIA GPU (Phold ProstT5 + Foldseek) |
| Storage | 100 GB | 250 GB (databases ~80 GB) |

### Software

- Linux (Ubuntu 20.04+)
- conda / mamba ≥ 23.x
- CUDA ≥ 11.8 (optional, for `phold --foldseek_gpu`)

---

## Installation

### Recommended: Bioconda

PhageFlow and all its dependencies install in one command from Bioconda:

```bash
# Into a fresh environment (recommended)
conda create -n phageflow -c conda-forge -c bioconda phageflow
conda activate phageflow

# Or with mamba (faster solver)
mamba create -n phageflow -c conda-forge -c bioconda phageflow
conda activate phageflow

# Verify
phageflow --version
```

This pulls in every tool the pipeline needs (fastp, bwa-mem2, SPAdes, MEGAHIT,
geNomad, CheckV, Pharokka, Phold, and the rest) at compatible versions.

### Alternative: from source (development)

```bash
git clone https://github.com/fcabezasmera/PhageFlow.git
cd PhageFlow

# Create environment with all tools
mamba env create -f environment.yml
conda activate phageflow

# Install PhageFlow in editable mode
pip install -e .

# Verify
phageflow --version
```

> **Note on dependencies**: RGI and bacphlip are incompatible with the phageflow
> environment (Python 3.11 + samtools ≥1.23 conflict). PhageFlow uses Phold's
> built-in CARD and VFDB databases instead.

---

## Database Setup

PhageFlow needs several reference databases. The easiest way to install them is the
built-in `download-databases` command, which fetches each one from its official source
and writes the resolved paths to `~/.config/phageflow/config.yaml` so the pipeline finds
them automatically.

### One-command setup (recommended)

```bash
# Download everything (~47 GB total) — CheckV, geNomad, Pharokka, Phold, Kraken2
phageflow download-databases --all -o ~/phageflow_databases
```

This writes `~/.config/phageflow/config.yaml` with the database locations, so subsequent
`phageflow run` calls require no extra configuration.

### Individual databases

Download only what you need:

```bash
phageflow download-databases --checkv --genomad -o ~/phageflow_databases
phageflow download-databases --pharokka --phold -o ~/phageflow_databases
phageflow download-databases --kraken2 -o ~/phageflow_databases
```

| Flag | Database | Size | Source |
|------|----------|------|--------|
| `--checkv` | CheckV | ~2 GB | CheckV portal |
| `--genomad` | geNomad | ~7 GB | geNomad release |
| `--pharokka` | Pharokka | ~2 GB | Pharokka databases |
| `--phold` | Phold | ~18 GB | Phold install |
| `--kraken2` | Kraken2 standard-16GB | ~16 GB | AWS genome-idx |
| `--all` | All of the above | ~47 GB | — |

### Configuration priority

Database paths are resolved in this order (highest first):

1. CLI flag: `phageflow -c /path/config.yaml`
2. Project config: `./config/config.yaml`
3. User config: `~/.config/phageflow/config.yaml` (written by `download-databases`)
4. Environment variable: `PHAGEFLOW_DB=/path/to/databases`
5. Bundled defaults

### Manual setup (alternative)

If you prefer to manage databases yourself, install each one and set its path in
`config/config.yaml` under `databases:`:

```bash
checkv download_database  databases/checkv_db
genomad download-database databases/genomad_db
install_databases.py -o   databases/pharokka_db
phold install -d          databases/phold_db
```

For Kraken2, download a pre-built index from
[https://benlangmead.github.io/aws-indexes/k2](https://benlangmead.github.io/aws-indexes/k2)
and set `databases.kraken2` in config.yaml.

---

## Quick Start

### Set up databases (first time only)

```bash
phageflow download-databases --all -o ~/phageflow_databases
```

### Run the full pipeline

```bash
# Purified phage — bwa-mem2 with NCBI accession (recommended)
phageflow run \
    --r1 raw/sample_R1.fastq.gz \
    --r2 raw/sample_R2.fastq.gz \
    --accessions GCF_000005845.2 \
    -o results/

# Multiple host accessions
phageflow run \
    --r1 raw/sample_R1.fastq.gz \
    --r2 raw/sample_R2.fastq.gz \
    --accessions GCF_000005845.2,GCF_000013425.1 \
    -o results/

# Virome / metagenome — Kraken2 mode (no accession needed)
phageflow run \
    --r1 raw/sample_R1.fastq.gz \
    --r2 raw/sample_R2.fastq.gz \
    -o results/

# Resume from a specific module
phageflow run \
    --r1 raw/sample_R1.fastq.gz \
    --r2 raw/sample_R2.fastq.gz \
    --accessions GCF_000005845.2 \
    -o results/ \
    --from-module quality
```

### Run individual modules

```bash
phageflow qc \
    --r1 raw/sample_R1.fastq.gz --r2 raw/sample_R2.fastq.gz \
    -o results/

phageflow host-removal \
    --r1 results/01_qc/sample_R1.fastq.gz \
    --r2 results/01_qc/sample_R2.fastq.gz \
    --accessions GCF_000005845.2 \
    -o results/

phageflow assembly \
    --r1 results/02_host_removal/sample_R1.fastq.gz \
    --r2 results/02_host_removal/sample_R2.fastq.gz \
    --s1 results/02_host_removal/sample_singletons.fastq.gz \
    -o results/

phageflow coverage \
    --r1 results/02_host_removal/sample_R1.fastq.gz \
    --r2 results/02_host_removal/sample_R2.fastq.gz \
    --contigs results/03_assembly/sample_contigs.fasta \
    -o results/

phageflow viral-id \
    --contigs results/03_assembly/sample_contigs.fasta \
    -o results/

phageflow quality \
    --virus-fna results/04_viral_id/sample_virus.fna \
    -o results/

phageflow annotate \
    --genome results/05_quality/sample/annotation_ready/phages/Ackermannviridae_candidate_001.fasta \
    -o results/

phageflow resistance \
    --candidate-id Ackermannviridae_candidate_001 \
    --sample-id sample \
    -o results/
```

---

## Module Reference

### Module 01: Quality Control

**Command**: `phageflow qc` | **Tools**: fastp · FastQC · MultiQC

Adapter trimming, quality filtering, PE overlap correction, poly-X removal, and low-complexity filtering. Automatically detects read length (PE150/250/300) and adjusts parameters accordingly.

#### Read-length adaptive parameters

| Parameter | PE150 | PE250 | PE300 | Rationale |
|-----------|-------|-------|-------|-----------|
| `length_required` | 50 bp | 75 bp | 100 bp | ~1/3 read length minimum post-trim |
| `overlap_len_require` | 15 bp | 30 bp | 50 bp | Proportional to expected overlap |
| `cut_right_window_size` | 4 bp | 6 bp | 8 bp | Smooths quality variation in longer reads |

#### Fixed parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `average_qual` | Q20 | MIUViG minimum (Roux et al. 2019) |
| `qualified_quality_phred` | Q25 | 99.7% per-base accuracy |
| `unqualified_percent_limit` | 15% | Lenient to retain divergent-phage reads |
| `low_complexity_filter` | enabled | Removes homopolymers (Roux et al. 2019 MIUViG §2.1) |
| `complexity_threshold` | 20% | Conservative filter |
| `trim_poly_x` | enabled | NovaSeq/NextSeq 2-color poly-G artefacts |

> **Duplication note**: high duplication (>50%) is expected in purified phage preparations
> at high coverage. Reads are NOT deduplicated (Head et al. 2014, Biotechniques 56:61).

---

### Module 02: Host Read Removal

**Command**: `phageflow host-removal` | **Tools**: bwa-mem2 · samtools ≥1.15 OR Kraken2 · seqtk

#### bwa-mem2 mode (recommended)

Permissive alignment parameters for divergent host strain detection:

```
bwa-mem2 mem -A 1 -B 2 -O 2,2 -M | samtools view -bF 2304 | samtools sort -n |
samtools fastq -f 4 -F 256 -1 R1_out -2 R2_out -s singletons
```

`-A 1 -B 2 -O 2,2`: lower mismatch/gap penalties capture divergent host strains
that default parameters would miss (Schmieder & Edwards 2011, Bioinformatics 27:863).

Singletons (one mate maps to host, one does not) are retained → passed to SPAdes `--s1`.
These bridge terminal repeat junctions (DTR/ITR) in circular viral genomes.

Host references: NCBI accession (`--accessions`), local FASTA (`--host-file`), or both.

#### Kraken2 mode (alternative)

Used when no host genome is available. Retains unclassified + viral reads (taxid 10239).
Parameters: `--confidence 0.5 --minimum-hit-groups 3` (strict, prevents false-positive
host classification of divergent viral sequences; Wood et al. 2019, Genome Biol 20:257).

#### Level A contamination check

Post-bwa-mem2 diagnostic: 50,000 phage reads classified by Kraken2.
Warns if bacterial fraction > `contamination_warn_pct` (default 5%). No reads removed.

---

### Module 03: De Novo Assembly

**Command**: `phageflow assembly` | **Tools**: SPAdes · MEGAHIT · cd-hit-est

Dual-assembler strategy. SPAdes and MEGAHIT use different graph algorithms and produce
complementary contig sets. Their union is deduplicated with cd-hit-est.

#### SPAdes

Standard mode (no `--isolate`, no `--meta`) — balanced for both pure and mixed samples.
`--isolate` is not used because it discards low-coverage contigs, which would miss
minority phages in mixed preparations.
Singletons passed via `--s1` (DTR/ITR boundary recovery).

#### MEGAHIT

`--min-count 2` filters sequencing errors. `--no-mercy` is NOT used — mercy k-mers
retain divergent viral sequences at low coverage (Li et al. 2015).

#### k-mer range auto-selection

| Read length | SPAdes k-max | MEGAHIT k-max |
|-------------|-------------|---------------|
| PE150 | 127 (v4.x limit) | 127 |
| PE250 | 127 | 241 |
| PE300 | 127 | 281 |

SPAdes v4.x hard limit: k ≤ 127. MEGAHIT supports k up to 255.
k_max < read_length (Bankevich et al. 2012; Li et al. 2015).

---

### Module 03b: Coverage Profiling

**Command**: `phageflow coverage` | **Tools**: CoverM

Per-contig coverage metrics before viral identification. Diagnostic only — no contigs removed.

| Metric | Description |
|--------|-------------|
| `mean` | Average depth across all positions |
| `trimmed_mean` | 5–95% trimmed mean (robust to terminal repeat depth artefacts) |
| `covered_bases` | Positions covered at ≥1x |
| `variance` | Depth variance (CV = std/mean) |

| Threshold | Value | Action |
|-----------|-------|--------|
| Low coverage | mean < 5x | WARN — possible assembly artefact |
| Ultra-high coverage | mean > 1000x | WARN — possible concatemer |
| Uneven depth | CV > 2.0 AND mean > 10x | WARN — possible chimera |

---

### Module 04: Viral Identification

**Command**: `phageflow viral-id` | **Tools**: geNomad

#### Score thresholds

| Tier | Score | Length | Action |
|------|-------|--------|--------|
| main | ≥ 0.7 | any | Viral (97% precision at this threshold) |
| rescued | 0.4–0.7 | ≥ 3 kb | Novel lineage rescue (captures Microviridae/Inoviridae) |
| discarded | < 0.4 | — | Not viral |

Rescue threshold reduced to 3 kb (from 10 kb) to capture complete Microviridae
(3–6 kb) and Inoviridae (6–9 kb) genomes (Roux et al. 2019 MIUViG).

#### geNomad parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `composition` | auto | Auto-detects sample composition (Camargo et al. 2023) |
| `sensitivity` | 4.2 | Maximum MMseqs2 sensitivity for hallmark detection |
| `lenient_taxonomy` | true | Resolves classification below ICTV species boundary |
| `disable_find_proviruses` | true | Post-host-removal contigs: provirus detection → false positives |
| `min_score` | 0.7 | ~97% precision (Camargo et al. 2023, Suppl Fig 2) |

---

### Module 05: Genome Quality Assessment

**Command**: `phageflow quality` | **Tools**: CheckV · blastn · minimap2

#### Quality tiers

| Tier | Criterion | Output |
|------|-----------|--------|
| Complete | CheckV DTR/ITR detected | annotation_ready/ |
| High-quality | completeness ≥ 90% | annotation_ready/ |
| Medium-quality | completeness ≥ 50% (MIUViG standard) | annotation_ready/ |
| large-ND | ND, length ≥ 20 kb, ≥1 viral gene | annotation_ready/ |
| lq-draft | length ≥ 10 kb, ≥1 viral gene | drafts/ or bin-rescue |
| bin-rescue | LQ drafts, same taxon, combined ≥ 30 kb | annotation_ready/ |
| discarded | length < 1.5 kb or no viral genes | — |

#### Dereplication strategy

- **Large genomes (≥ 20 kb)**: blastn all-vs-all (ANI ≥ 95%, coverage ≥ 80%).
  Detects non-overlapping fragments of the same genome that mash misses
  (k-mer content is disjoint between non-overlapping regions).
- **Circular rotation detection**: cumulative multi-hit blastn coverage ≥ 90% at
  pident ≥ 95%. Catches SPAdes/MEGAHIT assemblies of the same circular genome
  with different start positions.
- **Small genomes (< 20 kb)**: minimap2 asm5 all-vs-all (identity ≥ 95%, coverage ≥ 80%).
  mash is unreliable for small genomes (sparse sketch).

Threshold: 95% ANI = ICTV species boundary (Turner et al. 2021, Arch Virol 166:2633).

#### Topology resolution

CheckV termini_type takes precedence over geNomad topology:
CheckV confirms DTR/ITR in the actual assembled sequence (≥20 bp repeat), whereas
geNomad predicts topology from sequence composition.

---

### Module 06: Structural Annotation

**Command**: `phageflow annotate` | **Tools**: Pharokka · Phold · dnaapler

#### Two-tier cascade

**Tier 1 — Pharokka** (sequence-based):
- Gene calling: PHANOTATE (single-contig) or prodigal-gv (multi-contig)
- Database search: MMseqs2 + PyHMMER vs PHROGs
- Reorientation: dnaapler (single-contig only — finds terL, sets canonical origin)
- Genetic code: read from genome FASTA header `[genetic_code=N]` (set by Module 04)
  CrAss-like phages use code 15 (TGA=Trp); passed as `--coding_table` to Pharokka

**Tier 2 — Phold** (structure-based):
- `--hyps`: upgrades hypothetical proteins from Pharokka only (conservative)
- ProstT5 generates 3Di structural tokens from amino acid sequence
- Foldseek searches PHROGs-3D, CARD, VFDB, DefenseFinder, NetFlax
- `--finetune`: phage-finetuned ProstT5 (better 3Di quality for phage proteins)

**Why `--hyps` outperforms full re-annotation:**
Phold full mode can downgrade high-confidence Pharokka MMseqs2 hits when structural
evidence is ambiguous. `--hyps` preserves confident sequence-based annotations and
applies structure-based annotation only where Pharokka found no evidence.

#### Visualisation

| Genome type | Plot tool | Output |
|-------------|-----------|--------|
| Single-contig (circular/complete) | phold plot | Circular map, PNG + SVG |
| Multi-contig (fragmented) | pyGenomeViz | Linear multi-track, PNG + SVG |

---

### Module 07: Resistance and Biosafety Screening

**Command**: `phageflow resistance` | **Tools**: none (parses Pharokka + Phold outputs)

No additional tools required. Parses existing annotation outputs from Module 06.

#### Data sources

| Category | Source | Method |
|----------|--------|--------|
| AMR (CARD) | Pharokka merged TSV | MMseqs2 sequence search |
| AMR (CARD) | Phold sub_db_tophits | Foldseek structural search |
| Virulence (VFDB) | Pharokka merged TSV | MMseqs2 sequence search |
| Virulence (VFDB) | Phold sub_db_tophits | Foldseek structural search |
| Anti-CRISPR (ACR) | Phold sub_db_tophits | Foldseek (fident ≥ 20%) |
| Defense systems | Phold sub_db_tophits | Foldseek structural search |
| Toxin-antitoxin | Phold sub_db_tophits | Foldseek structural search |

#### Confidence levels

| Level | Criterion |
|-------|-----------|
| HIGH | Detected by BOTH Pharokka (sequence) AND Phold (structure) |
| MEDIUM | Pharokka sequence search only |
| LOW | Phold structure search only |

`biosafety_flag: YES` is raised only for AMR or virulence hits. ACR, defense, and
toxin-antitoxin are informational (relevant for phage-host interaction interpretation).

ACR-specific thresholds (fident ≥ 20%, qcov ≥ 35%): anti-CRISPR proteins are among
the most sequence-divergent phage proteins; structural homology at low sequence
identity is biologically meaningful (Pawluk et al. 2016, Science 351:aad8405).

---

## Configuration Reference

### QC (`qc:`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `average_qual` | 20 | Read-level mean Q floor (MIUViG minimum) |
| `qualified_quality_phred` | 25 | Per-base Q threshold |
| `unqualified_percent_limit` | 15 | Max % low-Q bases per read |
| `length_required` | 50 | Minimum post-trim length (PE150 default; auto-adjusted) |
| `cut_right_window_size` | 4 | Sliding-window size (auto-adjusted by read length) |
| `cut_right_mean_quality` | 25 | Sliding-window quality floor |
| `correction` | true | PE overlap correction |
| `overlap_len_require` | 15 | Min overlap for PE correction (auto-adjusted) |
| `low_complexity_filter` | true | Homopolymer filter |
| `complexity_threshold` | 20 | Minimum complexity (%) |
| `trim_poly_x` | true | Poly-X tail removal |
| `poly_x_min_len` | 10 | Min poly-X length to trigger |
| `n_base_limit` | 5 | Max N calls per read |

### Assembly (`assembly:`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `kmers` | `"21,33,55,77,99,127"` | Base k-mer list (SPAdes max=127; MEGAHIT extended automatically) |
| `min_length` | 200 | Minimum contig length (bp) for NR pool |

### Viral Identification (`genomad:`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_score` | 0.7 | Primary viral score threshold (~97% precision) |
| `rescue_min_score` | 0.4 | Lower threshold for rescue tier |
| `rescue_min_length_bp` | 3000 | Minimum length for rescue (captures Microviridae) |
| `sensitivity` | 4.2 | MMseqs2 sensitivity (maximum) |
| `enable_score_calibration` | true | Score calibration by sample composition |
| `composition` | auto | Auto-detects composition (Camargo et al. 2023) |
| `lenient_taxonomy` | true | Resolve below family rank |
| `disable_find_proviruses` | true | Skip provirus detection on post-host-removal contigs |
| `min_virus_hallmarks` | 0 | Primary filter is score + rescue_length |

### Quality Assessment (`checkv:`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `min_completeness` | 50 | MQ threshold % (MIUViG standard) |
| `min_viral_genes` | 1 | Minimum viral genes for draft rescue |
| `length_rescue` | 10000 | Minimum length for LQ draft (bp) |
| `min_contig_bp` | 1500 | Global minimum contig length (bp) |
| `large_nd_rescue_bp` | 20000 | Length threshold for large-ND rescue (bp) |
| `min_bin_rescue_bp` | 30000 | Combined length for bin rescue (bp) |
| `max_kmer_freq` | 1.5 | Concatemer warning threshold |
| `max_genome_copies` | 1.5 | Concatemer warning threshold |

### Annotation (`annotate:`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `phold_gpu` | true | Use `--foldseek_gpu` (set false for CPU-only) |
| `phold_finetune` | true | Phage-finetuned ProstT5 model |
| `phold_batch_size` | 1 | ProstT5 batch size (increase on GPU: 4–8) |
| `phold_sensitivity` | 9.5 | Foldseek sensitivity |
| `genetic_code` | 11 | Default genetic code (overridden per-genome from geNomad) |

---

## Output Structure

```
<output_dir>/
├── results/
│   ├── 01_qc/
│   │   ├── {sample}_R1.fastq.gz
│   │   └── {sample}_R2.fastq.gz
│   ├── 02_host_removal/
│   │   ├── {sample}_R1.fastq.gz
│   │   ├── {sample}_R2.fastq.gz
│   │   └── {sample}_singletons.fastq.gz
│   ├── 03_assembly/
│   │   ├── {sample}_spades/contigs.fasta
│   │   ├── {sample}_megahit/{sample}.contigs.fa
│   │   └── {sample}_contigs.fasta          ← NR combined contigs
│   ├── 03b_coverage/{sample}/
│   │   └── {sample}_coverage.tsv           ← per-contig metrics
│   ├── 04_viral_id/
│   │   ├── {sample}_virus.fna              ← viral contigs
│   │   └── {sample}_metadata.tsv           ← topology + genetic_code
│   ├── 05_quality/{sample}/
│   │   ├── annotation_ready/phages/
│   │   │   └── {naming_level}_candidate_NNN.fasta
│   │   ├── drafts/
│   │   └── checkv/
│   ├── 06_annotation/{candidate_id}/
│   │   ├── pharokka/                       ← Tier 1 (GBK + TSVs)
│   │   ├── phold/                          ← Tier 2 (canonical GBK)
│   │   └── plots/                          ← genome map (PNG + SVG)
│   └── 07_resistance/{candidate_id}/
│       ├── {candidate_id}_amr.tsv
│       ├── {candidate_id}_virulence.tsv
│       ├── {candidate_id}_acr.tsv
│       ├── {candidate_id}_defense.tsv
│       └── {candidate_id}_toxin_antitoxin.tsv
└── reports/
    ├── 01_qc/{sample}/qc_summary.tsv
    ├── 02_host_removal/{sample}/host_removal_summary.tsv
    ├── 03_assembly/{sample}/assembly_summary.tsv
    ├── 03b_coverage/{sample}/coverage_summary.tsv
    ├── 04_viral_id/{sample}/viral_id_summary.tsv
    ├── 05_quality/{sample}/
    │   └── rename_map.tsv                  ← critical metadata link
    ├── 06_annotation/
    │   ├── annotation_summary.tsv          ← aggregate (all candidates)
    │   └── {candidate_id}/
    │       └── {candidate_id}_delta_report.tsv
    └── 07_resistance/
        ├── resistance_aggregate.tsv        ← all candidates
        └── {candidate_id}/{candidate_id}_resistance_summary.tsv
```

### Key files for downstream use

| File | Description |
|------|-------------|
| `results/05_quality/{sample}/annotation_ready/phages/*.fasta` | Annotation-ready genomes |
| `reports/05_quality/{sample}/rename_map.tsv` | Full metadata per contig |
| `results/06_annotation/{candidate}/phold/{candidate}.gbk` | Canonical annotated GBK |
| `reports/06_annotation/annotation_summary.tsv` | Annotation statistics |
| `reports/06_annotation/{candidate}/{candidate}_delta_report.tsv` | Per-CDS provenance |
| `reports/07_resistance/resistance_aggregate.tsv` | Biosafety summary all candidates |

---

## CLI Reference

```
phageflow [OPTIONS] COMMAND [ARGS]

Commands:
  qc            Quality control and trimming (fastp + FastQC + MultiQC)
  host-removal  Remove host reads (bwa-mem2 or Kraken2)
  assembly      De novo assembly (SPAdes + MEGAHIT + cd-hit-est)
  coverage      Coverage profiling (CoverM)
  viral-id      Viral identification (geNomad)
  quality       Genome quality and selection (CheckV)
  annotate      Structural annotation (Pharokka → Phold)
  resistance    AMR + virulence + ACR + defense screening
  run           Execute full pipeline for one or more samples
  download-databases  Download all reference databases + write user config
  config        Copy default config.yaml template
  init          Initialise a new project directory

Options (all commands):
  -c, --config       Path to config.yaml [default: config/config.yaml]
  -w, --workdir      Pipeline working directory
  -o, --output-dir   Base output directory
  --reports-dir      Override reports directory
  -t, --threads      CPU threads [default: auto 90% of logical CPUs]
  --force            Force re-run even if outputs exist

Options (run only):
  --r1               R1 FASTQ path (repeatable for multiple samples)
  --r2               R2 FASTQ path (repeatable for multiple samples)
  --raw-dir          Directory to scan for paired FASTQ files
  --accessions       NCBI genome accession(s) for host removal (comma-separated)
  --host-file        Local host FASTA or folder of FASTAs
  --from-module      Resume from module: qc|host-removal|assembly|coverage|
                     viral-id|quality|annotate|resistance
```

---

## Scientific Basis

**No composition assumptions before Module 04**
Modules 01–03 treat reads as unknown composition. The non-host fraction may be viral,
plasmid, bacterial (residual), or chimeric. Viral identification is performed after
assembly, not before — this ensures that even highly divergent phages with no database
reference are assembled and passed to geNomad for classification.

**bwa-mem2 permissive parameters for host removal**
Default BWA parameters are optimised for variant calling (high specificity). Host removal
requires high sensitivity — divergent host strains that default parameters miss will
appear as false non-host signal. Parameters `-A 1 -B 2 -O 2,2` lower the mismatch/gap
penalties to capture divergent host sequences (Schmieder & Edwards 2011).

**SPAdes standard mode over --isolate**
`--isolate` is optimised for single high-coverage isolates and discards low-coverage
contigs. Mixed phage preparations or metagenomes may contain multiple viral genomes at
variable coverage. Standard mode preserves low-coverage contigs and performs comparably
on dominant high-coverage genomes (Bankevich et al. 2012).

**blastn for dereplication of large viral genomes**
mash k-mer sketches compare k-mer content between sequences. Non-overlapping fragments
of the same circular genome have disjoint k-mer content → mash reports high distance
even at 100% identity. blastn all-vs-all with cumulative multi-hit coverage correctly
identifies these as the same genome.

**Circular rotation detection**
SPAdes and MEGAHIT may assemble the same circular genome with different start positions.
Both assemblies are 100% identical but produce two BLAST hits (each covering ~50% of
the query) rather than one hit at 100% coverage. PhageFlow detects this by summing
all aligned bases between a pair and checking if cumulative coverage ≥ 90% at ≥ 95%
identity.

**--hyps strategy in Phold**
Phold full re-annotation can downgrade high-confidence Pharokka MMseqs2 hits when
structural evidence is ambiguous. --hyps preserves all confident sequence-based
annotations and applies structure-based annotation only to hypothetical proteins.
Empirically: --hyps achieves +25 upgrades vs +22 for full re-annotation on a
157 kb Ackermannviridae genome.

**Dual-evidence AMR/virulence screening**
AMR gene detection in phages requires both sensitivity (divergent sequences missed by
sequence search) and specificity (avoiding false structural homologs). PhageFlow
combines Pharokka MMseqs2 (sequence) with Phold Foldseek (structure) and reports
confidence levels (HIGH = both, MEDIUM = sequence only, LOW = structure only).

---

## References

- Al-Shayeb B et al. (2020) Clades of huge phages from across Earth's ecosystems. *Nature* 578:425
- Antipov D et al. (2016) plasmidSPAdes: assembling plasmids from whole genome sequencing data. *Bioinformatics* 32:i60
- Bankevich A et al. (2012) SPAdes: a new genome assembly algorithm. *J Comput Biol* 19:455
- Bondy-Denomy J et al. (2013) Bacteriophage genes that inactivate the CRISPR/Cas bacterial immune system. *Nature* 493:429
- Bouras G et al. (2023) Pharokka: a fast scalable bacteriophage annotation tool. *Bioinformatics* 39:btac776
- Bouras G et al. (2025) Phold: structure-based functional annotation of phage proteins. *bioRxiv* 2025.08.05.668817
- Camargo AP et al. (2023) Identification of mobile genetic elements with geNomad. *Nat Biotechnol* 41:1783
- Chen S et al. (2018) fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics* 34:i884
- Colomer-Lluch M et al. (2011) Bacteriophages carrying antibiotic resistance genes in human feces. *Antimicrob Agents Chemother* 55:4247
- Fu L et al. (2012) CD-HIT: accelerated for clustering next-generation sequencing data. *Bioinformatics* 28:3150
- Head SR et al. (2014) Library construction for next-generation sequencing. *Biotechniques* 56:61
- Heinzinger M et al. (2023) Bilingual language model for protein sequence and structure. *Bioinformatics* 39:btad436
- Holtgrewe M et al. (2013) Mason: a read simulator for second generation sequencing data. *PLoS ONE* 8:e61458
- Li D et al. (2015) MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly. *Bioinformatics* 31:1674
- McNair K et al. (2019) PHANOTATE: a novel approach to gene identification in phage genomes. *Bioinformatics* 35:4537
- Nayfach S et al. (2021) CheckV assesses the quality and completeness of metagenome-assembled viral genomes. *Nat Biotechnol* 39:578
- Pawluk A et al. (2016) Naturally occurring off-switches for CRISPR-Cas that silence microbial immune systems. *Science* 351:aad8405
- Prjibelski A et al. (2020) Using SPAdes de novo assembler. *Curr Protoc Bioinf* 70:e102
- Roux S et al. (2019) Minimum information about an uncultivated virus genome (MIUViG). *Nat Biotechnol* 37:29
- Schmieder R & Edwards R (2011) Fast identification and removal of sequence contamination from genomic and metagenomic datasets. *PLoS ONE* 6:e17288
- Shimoyama Y (2022) pyGenomeViz: a genome visualization python package for comparative genomics. *bioRxiv* 2022.11.24.517870
- Stanley SY & Maxwell KL (2018) Phage-encoded anti-CRISPR defenses. *Annu Rev Genet* 52:445
- Turner D et al. (2021) Abolishment of morphology-based taxa and change to binomial species names. *Arch Virol* 166:2633
- van Kempen M et al. (2024) Fast and accurate protein structure search with Foldseek. *Nat Biotechnol* 42:243
- Wick R et al. (2021) Assembling the perfect bacterial genome using Oxford Nanopore and Illumina sequencing. *Genome Biol* 22:241
- Wood DE et al. (2019) Improved metagenomic analysis with Kraken 2. *Genome Biol* 20:257
- Yutin N et al. (2018) Discovery of an expansive bacteriophage family that includes the most abundant viruses from the human gut. *Nat Microbiol* 3:1145

---

## License

MIT © 2025 Cabezas-Mera et al. — see [LICENSE](LICENSE).

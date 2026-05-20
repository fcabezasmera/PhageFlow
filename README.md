# PhageFlow

**Modular bacteriophage genomics pipeline for Illumina paired-end sequencing.**

PhageFlow processes raw reads from purified phage preparations through quality
control, host removal, assembly, viral identification, quality assessment,
annotation, biosafety screening, and lifecycle prediction — one sample at a
time, with full control at each step.

```
reads → QC → host removal → assembly → viral ID → quality →
annotation → biosafety → lifecycle → final genomes
```

> **Mode:** `purified_phage` — optimised for preparations where the dominant
> nucleic acid is phage DNA. Parameters and rescue thresholds reflect this
> assumption throughout.

---

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Configuration](#configuration)
- [Samples file](#samples-file)
- [Environment setup](#environment-setup)
- [Modules](#modules)
  - [qc](#qc)
  - [host-removal](#host-removal)
  - [assembly](#assembly)
  - [viral-id](#viral-id)
  - [quality](#quality)
  - [annotate](#annotate)
  - [safety](#safety)
  - [lifecycle](#lifecycle)
- [Output structure](#output-structure)
- [Running all samples in a loop](#running-all-samples-in-a-loop)
- [Databases](#databases)
- [Citation](#citation)
- [Authors](#authors)

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/fcabezasmera/PhageFlow.git
cd PhageFlow
```

### 2. Create the main conda environment

```bash
conda create -n phageflow python=3.11 -y
conda activate phageflow
pip install -e .
```

### 3. Install external tools

```bash
# Main environment
conda install -c bioconda -c conda-forge \
    fastp fastqc multiqc \
    bwa-mem2 samtools seqtk \
    spades megahit cd-hit \
    checkv pharokka seqkit \
    abricate \
    ncbi-datasets-cli kraken2 -y

# geNomad (separate environment)
conda create -n genomad -c conda-forge -c bioconda genomad -y

# Phold (separate environment)
conda create -n phold -c conda-forge -c bioconda phold -y

# BACPHLIP (separate environment)
conda create -n bacphlip -c conda-forge -c bioconda bacphlip -y
```

### 4. Register the PhageFlow CLI in every environment

Modules that require a separate environment still need access to the
`phageflow` command. Register it once per environment without touching
its dependencies:

```bash
conda activate genomad  && pip install -e . --no-deps -q
conda activate phold    && pip install -e . --no-deps -q
conda activate bacphlip && pip install -e . --no-deps -q
conda activate phageflow
```

### 5. Set up databases

```bash
# CheckV
checkv download_database databases/checkv_db

# Pharokka
install_databases.py -o databases/pharokka_db

# geNomad
conda activate genomad
genomad download-database databases/genomad_db

# Phold
conda activate phold
phold install -d databases/phold_db

# (Optional) Kraken2 — only for auto host detection
conda activate phageflow
kraken2-build --download-library bacteria --db databases/k2_db
kraken2-build --build --db databases/k2_db

conda activate phageflow
```

### 6. Verify installation

```bash
phageflow check-tools
```

---

## Quick Start

```bash
# Copy config templates
phageflow config   # creates config/config.yaml
phageflow samples  # creates config/samples.tsv

# Edit samples.tsv with your data, then run each module in order.
# Stay in the phageflow environment — use conda run for tool-specific envs.

conda activate phageflow

phageflow qc --sample-id s1 \
  --r1 data/s1_R1.fastq.gz --r2 data/s1_R2.fastq.gz

phageflow host-removal --sample-id s1 \
  --r1 results/01_qc/s1_R1.fastq.gz \
  --r2 results/01_qc/s1_R2.fastq.gz \
  --host-file /path/to/host.fasta

phageflow assembly --sample-id s1 \
  --r1 results/02_host_removal/s1_R1.fastq.gz \
  --r2 results/02_host_removal/s1_R2.fastq.gz

conda run -n genomad phageflow viral-id --sample-id s1 \
  --contigs results/03_assembly/combined/s1_contigs_nr.fasta

phageflow quality --sample-id s1 \
  --virus-fna results/04_viral_id/s1_virus.fna

# annotation-ready genomes are written per-candidate under:
# results/05_quality/annotation_ready/phages/<candidate>.fasta
# results/05_quality/annotation_ready/proviruses/<candidate>.fasta

phageflow annotate --sample-id s1 \
  --genome results/05_quality/annotation_ready/phages/<candidate>.fasta

phageflow safety --sample-id s1 \
  --genome results/05_quality/annotation_ready/phages/<candidate>.fasta

conda run -n bacphlip phageflow lifecycle --sample-id s1 \
  --genome results/05_quality/annotation_ready/phages/<candidate>.fasta
```

> **Tip:** you never need to leave the `phageflow` environment.
> Use `conda run -n <env>` for modules that require geNomad, Phold, or BACPHLIP.

---

## Configuration

Edit `config/config.yaml` before running. All paths are relative to the
pipeline working directory (project root by default).

```yaml
project: PhageFlow
version: "3.2.0"
mode: purified_phage

samples_file: config/samples.tsv

threads:   24
memory_gb: 32

dirs:
  databases: databases
  results:   results
  reports:   reports

databases:
  checkv:   databases/checkv_db/checkv-db-v1.5
  pharokka: databases/pharokka_db
  genomad:  databases/genomad_db
  phold:    databases/phold_db
  kraken2:  databases/k2_db    # only needed for Kraken2 host-removal mode

envs:
  main:      phageflow
  genomad:   genomad
  phold:     phold
  lifecycle: bacphlip

assembly:
  kmers:      "21,33,55,77,99,127"
  min_length: 200              # contigs < N bp discarded post-assembly

genomad:
  min_score:           0.7    # virus score threshold (Camargo 2023: 0.7 ≈ 97% precision)
  min_virus_hallmarks: 0      # 0 = score-only mode; increase to require hallmark genes

checkv:
  min_completeness:   50      # MQ threshold for annotation_ready/
  min_viral_genes:    1       # LQ rescue: min viral genes → drafts/
  length_rescue:      10000   # ND moderate rescue (≥10 kb + ≥1 gene) → drafts/
  min_contig_bp:      1500    # global filter: contigs < N bp discarded
  large_nd_rescue_bp: 30000   # ND large rescue (≥30 kb + ≥1 gene) → annotation_ready/
  min_bin_rescue_bp:  30000   # draft co-bin rescue: bin total ≥ N bp → annotation_ready/
  min_gene_density:   0.5     # viral genes/kb for density-based ND rescue → drafts/
```

Override threads at runtime:

```bash
phageflow qc --sample-id s1 --r1 R1.fq --r2 R2.fq -t 8
```

---

## Samples file

`config/samples.tsv` — tab-separated, one sample per line:

```tsv
sample_id	r1	r2
s1	/path/to/s1_R1.fastq.gz	/path/to/s1_R2.fastq.gz
s2	/path/to/s2_R1.fastq.gz	/path/to/s2_R2.fastq.gz
```

- Supports `.fastq` and `.fastq.gz`
- `sample_id` is the prefix for all output files
- Lines starting with `#` are ignored

---

## Environment setup

| Environment | Tools    | Used by modules         |
|-------------|----------|-------------------------|
| `phageflow` | all core tools | qc, host-removal, assembly, quality, annotate, safety |
| `genomad`   | geNomad  | viral-id                |
| `phold`     | Phold    | annotate (step 2/2)     |
| `bacphlip`  | BACPHLIP | lifecycle               |

Register the CLI once in each environment (done at install, not per run):

```bash
conda activate genomad  && pip install -e . --no-deps -q
conda activate phold    && pip install -e . --no-deps -q
conda activate bacphlip && pip install -e . --no-deps -q
conda activate phageflow
```

---

## Modules

### qc

**Quality control and trimming** using fastp + FastQC + MultiQC.

```bash
phageflow qc \
  --sample-id s1 \
  --r1 data/s1_R1.fastq.gz \
  --r2 data/s1_R2.fastq.gz
```

**fastp parameters (PE150, purified phage):**

| Parameter | Value | Rationale |
|---|---|---|
| Quality threshold | Q ≥ 20 | 99% base call accuracy |
| Minimum read length | 75 bp | reliable k-mer coverage for assembly |
| Sliding window | 4 bp / Q20 | 3′ quality trimming (Bolger et al. 2014) |
| Low-complexity filter | ≥ 30% | removes homopolymers (Roux et al. 2019) |
| N base limit | 5 | removes low-confidence reads |
| Adapter detection | auto | paired-end mode |

**Warnings emitted when:**
- Pass rate < 80% — possible contamination or low-quality run
- Q30 < 75% — reduced assembly contiguity expected
- Duplication > 30% — possible PCR over-amplification

**Outputs:**
```
results/01_qc/{sample}_R1.fastq.gz
results/01_qc/{sample}_R2.fastq.gz
reports/01_qc/{sample}_fastp.json
reports/01_qc/multiqc/multiqc_qc.html
reports/01_qc/qc_summary.tsv
```

---

### host-removal

**Host read removal** using bwa-mem2 alignment or Kraken2 classification.

```bash
# Mode 1 — local FASTA (single file, folder, or path list)
phageflow host-removal --sample-id s1 --r1 R1.fq.gz --r2 R2.fq.gz \
  --host-file /path/to/host.fasta

# Mode 2 — NCBI accessions (auto-download via datasets CLI)
phageflow host-removal --sample-id s1 --r1 R1.fq.gz --r2 R2.fq.gz \
  --accessions GCF_000013465.1,GCF_000007785.1

# Mode 3 — accessions from file (one per line)
phageflow host-removal --sample-id s1 --r1 R1.fq.gz --r2 R2.fq.gz \
  --accessions-file hosts.txt

# Mode 4 — Kraken2 auto-detection (no reference needed)
phageflow host-removal --sample-id s1 --r1 R1.fq.gz --r2 R2.fq.gz
```

Mode priority: `--host-file` > `--accessions` > `--accessions-file` > Kraken2.

**bwa-mem2 mode:** aligns all reads to the combined host reference; retains
unmapped pairs (flag `-f 12 -F 256`).

**Kraken2 mode:** classifies reads with `--confidence 0.2 --minimum-hit-groups 2`;
retains unclassified reads + Viruses (taxid 10239). Particularly useful when
no host reference is available or multiple unknown hosts are present.

**Warnings emitted when:**
- Phage fraction < 30% — likely wrong reference or heavy contamination
- Surviving reads < 10,000 — assembly failure risk
- Kraken2: < 30% unclassified — phage may already be in the database

**Outputs:**
```
results/02_host_removal/{sample}_R1.fastq.gz
results/02_host_removal/{sample}_R2.fastq.gz
reports/02_host_removal/host_removal_summary.tsv
```

---

### assembly

**De novo assembly** with metaSPAdes + MEGAHIT, followed by cd-hit-est
dereplication.

```bash
phageflow assembly \
  --sample-id s1 \
  --r1 results/02_host_removal/s1_R1.fastq.gz \
  --r2 results/02_host_removal/s1_R2.fastq.gz
```

| Tool | Key flags | Rationale |
|---|---|---|
| metaSPAdes | `--meta --only-assembler` | handles variable coverage; skips error correction that collapses tail fiber variants at >200× (Roux et al. 2019) |
| MEGAHIT | `--no-mercy --min-count 2` | disables mercy k-mers; optimal for high-coverage phage (Li et al. 2015) |
| cd-hit-est | `-c 1.00` | removes exact duplicates only; preserves biological SNP-level variants |

Both assemblers use k-mer range `21–127` (configurable). Post-assembly
length filter: contigs < `min_length` (default 200 bp) discarded.
`contigs.fasta` is used rather than `scaffolds.fasta` to avoid artificial
N-gaps that break ORF prediction.

**Warnings emitted when:**
- NR N50 < 5 kb — fragmented assembly
- NR contigs > 100 — unusual for purified phage; possible contamination

**Outputs:**
```
results/03_assembly/spades/{sample}/
results/03_assembly/megahit/{sample}/
results/03_assembly/combined/{sample}_contigs_nr.fasta   ← input to viral-id
reports/03_assembly/assembly_summary.tsv
```

---

### viral-id

**Viral identification** with geNomad.

```bash
conda run -n genomad phageflow viral-id \
  --sample-id s1 \
  --contigs results/03_assembly/combined/s1_contigs_nr.fasta
```

geNomad classifies each contig as virus / plasmid / chromosome and detects
integrated proviruses. Taxonomy is assigned at ICTV family/genus level; the
most specific resolved rank is reported as `best_taxon`.

Key parameters (configurable in `config.yaml`):

| Parameter | Default | Effect |
|---|---|---|
| `min_score` | 0.7 | virus probability threshold; 0.7 ≈ 97% precision (Camargo 2023) |
| `min_virus_hallmarks` | 0 | 0 = score-only; increase to require hallmark gene hits |
| `--splits` | 8 | parallel processing splits |
| `--cleanup` | — | removes large intermediates, keeps summary + FASTA |

**Warnings emitted when:**
- 0 viral contigs — lower `genomad.min_score` in config
- < 5% of contigs classified as viral — check host removal
- Provirus detected — review lifecycle module output carefully
- Best score < 0.7 — borderline classification; inspect contigs manually

**Outputs:**
```
results/04_viral_id/{sample}/                          geNomad working dir
results/04_viral_id/{sample}_virus.fna                 ← input to quality
reports/04_viral_id/genomad_summary.tsv
```

---

### quality

**Genome quality assessment, rescue, dereplication, co-binning, and
selection** using CheckV + cd-hit-est + seqkit.

```bash
phageflow quality \
  --sample-id s1 \
  --virus-fna results/04_viral_id/s1_virus.fna
```

**Pipeline (v3.2):**

```
CheckV end-to-end
    ↓
Tier selection + multi-path rescue
    ↓
Draft co-bin rescue (NEW v3.2)
    ↓
cd-hit-est 98% ANI dereplication (single-contigs only)
    ↓
HQ co-binning by geNomad taxonomy
    ↓
Taxonomy rename + seqkit -w 60 formatting
    ↓
One FASTA per candidate genome
```

**Tier logic and rescue paths:**

| Condition | Destination | Reference |
|---|---|---|
| Complete or High-quality | `annotation_ready/` always | Nayfach 2021 |
| MQ ≥ `min_completeness` (50%) | `annotation_ready/` | Nayfach 2021 |
| MQ < 50% | `drafts/` | — |
| LQ + ≥ `min_viral_genes` (1) OR ≥ `min_gene_density` (0.5 genes/kb) | `drafts/` | Roux 2019; Nayfach 2021 |
| ND + ≥ `large_nd_rescue_bp` (30 kb) + ≥ 1 viral gene | `annotation_ready/` ⬆ | Camargo 2023; Adriaenssens 2020 |
| ND + ≥ `length_rescue` (10 kb) + ≥ 1 viral gene | `drafts/` | Nayfach 2021 |
| ND + gene density ≥ 0.5 genes/kb | `drafts/` | Roux 2019 |
| < `min_contig_bp` (1,500 bp) | discarded | Roux 2019; Camargo 2023 |

**Draft co-bin rescue** (`_cobin_draft_rescue`, new in v3.2):
Contigs that individually fall below the completeness threshold but share
the same geNomad taxonomy are grouped. If the combined bin reaches
`min_bin_rescue_bp` (30 kb), it is promoted to `annotation_ready/` as a
multi-FASTA for Pharokka `--meta`. This is the primary mechanism for
recovering large fragmented phage genomes (Herelleviridae, Ackermannviridae)
where each assembly fragment is legitimately MQ at low individual completeness.

> **Why ND ≥ 30 kb goes directly to `annotation_ready/`:** CheckV's
> completeness estimation relies on HMM profiles derived from database
> sequences. Large myoviruses (Herelleviridae >100 kb, Ackermannviridae
> >150 kb) without close database relatives are consistently classified as
> ND regardless of actual completeness. In a purified phage preparation,
> a contig ≥ 30 kb with at least one viral gene is almost certainly phage.

**Dereplication:** cd-hit-est at 0.98 ANI intra-sample (Turner et al. 2021
ICTV species threshold = 0.95; 0.98 targets assembly artefacts, not distinct
strains). Multi-contig bins bypass dereplication — they are unique by
construction.

**Co-binning by taxonomy:** Single-contig HQ sequences sharing a geNomad
taxon are grouped into a single multi-FASTA. Multi-contig bins already
present (from draft rescue) are passed through unchanged to avoid incorrect
re-grouping.

**Rename convention:**
- Single contig: `{Family}_candidate_{NNN}.fasta`
- Multi-contig bin: `{Family}_multicontig_{NNN}.fasta`

Proviruses use CheckV's trimmed sequences from `proviruses.fna` (host
flanks removed) and are placed under `annotation_ready/proviruses/`.

**Warnings emitted when:**
- 0 HQ genomes + 0 draft bins promoted — lower `min_completeness` or check geNomad output
- Best completeness < 50% — fragmented genome; lifecycle confidence reduced
- Host genes > 3 — check host-removal step
- Max contamination > 5% — inspect `quality_summary.tsv`

**Outputs:**
```
results/05_quality/annotation_ready/phages/{candidate}.fasta   ← annotate input
results/05_quality/annotation_ready/proviruses/{candidate}.fasta
results/05_quality/drafts/{sample}_draft.fasta
results/05_quality/{sample}/quality_summary.tsv
reports/05_quality/checkv_summary.tsv
reports/05_quality/{sample}_rename_map.tsv
```

---

### annotate

**Structural and functional annotation** using Pharokka + Phold.

```bash
phageflow annotate \
  --sample-id s1 \
  --genome results/05_quality/annotation_ready/phages/<candidate>.fasta
```

**Step 1 — Pharokka** (Bouras et al. 2023, *Bioinformatics*):
- Gene calling with PHANOTATE
- Functional annotation against PHROG database
- `--dnaapler all`: genome reorientation to terminase large subunit
- Mode selected automatically: `--single` for one contig, `--meta` for
  multi-contig genomes (fragmented assemblies, co-binned candidates)

**Step 2 — Phold** (Bouras et al. 2024, *Bioinformatics*):
- Structure-based annotation via ProstT5 protein language model +
  Foldseek alignment to PHROGs
- Upgrades hypothetical proteins that lack sequence homology
- Invoked via `conda run -n phold` transparently

**Outputs:**
```
results/06_annotation/{sample}/pharokka/{sample}.gbk
results/06_annotation/{sample}/pharokka/{sample}.gff
results/06_annotation/{sample}/pharokka/{sample}_cds_functions.tsv
results/06_annotation/{sample}/phold/{sample}_phold.gbk   ← final GBK
reports/06_annotation/annotation_summary.tsv
```

> Run `annotate` before `safety` — the safety module reads Pharokka's
> CDS functions TSV for integrase detection.

---

### safety

**Biosafety screening** for phage therapy candidacy.

```bash
phageflow safety \
  --sample-id s1 \
  --genome results/05_quality/annotation_ready/phages/<candidate>.fasta
```

Three independent screens:

| Screen | Tool | Database | Threshold |
|---|---|---|---|
| Antimicrobial resistance genes (ARGs) | abricate | CARD | 80% id / 80% cov |
| Virulence factors (VFs) | abricate | VFDB | 80% id / 80% cov |
| Integrase / lysogeny genes | Pharokka CDS TSV | PHROG categories | — |

PHROG categories screened for lysogeny: `integration and excision`,
`lysogeny`, `transcription regulation` (CI repressor, Cro).

**Verdicts:**

| Verdict | Condition |
|---|---|
| ✓ PASS | No ARG, no VF, no integrase detected |
| ⚠ CAUTION | Integrase or VF detected — expert review required |
| ✗ FAIL | ARG detected — exclude from therapeutic use |

**Outputs:**
```
reports/07_safety/safety_summary.tsv
reports/07_safety/{sample}_safety_details.tsv
results/07_safety/{sample}_CARD.tsv
results/07_safety/{sample}_VFDB.tsv
```

---

### lifecycle

**Lifecycle prediction** (virulent / temperate) with BACPHLIP.

```bash
conda run -n bacphlip phageflow lifecycle \
  --sample-id s1 \
  --genome results/05_quality/annotation_ready/phages/<candidate>.fasta
```

BACPHLIP uses a random forest classifier trained on 62 lifestyle-associated
HMM profiles derived from >1,000 manually curated phage genomes (Hockenberry
& Wilke 2021, *PeerJ*). Profiles include integrases, CI repressors, Cro
proteins, excisionases (temperate markers) and holins, endolysins, spanins
(lytic markers). Reported accuracy: 98% on complete genomes; performance
decreases for fragmented or highly novel sequences.

| Score | Prediction | Implication |
|---|---|---|
| Virulent ≥ 0.9 | Confident lytic | Preferred for phage therapy |
| Temperate ≥ 0.9 | Confident lysogen | Expert review required before therapeutic use |
| Neither ≥ 0.9 | Ambiguous | Manual annotation review recommended; common for novel or fragmented genomes |

PhageFlow cross-checks BACPHLIP output against integrase genes detected by
the safety module. A conflict (BACPHLIP = Virulent + integrase detected)
triggers a warning.

**Outputs:**
```
results/08_lifecycle/{sample}.fasta.bacphlip
reports/08_lifecycle/lifecycle_summary.tsv
```

---

## Output structure

```
results/
├── 01_qc/
│   ├── {sample}_R1.fastq.gz          trimmed reads
│   └── {sample}_R2.fastq.gz
├── 02_host_removal/
│   ├── {sample}_R1.fastq.gz          host-filtered reads
│   └── {sample}_R2.fastq.gz
├── 03_assembly/
│   ├── spades/{sample}/              SPAdes working directory
│   ├── megahit/{sample}/             MEGAHIT working directory
│   └── combined/{sample}_contigs_nr.fasta   ← geNomad input
├── 04_viral_id/
│   ├── {sample}/                     geNomad working directory
│   └── {sample}_virus.fna            ← CheckV input
├── 05_quality/
│   ├── annotation_ready/
│   │   ├── phages/{candidate}.fasta  ← annotate / safety / lifecycle input
│   │   └── proviruses/{candidate}.fasta
│   ├── drafts/{sample}_draft.fasta   rescued contigs below HQ threshold
│   └── {sample}/quality_summary.tsv
├── 06_annotation/{sample}/
│   ├── pharokka/                     .gbk .gff .tsv
│   └── phold/                        *_phold.gbk  ← final annotation
├── 07_safety/
│   ├── {sample}_CARD.tsv
│   └── {sample}_VFDB.tsv
└── 08_lifecycle/
    └── {sample}.fasta.bacphlip

reports/
├── 01_qc/
│   ├── {sample}_fastp.json
│   ├── {sample}_fastp.html
│   ├── *_fastqc.html
│   ├── multiqc/multiqc_qc.html
│   └── qc_summary.tsv
├── 02_host_removal/
│   └── host_removal_summary.tsv
├── 03_assembly/
│   └── assembly_summary.tsv
├── 04_viral_id/
│   └── genomad_summary.tsv
├── 05_quality/
│   ├── checkv_summary.tsv
│   └── {sample}_rename_map.tsv
├── 06_annotation/
│   └── annotation_summary.tsv
├── 07_safety/
│   ├── safety_summary.tsv
│   └── {sample}_safety_details.tsv
└── 08_lifecycle/
    └── lifecycle_summary.tsv
```

---

## Running all samples in a loop

```bash
#!/usr/bin/env bash
set -euo pipefail

conda activate phageflow

CONFIG="config/config.yaml"
SAMPLES=(s1 s2 s3 s4 uce01 uce02 uce03 uce04)
HOST="/path/to/host.fasta"

for SAMPLE in "${SAMPLES[@]}"; do
    echo "=== Processing: $SAMPLE ==="

    R1=$(awk -v s="$SAMPLE" '$1==s{print $2}' config/samples.tsv)
    R2=$(awk -v s="$SAMPLE" '$1==s{print $3}' config/samples.tsv)

    phageflow qc -c "$CONFIG" --sample-id "$SAMPLE" --r1 "$R1" --r2 "$R2"

    phageflow host-removal -c "$CONFIG" --sample-id "$SAMPLE" \
        --r1 "results/01_qc/${SAMPLE}_R1.fastq.gz" \
        --r2 "results/01_qc/${SAMPLE}_R2.fastq.gz" \
        --host-file "$HOST"

    phageflow assembly -c "$CONFIG" --sample-id "$SAMPLE" \
        --r1 "results/02_host_removal/${SAMPLE}_R1.fastq.gz" \
        --r2 "results/02_host_removal/${SAMPLE}_R2.fastq.gz"

    conda run -n genomad phageflow viral-id -c "$CONFIG" \
        --sample-id "$SAMPLE" \
        --contigs "results/03_assembly/combined/${SAMPLE}_contigs_nr.fasta"

    phageflow quality -c "$CONFIG" --sample-id "$SAMPLE" \
        --virus-fna "results/04_viral_id/${SAMPLE}_virus.fna"

    # Iterate over all annotation-ready candidates for this sample
    for GENOME in results/05_quality/annotation_ready/phages/*.fasta \
                  results/05_quality/annotation_ready/proviruses/*.fasta; do
        [[ -f "$GENOME" ]] || continue

        phageflow annotate -c "$CONFIG" --sample-id "$SAMPLE" --genome "$GENOME"
        phageflow safety  -c "$CONFIG" --sample-id "$SAMPLE" --genome "$GENOME"

        conda run -n bacphlip phageflow lifecycle -c "$CONFIG" \
            --sample-id "$SAMPLE" --genome "$GENOME"
    done
done

echo "=== All samples done ==="
```

> **Note on multi-contig candidates:** genomes produced by co-binning
> (named `*_multicontig_*.fasta`) contain multiple sequences and are
> automatically processed by Pharokka in `--meta` mode.

---

## Databases

| Database | Version | Used by | Download command |
|---|---|---|---|
| CheckV | v1.5 | quality | `checkv download_database databases/checkv_db` |
| Pharokka DB | current | annotate | `install_databases.py -o databases/pharokka_db` |
| geNomad DB | v1.9+ | viral-id | `genomad download-database databases/genomad_db` |
| Phold DB | current | annotate | `phold install -d databases/phold_db` |
| CARD | bundled | safety | `abricate --setupdb` |
| VFDB | bundled | safety | `abricate --setupdb` |
| Kraken2 | optional | host-removal | see Kraken2 documentation |

---

## Citation

If you use PhageFlow, please cite the underlying tools:

- **fastp** — Chen et al. (2018) *Genome Biology* 19:274
- **SPAdes** — Bankevich et al. (2012) *J Comp Biol* 19:455–477
- **MEGAHIT** — Li et al. (2015) *Bioinformatics* 31:1674–1676
- **cd-hit** — Fu et al. (2012) *Bioinformatics* 28:3150–3152
- **bwa-mem2** — Vasimuddin et al. (2019) *IPDPS*
- **geNomad** — Camargo et al. (2023) *Nature Biotechnology*
- **CheckV** — Nayfach et al. (2021) *Nature Biotechnology* 39:578–585
- **Pharokka** — Bouras et al. (2023) *Bioinformatics* 39:btac776
- **Phold** — Bouras et al. (2024) *Bioinformatics*
- **BACPHLIP** — Hockenberry & Wilke (2021) *PeerJ* 9:e11396
- **CARD** — Alcock et al. (2023) *Nucleic Acids Research*
- **VFDB** — Liu et al. (2022) *Nucleic Acids Research*

---

## Authors

Fausto Cabezas-Mera · fcabezasmera@utem.cl  
Estefania Tisalema Guanopatin · etisalemag@correo.uss.cl  
Antonella Nole  
Dayra Valle

PhageFlow is released under the MIT License.

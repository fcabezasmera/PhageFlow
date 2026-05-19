# PhageFlow

**Modular bacteriophage genomics pipeline for Illumina paired-end sequencing.**

PhageFlow processes raw reads from purified phage preparations through quality control, host removal, assembly, viral identification, quality assessment, annotation, biosafety screening, and lifecycle prediction — one sample at a time, with full control at each step.

```
reads → QC → host removal → assembly → viral ID → quality →
annotation → biosafety → lifecycle → final genomes
```

---

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Configuration](#configuration)
- [Samples file](#samples-file)
- [Environment setup](#environment-setup)
- [Modules](#modules)
  - [01 qc](#01-qc)
  - [02 host-removal](#02-host-removal)
  - [03 assembly](#03-assembly)
  - [04 viral-id](#04-viral-id)
  - [05 quality](#05-quality)
  - [06 annotate](#06-annotate)
  - [07 safety](#07-safety)
  - [08 lifecycle](#08-lifecycle)
- [Output structure](#output-structure)
- [Running all samples in a loop](#running-all-samples-in-a-loop)
- [Databases](#databases)
- [Citation](#citation)

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

Most tools are installed via conda/bioconda:

```bash
# Main environment
conda install -c bioconda -c conda-forge \
    fastp fastqc multiqc \
    bwa-mem2 samtools seqtk \
    spades megahit cd-hit \
    checkv pharokka abricate \
    ncbi-datasets-cli kraken2 -y

# geNomad (separate environment recommended)
conda create -n genomad -c conda-forge -c bioconda genomad -y

# Phold (separate environment)
conda create -n pholdENV -c conda-forge -c bioconda phold -y

# BACPHLIP (separate environment)
conda create -n bacphlip_env -c conda-forge -c bioconda bacphlip -y
```

### 4. Register the PhageFlow CLI in every environment

Modules that require a separate environment (geNomad, Phold, BACPHLIP) still
need access to the `phageflow` command. Register it once per environment
without touching its dependencies:

```bash
conda activate genomad    && pip install -e . --no-deps -q
conda activate pholdENV   && pip install -e . --no-deps -q
conda activate bacphlip_env && pip install -e . --no-deps -q
conda activate phageflow  # return to main env
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
conda activate pholdENV
phold install -d databases/phold_db

# (Optional) Kraken2 — only for auto host detection
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

# Edit samples.tsv with your data, then run each module.
# Stay in the phageflow environment throughout — use conda run for
# modules that require a separate environment (viral-id, annotate, lifecycle).

conda activate phageflow

phageflow qc           --sample-id s1 --r1 data/s1_R1.fastq.gz --r2 data/s1_R2.fastq.gz

phageflow host-removal --sample-id s1 \
    --r1 results/01_qc/s1_R1.fastq.gz \
    --r2 results/01_qc/s1_R2.fastq.gz \
    --host-file /path/to/host.fasta

phageflow assembly     --sample-id s1 \
    --r1 results/02_host_removal/s1_R1.fastq.gz \
    --r2 results/02_host_removal/s1_R2.fastq.gz

conda run -n genomad phageflow viral-id --sample-id s1 \
    --contigs results/03_assembly/combined/s1_contigs_nr.fasta

phageflow quality   --sample-id s1 \
    --virus-fna results/04_viral_id/s1_virus.fna

phageflow annotate  --sample-id s1 \
    --genome results/05_quality/annotation_ready/s1_HQ.fasta

phageflow safety    --sample-id s1 \
    --genome results/05_quality/annotation_ready/s1_HQ.fasta

conda run -n bacphlip_env phageflow lifecycle --sample-id s1 \
    --genome results/05_quality/annotation_ready/s1_HQ.fasta
```

> **Tip:** you never need to leave the `phageflow` environment.
> Use `conda run -n <env>` for modules that require geNomad, Phold, or BACPHLIP.

---

## Configuration

Edit `config/config.yaml` before running:

```yaml
project: PhageFlow
mode: purified_phage      # only mode currently supported

samples_file: config/samples.tsv

threads:   22             # CPU threads
memory_gb: 64             # RAM limit (GB)

databases:
  checkv:   databases/checkv_db/checkv-db-v1.5
  pharokka: databases/pharokka_db
  genomad:  databases/genomad_db
  phold:    databases/phold_db
  kraken2:  databases/k2_db    # only needed for auto host detection

assembly:
  kmers: "21,33,55,77,99,127"   # SPAdes k-mer list

genomad:
  min_score: 0.7                 # geNomad virus score threshold

checkv:
  min_completeness: 50           # minimum % completeness for HQ selection
```

Override threads at runtime with `-t / --threads`:

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

- Supports both `.fastq` and `.fastq.gz`
- `sample_id` is used as prefix for all output files
- Lines starting with `#` are ignored

---

## Environment setup

PhageFlow uses a main environment (`phageflow`) plus three tool-specific
environments. You only need to activate `phageflow`; all other environments
are invoked via `conda run`.

| Environment    | Tools           | Used by modules    |
|----------------|-----------------|--------------------|
| `phageflow`    | all core tools  | 01–03, 05–07       |
| `genomad`      | geNomad         | 04 viral-id        |
| `pholdENV`     | Phold           | 06 annotate        |
| `bacphlip_env` | BACPHLIP        | 08 lifecycle       |

```bash
# One-time setup: register CLI in each env
conda activate genomad      && pip install -e . --no-deps -q
conda activate pholdENV     && pip install -e . --no-deps -q
conda activate bacphlip_env && pip install -e . --no-deps -q
conda activate phageflow
```

---

## Modules

### 01 qc

**Quality control and trimming** using fastp + FastQC + MultiQC.

```bash
phageflow qc \
  --sample-id s1 \
  --r1 data/s1_R1.fastq.gz \
  --r2 data/s1_R2.fastq.gz
```

**Parameters (PE150 phage-optimised):**

| Parameter | Value | Rationale |
|---|---|---|
| Quality threshold | Q ≥ 20 | 99% base accuracy |
| Min read length | 75 bp | reliable k-mer coverage |
| Sliding window | 4 bp / Q20 | 3' quality trimming |
| Low-complexity filter | ≥ 30% | removes homopolymers |
| Adapter detection | auto | paired-end mode |

**Outputs:** `results/01_qc/`, `reports/01_qc/multiqc/multiqc_qc.html`

---

### 02 host-removal

**Host read removal** using bwa-mem2 (when a reference is available) or Kraken2 (auto-detection mode).

```bash
# Mode 1: local FASTA
phageflow host-removal --sample-id s1 --r1 R1.fq.gz --r2 R2.fq.gz \
  --host-file /path/to/host.fasta

# Mode 2: NCBI accessions (auto-download)
phageflow host-removal --sample-id s1 --r1 R1.fq.gz --r2 R2.fq.gz \
  --accessions GCF_000013465.1,GCF_000007785.1

# Mode 3: accessions from file
phageflow host-removal --sample-id s1 --r1 R1.fq.gz --r2 R2.fq.gz \
  --accessions-file hosts.txt

# Mode 4: Kraken2 auto-detection (no reference needed)
phageflow host-removal --sample-id s1 --r1 R1.fq.gz --r2 R2.fq.gz
```

**Kraken2 filter parameters:** `--confidence 0.2 --minimum-hit-groups 2`  
Keeps unclassified reads (novel phage) + Viruses (taxid 10239).

**Outputs:** `results/02_host_removal/{sample}_R1/R2.fastq.gz`

---

### 03 assembly

**De novo assembly** with metaSPAdes + MEGAHIT, followed by 100% identity dereplication with cd-hit-est.

```bash
phageflow assembly \
  --sample-id s1 \
  --r1 results/02_host_removal/s1_R1.fastq.gz \
  --r2 results/02_host_removal/s1_R2.fastq.gz
```

**Tool configuration:**

| Tool | Key flags | Rationale |
|---|---|---|
| SPAdes | `--meta --only-assembler` | avoids error correction collapsing tail fiber variants |
| MEGAHIT | `--no-mercy --min-count 2` | optimal for high-coverage phage data |
| cd-hit-est | `-c 1.00` | removes exact duplicates only |

**Outputs:** `results/03_assembly/combined/{sample}_contigs_nr.fasta`

---

### 04 viral-id

**Viral identification** with geNomad.

> **Requires geNomad environment** — use `conda run` (no need to switch envs):

```bash
conda run -n genomad phageflow viral-id \
  --sample-id s1 \
  --contigs results/03_assembly/combined/s1_contigs_nr.fasta
```

geNomad classifies contigs as virus / plasmid / chromosome and detects integrated proviruses.
Default `--min-score 0.7` (configurable in `config.yaml`).

**Outputs:** `results/04_viral_id/{sample}_virus.fna`, `reports/04_viral_id/genomad_summary.tsv`

---

### 05 quality

**Genome quality assessment** with CheckV.

```bash
phageflow quality \
  --sample-id s1 \
  --virus-fna results/04_viral_id/s1_virus.fna
```

**Quality tiers (MIUVIG standard):**

| Tier | Completeness | Destination |
|---|---|---|
| Complete / High-quality | ≥ configured threshold | `annotation_ready/` |
| Medium-quality | below threshold | `drafts/` |
| Low-quality / Not-determined | — | discarded |

**Outputs:** `results/05_quality/annotation_ready/{sample}_HQ.fasta`,
`reports/05_quality/checkv_summary.tsv`

---

### 06 annotate

**Structural and functional annotation** using Pharokka (sequence homology) + Phold (structure-based).

> **Phold requires its own environment** — invoked automatically via `conda run`:

```bash
phageflow annotate \
  --sample-id s1 \
  --genome results/05_quality/annotation_ready/s1_HQ.fasta
```

**Workflow:**

1. **Pharokka** — PHANOTATE gene calling + PHROG database annotation + genome reorientation (`--dnaapler all`).
2. **Phold** — ProstT5 protein language model + Foldseek structural alignment to PHROGs; upgrades hypothetical proteins with functional predictions.

**Outputs:**
- `results/06_annotation/{sample}/pharokka/{sample}.gbk` — GenBank
- `results/06_annotation/{sample}/pharokka/{sample}.gff` — GFF3
- `results/06_annotation/{sample}/phold/{sample}_phold.gbk` — enhanced GenBank

---

### 07 safety

**Biosafety screening** for phage therapy candidacy assessment.

```bash
phageflow safety \
  --sample-id s1 \
  --genome results/05_quality/annotation_ready/s1_HQ.fasta
```

**Checks:**

| Screen | Tool | Database | Threshold |
|---|---|---|---|
| Antimicrobial resistance | abricate | CARD | 80% id / 80% cov |
| Virulence factors | abricate | VFDB | 80% id / 80% cov |
| Integrase / lysogeny | Pharokka annotation | PHROG categories | — |

**Safety verdict:**

| Verdict | Condition |
|---|---|
| ✓ PASS | No ARG, no VF, no integrase |
| ⚠ CAUTION | Integrase or VF detected — expert review required |
| ✗ FAIL | ARG detected — exclude from therapeutic use |

> Run `phageflow annotate` first for integrase detection (requires Pharokka output).

**Outputs:** `reports/07_safety/safety_summary.tsv`, `{sample}_safety_details.tsv`

---

### 08 lifecycle

**Lifecycle prediction** (virulent / temperate) with BACPHLIP.

> **Requires BACPHLIP environment** — use `conda run` (no need to switch envs):

```bash
conda run -n bacphlip_env phageflow lifecycle \
  --sample-id s1 \
  --genome results/05_quality/annotation_ready/s1_HQ.fasta
```

BACPHLIP uses a random forest classifier trained on 62 lifestyle-associated HMM profiles (integrases, CI repressors, holins, endolysins).

**Interpretation:**

| Score | Prediction |
|---|---|
| Virulent ≥ 0.9 | Confident lytic — preferred for phage therapy |
| Temperate ≥ 0.9 | Confident lysogen — requires expert review |
| Neither ≥ 0.9 | Ambiguous — manual annotation review recommended |

PhageFlow automatically cross-checks the result against integrase genes detected by the safety module and warns on conflicts.

**Outputs:** `results/08_lifecycle/{sample}.fasta.bacphlip`, `reports/08_lifecycle/lifecycle_summary.tsv`

---

## Output structure

```
results/
├── 01_qc/                  trimmed reads
├── 02_host_removal/        host-filtered reads
├── 03_assembly/
│   ├── spades/{sample}/    SPAdes output
│   ├── megahit/{sample}/   MEGAHIT output
│   └── combined/           *_contigs_nr.fasta  ← geNomad input
├── 04_viral_id/
│   ├── {sample}/           geNomad working directory
│   └── {sample}_virus.fna  ← CheckV input
├── 05_quality/
│   ├── annotation_ready/   *_HQ.fasta  ← annotation / safety / lifecycle input
│   └── drafts/             *_draft.fasta
├── 06_annotation/{sample}/
│   ├── pharokka/           .gbk .gff .tsv
│   └── phold/              *_phold.gbk
├── 07_safety/              CARD / VFDB abricate TSVs
└── 08_lifecycle/           .bacphlip files

reports/
├── 01_qc/                  fastp JSON, FastQC, MultiQC HTML
├── 02_host_removal/        host_removal_summary.tsv
├── 03_assembly/            assembly_summary.tsv, assembler logs
├── 04_viral_id/            genomad_summary.tsv
├── 05_quality/             checkv_summary.tsv, genome_selection.tsv
├── 06_annotation/          annotation_summary.tsv
├── 07_safety/              safety_summary.tsv, *_safety_details.tsv
└── 08_lifecycle/           lifecycle_summary.tsv
```

---

## Running all samples in a loop

PhageFlow processes one sample at a time by design. Use a shell loop for batch processing.
Stay in `phageflow` throughout — `conda run` handles tool-specific environments.

```bash
#!/usr/bin/env bash
set -euo pipefail

conda activate phageflow

CONFIG="config/config.yaml"
SAMPLES=(s1 s2 s3 s4 uce01 uce02 uce03 uce04)
HOST="/path/to/host.fasta"

for sample in "${SAMPLES[@]}"; do
    echo "=== Processing: $sample ==="

    phageflow qc -c "$CONFIG" --sample-id "$sample" \
        --r1 "$(grep "^${sample}" config/samples.tsv | cut -f2)" \
        --r2 "$(grep "^${sample}" config/samples.tsv | cut -f3)"

    phageflow host-removal -c "$CONFIG" --sample-id "$sample" \
        --r1 "results/01_qc/${sample}_R1.fastq.gz" \
        --r2 "results/01_qc/${sample}_R2.fastq.gz" \
        --host-file "$HOST"

    phageflow assembly -c "$CONFIG" --sample-id "$sample" \
        --r1 "results/02_host_removal/${sample}_R1.fastq.gz" \
        --r2 "results/02_host_removal/${sample}_R2.fastq.gz"

    conda run -n genomad phageflow viral-id -c "$CONFIG" --sample-id "$sample" \
        --contigs "results/03_assembly/combined/${sample}_contigs_nr.fasta"

    phageflow quality -c "$CONFIG" --sample-id "$sample" \
        --virus-fna "results/04_viral_id/${sample}_virus.fna"

    GENOME="results/05_quality/annotation_ready/${sample}_HQ.fasta"
    if [[ ! -f "$GENOME" ]]; then
        echo "  No HQ genome for $sample — skipping annotation"
        continue
    fi

    phageflow annotate -c "$CONFIG" --sample-id "$sample" --genome "$GENOME"

    phageflow safety -c "$CONFIG" --sample-id "$sample" --genome "$GENOME"

    conda run -n bacphlip_env phageflow lifecycle -c "$CONFIG" \
        --sample-id "$sample" --genome "$GENOME"

done

echo "=== All samples done ==="
```

---

## Databases

| Database | Used by | Download command |
|---|---|---|
| CheckV v1.5 | quality | `checkv download_database databases/checkv_db` |
| Pharokka DB | annotate | `install_databases.py -o databases/pharokka_db` |
| geNomad DB | viral-id | `genomad download-database databases/genomad_db` |
| Phold DB | annotate | `phold install -d databases/phold_db` |
| CARD | safety | bundled with abricate (`abricate --setupdb`) |
| VFDB | safety | bundled with abricate |
| Kraken2 (optional) | host-removal | see Kraken2 documentation |

---

## Citation

If you use PhageFlow, please cite the underlying tools:

- **fastp** — Chen et al. (2018) *Genome Biology* 19:274
- **SPAdes** — Bankevich et al. (2012) *J Comp Biol* 19:455–477
- **MEGAHIT** — Li et al. (2015) *Bioinformatics* 31:1674–1676
- **cd-hit** — Fu et al. (2012) *Bioinformatics* 28:3150–3152
- **geNomad** — Camargo et al. (2023) *Nature Biotechnology*
- **CheckV** — Nayfach et al. (2021) *Nature Biotechnology* 39:578–585
- **Pharokka** — Bouras et al. (2023) *Bioinformatics* 39:btac776
- **Phold** — Bouras et al. (2024) *Bioinformatics*
- **BACPHLIP** — Hockenberry & Wilke (2021) *PeerJ* 9:e11396
- **CARD** — Alcock et al. (2023) *Nucleic Acids Research*
- **VFDB** — Liu et al. (2022) *Nucleic Acids Research*
- **bwa-mem2** — Vasimuddin et al. (2019) *IPDPS*

---

## Author

Fausto Cabezas-Mera · fcabezasmera@utem.cl  
Estefania Tisalema Guanopatin · etisalemag@correo.uss.cl  
Antonella Nole  
Dayra Valle

PhageFlow is released under the MIT License.

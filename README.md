# PhageFlow

**Modular bacteriophage genomics pipeline for Illumina paired-end sequencing.**

PhageFlow processes raw reads from purified phage preparations through quality
control, host removal, assembly, viral identification, and quality assessment —
one sample at a time, with full control at each step.

```
reads → QC → host removal → assembly → viral ID → quality → [annotation → safety]
```

> **Mode:** `purified_phage` — optimised for preparations where the dominant
> nucleic acid is phage DNA. Parameters and rescue thresholds reflect this
> assumption throughout.

> **Status:** modules 01–05 fully implemented and validated.
> Annotation (Pharokka + Phold) and safety screening (CARD + VFDB) are
> available in the codebase but under active development.

---

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Configuration](#configuration)
- [Samples file](#samples-file)
- [Modules](#modules)
  - [qc](#qc)
  - [host-removal](#host-removal)
  - [assembly](#assembly)
  - [viral-id](#viral-id)
  - [quality](#quality)
- [Output structure](#output-structure)
- [Running all samples](#running-all-samples-in-a-loop)
- [Databases](#databases)
- [Citation](#citation)
- [Authors](#authors)

---

## Installation

### Requirements

- Linux (Ubuntu 20.04+)
- [Miniforge](https://github.com/conda-forge/miniforge) or Miniconda with **mamba**
- ~60 GB free disk space for databases
- NVIDIA GPU recommended for Phold (CPU fallback available)

### 1. Clone the repository

```bash
git clone https://github.com/fcabezasmera/PhageFlow.git
cd PhageFlow
```

### 2. Create the conda environment

All tools including geNomad and Phold live in a single `phageflow` environment.

```bash
mamba create -n phageflow -c conda-forge -c bioconda -y \
    python=3.11 \
    fastp fastqc multiqc \
    "bwa-mem2>=2.2" "samtools>=1.15" seqtk \
    spades megahit cd-hit \
    checkv pharokka \
    seqkit abricate \
    ncbi-datasets-cli kraken2

conda activate phageflow
```

Install GPU-optimised tools separately to resolve their CUDA dependencies:

```bash
# Phold (structure-based annotation, GPU-accelerated)
mamba install -y bioconda::phold "pytorch=*=cuda*"

# geNomad (viral identification)
mamba install -y bioconda::genomad

# Phynteny transformer (phage lifestyle prediction — future module)
mamba install -y bioconda::phynteny_transformer
```

Register PhageFlow in the environment:

```bash
pip install -e .
```

### 3. Set up databases

Run all downloads with `conda activate phageflow` from the project root.

```bash
# CheckV (Nayfach et al. 2021)
checkv download_database databases/checkv_db

# Pharokka / PHROG (Bouras et al. 2023)
install_databases.py -o databases/pharokka_db

# geNomad — v1.9+ required (Camargo et al. 2023)
genomad download-database databases/genomad_db

# Phold (Bouras et al. 2024)
phold install -d databases/phold_db

# ABRicate databases — CARD + VFDB
abricate --setupdb

# (Optional) Kraken2 — only for auto host-detection mode
kraken2-build --download-library bacteria --db databases/k2_db
kraken2-build --build --db databases/k2_db
```

### 4. Verify installation

```bash
phageflow check-tools
```

### 5. Copy configuration templates

```bash
phageflow config   # creates config/config.yaml
phageflow samples  # creates config/samples.tsv
```

---

## Quick Start

All commands run from the project root with `phageflow` active.
No secondary environments are needed.

```bash
conda activate phageflow

phageflow qc --sample-id s1 \
  --r1 data/s1_R1.fastq.gz \
  --r2 data/s1_R2.fastq.gz

phageflow host-removal --sample-id s1 \
  --r1 results/01_qc/s1_R1.fastq.gz \
  --r2 results/01_qc/s1_R2.fastq.gz \
  --host-file /path/to/host.fasta

phageflow assembly --sample-id s1 \
  --r1 results/02_host_removal/s1_R1.fastq.gz \
  --r2 results/02_host_removal/s1_R2.fastq.gz

phageflow viral-id --sample-id s1 \
  --contigs results/03_assembly/combined/s1_contigs_nr.fasta

phageflow quality --sample-id s1 \
  --virus-fna results/04_viral_id/s1_virus.fna
```

Annotation-ready genomes are written to:
```
results/05_quality/annotation_ready/phages/
results/05_quality/annotation_ready/proviruses/
```

---

## Configuration

Edit `config/config.yaml` before running.

```yaml
project: PhageFlow
version: "0.1.0"
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
  kraken2:  databases/k2_db

assembly:
  kmers:      "21,33,55,77,99,127"
  min_length: 200

genomad:
  min_score:           0.7
  min_virus_hallmarks: 0

checkv:
  min_completeness:   50
  min_viral_genes:    1
  length_rescue:      10000
  min_contig_bp:      1500
  large_nd_rescue_bp: 30000
  min_bin_rescue_bp:  30000
  min_gene_density:   0.5
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
| Minimum read length | 75 bp | reliable k-mer coverage |
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

> **Requirement:** bwa-mem2 mode requires `samtools ≥ 1.15`.
> Verify with `samtools --version`. If your environment has an older version,
> update with `mamba install -n phageflow "samtools>=1.15"`.

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

**Warnings emitted when:**
- Phage fraction < 30% — likely wrong reference or heavy contamination
- Surviving reads < 10,000 — assembly failure risk

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
| metaSPAdes | `--meta --only-assembler` | handles variable coverage; skips error correction that collapses tail fiber variants (Roux et al. 2019) |
| MEGAHIT | `--no-mercy --min-count 2` | disables mercy k-mers; optimal for high-coverage phage (Li et al. 2015) |
| cd-hit-est | `-c 1.00` | removes exact duplicates only; preserves biological SNP-level variants |

**Warnings emitted when:**
- NR N50 < 5 kb — fragmented assembly
- NR contigs > 100 — unusual for purified phage

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
phageflow viral-id --sample-id s1 \
  --contigs results/03_assembly/combined/s1_contigs_nr.fasta
```

geNomad classifies each contig as virus / plasmid / chromosome and detects
integrated proviruses. Taxonomy is assigned at ICTV family/genus level.

Key parameters:

| Parameter | Default | Effect |
|---|---|---|
| `min_score` | 0.7 | virus probability threshold; ~97% precision (Camargo 2023) |
| `min_virus_hallmarks` | 0 | 0 = score-only; increase to require hallmark gene hits |

**Warnings emitted when:**
- 0 viral contigs — lower `genomad.min_score` in config
- < 5% of contigs classified as viral — check host removal
- Provirus detected — review annotation output carefully

**Outputs:**
```
results/04_viral_id/{sample}/                          geNomad working dir
results/04_viral_id/{sample}_virus.fna                 ← input to quality
reports/04_viral_id/genomad_summary.tsv
```

---

### quality

**Genome quality assessment, rescue, dereplication, co-binning, and selection**
using CheckV + cd-hit-est + seqkit.

```bash
phageflow quality --sample-id s1 \
  --virus-fna results/04_viral_id/s1_virus.fna
```

**Pipeline:**

```
CheckV end-to-end
    ↓
Tier selection + multi-path rescue
    ↓
Draft co-bin rescue
    ↓
cd-hit-est 98% ANI dereplication (single-contig only)
    ↓
HQ co-binning by geNomad taxonomy
    ↓
Rename → seqkit -w 60
    ↓
One FASTA per candidate genome
```

**Tier logic and rescue paths:**

| Condition | Destination | Reference |
|---|---|---|
| Complete or High-quality | `annotation_ready/` always | Nayfach 2021 |
| MQ ≥ `min_completeness` (50%) | `annotation_ready/` | Nayfach 2021 |
| MQ < 50% | `drafts/` | — |
| LQ + ≥ `min_viral_genes` (1) OR density ≥ 0.5 genes/kb | `drafts/` | Roux 2019; Nayfach 2021 |
| ND + ≥ `large_nd_rescue_bp` (30 kb) + ≥ 1 viral gene | `annotation_ready/` ↑ | Camargo 2023; Adriaenssens 2020 |
| ND + ≥ `length_rescue` (10 kb) + ≥ 1 viral gene | `drafts/` | Nayfach 2021 |
| ND + density ≥ 0.5 genes/kb | `drafts/` | Roux 2019 |
| < `min_contig_bp` (1,500 bp) | discarded | Roux 2019; Camargo 2023 |

**Why ND ≥ 30 kb goes directly to `annotation_ready/`:** CheckV completeness
relies on HMM profiles from database sequences. Large myoviruses
(Herelleviridae >100 kb, Ackermannviridae >150 kb) without close database
relatives are consistently ND regardless of actual completeness. In a purified
phage preparation, a contig ≥ 30 kb with ≥ 1 viral gene is almost certainly
phage. (Camargo et al. 2023; Adriaenssens & Brister 2017)

**Draft co-bin rescue:** Contigs that individually fall below the completeness
threshold but share the same geNomad taxonomy are grouped. If the combined bin
reaches `min_bin_rescue_bp` (30 kb), it is promoted to `annotation_ready/` as
a multi-FASTA for Pharokka `--meta`. Confidence criteria for considering a
bin a single phage genome:

1. Each contig independently classified as viral by geNomad ≥ 0.7 (~97% precision)
2. All contigs share the same family-level taxonomy (common viral lineage)
3. Combined bin length ≥ 30 kb (biologically meaningful for large myoviruses)
4. At least 1 viral gene across the bin (viral marker confirmed)

*Limitation:* two co-purified phages of the same family would be incorrectly
merged. This is rare in purified preparations and is detectable downstream via
Pharokka annotation (unexpected gene content or unusual genome size).

**Dereplication:** cd-hit-est at 0.98 ANI intra-sample (Turner et al. 2021;
ICTV species threshold = 0.95, so 0.98 targets assembly artefacts).

**Naming convention:**

All candidates — single-contig and multi-contig — follow the same pattern:

```
{Family}_candidate_{NNN}.fasta
```

FASTA headers:
```
Single  : >Herelleviridae_candidate_001
Multi   : >Herelleviridae_candidate_002_ctg001
          >Herelleviridae_candidate_002_ctg002
```

The `rename_map.tsv` records `n_contigs` for each candidate. Downstream
annotation (Pharokka) auto-selects `--single` or `--meta` mode by counting
sequences in the FASTA file.

**Warnings emitted when:**
- 0 HQ genomes + 0 draft bins promoted — lower `min_completeness` or check geNomad output
- Best completeness < 50% — fragmented genome; annotation confidence reduced
- Host genes > 3 — check host-removal step
- Max contamination > 5% — inspect `quality_summary.tsv`

**Outputs:**
```
results/05_quality/annotation_ready/phages/{candidate}.fasta
results/05_quality/annotation_ready/proviruses/{candidate}.fasta
results/05_quality/drafts/{sample}_draft.fasta
results/05_quality/{sample}/quality_summary.tsv
reports/05_quality/checkv_summary.tsv
reports/05_quality/{sample}_rename_map.tsv
```

---

## Output structure

```
results/
├── 01_qc/
│   ├── {sample}_R1.fastq.gz
│   └── {sample}_R2.fastq.gz
├── 02_host_removal/
│   ├── {sample}_R1.fastq.gz
│   └── {sample}_R2.fastq.gz
├── 03_assembly/
│   ├── spades/{sample}/
│   ├── megahit/{sample}/
│   └── combined/{sample}_contigs_nr.fasta
├── 04_viral_id/
│   ├── {sample}/
│   └── {sample}_virus.fna
└── 05_quality/
    ├── annotation_ready/
    │   ├── phages/{Family}_candidate_{NNN}.fasta
    │   └── proviruses/{Family}_candidate_{NNN}.fasta
    ├── drafts/{sample}_draft.fasta
    └── {sample}/quality_summary.tsv

reports/
├── 01_qc/
│   ├── {sample}_fastp.json / .html
│   ├── *_fastqc.html
│   ├── multiqc/multiqc_qc.html
│   └── qc_summary.tsv
├── 02_host_removal/
│   └── host_removal_summary.tsv
├── 03_assembly/
│   └── assembly_summary.tsv
├── 04_viral_id/
│   └── genomad_summary.tsv
└── 05_quality/
    ├── checkv_summary.tsv
    └── {sample}_rename_map.tsv
```

---

## Running all samples in a loop

```bash
#!/usr/bin/env bash
set -euo pipefail

conda activate phageflow

CONFIG="config/config.yaml"
SAMPLES=(s1 s2 s3 s4)
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

    phageflow viral-id -c "$CONFIG" --sample-id "$SAMPLE" \
        --contigs "results/03_assembly/combined/${SAMPLE}_contigs_nr.fasta"

    phageflow quality -c "$CONFIG" --sample-id "$SAMPLE" \
        --virus-fna "results/04_viral_id/${SAMPLE}_virus.fna"

done

echo "=== All samples done ==="
```

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
| Kraken2 | optional | host-removal | `kraken2-build --download-library bacteria --db databases/k2_db` |

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
- **CARD** — Alcock et al. (2023) *Nucleic Acids Research*
- **VFDB** — Liu et al. (2022) *Nucleic Acids Research*

---

## Authors

Fausto Cabezas-Mera · fcabezasmera@utem.cl
Estefania Tisalema Guanopatin · etisalemag@correo.uss.cl
Antonella Nole
Dayra Valle

PhageFlow is released under the MIT License.

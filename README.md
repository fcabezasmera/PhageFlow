# PhageFlow

> **Modular bacteriophage genomics pipeline** for Illumina paired-end sequencing.
> One sample at a time. Full control at every step.

```
reads → QC → host removal → assembly → viral ID → quality → [annotation → safety]
```

**Mode:** `purified_phage` — optimised for preparations where the dominant nucleic acid is phage DNA.
**Status:** modules 01–05 fully implemented and validated on real datasets.

---

## Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Pipeline overview](#pipeline-overview)
- [Module reference](#module-reference)
- [Configuration](#configuration)
- [Output structure](#output-structure)
- [Running all samples](#running-all-samples)
- [Databases](#databases)
- [Citation](#citation)

---

## Installation

**Requirements:** Linux (Ubuntu 20.04+) · Miniforge/Miniconda with mamba · ~60 GB disk · NVIDIA GPU recommended for Phold

```bash
# 1. Clone
git clone https://github.com/fcabezasmera/PhageFlow.git
cd PhageFlow

# 2. Create environment
mamba create -n phageflow -c conda-forge -c bioconda -y \
    python=3.11 fastp fastqc multiqc \
    "bwa-mem2>=2.2" "samtools>=1.15" seqtk \
    spades megahit cd-hit mash \
    checkv pharokka seqkit abricate \
    ncbi-datasets-cli kraken2

conda activate phageflow

# GPU tools (optional but recommended)
mamba install -y bioconda::genomad bioconda::phold "pytorch=*=cuda*"

# Register
pip install -e .

# 3. Download databases
checkv download_database databases/checkv_db
genomad download-database databases/genomad_db
install_databases.py -o databases/pharokka_db
phold install -d databases/phold_db
abricate --setupdb

# 4. Verify
phageflow check-tools
```

> **Note:** `samtools >= 1.15` is required for singleton extraction. Verify with `samtools --version`.
> Update if needed: `mamba install -n phageflow "samtools>=1.15"`

---

## Quick Start

```bash
conda activate phageflow

phageflow qc           --sample-id s1 --r1 data/s1_R1.fastq.gz --r2 data/s1_R2.fastq.gz

phageflow host-removal --sample-id s1 \
    --r1 results/01_qc/s1_R1.fastq.gz --r2 results/01_qc/s1_R2.fastq.gz \
    --host-file /path/to/host.fasta

phageflow assembly     --sample-id s1 \
    --r1 results/02_host_removal/s1_R1.fastq.gz \
    --r2 results/02_host_removal/s1_R2.fastq.gz \
    --s1 results/02_host_removal/s1_singletons.fastq.gz

phageflow viral-id     --sample-id s1 \
    --contigs results/03_assembly/combined/s1_contigs_nr.fasta

phageflow quality      --sample-id s1 \
    --virus-fna results/04_viral_id/s1_virus.fna
```

Annotation-ready genomes → `results/05_quality/s1/annotation_ready/`

---

## Pipeline overview

| # | Module | Tools | Key output |
|---|--------|-------|-----------|
| 01 | `qc` | fastp · FastQC · MultiQC | trimmed reads · QC report |
| 02 | `host-removal` | bwa-mem2 · samtools · Kraken2 | phage reads · singletons |
| 03 | `assembly` | SPAdes · MEGAHIT · cd-hit-est | NR contigs FASTA |
| 04 | `viral-id` | geNomad | viral contigs · taxonomy |
| 05 | `quality` | CheckV · mash · seqkit | HQ candidate FASTAs |
| 06 | `annotate` ⚗️ | Pharokka · Phold · Phynteny | GBK · GFF · genome plot |
| 07 | `safety` ⚗️ | ABRicate · CARD · VFDB | AMR/VF report |

⚗️ Available in codebase, under active development.

---

## Module reference

### `qc` — Quality control

fastp parameters tuned for purified phage PE150:

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| Quality threshold | Q ≥ 20 (per-base) + mean Q ≥ 25 | 99% accuracy + read-level floor |
| Min read length | 75 bp | reliable k-mer coverage |
| Unqualified limit | ≤ 10% | tighter than fastp default 40% |
| Sliding window | 4 bp / Q20 | 3′ quality trimming |
| PE overlap correction | on | reduces base errors before assembly |
| Low-complexity filter | ≥ 30% | removes homopolymers |

Warnings emitted when pass rate < 75%, Q30 < 75%, or duplication > 70%
(fastp duplication is k-mer based — rates of 50–70% are **expected** at high phage coverage).

```bash
phageflow qc --sample-id s1 --r1 R1.fq.gz --r2 R2.fq.gz [-t 16]
```

---

### `host-removal` — Host read removal

Four modes, applied in priority order:

```bash
# Local FASTA (single file, folder, or path list)
phageflow host-removal --sample-id s1 --r1 R1.fq --r2 R2.fq \
    --host-file /path/to/host.fasta

# NCBI accessions (auto-download)
phageflow host-removal --sample-id s1 --r1 R1.fq --r2 R2.fq \
    --accessions GCF_000013465.1,GCF_000007785.1

# Accession file
phageflow host-removal --sample-id s1 --r1 R1.fq --r2 R2.fq \
    --accessions-file hosts.txt

# Kraken2 auto-detection (no reference needed)
phageflow host-removal --sample-id s1 --r1 R1.fq --r2 R2.fq
```

Singleton reads (mates of host-aligned reads that are themselves unmapped) are
written to `{sample}_singletons.fastq.gz`. Pass them to `assembly --s1` to improve
DTR/ITR boundary coverage and raise CheckV *Complete* classification rates
(Nayfach et al. 2021).

---

### `assembly` — De novo assembly

```bash
phageflow assembly --sample-id s1 \
    --r1 results/02_host_removal/s1_R1.fastq.gz \
    --r2 results/02_host_removal/s1_R2.fastq.gz \
    --s1 results/02_host_removal/s1_singletons.fastq.gz  # recommended
```

| Tool | Key flags | Why |
|------|-----------|-----|
| SPAdes | `--isolate --only-assembler` | isolate mode for uniform high-coverage phage; skip error correction at >200× |
| MEGAHIT | `--no-mercy --min-count 2` | optimal for high-coverage purified phage |
| cd-hit-est | `-c 1.00 -aS 0.85` | exact-duplicate removal; preserves biological variants |

k-mers: `21,33,55,77,99,127` (configurable). Output: NR contigs FASTA.

---

### `viral-id` — Viral identification

```bash
phageflow viral-id --sample-id s1 \
    --contigs results/03_assembly/combined/s1_contigs_nr.fasta
```

geNomad classifies each contig as virus / plasmid / chromosome and detects
integrated proviruses. Default `min_score = 0.7` (~97% precision, Camargo et al. 2023).

**Borderline rescue** (`purified_phage` mode): contigs with score in
[`rescue_min_score`, `min_score`) AND length ≥ `rescue_min_length_bp` are
added with a WARNING. Novel phage lineages without close database relatives
systematically score 0.4–0.6 despite being genuine phage.

---

### `quality` — Genome quality assessment

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
Draft co-bin rescue (shared-taxonomy bins ≥ 30 kb → annotation_ready/)
    ↓
Dereplication — mash for all sizes (circular-safe)
    ↓
HQ co-binning by geNomad taxonomy
    ↓
Rename + seqkit -w 60 → one FASTA per candidate
```

**Why mash for dereplication?** Phage genomes are circular. SPAdes and MEGAHIT
linearize the same circular sequence at different positions (circular permutation).
cd-hit-est uses linear alignment — individual HSPs each cover only a fraction of the
sequence, failing the `-aS` threshold even when combined coverage is 100%.
mash computes distances from k-mer composition, which is rotation-invariant:
permutations of the same circle have distance ≈ 0 and are correctly clustered.
This is especially important for ssDNA phages (e.g. Inoviridae ~6 kb).

**Tier rescue logic:**

| Condition | Destination |
|-----------|-------------|
| Complete or High-quality | `annotation_ready/` always |
| MQ ≥ 50% completeness | `annotation_ready/` |
| MQ < 50% | `drafts/` |
| LQ + ≥ 1 viral gene OR density ≥ 0.5 genes/kb | `drafts/` |
| ND + ≥ 30 kb + ≥ 1 viral gene | `annotation_ready/` (large myovirus rule) |
| ND + ≥ 10 kb + ≥ 1 viral gene | `drafts/` |
| < 1,500 bp | discarded |

**Naming:** `{Family}_candidate_{NNN}.fasta` — single and multi-contig genomes
follow the same convention. `rename_map.tsv` records `n_contigs` so downstream
annotation auto-selects `--single` or `--meta`.

---

## Configuration

```yaml
# config/config.yaml
mode: purified_phage
threads: 24
memory_gb: 32

databases:
  checkv:   databases/checkv_db/checkv-db-v1.5
  genomad:  databases/genomad_db
  pharokka: databases/pharokka_db
  phold:    databases/phold_db

assembly:
  kmers:      "21,33,55,77,99,127"
  min_length: 200

genomad:
  min_score:            0.7
  rescue_min_score:     0.4     # borderline rescue floor
  rescue_min_length_bp: 10000

checkv:
  min_completeness:   50    # MQ threshold for annotation_ready/
  min_contig_bp:      1500
  large_nd_rescue_bp: 30000  # ND large myovirus rule
  min_bin_rescue_bp:  30000  # draft co-bin rescue threshold
  min_gene_density:   0.5
```

Override threads at runtime: `phageflow quality --sample-id s1 ... -t 8`

Generate config and samples templates:

```bash
phageflow config    # → config/config.yaml
phageflow samples   # → config/samples.tsv
```

---

## Output structure

```
results/
├── 01_qc/{sample}_R1/R2.fastq.gz
├── 02_host_removal/{sample}_R1/R2.fastq.gz  ·  {sample}_singletons.fastq.gz
├── 03_assembly/combined/{sample}_contigs_nr.fasta
├── 04_viral_id/{sample}_virus.fna
└── 05_quality/{sample}/
    ├── annotation_ready/
    │   ├── phages/          {Family}_candidate_{NNN}.fasta
    │   └── proviruses/      {Family}_candidate_{NNN}.fasta
    ├── checkv/              quality_summary.tsv  proviruses.fna  …
    └── drafts/              {sample}_draft.fasta

reports/
├── 01_qc/         fastp.json · fastqc.html · multiqc_qc.html · qc_summary.tsv
├── 02_host_removal/  host_removal_summary.tsv
├── 03_assembly/      assembly_summary.tsv
├── 04_viral_id/      genomad_summary.tsv
└── 05_quality/       checkv_summary.tsv · {sample}_rename_map.tsv
```

Each sample's results are fully self-contained under `results/05_quality/{sample}/`.
Reprocessing or deleting one sample never touches another's files.

---

## Running all samples

```bash
#!/usr/bin/env bash
set -euo pipefail
conda activate phageflow

HOST="/path/to/host.fasta"
SAMPLES=(s1 s2 s3 s4)

for S in "${SAMPLES[@]}"; do
    R1=$(awk -v s="$S" '$1==s{print $2}' config/samples.tsv)
    R2=$(awk -v s="$S" '$1==s{print $3}' config/samples.tsv)

    phageflow qc           --sample-id "$S" --r1 "$R1" --r2 "$R2"
    phageflow host-removal --sample-id "$S" \
        --r1 "results/01_qc/${S}_R1.fastq.gz" \
        --r2 "results/01_qc/${S}_R2.fastq.gz" \
        --host-file "$HOST"
    phageflow assembly     --sample-id "$S" \
        --r1 "results/02_host_removal/${S}_R1.fastq.gz" \
        --r2 "results/02_host_removal/${S}_R2.fastq.gz" \
        --s1 "results/02_host_removal/${S}_singletons.fastq.gz"
    phageflow viral-id     --sample-id "$S" \
        --contigs "results/03_assembly/combined/${S}_contigs_nr.fasta"
    phageflow quality      --sample-id "$S" \
        --virus-fna "results/04_viral_id/${S}_virus.fna"
done
```

---

## Databases

| Database | Used by | Download |
|----------|---------|----------|
| CheckV v1.5 | quality | `checkv download_database databases/checkv_db` |
| geNomad DB v1.9+ | viral-id | `genomad download-database databases/genomad_db` |
| Pharokka DB | annotate | `install_databases.py -o databases/pharokka_db` |
| Phold DB | annotate | `phold install -d databases/phold_db` |
| CARD + VFDB | safety | `abricate --setupdb` |
| Kraken2 | host-removal (auto mode) | `kraken2-build --download-library bacteria --db databases/k2_db && kraken2-build --build --db databases/k2_db` |

---

## Citation

If you use PhageFlow, please cite the underlying tools:

- **fastp** — Chen et al. (2018) *Genome Biology* 19:274
- **SPAdes** — Bankevich et al. (2012) *J Comp Biol* 19:455–477
- **MEGAHIT** — Li et al. (2015) *Bioinformatics* 31:1674–1676
- **bwa-mem2** — Vasimuddin et al. (2019) *IPDPS*
- **geNomad** — Camargo et al. (2023) *Nature Biotechnology*
- **CheckV** — Nayfach et al. (2021) *Nature Biotechnology* 39:578–585
- **mash** — Ondov et al. (2016) *Genome Biology* 17:132
- **cd-hit** — Fu et al. (2012) *Bioinformatics* 28:3150–3152
- **Pharokka** — Bouras et al. (2023) *Bioinformatics* 39:btac776
- **Phold** — Bouras et al. (2024) *Bioinformatics*
- **CARD** — Alcock et al. (2023) *Nucleic Acids Research*
- **VFDB** — Liu et al. (2022) *Nucleic Acids Research*

---

## Authors

Fausto Cabezas-Mera · fcabezasmera@utem.cl
Estefania Tisalema Guanopatin · etisalemag@correo.uss.cl
Antonella Nole · Dayra Valle

PhageFlow is released under the [MIT License](LICENSE).

# PhageFlow

**Modular bacteriophage genomics pipeline for complete genome recovery from Illumina paired-end sequencing of purified phage preparations.**

---

## Authors

| Name | Affiliation |
|------|-------------|
| Fausto Cabezas-Mera | Universidad Tecnológica Metropolitana, Santiago, Chile |
| Estefania Tisalema-Guanopatin | Universidad Tecnológica Metropolitana, Santiago, Chile |
| Dayra Valle | Universidad Tecnológica Metropolitana, Santiago, Chile |
| Antonella Nole | Universidad Tecnológica Metropolitana, Santiago, Chile |

---

## Overview

PhageFlow implements a sequential, module-based workflow that takes raw Illumina reads from purified phage preparations through quality control, host depletion, de novo assembly, viral identification, genome quality assessment, and multi-tier structural annotation. Each module is independently re-runnable, produces structured outputs, and is parametrised by a single `config.yaml` file.

**Target output:** complete or high-quality phage genome(s) per sample, fully annotated at the functional level.

```
Raw reads (PE150 Illumina)
    │
    ▼  phageflow qc
01  Quality control + trimming          fastp · FastQC · MultiQC
    │
    ▼  phageflow host-removal
02  Host read removal                   bwa-mem2 · Kraken2
    │  ↳ singletons retained for DTR/ITR boundary coverage
    │
    ▼  phageflow assembly
03  De novo assembly                    SPAdes --isolate · MEGAHIT · cd-hit-est
    │
    ▼  phageflow viral-id
04  Viral identification                geNomad
    │
    ▼  phageflow quality
05  Genome quality + selection          CheckV · mash
    │
    ▼  phageflow annotate
06  Structural annotation               Pharokka → Phold → Phynteny
```

---

## Installation

```bash
git clone https://github.com/fcabezasmera/PhageFlow.git
cd PhageFlow

mamba env create -f environment.yml
conda activate phageflow
pip install -e .
```

> **Note:** if `abricate` or `blast-legacy` are installed in the same environment, force samtools ≥ 1.15:
> ```bash
> mamba install -n phageflow "samtools>=1.15"
> ```

---

## Quick Start

### Initialize a project

```bash
phageflow init /path/to/myproject
# Copy FASTQ files to /path/to/myproject/raw/
# Edit /path/to/myproject/config/config.yaml
```

### Run all modules (auto-detects samples from `raw/`)

```bash
phageflow run \
    --project /path/to/myproject \
    --accessions GCF_000005845.2   # propagation host accession(s)
```

### Run individual modules

```bash
# --sample-id is inferred from the R1 filename when omitted
phageflow qc --r1 raw/s1_R1.fastq.gz --r2 raw/s1_R2.fastq.gz

phageflow host-removal \
    --r1 results/01_qc/s1_R1.fastq.gz \
    --r2 results/01_qc/s1_R2.fastq.gz \
    --accessions GCF_000005845.2

phageflow assembly \
    --r1 results/02_host_removal/s1_R1.fastq.gz \
    --r2 results/02_host_removal/s1_R2.fastq.gz \
    --s1 results/02_host_removal/s1_singletons.fastq.gz

phageflow viral-id --contigs results/03_assembly/s1_contigs.fasta

phageflow quality --virus-fna results/04_viral_id/s1_virus.fna

phageflow annotate \
    --genome results/05_quality/s1/annotation_ready/phages/Podoviridae_candidate_001.fasta
```

### Resume an interrupted run

```bash
phageflow run --project /path/to/myproject --from-module quality
```

---

## Configuration

```bash
phageflow config          # write config/config.yaml template
phageflow check-tools     # verify all required tools and versions
phageflow status --project /path/to/myproject   # pipeline progress
```

Key parameters (all in `config/config.yaml`):

| Section | Parameter | Default | Notes |
|---------|-----------|---------|-------|
| `qc` | `average_qual` | 20 | per-read mean quality floor |
| `host_removal` | `always_include_accessions` | `[]` | propagation host(s), always removed |
| `genomad` | `min_score` | 0.7 | ~97% precision (Camargo et al. 2023) |
| `checkv` | `min_completeness` | 50 | MQ threshold for annotation_ready |
| `annotate` | `phold_gpu` | true | disable for CPU-only nodes |

---

## Output Structure

```
<project>/
├── results/
│   ├── 01_qc/                     trimmed reads
│   ├── 02_host_removal/           phage reads + singletons
│   ├── 03_assembly/               NR contigs
│   ├── 04_viral_id/               virus.fna + metadata
│   ├── 05_quality/
│   │   └── {sample}/
│   │       ├── annotation_ready/phages/   ← candidate FASTAs
│   │       └── checkv/
│   └── 06_annotation/
│       └── {sample}/{candidate}/
│           ├── {candidate}_annotated.gbk  ← canonical output
│           └── plots/
└── reports/                       logs, TSVs, MultiQC HTML
```

---

## Scientific Basis

Methodological decisions are documented inline in each module. Key references:

- Camargo et al. 2023, *Nat Biotechnol* 41:1783 — geNomad viral classification
- Nayfach et al. 2021, *Nat Biotechnol* 39:578 — CheckV genome quality, DTR/ITR
- Prjibelski et al. 2020, *Curr Protoc Bioinf* 70:e102 — SPAdes `--isolate`
- Bouras et al. 2023, *Bioinformatics* 39:btac776 — Pharokka annotation
- Bouras et al. 2025, bioRxiv — Phold structure-based annotation
- Grigson et al. 2025, bioRxiv — Phynteny synteny-aware annotation
- Turner et al. 2021, *Arch Virol* 166:2633 — ICTV 95% ANI species boundary

---

## License

MIT © 2025 Cabezas-Mera et al. — see [LICENSE](LICENSE).

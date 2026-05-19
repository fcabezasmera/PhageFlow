# PhageFlow v3.0

**Modular end-to-end pipeline for bacteriophage genomics** — from raw Illumina PE reads
to annotated, classified, and phylogenetically placed phage genomes. Designed for
reproducible research and publication-ready outputs.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
![Version](https://img.shields.io/badge/version-3.0.0-blue)
![Branch](https://img.shields.io/badge/branch-v3.0-green)
![Status](https://img.shields.io/badge/status-active%20development-orange)
![Python](https://img.shields.io/badge/python-3.11-blue)
![Platform](https://img.shields.io/badge/platform-Ubuntu%2024.04-lightgrey)

---

## Overview

PhageFlow processes purified phage shotgun sequencing (Illumina PE150) through 11
modular steps. Each module is an independent bash script sourcing a shared `utils.sh`,
with its own logging and QC checkpoints. The pipeline supports three execution modes
depending on the input type.

```
raw reads → QC → host removal → assembly (SPAdes + MEGAHIT) → dereplication (cd-hit-est)
         → geNomad (viral ID) → CheckV (quality gate) → Pharokka → Phold
         → biosafety → lifecycle → taxonomy → comparative genomics
```

### Three execution modes

| Mode | Input | Host reference | Use case |
|---|---|---|---|
| `purified_phage` | Illumina PE reads | Known (config.yaml) | Propagated phage stocks |
| `virome_reads` | Illumina PE reads | Auto (Kraken2) | Environmental viromes |
| `virome_contigs` | Assembled contigs | None | Pre-assembled datasets |

---

## Biological findings (v3.0 development dataset)

### UISEK cohort — *S. aureus* MRSA and *E. faecalis*

| ID | Phage | Family (geNomad) | Genome | CheckV | BACPHLIP |
|----|-------|------------------|--------|--------|----------|
| s1 | S1-Sat1 | Herelleviridae satellite | 17kb Complete | PICI/satellite element | LYTIC p=0.89 |
| s2 | S2-Ino1 | Inoviridae (*E. coli*) | 6.5kb Complete | Accidental isolation | LYTIC p=0.98 |
| s3 | **S3-Auto1** | **Autographiviridae** | **39kb Complete ★** | Complete | LYTIC p=1.00 |
| s4 | S4-Efae1 | Unclassified | 148kb draft | MQ fragmentado | LYTIC p=N/A |

### UCE cohort — *K. pneumoniae* ST258

| ID | Phage | Family | Genome | CheckV | Note |
|----|-------|--------|--------|--------|------|
| uce01 | **UCE-Ack1** | **Ackermannviridae** | **157.5kb Complete ★** | Complete | Tipo representante |
| uce02 | UCE-Ack1 | Ackermannviridae | 157.5kb Complete | Complete | Réplica (ANI=100%) |
| uce03 | UCE-Ack1 | Ackermannviridae | 157.5kb HQ | High-quality | Réplica (ANI=99.99%) |
| uce04 | UCE-Ack1 | Ackermannviridae | 157.3kb draft | LQ fragmentado | Réplica (ANI=100%) |

**Key finding:** uce01–04 represent the same Ackermannviridae phage recovered independently
from four samples at three geographically distinct urban sites, suggesting widespread
clonal distribution in Quito's urban wastewater system.

---

## Repository structure

```
PhageFlow/
├── README.md · LICENSE · CHANGELOG.md · .gitignore
├── config/
│   ├── config.yaml          # mode, threads, DB paths, host references
│   └── samples.tsv          # sample manifest with raw read paths
├── envs/                    # 3 conda environment definitions
│   ├── phageflow.yml        # mother env — 90% of tools
│   ├── genomad.yml          # geNomad dedicated
│   └── phold.yml            # Phold dedicated (GPU-enabled)
├── pipeline/
│   ├── utils.sh             # shared logging + helper functions
│   ├── 01_qc.sh             # fastp + FastQC + MultiQC
│   ├── 02_host_removal.sh   # bwa-mem2 — 10-genome combined reference
│   ├── 03_assembly.sh       # metaSPAdes + MEGAHIT + cd-hit-est NR
│   ├── 04_genomad.sh        # geNomad viral classification
│   ├── 05_checkv.sh         # CheckV quality + genome selection
│   ├── 06_pharokka.sh       # Pharokka structural annotation
│   ├── 07_phold.sh          # Phold structural annotation (hypotheticals)
│   ├── 08_safety.sh         # CARD + VFDB + integrase check
│   ├── 09_lifecycle.sh      # BACPHLIP lifecycle prediction
│   ├── 10_taxonomy.sh       # mash + IQ-TREE2 phylogeny (TerL)
│   └── 11_compare.sh        # Clinker + pyGenomeViz + dnadiff
├── data/
│   ├── raw/                 # reads — NOT tracked, symlinks to originals
│   ├── host_genomes/        # 10 reference FASTAs (see below)
│   └── databases/           # CheckV · Pharokka · geNomad · Phold — NOT tracked
├── reports/                 # QC summaries and logs — git-tracked
└── results/                 # pipeline outputs — NOT tracked
    ├── 01_qc/
    ├── 02_host_removal/
    ├── 03_assembly/         # spades/ · megahit/ · combined/ (NR)
    ├── 04_genomad/
    ├── 05_checkv/           # final_genomes/ · annotation_ready/
    ├── 06_pharokka/
    ├── 07_phold/
    ├── 08_safety/
    ├── 09_lifecycle/
    ├── 10_taxonomy/         # trees/ × 3 families · references/
    └── 11_compare/          # clinker/ · genomeviz/ · dotplots/
```

---

## Quick start

### 1. Clone and configure

```bash
git clone https://github.com/fcabezasmera/PhageFlow.git
cd PhageFlow
git checkout v3.0

# Edit sample paths
nano config/samples.tsv

# Review parameters
nano config/config.yaml
```

### 2. Create conda environments

```bash
# Mother environment (90% of tools)
conda env create -f envs/phageflow.yml

# geNomad dedicated
conda env create -f envs/genomad.yml

# Phold dedicated (GPU-enabled)
conda env create -f envs/phold.yml

# Lifecycle prediction (BACPHLIP — requires python 3.7)
conda create -n bacphlip_env python=3.7 numpy=1.19.2 \
  pandas=1.2.5 scikit-learn=0.23.2 -c conda-forge
conda activate bacphlip_env && conda install -c bioconda hmmer && pip install bacphlip
```

### 3. Download host reference genomes

```bash
bash pipeline/00_setup_hosts.sh
```

This downloads 10 reference genomes covering the target host species:

| Host | References | Accessions |
|---|---|---|
| *S. aureus* MRSA | USA300 · MW2 · MRSA252 | GCF_000013465.1 · GCF_000011265.1 · GCF_000013425.1 |
| *S. aureus* MSSA | NCTC8325 | GCF_000009645.1 |
| *E. faecalis* | V583 · OG1RF | GCF_000007785.1 · GCF_000172575.2 |
| *K. pneumoniae* | KPNIH1 · HS11286 · ST307 | GCF_000281535.2 · GCF_000240185.1 · GCF_002813595.1 |
| *E. coli* | K-12 MG1655 | GCF_000005845.2 |

### 4. Download databases

```bash
conda activate phageflow
checkv download_database data/databases/checkv_db

conda activate phage_annot
install_databases.py -o data/databases/pharokka_db

conda activate genomad
genomad download-database data/databases/genomad_db

conda activate pholdENV
phold install -d data/databases/phold_db
```

### 5. Run modules sequentially

```bash
conda activate phageflow
bash pipeline/01_qc.sh
bash pipeline/02_host_removal.sh
bash pipeline/03_assembly.sh

conda activate genomad
bash pipeline/04_genomad.sh

conda activate phageflow
bash pipeline/05_checkv.sh
bash pipeline/06_pharokka.sh

conda activate pholdENV
bash pipeline/07_phold.sh

conda activate phageflow
bash pipeline/08_safety.sh

conda activate bacphlip_env
bash pipeline/09_lifecycle.sh

conda activate phageflow
bash pipeline/10_taxonomy.sh
bash pipeline/11_compare.sh
```

---

## Key methodological decisions

### Assembly strategy
- **metaSPAdes `--meta --only-assembler`** — BayesHammer error correction disabled for
  high-coverage viral data (>200×), following Roux et al. 2019 (*eLife*)
- **MEGAHIT** — complementary assembler; different graph traversal recovers alternative contigs
- **cd-hit-est 100% identity** — removes exact duplicates between assemblers while
  preserving biological diversity (95% collapses true variants)

### Host removal — 10 combined references
Multi-strain reference panel compensates for unavailability of original host strain genomes.
*S. aureus* references include both MRSA (ST8/ST1/ST36) and MSSA (NCTC8325) since
phages do not distinguish resistance phenotype.

### Viral identification
geNomad v1.12 (score threshold 0.7) applied immediately post-assembly, replacing manual
coverage-based filtering. Coverage analysis was used in v1 (PhageFlow) as a diagnostic
tool but geNomad's gene content classification is more robust.

### Genome quality gate
Only Complete and High-quality genomes (CheckV) proceed to full annotation. Medium-quality
genomes are retained as drafts and reported in supplementary material.

### Replicate isolates
When ANI ≥ 99% between genomes from different samples, a single representative is
designated for annotation and phylogenetics. Replicate isolations are reported as
evidence of phage distribution/abundance.

---

## Conda environments

| Env | Key tools | Modules |
|-----|-----------|---------|
| `phageflow` | fastp · FastQC · MultiQC · bwa-mem2 · samtools · SPAdes · MEGAHIT · cd-hit-est · CheckV · Pharokka · DIAMOND · HMMER · MAFFT · IQ-TREE2 · MASH · trimAl · MUMmer4 · abricate · Clinker · pyGenomeViz · Kraken2 | 01-06 · 08 · 10-11 |
| `genomad` | geNomad 1.12 | 04 |
| `pholdENV` | Phold 1.2.5 · FoldSeek · PyTorch CUDA | 07 |
| `bacphlip_env` | BACPHLIP 0.9.6 · HMMER · python 3.7 | 09 |

### Known installation issues

**BACPHLIP** requires `scikit-learn==0.23.1` (python ≤3.8) and patches for `np.float`
deprecated in numpy ≥1.20:
```bash
BDIR=$(python3 -c "import bacphlip,os; print(os.path.dirname(bacphlip.__file__))")
sed -i 's/np\.float\b/float/g' "$BDIR/bacphlip.py"
find "$BDIR" -name "*.pyc" -delete
```

**Phold** truncates GenBank record IDs >18 chars. Fix via BioPython before running:
```python
from Bio import SeqIO
records = list(SeqIO.parse("pharokka.gbk", "genbank"))
for rec in records: rec.id = "shortID"; rec.name = "shortID"
SeqIO.write(records, "pharokka_fixed.gbk", "genbank")
```

---

## Two-paper strategy

### Paper 1 — Biology (current dataset)
> *"Isolation and genomic characterization of novel bacteriophages from urban wastewater
> targeting multidrug-resistant ESKAPE pathogens in Quito, Ecuador"*

**Target:** PHAGE Journal · Frontiers in Microbiology · Viruses
**Timeline:** 6–8 weeks

**Phages characterized:**
1. **UCE-Ack1** — Ackermannviridae, *K. pneumoniae* ST258, 157.5kb Complete
2. **S3-Auto1** — Autographiviridae (T7-like), *S. aureus* MRSA, 39kb Complete
3. **S1-Sat1** — Herelleviridae satellite/PICI element, 17kb (supplementary)
4. **S2-Ino1** — Inoviridae, *E. coli* contaminant, 6.5kb (supplementary)
5. **S4-Efae1** — Unclassified, *E. faecalis*, 148kb draft (supplementary)

### Paper 2 — Pipeline tool
> *"PhageFlow: an end-to-end modular pipeline for bacteriophage genomics
> from short-read sequencing"*

**Target:** Bioinformatics Advances · Briefings in Bioinformatics
**Timeline:** 5–6 months (requires benchmarking + synthetic dataset)

**Additional modules for Paper 2:**
- Synthetic dataset (ART-Illumina, 0–90% contamination levels)
- Benchmarking vs VirSorter2 · VIBRANT · Sphae
- `virome_reads` and `virome_contigs` mode validation
- iPHoP host prediction module

---

## Citations

Please cite the following tools:

| Tool | Reference |
|------|-----------|
| **geNomad** | Camargo et al., *Nat Biotechnol* 2023. doi:10.1038/s41587-023-01953-y |
| **Pharokka** | Bouras et al., *Bioinformatics* 2023. doi:10.1093/bioinformatics/btac821 |
| **Phold** | Bouras et al., *Bioinformatics Advances* 2024. doi:10.1093/bioadv/vbae074 |
| **CheckV** | Nayfach et al., *Nat Biotechnol* 2021. doi:10.1038/s41587-020-00774-7 |
| **SPAdes** | Prjibelski et al., *Curr Protoc* 2020. doi:10.1002/cpbi.102 |
| **MEGAHIT** | Li et al., *Bioinformatics* 2015. doi:10.1093/bioinformatics/btv033 |
| **BACPHLIP** | Hockenberry & Wilke, *PeerJ* 2021. doi:10.7717/peerj.11396 |
| **Clinker** | Gilchrist & Chooi, *Bioinformatics* 2021. doi:10.1093/bioinformatics/btab007 |
| **IQ-TREE2** | Minh et al., *Mol Biol Evol* 2020. doi:10.1093/molbev/msaa015 |
| **cd-hit** | Fu et al., *Bioinformatics* 2012. doi:10.1093/bioinformatics/bts565 |

---

## License

MIT — see [LICENSE](LICENSE)

## Contact

Fausto Cabezas-Mera — fcabezasmera@utem.cl
Estefanía Tisalema-Guanopatín — etisalemag@correo.uss.cl

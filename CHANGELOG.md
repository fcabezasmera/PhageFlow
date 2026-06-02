# Changelog

All notable changes to PhageFlow are documented here.
This project adheres to [Semantic Versioning](https://semver.org/).

## [1.0.0] - 2026-06-01

First production release. Complete 7-module pipeline validated end-to-end on
purified phage, virome, and metagenome datasets, including a multi-phage
recovery (4 distinct genomes from a single K. pneumoniae sample) on an HPC cluster.

### Pipeline
- **01 qc**: fastp + FastQC + MultiQC, read-length auto-detection (PE150/250/300)
- **02 host-removal**: bwa-mem2, Kraken2, and pass-through modes
- **03 assembly**: SPAdes + MEGAHIT + cd-hit-est, k-mer range by read length
- **03b coverage**: CoverM per-contig profiling
- **04 viral-id**: geNomad with 3 kb rescue tier, CrAss-like genetic-code detection
- **05 quality**: CheckV tiers, blastn/minimap2 dereplication, circular-rotation detection
- **06 annotate**: Pharokka + Phold two-tier cascade (--hyps strategy)
- **07 resistance**: dual-evidence AMR/VFDB + ACR + DefenseFinder + NetFlax
- **runner**: full pipeline with --from-module resume and per-sample error handling

### Infrastructure
- `phageflow download-databases`: one-command setup for CheckV, geNomad, Pharokka,
  Phold (incl. GPU structure DB), and Kraken2 standard-16GB; writes user config
- Config hierarchy: CLI > project > ~/.config/phageflow/ > bundled; `PHAGEFLOW_DB` env var
- Available on Bioconda: `conda install -c bioconda phageflow`

### Robustness fixes since initial bioconda submission
- CLI now loads user config from ~/.config/phageflow/config.yaml
- geNomad nested database path (genomad_db/genomad_db) resolved automatically
- Phold GPU database built during download-databases
- Phold auto-retries on CPU when Foldseek GPU search fails (no-GPU clusters)
- Phold falls back to Pharokka GBK for tiny/divergent genomes with no structural hits

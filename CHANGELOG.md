# Changelog — PhageFlow

## [0.1.0] — 2025-05

### Added
- Module 00: conda envs (9), host genomes (4), read symlinks + pigz compression
- Module 01: FastQC + MultiQC raw reads (8 samples, 2 cohorts)
- Module 02: fastp trimming — all samples >85% pass, s1 lowest at 87.3%
- Module 03: multi-host removal (bwa-mem2)
  - 03b: per-host removal (MRSA/Efae/Kpn)
  - 03c: combined 3-host removal
  - 03d: cross-mapping diagnosis — detected E. coli in s2/s3/s4
  - 03e: 4-host removal adding E. coli K-12 (final canonical reads)
- Module 05: metaSPAdes assembly (--meta, 8 samples)
- Module 06: Phables v1.5 GFA-based resolution
  - No complex components found — consistent with purified phage at high coverage
  - Gurobi 13 academic license required; internal YAML patches documented
- Module 07: CheckV v1.0.1 quality assessment
  - geNomad v1.12 added as primary viral classifier (DB v1.9)
  - 4 complete genomes: s3 (Autographiviridae), uce01/02/03 (Ackermannviridae)
- Module 09: Pharokka v1.9.1 structural annotation (in progress)

### Key biological findings
- s3: Autographiviridae (T7-like Podovirus) — 39kb, DTR, Complete, score=0.983
- s1: Herelleviridae (Myovirus) — 130kb, fragmentado
- s2: E. faecalis phage — 98kb, fragmentado (E. coli contamination)
- s4: E. faecalis phage — 148kb, fragmentado
- uce01-03: Ackermannviridae (Myovirus) — 157.5kb, Complete, DTR, 31 hallmarks
- uce04: Ackermannviridae — 157.3kb, fragmentado (3 contigs)
- E. coli contamination in UISEK s2/s3/s4: confirmed by genome coverage >89%

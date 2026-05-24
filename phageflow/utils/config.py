"""PhageFlow configuration and sample loading."""

from __future__ import annotations
import math
import multiprocessing
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

import yaml


# ---------------------------------------------------------------------------
# Thread auto-detection
# ---------------------------------------------------------------------------

def _default_threads() -> int:
    """90% of available logical CPUs, minimum 1.

    Reserves ~10% for OS scheduling and I/O daemons.
    Override via config.yaml ``threads:`` or CLI ``-t N``.
    """
    return max(1, math.floor(multiprocessing.cpu_count() * 0.9))


# ---------------------------------------------------------------------------
# Sample manifest
# ---------------------------------------------------------------------------

@dataclass
class Sample:
    sample_id: str
    r1:        Path
    r2:        Path

    def validate(self) -> List[str]:
        errors = []
        for attr, path in [("r1", self.r1), ("r2", self.r2)]:
            if not path.exists():
                errors.append(f"[{self.sample_id}] {attr} not found: {path}")
        return errors


def load_samples(tsv_path: Path) -> List[Sample]:
    """Load samples.tsv (columns: sample_id  r1  r2)."""
    samples: List[Sample] = []
    with open(tsv_path) as f:
        header: Optional[dict] = None
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if header is None:
                header = {h.strip(): i for i, h in enumerate(cols)}
                continue
            samples.append(Sample(
                sample_id = cols[header["sample_id"]].strip(),
                r1        = Path(cols[header["r1"]].strip()),
                r2        = Path(cols[header["r2"]].strip()),
            ))
    return samples


# ---------------------------------------------------------------------------
# Sub-configs
# ---------------------------------------------------------------------------

@dataclass
class AssemblyConfig:
    # k-mers shared between SPAdes and MEGAHIT for reproducibility.
    # k=141 is the MEGAHIT maximum; improves resolution of repetitive regions
    # in large genomes (Ackermannviridae, jumbo phages >150 kb).
    # Li et al. 2015, Bioinformatics 31:1674; Al-Shayeb et al. 2020, Nature 578:425.
    kmers:      str = "21,33,55,77,99,127,141"
    min_length: int = 200

    # Iterative refinement (opt-in). Requires quality.py to have run first.
    # See modules/assembly_refine.py — run AFTER quality.py.
    iterative_refinement: bool = False


@dataclass
class HostRemovalConfig:
    """Configuration for Module 02 — host read removal.

    always_include_accessions
        NCBI accessions merged with any --host-file / --accessions the user
        provides. Empty by default. Add lab-specific strains that are always
        present in your preparations (e.g. the propagation host).
        These are downloaded once per sample run and stored in
        reports/02_host_removal/{sample_id}/host_genomes/.
        Example:
            always_include_accessions:
              - GCF_000005845.2   # E. coli K-12 MG1655

    kraken2_postfilter  [Level B, opt-in]
        When True (Kraken2 mode only): after extracting phage reads with
        Kraken2, parse the classification report and download reference genomes
        for bacterial species exceeding postfilter_min_pct%. A second bwa-mem2
        pass removes residual contamination and captures singletons at
        DTR/ITR boundaries that Kraken2 discards (Nayfach et al. 2021).

    postfilter_min_pct
        Minimum percentage in the Kraken2 report to trigger per-species
        genome download for Level B post-filter. Default 1.0%.

    contamination_warn_pct  [Level A, always-on when Kraken2 DB available]
        After bwa-mem2 extraction, a subsample of phage reads is classified
        with Kraken2. If the bacterial fraction exceeds this threshold, a
        warning is emitted and the contaminating taxa are listed.
        Diagnostic only — no reads are removed.
        Only runs if databases.kraken2 is configured.
    """
    always_include_accessions: list = field(default_factory=list)
    kraken2_postfilter:         bool  = False
    postfilter_min_pct:         float = 1.0
    contamination_warn_pct:     float = 5.0


@dataclass
class GenomadConfig:
    # Classification score threshold. 0.7 ≈ 97% precision.
    # Camargo et al. 2023, Nat Biotechnol.
    min_score:    float = 0.7
    min_hallmarks: int  = 0

    # Borderline large-contig rescue (purified phage):
    # Novel lineages without close DB relatives score 0.4–0.6.
    rescue_min_score:     float = 0.4
    rescue_min_length_bp: int   = 10_000

    # MMseqs2 marker search sensitivity (default 4.2).
    # Increase to 5.7–7.5 to detect distant relatives at speed cost.
    sensitivity: float = 4.2

    # --splits: 0 = auto memory management (geNomad default).
    # Only increase if MMseqs2 OOMs. NOT a parallelism parameter.
    # Camargo et al. 2023 — splits controls DB batching, not speed.
    splits: int = 0

    # Provirus detection: disabled by default for purified phage.
    # In purified preparations, bacterial chromosomes should be absent.
    # If proviruses appear, it indicates host-removal failure.
    disable_find_proviruses: bool = True

    # Score calibration adjusts scores for sample composition.
    # --composition virome: most sequences are assumed viral.
    # Reduces false negatives in novel lineages. Camargo et al. 2023.
    enable_score_calibration: bool = True
    composition:              str  = "virome"

    # Taxonomy resolution below family rank (subfamily, genus, species).
    # Required for biologically meaningful co-binning in quality module.
    # Without this, all members of a family co-bin regardless of genus.
    lenient_taxonomy: bool = True


@dataclass
class CheckvConfig:
    # Nayfach et al. 2021, Nat Biotechnol 39:578.
    min_completeness:   float = 50.0   # MQ threshold for annotation_ready
    min_viral_genes:    int   = 1      # LQ rescue: ≥N viral genes → drafts
    length_rescue:      int   = 10_000 # ND ≥N bp + ≥1 gene → drafts
    min_contig_bp:      int   = 1_500  # contigs shorter than this discarded
    large_nd_rescue_bp: int   = 30_000 # ND ≥N bp + ≥1 gene → annotation_ready
    min_bin_rescue_bp:  int   = 30_000 # draft bin (same taxon) ≥N bp → annotation_ready
    min_gene_density:   float = 0.5    # genes/kb minimum for density rescue

    # kmer_freq > this threshold indicates concatemer / genome copies.
    # CheckV PyPI docs: values >1.5–1.7 are problematic.
    max_kmer_freq:    float = 1.5
    # genome_copies > this threshold flags possible concatemers.
    max_genome_copies: float = 1.5


@dataclass
class AnnotateConfig:
    """Settings for Module 06 (Pharokka → Phold → Phynteny → merge)."""
    # Phold
    phold_gpu:        bool  = True  # --foldseek_gpu; False → --cpu
    phold_finetune:   bool  = True  # phage-finetuned ProstT5 (recommended)
    phold_batch_size: int   = 1     # increase on GPU (4–8)
    phold_sensitivity: float = 9.5  # Foldseek sensitivity; 7.5 is 2× faster

    # Phynteny Transformer
    phynteny_confidence: float = 0.8  # default in tool source (not 0.7)

    # Genetic code default. Overridden at runtime from geNomad _virus_genes.tsv.
    # Code 11: standard phage/bacterial. Code 15: crAss-like (TAG→Gln).
    # If geNomad reports a different code per candidate, that value takes priority.
    genetic_code: int = 11


@dataclass
class DatabaseConfig:
    checkv:   Path           = Path("databases/checkv_db/checkv-db-v1.5")
    pharokka: Path           = Path("databases/pharokka_db")
    genomad:  Path           = Path("databases/genomad_db")
    phold:    Path           = Path("databases/phold_db")
    phynteny: Path           = Path("databases/phynteny_db")
    kraken2:  Optional[Path] = None


@dataclass
class EnvConfig:
    main: str = "phageflow"


# ---------------------------------------------------------------------------
# Main Config
# ---------------------------------------------------------------------------

@dataclass
class Config:
    project:      str  = "PhageFlow"
    version:      str  = "0.1.0"
    threads:      int  = field(default_factory=_default_threads)
    memory_gb:    int  = 32
    samples_file: Path = Path("config/samples.tsv")
    workdir:      Path = Path(".")

    # Optional output/report directory overrides (set by -o / --reports-dir CLI flags
    # or by dirs.results / dirs.reports in config.yaml).
    # When None, fall back to workdir/results and workdir/reports.
    _output_dir:  Optional[Path] = field(default=None, repr=False, compare=False)
    _reports_dir: Optional[Path] = field(default=None, repr=False, compare=False)

    databases:    DatabaseConfig    = field(default_factory=DatabaseConfig)
    envs:         EnvConfig         = field(default_factory=EnvConfig)
    assembly:     AssemblyConfig    = field(default_factory=AssemblyConfig)
    host_removal: HostRemovalConfig = field(default_factory=HostRemovalConfig)
    genomad:      GenomadConfig     = field(default_factory=GenomadConfig)
    checkv:       CheckvConfig      = field(default_factory=CheckvConfig)
    annotate:     AnnotateConfig    = field(default_factory=AnnotateConfig)

    # ------------------------------------------------------------------
    # Directory resolution
    # ------------------------------------------------------------------

    @property
    def results_dir(self) -> Path:
        """Base directory for pipeline results.

        Priority (highest to lowest):
          1. -o / --output-dir CLI flag
          2. dirs.results in config.yaml
          3. workdir/results  (legacy default)
        """
        return self._output_dir or (self.workdir / "results")

    @property
    def reports_dir(self) -> Path:
        """Base directory for reports (logs, TSVs, MultiQC).

        Priority (highest to lowest):
          1. --reports-dir CLI flag  or  set_reports_dir()
          2. dirs.reports in config.yaml
          3. workdir/reports  (default)

        When -o is used, the CLI explicitly sets both results and reports
        under the project base directory, so no special logic is needed here.
        """
        return self._reports_dir or (self.workdir / "reports")

    def set_output_dir(self, path: Optional[Path]) -> None:
        if path is not None:
            self._output_dir = Path(path)

    def set_reports_dir(self, path: Optional[Path]) -> None:
        if path is not None:
            self._reports_dir = Path(path)

    def results(self, step: str) -> Path:
        return self.results_dir / step

    def reports(self, step: str) -> Path:
        return self.reports_dir / step


# ---------------------------------------------------------------------------
# YAML loader
# ---------------------------------------------------------------------------

def load_config(config_path: Path, workdir: Optional[Path] = None) -> Config:
    """Load config.yaml and return a fully resolved Config object."""
    with open(config_path) as f:
        raw = yaml.safe_load(f)

    cfg         = Config()
    cfg.workdir = workdir or config_path.parent.parent
    cfg.project = raw.get("project", cfg.project)

    # Threads: YAML overrides auto-detection only when explicitly set.
    if "threads" in raw:
        cfg.threads = int(raw["threads"])
    if "memory_gb" in raw:
        cfg.memory_gb = int(raw["memory_gb"])

    # Output / report directories from YAML dirs block.
    if "dirs" in raw:
        d = raw["dirs"]
        if "results" in d:
            p = d["results"]
            cfg._output_dir = Path(p) if Path(p).is_absolute() else cfg.workdir / p
        if "reports" in d:
            p = d["reports"]
            cfg._reports_dir = Path(p) if Path(p).is_absolute() else cfg.workdir / p

    # samples_file
    sf = raw.get("samples_file", str(cfg.samples_file))
    cfg.samples_file = Path(sf) if Path(sf).is_absolute() else cfg.workdir / sf

    # databases
    if "databases" in raw:
        db = raw["databases"]
        def _db(key: str, default: Path) -> Path:
            v = db.get(key, str(default))
            return Path(v) if Path(v).is_absolute() else cfg.workdir / v
        cfg.databases.checkv   = _db("checkv",   cfg.databases.checkv)
        cfg.databases.pharokka = _db("pharokka", cfg.databases.pharokka)
        cfg.databases.genomad  = _db("genomad",  cfg.databases.genomad)
        cfg.databases.phold    = _db("phold",    cfg.databases.phold)
        cfg.databases.phynteny = _db("phynteny", cfg.databases.phynteny)
        if "kraken2" in db:
            cfg.databases.kraken2 = cfg.workdir / db["kraken2"]

    if "envs" in raw:
        cfg.envs.main = raw["envs"].get("main", cfg.envs.main)

    # assembly
    if "assembly" in raw:
        a = raw["assembly"]
        cfg.assembly.kmers      = a.get("kmers",      cfg.assembly.kmers)
        cfg.assembly.min_length = int(a.get("min_length", cfg.assembly.min_length))

    # assembly_refine (optional iterative refinement flag)
    if "assembly" in raw:
        a = raw.get("assembly", {})
        cfg.assembly.iterative_refinement = bool(
            a.get("iterative_refinement", cfg.assembly.iterative_refinement)
        )

    # host_removal
    if "host_removal" in raw:
        hr = raw["host_removal"]
        cfg.host_removal.always_include_accessions = list(
            hr.get("always_include_accessions", [])
        )
        cfg.host_removal.kraken2_postfilter   = bool( hr.get("kraken2_postfilter",   cfg.host_removal.kraken2_postfilter))
        cfg.host_removal.postfilter_min_pct   = float(hr.get("postfilter_min_pct",   cfg.host_removal.postfilter_min_pct))
        cfg.host_removal.contamination_warn_pct = float(hr.get("contamination_warn_pct", cfg.host_removal.contamination_warn_pct))

    # genomad
    if "genomad" in raw:
        g = raw["genomad"]
        cfg.genomad.min_score    = float(g.get("min_score",    cfg.genomad.min_score))
        cfg.genomad.min_hallmarks = int(g.get(
            "min_virus_hallmarks",
            g.get("min_hallmarks", cfg.genomad.min_hallmarks),
        ))
        cfg.genomad.rescue_min_score     = float(g.get("rescue_min_score",     cfg.genomad.rescue_min_score))
        cfg.genomad.rescue_min_length_bp = int(  g.get("rescue_min_length_bp", cfg.genomad.rescue_min_length_bp))
        cfg.genomad.sensitivity          = float(g.get("sensitivity",          cfg.genomad.sensitivity))
        cfg.genomad.splits               = int(  g.get("splits",               cfg.genomad.splits))
        cfg.genomad.disable_find_proviruses  = bool(g.get("disable_find_proviruses",  cfg.genomad.disable_find_proviruses))
        cfg.genomad.enable_score_calibration = bool(g.get("enable_score_calibration", cfg.genomad.enable_score_calibration))
        cfg.genomad.composition   = g.get("composition",   cfg.genomad.composition)
        cfg.genomad.lenient_taxonomy = bool(g.get("lenient_taxonomy", cfg.genomad.lenient_taxonomy))

    # checkv
    if "checkv" in raw:
        c = raw["checkv"]
        cfg.checkv.min_completeness   = float(c.get("min_completeness",   cfg.checkv.min_completeness))
        cfg.checkv.min_viral_genes    = int(  c.get("min_viral_genes",    cfg.checkv.min_viral_genes))
        cfg.checkv.length_rescue      = int(  c.get("length_rescue",      cfg.checkv.length_rescue))
        cfg.checkv.min_contig_bp      = int(  c.get("min_contig_bp",      cfg.checkv.min_contig_bp))
        cfg.checkv.large_nd_rescue_bp = int(  c.get("large_nd_rescue_bp", cfg.checkv.large_nd_rescue_bp))
        cfg.checkv.min_bin_rescue_bp  = int(  c.get("min_bin_rescue_bp",  cfg.checkv.min_bin_rescue_bp))
        cfg.checkv.min_gene_density   = float(c.get("min_gene_density",   cfg.checkv.min_gene_density))
        cfg.checkv.max_kmer_freq      = float(c.get("max_kmer_freq",      cfg.checkv.max_kmer_freq))
        cfg.checkv.max_genome_copies  = float(c.get("max_genome_copies",  cfg.checkv.max_genome_copies))

    # annotate
    if "annotate" in raw:
        a = raw["annotate"]
        cfg.annotate.phold_gpu           = bool( a.get("phold_gpu",           cfg.annotate.phold_gpu))
        cfg.annotate.phold_finetune      = bool( a.get("phold_finetune",      cfg.annotate.phold_finetune))
        cfg.annotate.phold_batch_size    = int(  a.get("phold_batch_size",    cfg.annotate.phold_batch_size))
        cfg.annotate.phold_sensitivity   = float(a.get("phold_sensitivity",   cfg.annotate.phold_sensitivity))
        cfg.annotate.phynteny_confidence = float(a.get("phynteny_confidence", cfg.annotate.phynteny_confidence))
        cfg.annotate.genetic_code        = int(  a.get("genetic_code",        cfg.annotate.genetic_code))

    return cfg

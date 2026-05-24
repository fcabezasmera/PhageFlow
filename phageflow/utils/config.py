"""PhageFlow configuration loading."""

from __future__ import annotations
import math
import multiprocessing
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import yaml


# ---------------------------------------------------------------------------
# Thread auto-detection
# ---------------------------------------------------------------------------

def _default_threads() -> int:
    """90% of available logical CPUs, minimum 1."""
    return max(1, math.floor(multiprocessing.cpu_count() * 0.9))


# ---------------------------------------------------------------------------
# Sub-configs
# ---------------------------------------------------------------------------

@dataclass
class QcConfig:
    # fastp parameters — purified phage PE150
    # Chen 2018; Wick & Holt 2022; Roux 2019 (MIUViG)
    average_qual:               int  = 22
    qualified_quality_phred:    int  = 25
    unqualified_percent_limit:  int  = 10
    length_required:            int  = 80
    cut_right_window_size:      int  = 4
    cut_right_mean_quality:     int  = 25
    correction:                 bool = True
    overlap_len_require:        int  = 15
    overlap_diff_percent_limit: int  = 5
    low_complexity_filter:      bool = True
    complexity_threshold:       int  = 20
    trim_poly_x:                bool = True
    poly_x_min_len:             int  = 10
    n_base_limit:               int  = 5


@dataclass
class AssemblyConfig:
    kmers:      str  = "21,33,55,77,99,127,141"
    min_length: int  = 200
    iterative_refinement: bool = False


@dataclass
class HostRemovalConfig:
    always_include_accessions: list  = field(default_factory=list)
    kraken2_postfilter:        bool  = False
    postfilter_min_pct:        float = 1.0
    contamination_warn_pct:    float = 5.0


@dataclass
class GenomadConfig:
    min_score:            float = 0.7
    min_hallmarks:        int   = 0
    rescue_min_score:     float = 0.4
    rescue_min_length_bp: int   = 10_000
    sensitivity:          float = 4.2
    splits:               int   = 0
    disable_find_proviruses:    bool = True
    enable_score_calibration:   bool = True
    composition:          str   = "virome"
    lenient_taxonomy:     bool  = True


@dataclass
class CheckvConfig:
    min_completeness:    float = 50.0
    min_viral_genes:     int   = 1
    length_rescue:       int   = 10_000
    min_contig_bp:       int   = 1_500
    large_nd_rescue_bp:  int   = 30_000
    min_bin_rescue_bp:   int   = 30_000
    min_gene_density:    float = 0.5
    max_kmer_freq:       float = 1.5
    max_genome_copies:   float = 1.5


@dataclass
class AnnotateConfig:
    phold_gpu:           bool  = True
    phold_finetune:      bool  = True
    phold_batch_size:    int   = 1
    phold_sensitivity:   float = 9.5
    phynteny_confidence: float = 0.8
    genetic_code:        int   = 11


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
    project:   str  = "PhageFlow"
    version:   str  = "0.1.0"
    threads:   int  = field(default_factory=_default_threads)
    memory_gb: int  = 32
    workdir:   Path = Path(".")

    _output_dir:  Optional[Path] = field(default=None, repr=False, compare=False)
    _reports_dir: Optional[Path] = field(default=None, repr=False, compare=False)

    databases:    DatabaseConfig    = field(default_factory=DatabaseConfig)
    envs:         EnvConfig         = field(default_factory=EnvConfig)
    qc:           QcConfig          = field(default_factory=QcConfig)
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
        return self._output_dir or (self.workdir / "results")

    @property
    def reports_dir(self) -> Path:
        return self._reports_dir or (self.workdir / "reports")

    def set_output_dir(self, path: Optional[Path]) -> None:
        if path is not None:
            self._output_dir = Path(path)

    def set_reports_dir(self, path: Optional[Path]) -> None:
        if path is not None:
            self._reports_dir = Path(path)

    def results(self, step: str, sample_id: str = "") -> Path:
        """Get results directory for a step and optionally a sample.
        
        Examples:
            cfg.results("01_qc", "sample_001")  → results_dir/sample_001/01_qc/
            cfg.results("01_qc")                 → results_dir/01_qc/  (backward compat)
        """
        base = self.results_dir
        if sample_id:
            return base / sample_id / step
        return base / step

    def logs(self, sample_id: str) -> Path:
        """Get logs directory for a sample.
        
        Example:
            cfg.logs("sample_001")  → results_dir/sample_001/logs/
        """
        return self.results_dir / sample_id / "logs"

    def reports(self, step: str) -> Path:
        """Get reports directory for a step (deprecated in v0.2.0)."""
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

    if "threads" in raw:
        cfg.threads = int(raw["threads"])
    if "memory_gb" in raw:
        cfg.memory_gb = int(raw["memory_gb"])

    # Output / report directories
    if "dirs" in raw:
        d = raw["dirs"]
        if "results" in d:
            p = d["results"]
            cfg._output_dir = Path(p) if Path(p).is_absolute() else cfg.workdir / p
        if "reports" in d:
            p = d["reports"]
            cfg._reports_dir = Path(p) if Path(p).is_absolute() else cfg.workdir / p

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

    # qc
    if "qc" in raw:
        q = raw["qc"]
        cfg.qc.average_qual               = int(  q.get("average_qual",               cfg.qc.average_qual))
        cfg.qc.qualified_quality_phred    = int(  q.get("qualified_quality_phred",    cfg.qc.qualified_quality_phred))
        cfg.qc.unqualified_percent_limit  = int(  q.get("unqualified_percent_limit",  cfg.qc.unqualified_percent_limit))
        cfg.qc.length_required            = int(  q.get("length_required",            cfg.qc.length_required))
        cfg.qc.cut_right_window_size      = int(  q.get("cut_right_window_size",      cfg.qc.cut_right_window_size))
        cfg.qc.cut_right_mean_quality     = int(  q.get("cut_right_mean_quality",     cfg.qc.cut_right_mean_quality))
        cfg.qc.correction                 = bool( q.get("correction",                 cfg.qc.correction))
        cfg.qc.overlap_len_require        = int(  q.get("overlap_len_require",        cfg.qc.overlap_len_require))
        cfg.qc.overlap_diff_percent_limit = int(  q.get("overlap_diff_percent_limit", cfg.qc.overlap_diff_percent_limit))
        cfg.qc.low_complexity_filter      = bool( q.get("low_complexity_filter",      cfg.qc.low_complexity_filter))
        cfg.qc.complexity_threshold       = int(  q.get("complexity_threshold",       cfg.qc.complexity_threshold))
        cfg.qc.trim_poly_x                = bool( q.get("trim_poly_x",               cfg.qc.trim_poly_x))
        cfg.qc.poly_x_min_len             = int(  q.get("poly_x_min_len",             cfg.qc.poly_x_min_len))
        cfg.qc.n_base_limit               = int(  q.get("n_base_limit",               cfg.qc.n_base_limit))

    # assembly
    if "assembly" in raw:
        a = raw["assembly"]
        cfg.assembly.kmers               = a.get("kmers",      cfg.assembly.kmers)
        cfg.assembly.min_length          = int(a.get("min_length", cfg.assembly.min_length))
        cfg.assembly.iterative_refinement = bool(a.get("iterative_refinement", cfg.assembly.iterative_refinement))

    # host_removal
    if "host_removal" in raw:
        hr = raw["host_removal"]
        cfg.host_removal.always_include_accessions = list(hr.get("always_include_accessions", []))
        cfg.host_removal.kraken2_postfilter        = bool(hr.get("kraken2_postfilter",   cfg.host_removal.kraken2_postfilter))
        cfg.host_removal.postfilter_min_pct        = float(hr.get("postfilter_min_pct",  cfg.host_removal.postfilter_min_pct))
        cfg.host_removal.contamination_warn_pct    = float(hr.get("contamination_warn_pct", cfg.host_removal.contamination_warn_pct))

    # genomad
    if "genomad" in raw:
        g = raw["genomad"]
        cfg.genomad.min_score            = float(g.get("min_score",    cfg.genomad.min_score))
        cfg.genomad.min_hallmarks        = int(g.get("min_virus_hallmarks", g.get("min_hallmarks", cfg.genomad.min_hallmarks)))
        cfg.genomad.rescue_min_score     = float(g.get("rescue_min_score",     cfg.genomad.rescue_min_score))
        cfg.genomad.rescue_min_length_bp = int(  g.get("rescue_min_length_bp", cfg.genomad.rescue_min_length_bp))
        cfg.genomad.sensitivity          = float(g.get("sensitivity",          cfg.genomad.sensitivity))
        cfg.genomad.splits               = int(  g.get("splits",               cfg.genomad.splits))
        cfg.genomad.disable_find_proviruses  = bool(g.get("disable_find_proviruses",  cfg.genomad.disable_find_proviruses))
        cfg.genomad.enable_score_calibration = bool(g.get("enable_score_calibration", cfg.genomad.enable_score_calibration))
        cfg.genomad.composition              = g.get("composition",   cfg.genomad.composition)
        cfg.genomad.lenient_taxonomy         = bool(g.get("lenient_taxonomy", cfg.genomad.lenient_taxonomy))

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

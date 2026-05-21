"""PhageFlow configuration and sample loading (purified phage mode)."""

from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional
import yaml


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

    def validate_reads(self) -> List[str]:
        return self.validate()


def load_samples(tsv_path: Path) -> List[Sample]:
    """Load samples.tsv (columns: sample_id, r1, r2)."""
    samples = []
    with open(tsv_path) as f:
        header = None
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


@dataclass
class AssemblyConfig:
    kmers:      str = "21,33,55,77,99,127"
    min_length: int = 200


@dataclass
class GenomadConfig:
    min_score:            float = 0.7
    min_hallmarks:        int   = 0
    rescue_min_score:     float = 0.4   # NEW: floor score for length-based rescue
    rescue_min_length_bp: int   = 10_000  # NEW: min length (bp) for borderline rescue



@dataclass
class CheckvConfig:
    min_completeness:   float = 50.0
    min_viral_genes:    int   = 1
    length_rescue:      int   = 10_000
    min_contig_bp:      int   = 1_500
    large_nd_rescue_bp: int   = 30_000
    min_bin_rescue_bp:  int   = 30_000
    min_gene_density:   float = 0.5


@dataclass
class AnnotateConfig:
    """Settings for Module 06 (Pharokka → Phold → Phynteny)."""
    phold_gpu:           bool  = True   # --foldseek_gpu (False → --cpu)
    phold_finetune:      bool  = True   # --finetune: phage-finetuned ProstT5
    phold_batch_size:    int   = 1      # ProstT5 batch size (increase on GPU)
    phynteny_confidence: float = 0.7   # confidence threshold (0.0–1.0)


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


@dataclass
class Config:
    project:      str   = "PhageFlow"
    version:      str   = "0.1.0"
    mode:         str   = "purified_phage"
    threads:      int   = 8
    memory_gb:    int   = 32
    samples_file: Path  = Path("config/samples.tsv")
    workdir:      Path  = Path(".")

    databases: DatabaseConfig = field(default_factory=DatabaseConfig)
    envs:      EnvConfig      = field(default_factory=EnvConfig)
    assembly:  AssemblyConfig = field(default_factory=AssemblyConfig)
    genomad:   GenomadConfig  = field(default_factory=GenomadConfig)
    checkv:    CheckvConfig   = field(default_factory=CheckvConfig)
    annotate:  AnnotateConfig = field(default_factory=AnnotateConfig)

    @property
    def results_dir(self) -> Path:
        return self.workdir / "results"

    @property
    def reports_dir(self) -> Path:
        return self.workdir / "reports"

    def results(self, step: str) -> Path:
        return self.results_dir / step

    def reports(self, step: str) -> Path:
        return self.reports_dir / step


def load_config(config_path: Path, workdir: Optional[Path] = None) -> Config:
    with open(config_path) as f:
        raw = yaml.safe_load(f)

    cfg         = Config()
    cfg.workdir = workdir or config_path.parent.parent
    cfg.project   = raw.get("project",   cfg.project)
    cfg.mode      = raw.get("mode",      cfg.mode)
    cfg.threads   = int(raw.get("threads",   cfg.threads))
    cfg.memory_gb = int(raw.get("memory_gb", cfg.memory_gb))

    sf = raw.get("samples_file", str(cfg.samples_file))
    cfg.samples_file = cfg.workdir / sf if not Path(sf).is_absolute() else Path(sf)

    if "databases" in raw:
        db = raw["databases"]
        cfg.databases.checkv   = cfg.workdir / db.get("checkv",   str(cfg.databases.checkv))
        cfg.databases.pharokka = cfg.workdir / db.get("pharokka", str(cfg.databases.pharokka))
        cfg.databases.genomad  = cfg.workdir / db.get("genomad",  str(cfg.databases.genomad))
        cfg.databases.phold    = cfg.workdir / db.get("phold",    str(cfg.databases.phold))
        cfg.databases.phynteny = cfg.workdir / db.get("phynteny", str(cfg.databases.phynteny))
        if "kraken2" in db:
            cfg.databases.kraken2 = cfg.workdir / db["kraken2"]

    if "envs" in raw:
        e = raw["envs"]
        cfg.envs.main = e.get("main", cfg.envs.main)

    if "assembly" in raw:
        a = raw["assembly"]
        cfg.assembly.kmers      = a.get("kmers",      cfg.assembly.kmers)
        cfg.assembly.min_length = int(a.get("min_length", cfg.assembly.min_length))

    if "genomad" in raw:
        g = raw["genomad"]
        cfg.genomad.min_score     = float(g.get("min_score", cfg.genomad.min_score))
        cfg.genomad.min_hallmarks = int(
            g.get("min_virus_hallmarks", g.get("min_hallmarks", cfg.genomad.min_hallmarks))
        )
        
        cfg.genomad.rescue_min_score     = float(
            g.get("rescue_min_score",     cfg.genomad.rescue_min_score)
        )
        cfg.genomad.rescue_min_length_bp = int(
            g.get("rescue_min_length_bp", cfg.genomad.rescue_min_length_bp)
        )
 
    if "checkv" in raw:
        c = raw["checkv"]
        cfg.checkv.min_completeness   = float(c.get("min_completeness",   cfg.checkv.min_completeness))
        cfg.checkv.min_viral_genes    = int(  c.get("min_viral_genes",    cfg.checkv.min_viral_genes))
        cfg.checkv.length_rescue      = int(  c.get("length_rescue",      cfg.checkv.length_rescue))
        cfg.checkv.min_contig_bp      = int(  c.get("min_contig_bp",      cfg.checkv.min_contig_bp))
        cfg.checkv.large_nd_rescue_bp = int(  c.get("large_nd_rescue_bp", cfg.checkv.large_nd_rescue_bp))
        cfg.checkv.min_bin_rescue_bp  = int(  c.get("min_bin_rescue_bp",  cfg.checkv.min_bin_rescue_bp))
        cfg.checkv.min_gene_density   = float(c.get("min_gene_density",   cfg.checkv.min_gene_density))

    if "annotate" in raw:
        a = raw["annotate"]
        cfg.annotate.phold_gpu           = bool( a.get("phold_gpu",           cfg.annotate.phold_gpu))
        cfg.annotate.phold_finetune      = bool( a.get("phold_finetune",      cfg.annotate.phold_finetune))
        cfg.annotate.phold_batch_size    = int(  a.get("phold_batch_size",    cfg.annotate.phold_batch_size))
        cfg.annotate.phynteny_confidence = float(a.get("phynteny_confidence", cfg.annotate.phynteny_confidence))

    return cfg

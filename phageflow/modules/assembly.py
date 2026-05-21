"""PhageFlow Module 03 — De novo assembly.

Goal: produce complete or near-complete phage genome contigs for downstream
CheckV Complete/High-quality classification (Nayfach et al. 2021).

Tools and parameters (literature-based for purified phage sequencing):

SPAdes v3.15+ (Bankevich et al. 2012, Journal of Computational Biology 19:455):
    --isolate
        Isolate mode: optimised for high-coverage single-organism preparations.
        Avoids the variable-coverage heuristics in --meta (designed for
        metagenomes spanning 0.1×–10 000×) that can fragment the De Bruijn
        graph in uniform high-coverage phage data.
        (Prjibelski et al. 2020, Current Protocols; Bankevich et al. 2012)

    --only-assembler
        Skip BayesHammer error correction — required at >200× coverage
        (BayesHammer collapses true variants into consensus, losing
        biological diversity in tail fiber / receptor binding proteins).
        Compatible with --isolate (SPAdes ≥3.15).
        (Roux et al. 2019, eLife 8:e42923)

    -k 21,33,55,77,99,127
        Multiple k-mers ensure coverage of both short repeats (k=21)
        and long unique regions (k=127). Standard for PE150 phage data
        (Bankevich et al. 2012; SPAdes manual).

    --s1 (singletons, optional)
        Unpaired reads from host-removal step (reads whose pair aligned to host
        but which are themselves unmapped). These represent phage fragments at
        DTR/ITR boundaries; including them improves terminal coverage and raises
        the probability of CheckV 'Complete' classification
        (Nayfach et al. 2021, Nat Biotechnol 39:578).

    Post-assembly length filter: contigs < min_length (default 200 bp) discarded.

    Note: --cov-cutoff auto is INCOMPATIBLE with --isolate (SPAdes ≥3.15)
    Note: contigs.fasta used (NOT scaffolds.fasta) — scaffolds introduce
          artificial N-gaps that break ORF prediction in Pharokka/Phold.

MEGAHIT v1.2+ (Li et al. 2015, Bioinformatics 31:1674):
    --k-list {cfg.assembly.kmers}
        Reuses the same k-mer list from config.yaml (default 21,33,55,77,99,127)
        for reproducibility. Previously hardcoded; now synced with SPAdes.

    --min-count 2
        Minimum solid k-mer multiplicity — removes sequencing errors
        without losing low-coverage genuine phage k-mers.

    --no-mercy
        Disable mercy k-mers at high coverage (>200×).
        (Li et al. 2015)

    --min-contig-len {cfg.assembly.min_length}
        Consistent with SPAdes post-filter (default 200 bp).

cd-hit-est (Fu et al. 2012, Bioinformatics 28:3150):
    -c 1.00 : 100% sequence identity — removes exact duplicates only,
              preserving biological SNP-level variants.
    -aS 0.85 : ≥85% alignment coverage of the shorter sequence.
               Changed from 1.00: with -aS 1.00, a MEGAHIT contig that is
               95% the length of an equivalent SPAdes contig (same region,
               100% identity) would NOT be dereplicated, creating functional
               duplicates in downstream modules. -aS 0.85 eliminates these
               while preserving genuine biological variants.
               (Fu et al. 2012 — recommended for variable-length sequences)
"""

from __future__ import annotations
import shutil
from pathlib import Path
from typing import List, Optional

from rich.panel import Panel
from rich.text import Text
from rich.progress import (
    Progress, SpinnerColumn, TextColumn,
    BarColumn, MofNCompleteColumn, TimeElapsedColumn,
)

from phageflow.utils.config import Config
from phageflow.utils.logger import log_step, log_info, log_ok, log_warn, console
from phageflow.utils.tools import require_tools, run_silent, human_size, mkdirs, fasta_stats

STEP  = "03_assembly"
TOOLS = ["spades.py", "megahit", "cd-hit-est"]

# ── Thresholds ────────────────────────────────────────────────────────────────
# For purified phage the goal is 1–5 contigs covering the full genome.
_N50_GOOD         = 50_000  # N50 ≥ 50 kb: likely near-complete assembly
_N50_WARN         = 20_000  # N50 < 20 kb: fragmented; CheckV Complete unlikely
_N50_CRITICAL     = 5_000   # N50 < 5 kb:  severely fragmented; annotation unreliable
_MAX_CONTIGS_WARN = 20      # >20 NR contigs unusual for purified phage
_MAX_CONTIGS_CRIT = 100     # >100 NR contigs suggests assembly failure or contamination
_MIN_CONTIGS      = 1


# ── Public entry point ────────────────────────────────────────────────────────

def run(
    cfg:       Config,
    sample_id: str,
    r1:        Path,
    r2:        Path,
    s1:        Optional[Path] = None,
    force:     bool           = False,
) -> Path:
    """
    Assemble a single sample with SPAdes (--isolate) + MEGAHIT, then dereplicate.

    Parameters
    ----------
    s1 : singletons FASTQ from host-removal (passed to SPAdes --s1).
         Improves terminal DTR/ITR coverage → raises CheckV 'Complete' rate
         (Nayfach et al. 2021). Pass None to skip.

    Returns
    -------
    nr_fa : Path to non-redundant contigs FASTA (input for geNomad)
    """
    require_tools(*TOOLS)

    r1 = Path(r1); r2 = Path(r2)
    out_dir = cfg.results(STEP)
    rpt_dir = cfg.reports(STEP)
    mkdirs(out_dir / "spades", out_dir / "megahit", out_dir / "combined", rpt_dir)

    nr_fa = out_dir / "combined" / f"{sample_id}_contigs_nr.fasta"

    log_step(f"Module 03 — assembly [{sample_id}]")
    _print_input_table(sample_id, r1, r2, s1)
    log_info(f"  SPAdes   : --isolate --only-assembler  (Prjibelski et al. 2020; Roux et al. 2019)")
    log_info(f"  MEGAHIT  : --no-mercy --min-count 2    (Li et al. 2015)")
    log_info(f"  Min len  : {cfg.assembly.min_length} bp  |  k-mers: {cfg.assembly.kmers}")
    log_info(f"  cd-hit-est: -c 1.00 -aS 0.85  (Fu et al. 2012)")
    if s1 and Path(s1).exists() and Path(s1).stat().st_size > 0:
        log_info(f"  Singletons: {s1}  ({human_size(Path(s1))}) → SPAdes --s1")
    else:
        log_info("  Singletons: none provided  (pass --s1 to include DTR/ITR boundary reads)")

    # ── Pipeline steps ────────────────────────────────────────────────────────
    sp_fa: Optional[Path] = None
    mh_fa: Optional[Path] = None

    with Progress(
        SpinnerColumn(spinner_name="dots2"),
        TextColumn("[bold cyan]{task.description:<48}"),
        BarColumn(bar_width=20, complete_style="cyan", finished_style="green"),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        console=console,
        transient=True,
    ) as progress:
        task = progress.add_task("Initializing...", total=3)

        # 1/3 — SPAdes + length filter
        progress.update(task, description="[1/3] SPAdes    — de novo assembly (--isolate)")
        sp_fa = _run_spades(cfg, sample_id, r1, r2, s1, out_dir, rpt_dir, force)
        progress.advance(task)

        # 2/3 — MEGAHIT
        progress.update(task, description="[2/3] MEGAHIT   — de novo assembly")
        mh_fa = _run_megahit(cfg, sample_id, r1, r2, out_dir, rpt_dir, force)
        progress.advance(task)

        # 3/3 — cd-hit-est dereplication
        progress.update(task, description="[3/3] cd-hit-est — dereplication")
        nr_fa = _run_cdhit(cfg, sample_id, sp_fa, mh_fa, out_dir, rpt_dir, force)
        progress.advance(task)

    # ── Metrics, display, save ────────────────────────────────────────────────
    stats = _collect_stats(sp_fa, mh_fa, nr_fa)
    _validate_output(sample_id, nr_fa, stats)
    _save_tsv(sample_id, sp_fa, mh_fa, nr_fa, rpt_dir / "assembly_summary.tsv")
    _print_summary_table(sample_id, stats)
    _check_warnings(stats)
    _print_completion_panel(sample_id, nr_fa, rpt_dir, stats)

    log_step(f"Module 03 completed ✓  [{sample_id}]")
    log_info(
        f"  Next: phageflow viral-id --sample-id {sample_id} "
        f"--contigs {nr_fa}"
    )
    return nr_fa


# ── Step implementations ──────────────────────────────────────────────────────

def _run_spades(
    cfg:       Config,
    sample_id: str,
    r1:        Path,
    r2:        Path,
    s1:        Optional[Path],
    out_dir:   Path,
    rpt_dir:   Path,
    force:     bool,
) -> Path:
    """
    Run SPAdes in --isolate mode with optional singleton support.

    --isolate is preferred over --meta for purified phage because:
      - --meta applies variable-coverage heuristics (designed for 0.1×–10 000×
        metagenomes) that can fragment the De Bruijn graph at the uniformly
        high coverage typical of phage preparations (Prjibelski et al. 2020).
      - --isolate uses an optimised graph simplification for single-organism
        high-coverage data, producing longer contigs.

    --only-assembler skips BayesHammer error correction, which is
    counterproductive at >200× (collapses biological variants).
    Compatible with --isolate in SPAdes ≥3.15 (Roux et al. 2019).

    --s1 passes singleton reads from host removal as unpaired input.
    SPAdes incorporates them in the De Bruijn graph, improving coverage
    at DTR/ITR boundaries (Nayfach et al. 2021).
    """
    sp_dir = out_dir / "spades" / sample_id
    sp_raw = sp_dir / "contigs.fasta"
    sp_fa  = sp_dir / "contigs_filtered.fasta"
    mkdirs(sp_dir)

    if sp_fa.exists() and sp_fa.stat().st_size > 0 and not force:
        log_info("  [SPAdes] already exists — skipping")
        return sp_fa

    cmd = [
        "spades.py",
        "--pe1-1", str(r1), "--pe1-2", str(r2),
        "--isolate",          # ← purified phage high-coverage mode
        "--only-assembler",   # ← skip error correction at >200× (Roux et al. 2019)
        "-k",        cfg.assembly.kmers,
        "--threads", str(cfg.threads),
        "--memory",  str(cfg.memory_gb),
        "-o",        str(sp_dir),
    ]

    # Singletons from host-removal → improve DTR/ITR boundary coverage
    s1_path = Path(s1) if s1 else None
    if s1_path and s1_path.exists() and s1_path.stat().st_size > 0:
        cmd += ["--s1", str(s1_path)]
        log_info(
            f"  [SPAdes] --s1 {s1_path.name} included "
            f"({human_size(s1_path)}) — DTR/ITR boundary coverage "
            f"(Nayfach et al. 2021)"
        )
    else:
        log_info("  [SPAdes] --isolate --only-assembler | no singletons")

    try:
        run_silent(cmd, log_file=rpt_dir / f"{sample_id}_spades.log")
    except Exception as e:
        log_warn(f"  [SPAdes] warning: {e}")

    if sp_raw.exists() and sp_raw.stat().st_size > 0:
        n_raw, n_kept = _filter_fasta(sp_raw, sp_fa, cfg.assembly.min_length)
        log_ok(
            f"  [SPAdes] length filter ≥{cfg.assembly.min_length}bp: "
            f"{n_raw} → {n_kept} contigs  (discarded {n_raw - n_kept})"
        )
    else:
        log_warn("  [SPAdes] no assembly output — check log")

    return sp_fa


def _run_megahit(
    cfg:       Config,
    sample_id: str,
    r1:        Path,
    r2:        Path,
    out_dir:   Path,
    rpt_dir:   Path,
    force:     bool,
) -> Path:
    """
    Run MEGAHIT assembly.

    k-list is read from cfg.assembly.kmers (same as SPAdes) for reproducibility.
    Previously hardcoded as "21,29,39,59,79,99,119,141", which diverged from
    the config value and was not user-configurable.

    MEGAHIT and SPAdes use the same comma-separated k-mer format; the config
    default "21,33,55,77,99,127" is valid for both tools.
    (Li et al. 2015, Bioinformatics 31:1674)
    """
    mh_dir = out_dir / "megahit" / sample_id
    mh_fa  = mh_dir / "final.contigs.fa"

    if mh_fa.exists() and mh_fa.stat().st_size > 0 and not force:
        log_info("  [MEGAHIT] already exists — skipping")
        return mh_fa

    if mh_dir.exists():
        shutil.rmtree(mh_dir)
    mh_dir.parent.mkdir(parents=True, exist_ok=True)

    try:
        run_silent([
            "megahit",
            "-1", str(r1), "-2", str(r2),
            "-o",               str(mh_dir),
            "--k-list",         cfg.assembly.kmers,   # ← synced with config (was hardcoded)
            "--min-count",      "2",
            "--no-mercy",
            "--min-contig-len", str(cfg.assembly.min_length),
            "-t",               str(cfg.threads),
        ], log_file=rpt_dir / f"{sample_id}_megahit.log")
        log_ok("  [MEGAHIT] assembly complete")
    except Exception as e:
        log_warn(f"  [MEGAHIT] warning: {e}")

    return mh_fa


def _run_cdhit(
    cfg:       Config,
    sample_id: str,
    sp_fa:     Optional[Path],
    mh_fa:     Optional[Path],
    out_dir:   Path,
    rpt_dir:   Path,
    force:     bool,
) -> Path:
    """
    Merge SPAdes + MEGAHIT contigs and dereplicate.

    Parameters
    ----------
    -c 1.00  : 100% sequence identity — only exact duplicates collapsed.
    -aS 0.85 : ≥85% alignment coverage of the shorter sequence.

    Rationale for -aS 0.85 (changed from 1.00):
        With -aS 1.00, two contigs must be nearly identical in length to be
        considered duplicates. A MEGAHIT contig covering 95% of an equivalent
        SPAdes contig (same genomic region, 100% identity) would NOT be
        dereplicated, producing functional duplicate candidates in downstream
        CheckV/geNomad modules. -aS 0.85 removes these while preserving
        genuine biological variants (SNPs, indels).
        (Fu et al. 2012, Bioinformatics 28:3150 — recommended for
        variable-length sequence dereplication)
    """
    nr_fa  = out_dir / "combined" / f"{sample_id}_contigs_nr.fasta"
    tmp_fa = out_dir / "combined" / f"{sample_id}_all_tmp.fasta"

    if nr_fa.exists() and nr_fa.stat().st_size > 0 and not force:
        log_info("  \\[cd-hit-est] already exists — skipping")
        return nr_fa

    n_sp = n_mh = 0
    with open(tmp_fa, "w") as out:
        if sp_fa and sp_fa.exists():
            with open(sp_fa) as f:
                for line in f:
                    if line.startswith(">"):
                        out.write(f">{sample_id}_sp_{line[1:]}")
                        n_sp += 1
                    else:
                        out.write(line)
        if mh_fa and mh_fa.exists():
            with open(mh_fa) as f:
                for line in f:
                    if line.startswith(">"):
                        out.write(f">{sample_id}_mh_{line[1:]}")
                        n_mh += 1
                    else:
                        out.write(line)

    n_input = n_sp + n_mh
    if n_input == 0:
        log_warn("  [cd-hit-est] no contigs from either assembler")
        return nr_fa

    try:
        run_silent([
            "cd-hit-est",
            "-i", str(tmp_fa), "-o", str(nr_fa),
            "-c", "1.00",
            "-aS", "0.85",    # ← was 1.00; correct for variable-length dedup
            "-n", "10", "-d", "0",
            "-T", str(cfg.threads),
            "-M", str(cfg.memory_gb * 1000),
        ], log_file=rpt_dir / f"{sample_id}_cdhit.log")
        log_ok("  \\[cd-hit-est] dereplication complete")
    except Exception as e:
        log_warn(f"  \\[cd-hit-est] warning: {e}")
    finally:
        tmp_fa.unlink(missing_ok=True)
        Path(str(nr_fa) + ".clstr").unlink(missing_ok=True)

    return nr_fa


# ── Validation ────────────────────────────────────────────────────────────────

def _validate_output(sample_id: str, nr_fa: Path, stats: dict) -> None:
    """Warn if the NR FASTA is missing or empty after assembly."""
    if not nr_fa.exists() or nr_fa.stat().st_size == 0:
        log_warn(f"  validate · NR contigs file missing or empty: {nr_fa}")
        return

    n = stats.get("nr_n", 0)
    if n < _MIN_CONTIGS:
        log_warn(
            f"  validate · 0 NR contigs assembled — "
            "check read quality and host removal step."
        )
    else:
        log_ok(f"  validate · {n} NR contigs — OK")


# ── Stats collection ──────────────────────────────────────────────────────────

def _collect_stats(
    sp_fa: Optional[Path],
    mh_fa: Optional[Path],
    nr_fa: Optional[Path],
) -> dict:
    """Gather fasta_stats for each assembler output into a flat dict."""
    def _s(fa):
        if fa and fa.exists() and fa.stat().st_size > 0:
            return fasta_stats(fa)
        return {"n": 0, "total_bp": 0, "largest_bp": 0, "n50": 0}

    sp = _s(sp_fa)
    mh = _s(mh_fa)
    nr = _s(nr_fa)

    n_input = sp["n"] + mh["n"]
    redund  = n_input - nr["n"]
    pct_red = f"{redund / n_input * 100:.1f}%" if n_input else "N/A"

    return {
        "sp_n":          sp["n"],
        "sp_largest":    sp["largest_bp"],
        "sp_n50":        sp["n50"],
        "sp_total":      sp["total_bp"],
        "mh_n":          mh["n"],
        "mh_largest":    mh["largest_bp"],
        "mh_n50":        mh["n50"],
        "mh_total":      mh["total_bp"],
        "nr_n":          nr["n"],
        "nr_largest":    nr["largest_bp"],
        "nr_n50":        nr["n50"],
        "nr_total":      nr["total_bp"],
        "pct_redundant": pct_red,
    }


# ── Helpers ───────────────────────────────────────────────────────────────────

def _filter_fasta(src: Path, dst: Path, min_len: int) -> tuple[int, int]:
    """Write sequences from src to dst keeping only those >= min_len bp."""
    seqs: List[tuple[str, str]] = []
    header = seq = ""

    with open(src) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header:
                    seqs.append((header, seq))
                header, seq = line, ""
            else:
                seq += line
    if header:
        seqs.append((header, seq))

    kept = [(h, s) for h, s in seqs if len(s) >= min_len]

    with open(dst, "w") as f:
        for h, s in kept:
            f.write(h + "\n" + s + "\n")

    return len(seqs), len(kept)


# ── Rich display helpers ──────────────────────────────────────────────────────

def _fmt(val, unit: str = "") -> str:
    """Format metric value; dim 'N/A' or missing entries."""
    if val in ("N/A", "", None):
        return "[dim]N/A[/dim]"
    return f"{val}{unit}"


def _color_n50(val) -> str:
    """Color N50 value by completeness expectation for purified phage."""
    try:
        v = int(val)
        if v >= _N50_GOOD:     return f"[bold green]{v:,}bp[/bold green]"
        if v >= _N50_WARN:     return f"[bold yellow]{v:,}bp[/bold yellow]"
        return f"[bold red]{v:,}bp[/bold red]"
    except (ValueError, TypeError):
        return "[dim]N/A[/dim]"


def _print_input_table(
    sample_id: str, r1: Path, r2: Path, s1: Optional[Path]
) -> None:
    log_info(f"  R1 : {r1}  ({human_size(r1)})")
    log_info(f"  R2 : {r2}  ({human_size(r2)})")
    if s1 and Path(s1).exists():
        log_info(f"  S1 : {s1}  ({human_size(Path(s1))})  [singletons]")


def _print_summary_table(sample_id: str, stats: dict) -> None:
    """Print assembly metrics as compact log lines."""
    log_ok(
        f"  SPAdes   : {_fmt(stats['sp_n'])} contigs | "
        f"largest={_fmt(stats['sp_largest'], 'bp')} | "
        f"N50={_color_n50(stats['sp_n50'])} | "
        f"total={_fmt(stats['sp_total'], 'bp')}"
    )
    log_ok(
        f"  MEGAHIT  : {_fmt(stats['mh_n'])} contigs | "
        f"largest={_fmt(stats['mh_largest'], 'bp')} | "
        f"N50={_color_n50(stats['mh_n50'])} | "
        f"total={_fmt(stats['mh_total'], 'bp')}"
    )
    log_ok(
        f"  NR final : {_fmt(stats['nr_n'])} contigs | "
        f"largest={_fmt(stats['nr_largest'], 'bp')} | "
        f"N50={_color_n50(stats['nr_n50'])} | "
        f"−{stats['pct_redundant']} redundant"
    )


def _check_warnings(stats: dict) -> None:
    """Emit contextual warnings based on assembly metrics."""
    if stats["nr_n"] == 0:
        log_warn(
            "  No contigs assembled — check read quality and host removal step."
        )
        return

    try:
        n50 = int(stats["nr_n50"])
        if n50 < _N50_CRITICAL:
            log_warn(
                f"  NR N50 = {n50:,}bp < {_N50_CRITICAL:,}bp — "
                "severely fragmented assembly; Complete/HQ classification in "
                "CheckV is unlikely. Check host removal and read quality."
            )
        elif n50 < _N50_WARN:
            log_warn(
                f"  NR N50 = {n50:,}bp < {_N50_WARN:,}bp — "
                "fragmented assembly; CheckV 'Complete' classification requires "
                "contigs spanning the full genome with DTR/ITR termini. "
                "Consider increasing coverage or checking for contamination."
            )
    except (TypeError, ValueError):
        pass

    try:
        nr_n = int(stats["nr_n"])
        if nr_n > _MAX_CONTIGS_CRIT:
            log_warn(
                f"  {nr_n} NR contigs — assembly failure or heavy contamination; "
                "verify host-removal step and sample purity."
            )
        elif nr_n > _MAX_CONTIGS_WARN:
            log_warn(
                f"  {nr_n} NR contigs — unusually high for purified phage "
                "(expected 1–{_MAX_CONTIGS_WARN}); "
                "check for contamination or consider stricter host removal."
            )
    except (TypeError, ValueError):
        pass


def _print_completion_panel(
    sample_id: str, nr_fa: Path, rpt_dir: Path, stats: dict
) -> None:
    sp_n = stats.get("sp_n", "?")
    mh_n = stats.get("mh_n", "?")
    nr_n = stats.get("nr_n", "?")
    pct  = stats.get("pct_redundant", "?")
    n50  = stats.get("nr_n50", 0)

    text = Text()
    text.append("✓ ", style="bold green")
    text.append(f"SPAdes={sp_n} + MEGAHIT={mh_n} contigs  →  ", style="dim white")
    text.append(f"{nr_n}", style="bold green")
    text.append(f" NR contigs  (−{pct} redundant)\n", style="cyan")
    try:
        n50_int = int(n50)
        style = (
            "bold green"  if n50_int >= _N50_GOOD else
            "bold yellow" if n50_int >= _N50_WARN else
            "bold red"
        )
        text.append(f"  NR N50 = ", style="dim white")
        text.append(f"{n50_int:,} bp\n\n", style=style)
    except (TypeError, ValueError):
        text.append("\n")
    text.append("NR contigs : ", style="dim white")
    text.append(str(nr_fa) + "\n", style="white")
    text.append("Summary    : ", style="dim white")
    text.append(str(rpt_dir / "assembly_summary.tsv"), style="white")

    console.print(Panel(
        text,
        title=f"[bold cyan]Assembly complete — {sample_id}[/bold cyan]",
        border_style="cyan",
        padding=(0, 2),
        width=90,
    ))


# ── TSV summary ───────────────────────────────────────────────────────────────

_TSV_HEADERS = ["sample", "assembler", "contigs", "largest_bp", "n50_bp", "total_bp"]


def _save_tsv(
    sample_id: str,
    sp_fa:     Optional[Path],
    mh_fa:     Optional[Path],
    nr_fa:     Optional[Path],
    path:      Path,
) -> None:
    """Write / update the assembly summary TSV."""
    def _s(fa):
        if fa and fa.exists() and fa.stat().st_size > 0:
            s = fasta_stats(fa)
            return str(s["n"]), str(s["largest_bp"]), str(s["n50"]), str(s["total_bp"])
        return "0", "0", "0", "0"

    rows: dict[str, dict] = {}
    if path.exists():
        with open(path) as f:
            old_hdrs = f.readline().rstrip("\n").split("\t")
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if len(cols) >= 2:
                    rows[f"{cols[0]}_{cols[1]}"] = dict(zip(old_hdrs, cols))

    for asm, fa in [("spades", sp_fa), ("megahit", mh_fa), ("NR(cdhit)", nr_fa)]:
        n, lg, n50, tot = _s(fa)
        key = f"{sample_id}_{asm}"
        rows[key] = dict(zip(_TSV_HEADERS, [sample_id, asm, n, lg, n50, tot]))

    with open(path, "w") as f:
        f.write("\t".join(_TSV_HEADERS) + "\n")
        for row in rows.values():
            f.write("\t".join(row.get(h, "") for h in _TSV_HEADERS) + "\n")

"""PhageFlow external tool execution utilities."""

from __future__ import annotations
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Union

from .logger import log_info, log_ok, log_warn, log_error


# ---------------------------------------------------------------------------
# Tool validation
# ---------------------------------------------------------------------------

def check_tool(name: str, version_flag: str = "--version") -> bool:
    """Check if an external tool is available in PATH."""
    if shutil.which(name) is None:
        return False
    return True


def require_tools(*tools: str) -> None:
    """Exit with error if any required tool is missing."""
    missing = [t for t in tools if not check_tool(t)]
    if missing:
        log_error(f"Required tools not found in PATH: {', '.join(missing)}")
        log_error("Please activate the correct conda environment.")
        sys.exit(1)


def get_tool_version(name: str, flag: str = "--version") -> str:
    """Get version string of an external tool."""
    try:
        result = subprocess.run(
            [name, flag],
            capture_output=True, text=True, timeout=10
        )
        out = (result.stdout + result.stderr).strip().split("\n")[0]
        return out
    except Exception:
        return "unknown"


# ---------------------------------------------------------------------------
# Command execution
# ---------------------------------------------------------------------------

def run(
    cmd: Union[str, List[str]],
    log_file: Optional[Path] = None,
    check: bool = True,
    shell: bool = False,
    cwd: Optional[Path] = None,
    env: Optional[dict] = None,
) -> subprocess.CompletedProcess:
    """
    Run an external command, optionally logging stdout/stderr to a file.

    Parameters
    ----------
    cmd       : command as list (preferred) or string (shell=True)
    log_file  : path to append stdout+stderr
    check     : raise CalledProcessError on non-zero exit
    shell     : run via shell (avoid when possible)
    cwd       : working directory
    env       : environment variables override
    """
    if isinstance(cmd, list):
        cmd_str = " ".join(str(c) for c in cmd)
    else:
        cmd_str = cmd

    log_info(f"  $ {cmd_str}")

    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        with open(log_file, "a") as lf:
            lf.write(f"\n$ {cmd_str}\n")
            result = subprocess.run(
                cmd, stdout=lf, stderr=lf,
                shell=shell, cwd=cwd, env=env,
            )
    else:
        result = subprocess.run(
            cmd, capture_output=False,
            shell=shell, cwd=cwd, env=env,
        )

    if check and result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, cmd)

    return result


def run_silent(
    cmd: Union[str, List[str]],
    log_file: Optional[Path] = None,
    check: bool = True,
    cwd: Optional[Path] = None,
    shell: bool = False,
) -> subprocess.CompletedProcess:
    """Run command silently (capture all output)."""
    if isinstance(cmd, list):
        cmd_str = " ".join(str(c) for c in cmd)
    else:
        cmd_str = cmd

    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        with open(log_file, "a") as lf:
            lf.write(f"\n$ {cmd_str}\n")
            result = subprocess.run(
                cmd, stdout=lf, stderr=lf, cwd=cwd, shell=shell
            )
    else:
        result = subprocess.run(
            cmd, capture_output=True, text=True, cwd=cwd, shell=shell
        )

    if check and result.returncode != 0:
        if log_file:
            log_warn(f"  Non-zero exit ({result.returncode}) — see {log_file}")
        raise subprocess.CalledProcessError(result.returncode, cmd)

    return result


# ---------------------------------------------------------------------------
# File utilities
# ---------------------------------------------------------------------------

def count_fasta(path: Path) -> int:
    """Count sequences in a FASTA file."""
    if not path.exists():
        return 0
    count = 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


def fasta_stats(path: Path) -> dict:
    """Return basic stats for a FASTA file."""
    if not path.exists():
        return {"n": 0, "total_bp": 0, "largest_bp": 0, "n50": 0}

    lengths = []
    length = 0
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if length:
                    lengths.append(length)
                length = 0
            else:
                length += len(line)
    if length:
        lengths.append(length)

    if not lengths:
        return {"n": 0, "total_bp": 0, "largest_bp": 0, "n50": 0}

    lengths.sort(reverse=True)
    total = sum(lengths)
    cumsum = 0
    n50 = 0
    for ln in lengths:
        cumsum += ln
        if cumsum >= total / 2:
            n50 = ln
            break

    return {
        "n":          len(lengths),
        "total_bp":   total,
        "largest_bp": lengths[0],
        "n50":        n50,
    }


def human_size(path: Path) -> str:
    """Return human-readable file size."""
    if not path.exists():
        return "N/A"
    size = path.stat().st_size
    for unit in ["B", "KB", "MB", "GB"]:
        if size < 1024:
            return f"{size:.1f}{unit}"
        size /= 1024
    return f"{size:.1f}TB"


def mkdirs(*paths: Path) -> None:
    """Create directories, including parents."""
    for p in paths:
        Path(p).mkdir(parents=True, exist_ok=True)

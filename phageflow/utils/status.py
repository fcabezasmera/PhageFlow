"""Pipeline status tracking for PhageFlow.

Writes and reads pipeline_status.tsv in the project reports directory.
Allows the runner and `phageflow status` command to track module completion.
"""
from __future__ import annotations
import csv
from datetime import datetime
from pathlib import Path
from typing import Optional


_STATUS_HEADERS = [
    "sample_id", "module", "status",
    "started_at", "finished_at", "output_path", "notes",
]

_STATUS_VALUES = {"pending", "running", "completed", "failed", "skipped"}


def _status_file(reports_dir: Path) -> Path:
    return reports_dir / "pipeline_status.tsv"


def _load(reports_dir: Path) -> list[dict]:
    path = _status_file(reports_dir)
    if not path.exists():
        return []
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(dict(row))
    return rows


def _save(rows: list[dict], reports_dir: Path) -> None:
    path = _status_file(reports_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=_STATUS_HEADERS,
                                delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _key(sample_id: str, module: str) -> tuple[str, str]:
    return (sample_id.strip(), module.strip())


def mark(
    reports_dir: Path,
    sample_id:   str,
    module:      str,
    status:      str,
    output_path: str = "",
    notes:       str = "",
) -> None:
    """Write or update a status row for (sample_id, module)."""
    assert status in _STATUS_VALUES, f"Invalid status: {status}"
    rows = _load(reports_dir)
    now  = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    target_key = _key(sample_id, module)

    updated = False
    for row in rows:
        if _key(row["sample_id"], row["module"]) == target_key:
            if status == "running":
                row["started_at"]  = now
                row["finished_at"] = ""
            else:
                row["finished_at"] = now
            row["status"]      = status
            row["output_path"] = output_path or row.get("output_path", "")
            row["notes"]       = notes or row.get("notes", "")
            updated = True
            break

    if not updated:
        rows.append({
            "sample_id":   sample_id,
            "module":      module,
            "status":      status,
            "started_at":  now if status == "running" else "",
            "finished_at": now if status != "running" else "",
            "output_path": output_path,
            "notes":       notes,
        })

    _save(rows, reports_dir)


def is_completed(reports_dir: Path, sample_id: str, module: str) -> bool:
    for row in _load(reports_dir):
        if _key(row["sample_id"], row["module"]) == _key(sample_id, module):
            return row.get("status") == "completed"
    return False


def get_all(reports_dir: Path) -> list[dict]:
    return _load(reports_dir)

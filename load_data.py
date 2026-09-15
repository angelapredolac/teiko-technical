"""Validate the supplied CSV and atomically build the root SQLite database.

Run directly with ``python load_data.py``; no third-party packages are needed.
"""

from __future__ import annotations

import csv
import hashlib
import os
import sqlite3
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "cell-count.csv"
DB_PATH = ROOT / "cell_counts.db"
POPULATIONS = ("b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte")
SUBJECT_FIELDS = ("project", "condition", "age", "sex", "treatment", "response")
HEADERS = (
    "project",
    "subject",
    "condition",
    "age",
    "sex",
    "treatment",
    "response",
    "sample",
    "sample_type",
    "time_from_treatment_start",
    *POPULATIONS,
)


def read_rows(path: Path) -> list[dict]:
    """Reject invalid records rather than silently dropping or repairing them."""
    subjects: dict[str, tuple] = {}
    samples: set[str] = set()
    rows = []
    with path.open(newline="", encoding="utf-8-sig") as source:
        reader = csv.DictReader(source)
        if reader.fieldnames != list(HEADERS):
            raise ValueError(f"Expected CSV columns in this order: {HEADERS}")
        for line, raw in enumerate(reader, start=2):
            try:
                if None in raw or any(value is None for value in raw.values()):
                    raise ValueError("incorrect number of fields")
                row = {key: value.strip() for key, value in raw.items()}
                if any(not value for key, value in row.items() if key != "response"):
                    raise ValueError("missing required value")
                for key in ("age", "time_from_treatment_start", *POPULATIONS):
                    row[key] = int(row[key])
                if row["age"] < 0 or any(row[p] < 0 for p in POPULATIONS):
                    raise ValueError("age and cell counts must be nonnegative")
                if sum(row[p] for p in POPULATIONS) == 0:
                    raise ValueError("sample has zero total cells; percentages are undefined")
                if row["sex"] not in ("M", "F") or row["response"] not in ("yes", "no", ""):
                    raise ValueError("invalid sex or response category")
                row["response"] = row["response"] or None
                if row["sample"] in samples:
                    raise ValueError(f"duplicate sample {row['sample']}")
                samples.add(row["sample"])
                metadata = tuple(row[key] for key in SUBJECT_FIELDS)
                if subjects.setdefault(row["subject"], metadata) != metadata:
                    raise ValueError(f"conflicting metadata for subject {row['subject']}")
                rows.append(row)
            except (ValueError, TypeError) as exc:
                raise ValueError(f"{path.name}, line {line}: {exc}") from exc
    if not rows:
        raise ValueError("CSV contains no samples")
    return rows


def load_data(csv_path: Path = CSV_PATH, db_path: Path = DB_PATH) -> Path:
    """Build in a temporary file; preserve the existing database on failure."""
    rows = read_rows(csv_path)
    fd, name = tempfile.mkstemp(prefix=".cell_counts-", suffix=".db", dir=db_path.parent)
    os.close(fd)
    temporary = Path(name)
    try:
        connection = sqlite3.connect(temporary)
        try:
            connection.executescript((ROOT / "schema.sql").read_text())
            with connection:
                connection.executemany(
                    "INSERT INTO projects VALUES (?)", sorted({(r["project"],) for r in rows})
                )
                subjects = {
                    r["subject"]: (r["subject"], *(r[k] for k in SUBJECT_FIELDS)) for r in rows
                }
                connection.executemany(
                    "INSERT INTO subjects VALUES (?, ?, ?, ?, ?, ?, ?)", subjects.values()
                )
                connection.executemany(
                    "INSERT INTO populations VALUES (?)", [(p,) for p in POPULATIONS]
                )
                connection.executemany(
                    "INSERT INTO samples VALUES (?, ?, ?, ?)",
                    [
                        (
                            r["sample"],
                            r["subject"],
                            r["sample_type"],
                            r["time_from_treatment_start"],
                        )
                        for r in rows
                    ],
                )
                connection.executemany(
                    "INSERT INTO cell_counts VALUES (?, ?, ?)",
                    [(r["sample"], p, r[p]) for r in rows for p in POPULATIONS],
                )
                connection.execute(
                    "CREATE TABLE provenance (source TEXT, sha256 TEXT, sample_count INTEGER)"
                )
                connection.execute(
                    "INSERT INTO provenance VALUES (?, ?, ?)",
                    (csv_path.name, hashlib.sha256(csv_path.read_bytes()).hexdigest(), len(rows)),
                )
                if connection.execute("PRAGMA foreign_key_check").fetchall():
                    raise ValueError("Foreign key validation failed")
        finally:
            connection.close()
        os.replace(temporary, db_path)
    finally:
        temporary.unlink(missing_ok=True)
    return db_path


if __name__ == "__main__":
    print(f"Created {load_data()}")

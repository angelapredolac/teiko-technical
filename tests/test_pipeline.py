"""Independent numerical checks and failure-path regression tests."""

import csv
import sqlite3

import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

from analysis import baseline_tables, frequencies, holm_adjust, statistics, subject_values
from load_data import CSV_PATH, DB_PATH, HEADERS, POPULATIONS, ROOT, load_data, read_rows


@pytest.fixture
def source_rows():
    with CSV_PATH.open(newline="") as file:
        return list(csv.DictReader(file))


def write_csv(path, rows):
    with path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=HEADERS)
        writer.writeheader()
        writer.writerows(rows)
    return path


def test_full_import_and_percentages(tmp_path, source_rows):
    db = load_data(db_path=tmp_path / "test.db")
    frame = frequencies(db)
    assert len(frame) == len(source_rows) * 5
    assert np.allclose(frame.groupby("sample").percentage.sum(), 100)
    for row in (source_rows[0], source_rows[len(source_rows) // 2], source_rows[-1]):
        total = sum(int(row[p]) for p in POPULATIONS)
        for p in POPULATIONS:
            record = frame.loc[frame["sample"].eq(row["sample"]) & frame.population.eq(p)].iloc[0]
            assert record.total_count == total
            assert record["count"] == int(row[p])
            assert record.percentage == pytest.approx(100 * int(row[p]) / total)
    load_data(db_path=db)
    with sqlite3.connect(db) as connection:
        assert connection.execute("SELECT COUNT(*) FROM samples").fetchone()[0] == len(source_rows)
        assert connection.execute("PRAGMA integrity_check").fetchone()[0] == "ok"
        assert connection.execute("PRAGMA foreign_key_check").fetchall() == []
        assert (
            connection.execute("SELECT COUNT(*) FROM subjects WHERE response IS NULL").fetchone()[0]
            == 474
        )


@pytest.mark.parametrize(
    "case", ["duplicate", "negative", "zero", "conflict", "fractional", "missing"]
)
def test_invalid_input_preserves_database(tmp_path, source_rows, case):
    rows = [dict(r) for r in source_rows[:2]]
    if case == "duplicate":
        rows[1]["sample"] = rows[0]["sample"]
    elif case == "negative":
        rows[0]["b_cell"] = "-1"
    elif case == "zero":
        rows[0].update(dict.fromkeys(POPULATIONS, "0"))
    elif case == "conflict":
        rows[1]["age"] = str(int(rows[0]["age"]) + 1)
    elif case == "fractional":
        rows[0]["b_cell"] = "1.5"
    else:
        rows[0]["subject"] = ""
    db = tmp_path / "existing.db"
    db.write_bytes(b"previous database must survive")
    with pytest.raises(ValueError, match="line"):
        load_data(write_csv(tmp_path / "bad.csv", rows), db)
    assert db.read_bytes() == b"previous database must survive"


def test_baseline_counts_and_distinct_subjects(tmp_path, source_rows):
    # Add another baseline sample from the same subject: sample counts must rise,
    # while subject counts stay unchanged.
    baseline = next(
        r
        for r in source_rows
        if r["condition"] == "melanoma"
        and r["treatment"] == "miraclib"
        and r["sample_type"] == "PBMC"
        and r["time_from_treatment_start"] == "0"
    )
    duplicate_visit = dict(baseline, sample="extra_baseline")
    db = load_data(
        write_csv(tmp_path / "source.csv", source_rows + [duplicate_visit]), tmp_path / "data.db"
    )
    tables = baseline_tables(db)
    assert len(tables["samples"]) == 657
    assert tables["response"].set_index("response").subjects.to_dict() == {"no": 325, "yes": 331}
    assert tables["sex"].set_index("sex").subjects.to_dict() == {"M": 344, "F": 312}
    assert tables["projects"]["samples"].sum() == 657


def test_holm_known_example():
    np.testing.assert_allclose(
        holm_adjust([0.01, 0.04, 0.03, 0.2, 0.5]), [0.05, 0.12, 0.12, 0.4, 0.5]
    )


def test_subjects_have_equal_weight():
    frame = pd.DataFrame(
        {
            "subject": ["a", "a", "a", "b"],
            "response": ["yes"] * 4,
            "population": ["b_cell"] * 4,
            "percentage": [10, 20, 30, 90],
            "time_from_treatment_start": [0, 7, 14, 0],
        }
    )
    result = subject_values(frame).set_index("subject")
    assert result.loc["a", "percentage"] == 20
    assert result.percentage.mean() == 55
    assert subject_values(frame, True).set_index("subject").loc["a", "percentage"] == 10


def test_effect_direction_and_missing_group():
    values = pd.DataFrame(
        {
            "population": ["b_cell"] * 8,
            "response": ["yes"] * 4 + ["no"] * 4,
            "percentage": [60, 70, 80, 90, 10, 20, 30, 40],
        }
    )
    result = statistics(values).set_index("population")
    assert result.loc["b_cell", "rank_biserial"] == 1
    assert result.loc["b_cell", "median_difference_pp"] == 50
    assert result.loc["monocyte", "conclusion"] == "Insufficient subjects"
    assert np.isnan(result.loc["monocyte", "p_holm"])


def test_header_validation(tmp_path):
    path = tmp_path / "bad.csv"
    path.write_text("sample,count\na,1\n")
    with pytest.raises(ValueError, match="columns"):
        read_rows(path)


def test_dashboard_interactions():
    if not DB_PATH.exists():
        load_data()
    app = AppTest.from_file(str(ROOT / "app.py"), default_timeout=30).run()
    assert not app.exception
    app.radio[0].set_value("Baseline only").run()
    assert not app.exception
    app.radio[1].set_value("Subject means (statistical test unit)").run()
    assert not app.exception
    app.multiselect[0].set_value([]).run()
    assert not app.exception
    assert app.dataframe[0].value.empty


def test_cloud_entry_point():
    if not DB_PATH.exists():
        load_data()
    app = AppTest.from_file(str(ROOT / "cloud_app.py"), default_timeout=30).run()
    assert not app.exception

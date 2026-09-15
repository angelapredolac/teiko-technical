"""Database queries, subject-level inference, and reproducible report exports."""

from __future__ import annotations

import sqlite3
from contextlib import closing
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
from scipy.stats import mannwhitneyu

from load_data import DB_PATH, POPULATIONS, ROOT

COHORT_WHERE = "condition = 'melanoma' AND treatment = 'miraclib' AND sample_type = 'PBMC'"


def query(sql: str, db_path: Path = DB_PATH, params: tuple = ()) -> pd.DataFrame:
    """Open an existing database read-only; never create an empty one by accident."""
    with closing(sqlite3.connect(db_path.resolve().as_uri() + "?mode=ro", uri=True)) as db:
        return pd.read_sql_query(sql, db, params=params)


def frequencies(db_path: Path = DB_PATH) -> pd.DataFrame:
    return query("SELECT * FROM sample_frequencies ORDER BY sample, population", db_path)


def comparison_samples(db_path: Path = DB_PATH) -> pd.DataFrame:
    return query(
        f"SELECT m.*, f.population, f.count, f.total_count, f.percentage "
        f"FROM sample_metadata m JOIN sample_frequencies f USING(sample) "
        f"WHERE {COHORT_WHERE} AND response IN ('yes', 'no') ORDER BY sample, population",
        db_path,
    )


def baseline_tables(db_path: Path = DB_PATH) -> dict[str, pd.DataFrame]:
    """Count samples by project and DISTINCT subjects by response and sex in SQL."""
    return {
        "samples": query(
            "SELECT * FROM baseline_samples ORDER BY project, subject, sample", db_path
        ),
        "projects": query(
            "SELECT project, COUNT(*) AS samples FROM baseline_samples GROUP BY project", db_path
        ),
        "response": query(
            "SELECT COALESCE(response, 'unknown') AS response, COUNT(DISTINCT subject) AS subjects "
            "FROM baseline_samples GROUP BY response",
            db_path,
        ),
        "sex": query(
            "SELECT sex, COUNT(DISTINCT subject) AS subjects FROM baseline_samples GROUP BY sex",
            db_path,
        ),
    }


def subject_values(samples: pd.DataFrame, baseline_only: bool = False) -> pd.DataFrame:
    """One equally weighted observation per subject and population."""
    selected = samples.loc[samples.time_from_treatment_start.eq(0)] if baseline_only else samples
    return selected.groupby(["subject", "response", "population"], as_index=False).percentage.mean()


def holm_adjust(pvalues: list[float]) -> np.ndarray:
    """Holm step-down adjustment, valid under arbitrary dependence among tests."""
    p = np.asarray(pvalues, dtype=float)
    if not len(p):
        return p
    order = np.argsort(p)
    adjusted = np.empty_like(p)
    adjusted[order] = np.minimum(1, np.maximum.accumulate(p[order] * np.arange(len(p), 0, -1)))
    return adjusted


def statistics(values: pd.DataFrame) -> pd.DataFrame:
    """Two-sided Mann–Whitney U tests on independent subject summaries.

    Rank-biserial correlation is positive when responder values tend to be larger.
    Insufficient groups are explicitly untestable, not declared nonsignificant.
    """
    results = []
    for population in POPULATIONS:
        part = values.loc[values.population.eq(population)]
        yes = part.loc[part.response.eq("yes"), "percentage"].to_numpy()
        no = part.loc[part.response.eq("no"), "percentage"].to_numpy()
        testable = len(yes) >= 2 and len(no) >= 2
        u, p = (
            mannwhitneyu(yes, no, alternative="two-sided", method="asymptotic")
            if testable
            else (np.nan, np.nan)
        )
        results.append(
            {
                "population": population,
                "responders": len(yes),
                "nonresponders": len(no),
                "responder_median_pct": np.median(yes) if len(yes) else np.nan,
                "nonresponder_median_pct": np.median(no) if len(no) else np.nan,
                "median_difference_pp": np.median(yes) - np.median(no)
                if len(yes) and len(no)
                else np.nan,
                "u_statistic": u,
                "p_value": p,
                "rank_biserial": 2 * u / (len(yes) * len(no)) - 1 if testable else np.nan,
            }
        )
    result = pd.DataFrame(results)
    # Retain all five planned comparisons in the multiplicity adjustment.
    result["p_holm"] = holm_adjust(result.p_value.fillna(1).tolist())
    result.loc[result.p_value.isna(), "p_holm"] = np.nan
    result["conclusion"] = np.where(
        result.p_value.isna(),
        "Insufficient subjects",
        np.where(result.p_holm < 0.05, "Significant", "Not significant"),
    )
    return result


def boxplot(data: pd.DataFrame, title: str):
    """Show all five populations and explicit responder labels."""
    plot_data = data.assign(response=data.response.map({"yes": "Responder", "no": "Non-responder"}))
    fig = px.box(
        plot_data,
        x="population",
        y="percentage",
        color="response",
        points="outliers",
        category_orders={
            "population": list(POPULATIONS),
            "response": ["Responder", "Non-responder"],
        },
        color_discrete_map={"Responder": "#087e8b", "Non-responder": "#cb6843"},
        labels={
            "population": "Cell population",
            "percentage": "Relative frequency (%)",
            "response": "Response",
        },
        title=title,
        template="plotly_white",
    )
    fig.update_layout(legend_title_text="", boxmode="group", height=480)
    return fig


def run_analysis(db_path: Path = DB_PATH, output_dir: Path = ROOT / "outputs") -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    frequencies(db_path).to_csv(output_dir / "sample_frequencies.csv", index=False)
    samples = comparison_samples(db_path)
    reports = {}
    for mode, baseline_only in (("all_timepoints", False), ("baseline", True)):
        values = subject_values(samples, baseline_only)
        report = statistics(values)
        reports[mode] = report
        report.to_csv(output_dir / f"statistics_{mode}.csv", index=False)
        values.to_csv(output_dir / f"subject_frequencies_{mode}.csv", index=False)
    boxplot(samples, "All qualifying PBMC samples · descriptive comparison").write_html(
        output_dir / "response_boxplots.html", include_plotlyjs=True
    )
    for name, table in baseline_tables(db_path).items():
        table.to_csv(output_dir / f"baseline_{name}.csv", index=False)
    lines = [
        "# Analysis results",
        "",
        "Melanoma · miraclib · PBMC. Alpha = 0.05; Holm correction across five populations.",
        "",
    ]
    for mode, report in reports.items():
        lines += [f"## {mode}", "", "```", report.to_string(index=False), "```", ""]
    lines += [
        "All-timepoint tests use one mean percentage per subject. Baseline tests use day 0 only.",
        "These are exploratory associations, not evidence of causal drug effects "
        "or validated prediction.",
    ]
    (output_dir / "results.md").write_text("\n".join(lines) + "\n")
    print(f"Exported tables, statistical results and interactive boxplots to {output_dir}")


if __name__ == "__main__":
    run_analysis()

"""Interactive Parts 2–4 dashboard; all analytical data comes from SQLite."""

import streamlit as st

from analysis import (
    baseline_tables,
    boxplot,
    comparison_samples,
    frequencies,
    query,
    statistics,
    subject_values,
)
from load_data import DB_PATH

st.set_page_config(page_title="Loblaw Bio | Immune cell analysis", page_icon="🧬", layout="wide")
st.title("Immune cell response explorer")
st.caption("Loblaw Bio · Clinical trial analysis · Melanoma / miraclib / PBMC")
if not DB_PATH.exists():
    st.error(
        "Database not found. Run `make setup` followed by `make pipeline`, then reload this page."
    )
    st.stop()


def show_table(frame, filename):
    columns = {
        "percentage": st.column_config.NumberColumn("Percentage (%)", format="%.3f"),
        "p_value": st.column_config.NumberColumn("Raw p", format="%.4g"),
        "p_holm": st.column_config.NumberColumn("Adjusted p", format="%.4g"),
        "rank_biserial": st.column_config.NumberColumn("Rank-biserial", format="%.3f"),
        "median_difference_pp": st.column_config.NumberColumn("Median Δ (pp)", format="%.3f"),
        "responder_median_pct": st.column_config.NumberColumn(
            "Responder median (%)", format="%.3f"
        ),
        "nonresponder_median_pct": st.column_config.NumberColumn(
            "Non-responder median (%)", format="%.3f"
        ),
    }
    order = list(frame.columns)
    if "p_holm" in frame:
        leading = ["population", "conclusion", "p_holm", "p_value", "median_difference_pp"]
        order = leading + [column for column in order if column not in leading]
    st.dataframe(
        frame,
        use_container_width=True,
        hide_index=True,
        column_config=columns,
        column_order=order,
    )
    st.download_button(
        "Download CSV", frame.to_csv(index=False), file_name=filename, mime="text/csv", key=filename
    )


overview, comparison, baseline, methods = st.tabs(
    ["Sample overview", "Response comparison", "Baseline cohort", "Methods"]
)

with overview:
    st.header("Cell frequencies in every sample")
    summary = frequencies()
    metadata = query("SELECT * FROM sample_metadata ORDER BY sample")
    c1, c2, c3 = st.columns(3)
    c1.metric("Samples", metadata.shape[0])
    c2.metric("Subjects", metadata.subject.nunique())
    c3.metric("Cell populations", summary.population.nunique())
    projects = st.multiselect(
        "Projects", sorted(metadata.project.unique()), default=sorted(metadata.project.unique())
    )
    search = st.text_input("Sample ID contains", placeholder="e.g. sample00000")
    eligible = metadata.loc[metadata.project.isin(projects), "sample"]
    displayed = summary.loc[
        summary["sample"].isin(eligible) & summary["sample"].str.contains(search, regex=False)
    ]
    st.caption(
        "Percentage = 100 × population count ÷ sum of all five population counts in that sample."
    )
    show_table(displayed, "sample_frequencies.csv")

with comparison:
    st.header("Responders versus non-responders")
    st.caption(
        "Fixed cohort: melanoma patients receiving miraclib, PBMC samples, known yes/no response."
    )
    samples = comparison_samples()
    mode = st.radio("Analysis window", ["All timepoints", "Baseline only"], horizontal=True)
    baseline_only = mode == "Baseline only"
    selected = samples.loc[samples.time_from_treatment_start.eq(0)] if baseline_only else samples
    values = subject_values(samples, baseline_only)
    counts = values[["subject", "response"]].drop_duplicates().response.value_counts()
    c1, c2, c3 = st.columns(3)
    c1.metric("Responders", int(counts.get("yes", 0)))
    c2.metric("Non-responders", int(counts.get("no", 0)))
    c3.metric("Samples", selected["sample"].nunique())
    unit = st.radio(
        "Boxplot observations",
        ["Individual samples", "Subject means (statistical test unit)"],
        horizontal=True,
    )
    plotted = selected if unit == "Individual samples" else values
    if plotted.empty:
        st.info("No samples match this analysis window.")
    else:
        st.plotly_chart(boxplot(plotted, f"{mode} · {unit}"), use_container_width=True)
    st.caption(
        "Sample plots are descriptive. Tests use one mean frequency per subject, "
        "so repeat visits do not count as independent patients."
    )
    report = statistics(values)
    significant = report.loc[report.conclusion.eq("Significant"), "population"].tolist()
    if significant:
        st.success("Significant after Holm correction: " + ", ".join(significant))
    elif report.p_value.notna().any():
        st.info("No tested population meets the Holm-adjusted p < 0.05 threshold.")
    else:
        st.info("Insufficient subjects to perform the comparison.")
    st.subheader("Statistical evidence")
    st.caption(
        "Two-sided Mann–Whitney U; Holm-adjusted p-values across five tests. "
        "Median differences are percentage points; positive effect sizes indicate "
        "higher responder frequencies."
    )
    show_table(report, f"statistics_{'baseline' if baseline_only else 'all_timepoints'}.csv")
    st.info(
        "Post-treatment differences are associations with response. Baseline-only results "
        "explore potential predictors; predictive performance has not been validated."
    )

with baseline:
    st.header("Baseline melanoma cohort")
    st.caption("Database filter: melanoma · miraclib · PBMC · time_from_treatment_start = 0.")
    tables = baseline_tables()
    c1, c2 = st.columns(2)
    c1.metric("Baseline samples", len(tables["samples"]))
    c2.metric("Unique subjects", tables["samples"].subject.nunique())
    for col, key, title in zip(
        st.columns(3),
        ("projects", "response", "sex"),
        ("Samples by project", "Subjects by response", "Subjects by sex"),
        strict=True,
    ):
        with col:
            st.subheader(title)
            show_table(tables[key], f"baseline_{key}.csv")
    st.subheader("Matching samples")
    show_table(tables["samples"], "baseline_samples.csv")

with methods:
    st.header("How to interpret the results")
    st.markdown("""
    - **Denominator:** the sum of the five measured populations, not all possible immune cells.
    - **Repeated visits:** average each subject's sample percentages before the primary test;
      each subject receives equal weight. Baseline-only analysis uses day 0.
    - **Test:** two-sided Mann–Whitney U with asymptotic, tie-corrected p-values
      and continuity correction.
      It tests distributional differences; it is not generally a test of medians alone.
    - **Multiplicity:** Holm correction controls the family-wise error rate across the five
      populations, separately for each analysis window, at 0.05.
    - **Effect size:** rank-biserial correlation ranges from −1 to +1. Positive values mean
      responder frequencies tend to be higher. Median differences use percentage points.
    - **Missing response:** retained as SQL NULL in the database; excluded from Part 3 tests.
      Baseline counts would report these subjects as unknown.
    - **Limits:** cell percentages sum to 100% and are dependent. Project, age, sex and other
      confounders are not adjusted for. Averaging visits does not model time trends or prove
      a drug effect. Both analyses are exploratory; multiple windows add opportunities
      for discovery.
      Prediction would require a baseline model evaluated on held-out subjects.
    """)
    st.subheader("Source provenance")
    st.dataframe(query("SELECT * FROM provenance"), hide_index=True, use_container_width=True)

"""Statistical Inference page (OpenMS-Insight engine)."""
import polars as pl
import streamlit as st

from src.common.common import page_setup
from src.common.results_helpers import get_abundance_data, get_id_column, get_sample_group_map
from openms_insight.analysis.statistics import adjust_fdr_lazy, calculate_statistical_tests

params = page_setup(page="workflow")
st.title("Statistical Inference")

st.markdown("Run differential abundance analysis to identify significant features across groups.")

if "workspace" not in st.session_state:
    st.warning("Please initialize your workspace first.")
    st.stop()

workspace = st.session_state["workspace"]
result = get_abundance_data(workspace)
if result is None:
    st.info("No feature matrix found yet. Run **UmetaFlow TOPP Workflow** first.")
    st.page_link("content/umetaflow_run.py", label="Go to Run", icon="🚀")
    st.stop()

pivot_df, expr_df, sample_cols = result
id_col = get_id_column()
sample_group_map = get_sample_group_map(workspace, sample_cols)

if st.session_state.get("normalized_df") is not None:
    base_df = st.session_state["normalized_df"]
    st.info("🔄 Using data from the **Normalization** step.")
elif st.session_state.get("imputed_df") is not None:
    base_df = st.session_state["imputed_df"]
    st.warning("⚠️ Normalization skipped — using data from the **Imputation** step.")
elif st.session_state.get("filtered_df") is not None:
    base_df = st.session_state["filtered_df"]
    st.warning("⚠️ Preprocessing skipped — using data from the **Filtering** step.")
else:
    base_df = pivot_df
    st.warning("⚠️ No preprocessing history found — operating on the raw feature matrix.")

unique_groups = sorted({sample_group_map[s] for s in sample_cols if sample_group_map[s]})
group_count = len(unique_groups)

st.subheader("Input Table Overview")
st.markdown(
    f"**{base_df.shape[0]}** rows across **{len(sample_cols)}** samples "
    f"in **{group_count}** group(s) ({', '.join(unique_groups)})."
)
st.dataframe(base_df, use_container_width=True)

st.markdown("---")
st.subheader("Configure Statistical Test")

metadata_rows = [
    {"sample_id": s, "group": sample_group_map[s]} for s in sample_cols if sample_group_map[s]
]
metadata_pl = pl.DataFrame(metadata_rows, schema={"sample_id": pl.String, "group": pl.String})

col1, col2 = st.columns(2)
with col1:
    if group_count == 2:
        method_options = ["limma_like", "welch", "paired"]
    elif group_count >= 3:
        method_options = ["limma_like", "anova"]
    else:
        st.error("Statistical testing requires at least 2 sample groups. Assign groups on the Filtering page.")
        st.stop()
    selected_method = st.selectbox("Statistical test", options=method_options)
with col2:
    selected_fdr = st.selectbox("FDR adjustment", options=["BH", "Bonferroni", "None"])

if st.button("Run Statistical Analysis", type="primary"):
    stats_lazy = pl.from_pandas(base_df).lazy()
    try:
        stats_lazy = calculate_statistical_tests(
            quantification_data=stats_lazy, metadata=metadata_pl, method=selected_method
        )
        stats_lazy = adjust_fdr_lazy(quantification_data=stats_lazy, strategy=selected_fdr)
        statistics_df = stats_lazy.collect().to_pandas()
        st.session_state["statistics_df"] = statistics_df

        st.success(f"Calculated **{selected_method}** test with **{selected_fdr}** FDR correction.")

        extreme = (statistics_df["log2FC"].abs() > 100).sum()
        if extreme:
            st.warning(
                f"⚠️ {extreme} feature(s) have |log2FC| > 100. `log2FC` here is just a difference of "
                "the input values, so this usually means the intensities weren't log2-transformed "
                "first — go to **Normalization** and set Transformation = 'log2', then re-run this "
                "test. (Rendering extreme values on the Volcano plot can also crash the page.)"
            )

        st.subheader("Results")
        st.markdown(f"Columns: `{id_col}`, `log2FC`, `stat`, `p-value`, `p-adj`")
        st.dataframe(statistics_df, use_container_width=True)

        st.page_link("content/downstream_results.py", label="View Volcano / PCA / Heatmap", icon="📊")
    except ValueError as val_err:
        st.error(f"Validation error: {val_err}")
    except Exception as e:
        st.error(f"Unexpected error: {e}")

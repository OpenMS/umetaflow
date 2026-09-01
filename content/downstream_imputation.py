"""Missing Value Imputation page (OpenMS-Insight engine)."""
import polars as pl
import streamlit as st

from src.common.common import page_setup
from src.common.results_helpers import get_abundance_data, get_id_column, get_sample_group_map
from openms_insight.analysis.imputation import impute_mar, impute_smallest_value

params = page_setup(page="workflow")
st.title("Missing Value Imputation")

st.markdown(
    "Handle missing values (zero-intensity features) using biological group-aware (MAR) "
    "or lowest-detected-value (MNAR) strategies."
)

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

if "filtered_df" in st.session_state and st.session_state["filtered_df"] is not None:
    base_df = st.session_state["filtered_df"]
    st.info("🔄 Using data from the **Filtering** step.")
else:
    base_df = pivot_df
    st.warning("⚠️ No filtering history found — operating on the unfiltered feature matrix.")

st.subheader("Input Matrix Overview")
st.markdown(f"**{base_df.shape[0]}** features across **{len(sample_cols)}** samples.")
st.dataframe(base_df, use_container_width=True)

st.markdown("---")
st.subheader("Configure Imputation")

metadata_rows = [
    {"sample_id": s, "group": sample_group_map[s]} for s in sample_cols if sample_group_map[s]
]
if not metadata_rows:
    st.warning("Assign sample groups on the **Filtering** page first.")
    st.stop()
metadata_pl = pl.DataFrame(metadata_rows, schema={"sample_id": pl.String, "group": pl.String})

impute_category = st.selectbox(
    "Imputation class",
    options=["MAR (Missing At Random)", "MNAR (Missing Not At Random)"],
    help="MAR fills with a group mean/median. MNAR fills with the smallest observed value "
    "(assumes missingness is due to detection-limit dropout).",
)

if impute_category == "MAR (Missing At Random)":
    strategy_opt = st.radio("Strategy", options=["median", "mean"], horizontal=True)
else:
    scope_opt = st.radio(
        "Detection-minimum scope",
        options=["row", "global"],
        horizontal=True,
        help="'row' uses that feature's own minimum; 'global' uses the smallest value in the whole matrix.",
    )

if st.button("Apply Imputation", type="primary"):
    quant_lazy = pl.from_pandas(base_df).lazy()

    if impute_category == "MAR (Missing At Random)":
        imputed_lazy = impute_mar(
            quantification_data=quant_lazy,
            metadata=metadata_pl,
            group_column="group",
            strategy=strategy_opt,
        )
    else:
        imputed_lazy = impute_smallest_value(
            quantification_data=quant_lazy, metadata=metadata_pl, scope=scope_opt
        )

    imputed_df = imputed_lazy.collect().to_pandas()
    st.session_state["imputed_df"] = imputed_df

    st.success(f"Applied **{impute_category}** imputation.")
    st.subheader("Imputed Feature Matrix")
    st.dataframe(imputed_df, use_container_width=True)

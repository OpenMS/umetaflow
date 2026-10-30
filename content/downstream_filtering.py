"""Data Filtering page (OpenMS-Insight engine)."""
import pandas as pd
import polars as pl
import streamlit as st

from src.common.common import page_setup
from src.common.results_helpers import (
    get_abundance_data,
    get_id_column,
    get_sample_group_map,
    render_group_assignment,
)
from openms_insight.analysis.filter import (
    filter_low_abundance,
    filter_low_repeatability,
    filter_low_variance,
)

params = page_setup(page="workflow")
st.title("Data Filtering")

st.markdown(
    "Filter out low-quality features from your untargeted feature matrix based on "
    "abundance, repeatability, or variance thresholds."
)

if "workspace" not in st.session_state:
    st.warning("Please initialize your workspace first.")
    st.stop()

workspace = st.session_state["workspace"]
result = get_abundance_data(workspace)
if result is None:
    st.info(
        "No feature matrix found yet. Run **UmetaFlow TOPP Workflow** (Pre-Processing) first."
    )
    st.page_link("content/umetaflow_run.py", label="Go to Run", icon="🚀")
    st.stop()

pivot_df, expr_df, sample_cols = result
id_col = get_id_column()

sample_group_map = render_group_assignment(workspace, sample_cols)

st.markdown("---")

st.subheader("Original Feature Matrix")
st.markdown(
    f"Currently displaying **{pivot_df.shape[0]}** features across **{len(sample_cols)}** samples "
    "before filtering."
)
st.dataframe(pivot_df, use_container_width=True)

st.markdown("---")
st.subheader("Configure Filter")

metadata_rows = [
    {"sample_id": s, "group": sample_group_map[s]} for s in sample_cols if sample_group_map[s]
]
if not metadata_rows:
    st.warning("Assign at least one sample to a group above before filtering.")
    st.stop()

metadata_pl = pl.DataFrame(metadata_rows, schema={"sample_id": pl.String, "group": pl.String})

filter_method = st.selectbox(
    "Filtering method",
    options=["Low Abundance", "Low Repeatability", "Low Variance"],
    help="Choose the criteria used to prune unreliable feature rows.",
)

if filter_method == "Low Abundance":
    st.markdown(
        "Keeps rows where at least one group's median is above the selected percentile threshold."
    )
    threshold = st.slider("Threshold percentile (%)", 0.0, 100.0, 10.0, 5.0)
elif filter_method == "Low Repeatability":
    st.markdown(
        "Keeps rows where at least one group has a missing-value ratio within the allowed maximum."
    )
    threshold = st.slider(
        "Max missing ratio (%)", 0.0, 100.0, 50.0, 5.0,
        help="Allowed missing value (zero) ratio per group.",
    )
else:
    st.markdown(
        "Keeps rows where at least one group's variance is above the selected percentile threshold."
    )
    threshold = st.slider("Threshold percentile (%)", 0.0, 100.0, 10.0, 5.0)

if st.button("Apply Filter", type="primary"):
    quant_lazy = pl.from_pandas(pivot_df).lazy()

    if filter_method == "Low Abundance":
        filtered_lazy = filter_low_abundance(
            quantification_data=quant_lazy,
            metadata=metadata_pl,
            group_column="group",
            threshold_percentile=threshold,
        )
    elif filter_method == "Low Repeatability":
        filtered_lazy = filter_low_repeatability(
            quantification_data=quant_lazy,
            metadata=metadata_pl,
            group_column="group",
            max_missing_ratio=threshold / 100.0,
        )
    else:
        filtered_lazy = filter_low_variance(
            quantification_data=quant_lazy,
            metadata=metadata_pl,
            group_column="group",
            threshold_percentile=threshold,
        )

    filtered_df = filtered_lazy.collect().to_pandas()
    st.session_state["filtered_df"] = filtered_df

    st.success(f"Applied **{filter_method}** filter.")

    col1, col2, col3 = st.columns(3)
    col1.metric("Original features", pivot_df.shape[0])
    col2.metric("Filtered features", filtered_df.shape[0])
    col3.metric("Removed features", pivot_df.shape[0] - filtered_df.shape[0])

    st.subheader("Filtered Feature Matrix")
    if filtered_df.empty:
        st.warning("The filtered table is empty. Try relaxing the threshold.")
    else:
        st.dataframe(filtered_df, use_container_width=True)

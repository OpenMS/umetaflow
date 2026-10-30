"""Data Normalization & Scaling page (OpenMS-Insight engine)."""
import pandas as pd
import polars as pl
import streamlit as st

from src.common.common import page_setup
from src.common.results_helpers import get_abundance_data, get_id_column, get_sample_group_map
from openms_insight.analysis.normalization import normalize_samples, scale_data, transform_data

params = page_setup(page="workflow")
st.title("Data Normalization & Scaling")

st.markdown(
    "Correct for technical variation between samples and stabilize the intensity "
    "distribution before statistical testing."
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

filtered_df = st.session_state.get("filtered_df")
imputed_df = st.session_state.get("imputed_df")
normalized_df = st.session_state.get("normalized_df")

if imputed_df is not None:
    base_df = imputed_df
    st.info("🔄 Using data from the **Imputation** step.")
elif filtered_df is not None:
    base_df = filtered_df
    st.warning("⚠️ Imputation skipped — using data from the **Filtering** step.")
else:
    base_df = pivot_df
    st.warning("⚠️ No preprocessing history found — operating on the raw feature matrix.")

st.subheader("Pipeline Overview")
st.caption("Data flows: Filtering → Imputation → Normalization")
step_rows = [
    {"Step": "Filtering", "Status": "Done" if filtered_df is not None else "Not run"},
    {"Step": "Imputation", "Status": "Done" if imputed_df is not None else "Not run"},
    {"Step": "Normalization", "Status": "Done" if normalized_df is not None else "Not run"},
]
st.dataframe(pd.DataFrame(step_rows), hide_index=True, use_container_width=True)

st.markdown("---")
st.subheader("Configure Normalization")

metadata_rows = [
    {"sample_id": s, "group": sample_group_map[s]} for s in sample_cols if sample_group_map[s]
]
if not metadata_rows:
    st.warning("Assign sample groups on the **Filtering** page first.")
    st.stop()
metadata_pl = pl.DataFrame(metadata_rows, schema={"sample_id": pl.String, "group": pl.String})

col1, col2, col3 = st.columns(3)
with col1:
    st.markdown("**1. Transformation**")
    transform_strategy = st.selectbox(
        "Transformation", options=["None", "log2", "log10", "square_root", "cube_root"]
    )
with col2:
    st.markdown("**2. Sample Normalization**")
    norm_strategy = st.selectbox(
        "Normalization",
        options=["None", "sum", "median", "pqn", "reference_feature", "quantile"],
    )
    ref_feature_input = None
    if norm_strategy == "reference_feature":
        ref_feature_input = st.text_input(
            "Reference feature",
            placeholder="exact value from the metabolite column",
        )
with col3:
    st.markdown("**3. Row Scaling**")
    scaling_strategy = st.selectbox(
        "Scaling",
        options=["None", "mean_centering", "auto_scaling", "pareto_scaling", "range_scaling"],
    )

if st.button("Apply Normalization", type="primary"):
    if norm_strategy == "reference_feature" and not ref_feature_input:
        st.error("Provide a reference feature to use the 'reference_feature' strategy.")
        st.stop()

    processing_lazy = pl.from_pandas(base_df).lazy()
    try:
        processing_lazy = transform_data(
            quantification_data=processing_lazy, metadata=metadata_pl, strategy=transform_strategy
        )
        processing_lazy = normalize_samples(
            quantification_data=processing_lazy,
            metadata=metadata_pl,
            strategy=norm_strategy,
            id_col=id_col,
            reference_feature=ref_feature_input if norm_strategy == "reference_feature" else None,
        )
        processing_lazy = scale_data(
            quantification_data=processing_lazy, metadata=metadata_pl, strategy=scaling_strategy
        )

        normalized_df = processing_lazy.collect().to_pandas()
        st.session_state["normalized_df"] = normalized_df

        st.success("Normalization applied.")
        st.subheader("Normalized Feature Matrix")
        st.dataframe(normalized_df, use_container_width=True)
    except ValueError as val_err:
        st.error(f"Configuration error: {val_err}")
    except Exception as e:
        st.error(f"Unexpected error: {e}")

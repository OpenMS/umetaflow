"""Volcano / PCA / Clustered Heatmap results (OpenMS-Insight components)."""
import numpy as np
import polars as pl
import streamlit as st

from src.common.common import page_setup
from src.common.results_helpers import get_abundance_data, get_id_column, get_sample_group_map
from openms_insight import ClusteredHeatmap, PCAPlot, VolcanoPlot

params = page_setup(page="workflow")
st.title("Statistics Results")

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

state_manager = st.session_state.get("state")
GROUP_PALETTE = [
    "#00BFC4", "#F8766D", "#7CAE00", "#C77CFF", "#E7B800",
    "#619CFF", "#FF61C3", "#00BA38", "#FF8C42", "#00B0F6",
]

tab_volcano, tab_pca, tab_heatmap = st.tabs(["🌋 Volcano", "🧭 PCA", "🔥 Clustered Heatmap"])

# ---------------------------------------------------------------- Volcano ---
with tab_volcano:
    statistics_df = st.session_state.get("statistics_df")
    if statistics_df is None or statistics_df.empty:
        st.info("Run the Statistical Inference page first.")
        st.page_link("content/downstream_statistics.py", label="Go to Statistical Inference", icon="🔬")
    else:
        volcano_df = statistics_df.dropna(subset=["log2FC", "p-adj"]).copy()
        volcano_df = volcano_df[np.isfinite(volcano_df["log2FC"])]

        # calculate_statistical_tests() assumes its input is already log2-transformed and
        # just diffs group means; if Normalization's "log2" transform was skipped, log2FC
        # ends up on the raw-intensity scale (huge magnitudes) instead of a real fold change,
        # which can crash the volcano plot's axis rendering. Guard against that here.
        EXTREME_LOG2FC = 100
        n_extreme = (volcano_df["log2FC"].abs() > EXTREME_LOG2FC).sum()
        if n_extreme:
            st.warning(
                f"⚠️ {n_extreme} feature(s) have |log2FC| > {EXTREME_LOG2FC}, which usually means the "
                "input wasn't log2-transformed before running Statistical Inference "
                "(go to **Normalization** and set Transformation = 'log2' first). "
                "These features are excluded from the plot below to avoid rendering issues."
            )
            volcano_df = volcano_df[volcano_df["log2FC"].abs() <= EXTREME_LOG2FC]

        if volcano_df.empty:
            st.info("No features left to plot after filtering out non-finite/extreme log2FC values.")
        else:
            # VolcanoPlot has no built-in point cap (unlike Heatmap's min_points), and
            # rendering the full feature set (can be tens of thousands of rows) has been
            # observed to crash the whole Streamlit process, not just the browser tab.
            # Cap what's actually sent to the component, keeping the most significant
            # points first (lowest p-adj) since those are what a volcano plot is for.
            MAX_VOLCANO_POINTS = 5000
            if len(volcano_df) > MAX_VOLCANO_POINTS:
                st.info(
                    f"Showing the {MAX_VOLCANO_POINTS:,} most significant of {len(volcano_df):,} "
                    "features (by p-adj) to keep the plot responsive."
                )
                volcano_df = volcano_df.nsmallest(MAX_VOLCANO_POINTS, "p-adj")

            volcano_pl_lazy = pl.from_pandas(volcano_df).lazy()

            fc_thresh = st.slider("log2 Fold Change threshold", 0.5, 3.0, 1.0, 0.1, key="volcano_fc")
            p_thresh = st.slider("p-adj (FDR) threshold", 0.001, 0.1, 0.05, 0.001, key="volcano_p")

            volcano_component = VolcanoPlot(
                cache_id="umetaflow_volcano_plot",
                cache_path=str(workspace),
                data=volcano_pl_lazy,
                log2fc_column="log2FC",
                pvalue_column="p-adj",
                label_column=id_col,
                up_color="#E74C3C",
                down_color="#3498DB",
                ns_color="#95A5A6",
                show_threshold_lines=True,
                threshold_line_style="dash",
            )
            volcano_component(
                state_manager=state_manager,
                fc_threshold=fc_thresh,
                p_threshold=p_thresh,
                max_labels=10,
                height=600,
            )

# -------------------------------------------------------------------- PCA ---
with tab_pca:
    if st.session_state.get("normalized_df") is not None:
        base_df = st.session_state["normalized_df"]
    elif st.session_state.get("imputed_df") is not None:
        base_df = st.session_state["imputed_df"]
    elif st.session_state.get("filtered_df") is not None:
        base_df = st.session_state["filtered_df"]
    else:
        base_df = pivot_df

    unique_groups = sorted({sample_group_map[s] for s in sample_cols if sample_group_map[s]})
    if len(sample_cols) < 2:
        st.info("PCA requires at least 2 samples.")
    else:
        if len(unique_groups) < 2:
            st.warning("Only one biological group detected — points will plot without meaningful group coloring.")

        expr_df_wide = base_df.set_index(id_col)[sample_cols]
        max_available = expr_df_wide.shape[0]
        if max_available <= 20:
            top_n = max_available
            st.caption(f"Using all {top_n} features for PCA (dataset too small for variance filtering).")
        else:
            top_n = st.slider(
                "Number of features (highest variance)", 20, min(5000, max_available),
                min(500, max_available), 10, key="pca_top_n",
            )

        top_features = expr_df_wide.var(axis=1).sort_values(ascending=False).head(top_n).index
        expr_df_pca = expr_df_wide.loc[top_features].reset_index()

        if expr_df_pca.shape[0] < 2:
            st.info("Not enough features after variance filtering for PCA.")
        else:
            metadata_pl = pl.DataFrame(
                [{"sample_id": s, "group": sample_group_map[s]} for s in sample_cols if sample_group_map[s]],
                schema={"sample_id": pl.String, "group": pl.String},
            )
            pca_lazy = pl.from_pandas(expr_df_pca).lazy()
            try:
                pca_component = PCAPlot(
                    cache_id="umetaflow_pca_plot",
                    cache_path=str(workspace),
                    data=pca_lazy,
                    metadata=metadata_pl,
                    sample_id_field="sample_id",
                    group_field="group",
                    n_components=5,
                    title="Sample PCA",
                )
                variance_ratio = pca_component.get_variance_ratio()
                pc_columns = pca_component.get_pc_columns()

                col1, col2 = st.columns(2)
                with col1:
                    pc_x_label = st.selectbox("X-axis component", pc_columns, index=0, key="pca_pc_x")
                with col2:
                    default_y = 1 if len(pc_columns) > 1 else 0
                    pc_y_label = st.selectbox("Y-axis component", pc_columns, index=default_y, key="pca_pc_y")

                pc_x = int(pc_x_label.replace("PC", ""))
                pc_y = int(pc_y_label.replace("PC", ""))
                pca_component(state_manager=state_manager, pc_x=pc_x, pc_y=pc_y, height=600)

                st.markdown(
                    "**Explained variance:** "
                    + ", ".join(f"{c} {r * 100:.1f}%" for c, r in zip(pc_columns, variance_ratio))
                )
                st.markdown(f"**Features used:** {expr_df_pca.shape[0]} (top {top_n} by variance)")
            except ValueError as e:
                st.error(f"PCA computation failed: {e}")

# --------------------------------------------------------- Clustered Heatmap ---
with tab_heatmap:
    if expr_df.empty:
        st.info("No data available for heatmap.")
    else:
        top_n_hm = st.slider("Number of features (highest variance)", 10, 200, 30, key="clustered_heatmap_top_n")

        var_series = expr_df.var(axis=1)
        top_features = var_series.sort_values(ascending=False).head(top_n_hm).index
        heatmap_df = expr_df.loc[top_features]

        heatmap_z = heatmap_df.sub(heatmap_df.mean(axis=1), axis=0).div(heatmap_df.std(axis=1), axis=0)
        heatmap_z = heatmap_z.replace([np.inf, -np.inf], np.nan).dropna()

        if heatmap_z.empty:
            st.warning("Insufficient data to generate the heatmap.")
        else:
            heatmap_lazy = pl.from_pandas(heatmap_z.reset_index()).lazy()
            metadata_pl = pl.DataFrame(
                [{"sample_id": s, "group": sample_group_map[s]} for s in sample_cols if sample_group_map[s]],
                schema={"sample_id": pl.String, "group": pl.String},
            )

            heatmap_component = ClusteredHeatmap(
                cache_id="umetaflow_clustered_heatmap",
                cache_path=str(workspace),
                id_col=id_col,
                data=heatmap_lazy,
                metadata=metadata_pl,
                row_cluster=True,
                col_cluster=True,
                title="Feature Abundance Heatmap (Z-score, clustered)",
                x_label="Samples",
                y_label="Features",
                colorscale=[[0, "#6699E0"], [0.5, "#FFFFFF"], [1, "#E06666"]],
                reversescale=False,
                intensity_label="Z-score",
            )
            heatmap_component(state_manager=state_manager, height=700)

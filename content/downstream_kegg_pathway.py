"""KEGG Pathway mapping & enrichment for annotated metabolite features.

Resolves feature compound-name annotations (from SIRIUS/CSI:FingerID, GNPS or
in-house library matches) to KEGG compound IDs via the KEGG REST API, maps them to
KEGG pathways, and runs a Fisher's exact enrichment test comparing statistically
significant features (foreground) against all annotated features (background).
"""
import math

import pandas as pd
import plotly.express as px
import streamlit as st

from src.common.common import page_setup
from src.common.kegg import map_names_to_pathways, run_pathway_enrichment

params = page_setup(page="workflow")
st.title("KEGG Pathway Analysis")

st.markdown(
    "Map annotated features to KEGG compounds and pathways, then test which "
    "pathways are over-represented among your statistically significant features."
)

if "workspace" not in st.session_state:
    st.warning("Please initialize your workspace first.")
    st.stop()

statistics_df = st.session_state.get("statistics_df")
if statistics_df is None or statistics_df.empty:
    st.info("Run the Statistical Inference page first — KEGG enrichment needs log2FC/p-adj values.")
    st.page_link("content/downstream_statistics.py", label="Go to Statistical Inference", icon="🔬")
    st.stop()

# Any column whose name suggests it carries a human-readable compound name, e.g.
# SIRIUS "CSI:FingerID_<file>_name", GNPS/library "compound_name", or umetaflow's
# in-house-library annotation column "MS1 annotation" / "MS2 annotation".
NAME_HINTS = ("name", "annotation")
name_candidates = [
    c for c in statistics_df.columns if any(hint in c.lower() for hint in NAME_HINTS)
]
if not name_candidates:
    st.warning(
        "No compound-name column found in the statistics table. Run SIRIUS/CSI:FingerID "
        "or in-house library annotation in the UmetaFlow workflow first — KEGG lookups "
        "are done by compound name, not by m/z."
    )
    st.stop()

name_col = st.selectbox(
    "Compound name column",
    options=name_candidates,
    help="Column holding a human-readable compound name to look up against KEGG.",
)

col1, col2 = st.columns(2)
p_cutoff = col1.number_input("p-adj cutoff (foreground)", value=0.05, step=0.01, format="%.3f")
fc_cutoff = col2.number_input("|log2FC| cutoff (foreground)", value=1.0, step=0.1)

df = statistics_df.dropna(subset=[name_col, "log2FC", "p-adj"]).copy()
df[name_col] = df[name_col].astype(str).str.strip()
df = df[df[name_col].ne("") & df[name_col].ne("nan")]

background_names = set(df[name_col])
foreground_names = set(
    df[(df["p-adj"] < p_cutoff) & (df["log2FC"].abs() > fc_cutoff)][name_col]
)

st.markdown(
    f"**Background:** {len(background_names)} annotated feature(s) &nbsp;|&nbsp; "
    f"**Foreground (significant):** {len(foreground_names)} feature(s)"
)

if len(foreground_names) < 3:
    st.warning(f"Only {len(foreground_names)} significant feature(s) — need at least 3 to test enrichment.")
    st.stop()

if st.button("Run KEGG pathway mapping + enrichment", type="primary"):
    with st.spinner("Querying KEGG REST API..."):
        mapping_df = map_names_to_pathways(
            sorted(background_names), progress_label="Resolving compounds against KEGG"
        )

    if mapping_df.empty:
        st.warning("None of the compound names could be resolved against KEGG.")
        st.stop()

    st.session_state["kegg_mapping_df"] = mapping_df

    resolved = mapping_df["name"].nunique()
    st.success(f"Resolved {resolved}/{len(background_names)} compound name(s) to KEGG.")

    with st.expander("Feature → KEGG compound → pathway mapping table"):
        st.dataframe(mapping_df, use_container_width=True)

    enrichment_df = run_pathway_enrichment(mapping_df, foreground_names, background_names)
    st.session_state["kegg_enrichment_df"] = enrichment_df

    if enrichment_df.empty:
        st.info("No KEGG pathways found among the significant (foreground) features.")
    else:
        st.subheader("Pathway Enrichment Results")
        top = enrichment_df.head(25).copy()
        top["-log10(p)"] = -top["p_value"].apply(lambda p: 0 if p <= 0 else math.log10(p))
        fig = px.bar(
            top.sort_values("-log10(p)"),
            x="-log10(p)",
            y="pathway_name",
            orientation="h",
            hover_data=["fg_count", "bg_count", "p_value"],
            title="KEGG Pathway Enrichment",
        )
        st.plotly_chart(fig, use_container_width=True)
        st.dataframe(enrichment_df, use_container_width=True)

"""KEGG Pathway mapping & enrichment for annotated metabolite features.

Resolves feature compound annotations (from SIRIUS/CSI:FingerID, GNPS, MS2Query or
in-house library matches) to KEGG compound IDs - preferring a structure-based match
(SMILES -> PubChem -> KEGG) over a name-based one when a SMILES column is available -
maps them to KEGG pathways, and runs a Fisher's exact enrichment test comparing
statistically significant features (foreground) against all annotated features
(background).
"""
import math

import pandas as pd
import plotly.express as px
import streamlit as st

from src.common.common import page_setup
from src.common.kegg import map_compounds_to_pathways, run_pathway_enrichment

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

# Optional SMILES column (e.g. MS2Query's "MS2Query_smiles") for a structure-based match,
# tried before falling back to the name-based KEGG search above.
smiles_candidates = [c for c in statistics_df.columns if "smiles" in c.lower()]
smiles_col = None
if smiles_candidates:
    smiles_col = st.selectbox(
        "SMILES column (optional, tried first)",
        options=["(none)"] + smiles_candidates,
        index=1,
        help="If set, each compound is looked up by structure (SMILES -> PubChem -> KEGG) "
        "before falling back to the name-based search. More reliable for systematic/IUPAC "
        "names that don't exist verbatim in KEGG.",
    )
    if smiles_col == "(none)":
        smiles_col = None

col1, col2 = st.columns(2)
p_cutoff = col1.number_input("p-adj cutoff (foreground)", value=0.05, step=0.01, format="%.3f")
fc_cutoff = col2.number_input("|log2FC| cutoff (foreground)", value=1.0, step=0.1)

df = statistics_df.dropna(subset=[name_col, "log2FC", "p-adj"]).copy()
df[name_col] = df[name_col].astype(str).str.strip()
df = df[df[name_col].ne("") & df[name_col].ne("nan")]

# name -> smiles lookup (first non-empty SMILES seen for that name), used when calling
# map_compounds_to_pathways below.
if smiles_col:
    name_to_smiles = (
        df[[name_col, smiles_col]]
        .dropna(subset=[smiles_col])
        .drop_duplicates(subset=[name_col])
        .set_index(name_col)[smiles_col]
        .to_dict()
    )
else:
    name_to_smiles = {}

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

# Each compound resolution is a PubChem + KEGG round trip (~1-5s incl. retries), run
# 4-wide in parallel. A background set of a thousand-plus unique names (common with
# MS2Query analog search) can still take minutes; cap it so a single run has a bounded,
# predictable runtime instead of appearing to hang / drop the connection. Significant
# (foreground) compounds are always kept in full since dropping any of them would bias
# the enrichment test — only the non-significant background is truncated.
MAX_BACKGROUND_COMPOUNDS = 300
capped_background_names = background_names
if len(background_names) > MAX_BACKGROUND_COMPOUNDS:
    extra_slots = max(MAX_BACKGROUND_COMPOUNDS - len(foreground_names), 0)
    other_names = sorted(background_names - foreground_names)[:extra_slots]
    capped_background_names = foreground_names | set(other_names)
    st.warning(
        f"⚠️ Background has {len(background_names)} unique compound(s), which would take a long "
        f"time to resolve against KEGG. Limiting to {len(capped_background_names)} "
        f"(all {len(foreground_names)} significant + {len(other_names)} others) for this run."
    )

if st.button("Run KEGG pathway mapping + enrichment", type="primary"):
    entries = [(name, name_to_smiles.get(name)) for name in sorted(capped_background_names)]
    spinner_msg = (
        "Querying PubChem (structure) + KEGG..." if smiles_col else "Querying KEGG REST API..."
    )
    with st.spinner(spinner_msg):
        mapping_df = map_compounds_to_pathways(
            entries, progress_label="Resolving compounds against KEGG"
        )

    if mapping_df.empty:
        st.warning("None of the compounds could be resolved against KEGG.")
        st.stop()

    st.session_state["kegg_mapping_df"] = mapping_df

    resolved = mapping_df["name"].nunique()
    method_counts = mapping_df.drop_duplicates("name")["match_method"].value_counts().to_dict()
    method_breakdown = ", ".join(f"{v} by {k}" for k, v in method_counts.items())
    st.success(
        f"Resolved {resolved}/{len(capped_background_names)} compound(s) to KEGG "
        f"({method_breakdown})."
    )

    # One row per resolved compound: the SMILES it was looked up with (if any) and the
    # KEGG ID it resolved to, plus whether that came from the structure (SMILES) or the
    # name-based fallback search.
    compound_table = (
        mapping_df.drop_duplicates("name")[["name", "kegg_compound_id", "match_method"]]
        .rename(
            columns={
                "name": "compound_name",
                "kegg_compound_id": "kegg_id",
                "match_method": "matched_by",
            }
        )
        .reset_index(drop=True)
    )
    compound_table.insert(1, "smiles", compound_table["compound_name"].map(name_to_smiles))

    st.subheader("SMILES → KEGG ID Conversion")
    st.dataframe(compound_table, use_container_width=True)
    st.download_button(
        "Download SMILES → KEGG ID table (CSV)",
        data=compound_table.to_csv(index=False).encode("utf-8"),
        file_name="smiles_to_kegg_id.csv",
        mime="text/csv",
    )

    with st.expander("Feature → KEGG compound → pathway mapping table (long format)"):
        st.dataframe(mapping_df, use_container_width=True)
        st.download_button(
            "Download full compound-pathway mapping (CSV)",
            data=mapping_df.to_csv(index=False).encode("utf-8"),
            file_name="kegg_compound_pathway_mapping.csv",
            mime="text/csv",
        )

    enrichment_df = run_pathway_enrichment(mapping_df, foreground_names, capped_background_names)
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
        st.download_button(
            "Download pathway enrichment results (CSV)",
            data=enrichment_df.to_csv(index=False).encode("utf-8"),
            file_name="kegg_pathway_enrichment.csv",
            mime="text/csv",
        )

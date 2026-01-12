import streamlit as st
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Draw

import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from itertools import cycle

from src.common.common import show_fig, load_parquet

COLOR_SCALE = [
    (0.00, "rgba(233, 233, 233, 1.0)"),
    (0.01, "rgba(243, 236, 166, 1.0)"),
    (0.1, "rgba(255, 168, 0, 1.0)"),
    (0.2, "rgba(191, 0, 191, 1.0)"),
    (0.4, "rgba(68, 0, 206, 1.0)"),
    (1.0, "rgba(33, 0, 101, 1.0)"),
]


def add_color_column(df):
    color_cycle = cycle(px.colors.qualitative.Plotly)
    df["color"] = [next(color_cycle) for _ in range(len(df))]
    return df


@st.dialog("🔎 Filter Feature Matrix")
def filter_dialog(df):
    len_unfiltered = len(df)
    mz = st.slider(
        "*m/z* range",
        df["mz"].min(),
        df["mz"].max(),
        value=(df["mz"].min(), df["mz"].max()),
    )
    rt = st.slider(
        "RT range",
        df["RT"].min(),
        df["RT"].max(),
        value=(df["RT"].min(), df["RT"].max()),
    )
    c1, c2 = st.columns(2)
    charge = c1.selectbox(
        "charge state", ["all"] + sorted(df["charge"].unique().tolist())
    )
    adduct = "all"
    if "adduct" in df.columns:
        adduct = c2.selectbox(
            "adduct", ["all"] + sorted(df["adduct"].unique().tolist())
        )
    c1, c2 = st.columns(2)
    filter_annotation = c1.toggle("filter for annotation", False)
    filter_annotation_type = c2.selectbox("annotation type", ["all", "Spectral Matcher", "SIRIUS", "MS2Query", "GNPS"], 0)
    # filter text
    filter_text = ""
    if rt[0] > df["RT"].min():
        filter_text += f" **RT** min = {rt[0]};"
    if rt[1] < df["RT"].max():
        filter_text += f" **RT** max = {rt[1]};"
    if mz[0] > df["mz"].min():
        filter_text += f" ***m/z*** min = {mz[0]};"
    if mz[1] < df["mz"].max():
        filter_text += f" ***m/z*** max = {mz[1]};"
    if filter_annotation:
        filter_text += f" **Annotations:** {filter_annotation_type};"
        cols = ["SpectralMatch", "SIRIUS_", "CSI:FingerID", "CANOPUS", "MS2Query", "GNPS_"]
        if filter_annotation_type == "all":
            mask = cols
        elif filter_annotation_type == "Spectral Matcher":
            mask = [cols[0]]
        elif filter_annotation_type == "SIRIUS":
            mask = cols[1:4]
        elif filter_annotation_type == "MS2Query":
            mask = [cols[4]]
        elif filter_annotation_type == "GNPS":
            mask = [cols[5]]
        df_annotation = df[
            [c for c in df.columns if any([c.startswith(k) for k in mask])]
        ].replace('', pd.NA).dropna(how="all")
        df = df.loc[df_annotation.index, :]
    if charge != "all":
        filter_text += f" **charge** = {charge};"
        df = df[df["charge"] == int(charge)]
    if adduct != "all":
        filter_text += f" **adduct** = {adduct};"
        df = df[df["adduct"] == adduct]
    df = df[(df["mz"] >= mz[0]) & (df["mz"] <= mz[1])]
    df = df[(df["RT"] >= rt[0]) & (df["RT"] <= rt[1])]
    if df.empty:
        st.warning(
            "⚠️ Feature Matrix is empty after filtering. Filter will not be applied."
        )
    _, _, c1, c2 = st.columns(4)
    if c1.button("Cancel", use_container_width=True):
        st.rerun()

    if c2.button("Apply", type="primary", use_container_width=True):
        if len(df) != len_unfiltered and not df.empty:
            st.session_state["feature-matrix-filtered"] = df
            st.session_state["fm-filter-info"] = filter_text.rstrip(";")
        st.rerun()


def metabolite_selection():
    df = load_parquet(
        Path(st.session_state.results_dir, "consensus-dfs", "feature-matrix.parquet")
    )

    if df.empty:
        st.error("FeatureMatrix is empty.")
        return None

    sample_cols = sorted([col for col in df.columns if col.endswith(".mzML")])

    # Insert a column with normalized intensity values to display as barchart column in dataframe
    df.insert(
        1,
        "intensity",
        df.apply(lambda row: [int(row[col]) for col in sample_cols], axis=1),
    )
    df["intensity"] = df["intensity"].apply(
        lambda intensities: [i / max(intensities) for i in intensities]
    )
    c1, c2, c3 = st.columns([0.5, 0.25, 0.25])
    if "feature-matrix-filtered" in st.session_state:
        if c2.button("❌ Reset", use_container_width=True):
            del st.session_state["feature-matrix-filtered"]
            del st.session_state["fm-filter-info"]
            st.rerun()
        st.success(st.session_state["fm-filter-info"])
    if c3.button("🔎 Filter", use_container_width=True):
        filter_dialog(df)
    if "feature-matrix-filtered" in st.session_state:
        if not st.session_state["feature-matrix-filtered"].empty:
            df = st.session_state["feature-matrix-filtered"]
    c1.markdown(f"**Feature Matrix** ({df.shape[0]} metabolites)")

    tab1, tab2 = st.tabs(["✅ **Selection**", "👀 View"])
    with tab2:
        fig = plot_consensus_map(df)
        show_fig(fig, "consensus-map")
    with tab1:
        event = st.dataframe(
            df,
            column_order=["intensity", "RT", "mz", "charge", "adduct"],
            hide_index=False,
            column_config={
                "intensity": st.column_config.BarChartColumn(
                    width="small",
                    help=", ".join(
                        [
                            str(Path(col).stem)
                            for col in sorted(df.columns)
                            if col.endswith(".mzML")
                        ]
                    ),
                ),
            },
            height=300,
            use_container_width=True,
            on_select="rerun",
            selection_mode="single-row",
        )
        rows = event.selection.rows
        if rows:
            return df.iloc[rows[0], :]
        st.info(
            "💡 Select a row (metabolite) in the feature matrix for more information."
        )
        return None


def metabolite_metrics(metabolite):
    cols = st.columns(5)
    with cols[0]:
        st.metric("*m/z* (monoisotopic)", round(metabolite["mz"], 4))
    with cols[1]:
        st.metric("RT (seconds)", round(metabolite["RT"], 1))
    with cols[2]:
        st.metric("charge", metabolite["charge"])
    with cols[3]:
        st.metric("re-quantified", metabolite["re-quantified"])
    with cols[4]:
        if "adduct" in metabolite:
            st.metric("adduct", metabolite["adduct"])


@st.cache_data
def get_chroms_for_each_sample(metabolite):
    # Get index of row in df where "metabolite" is equal to metabolite
    all_samples = [
        i.replace(".mzML", "") for i in metabolite.index if i.endswith("mzML")
    ]
    dfs = []
    for sample in all_samples:
        # Get feature ID for sample
        fid = metabolite[sample + ".mzML_IDs"]
        path = Path(
            st.session_state.results_dir,
            "ffmid-df" if metabolite["re-quantified"] else "ffm-df",
            sample + ".parquet",
        )
        f_df = load_parquet(path)
        if fid in f_df.index:
            f_df = f_df.loc[[fid]]
            f_df["sample"] = [sample]
            dfs.append(f_df.loc[[fid]])
        else:
            dfs.append(pd.DataFrame({"sample": [sample], "chrom_RT": [None], "chrom_intensity": [None]}))

    if not dfs:
        return pd.DataFrame()
    df = pd.concat(dfs)
    return df


@st.cache_resource
def get_feature_chromatogram_plot(df):
    df = df.sort_values("sample")
    df = add_color_column(df)
    # Create an empty figure
    fig = go.Figure()
    # Loop through each row in the DataFrame and add a line trace for each
    for _, row in df.iterrows():
        if row["chrom_RT"] is not None:
            fig.add_trace(
                go.Scatter(
                    x=row["chrom_RT"],  # Assuming chrom_RT is a list of values
                    y=row[
                        "chrom_intensity"
                    ],  # Assuming chrom_intensity is a list of values
                    mode="lines",  # Line plot
                    name=row["sample"],  # Giving each line a name based on its index
                    marker=dict(color=row["color"]),
                )
            )
    # Update layout of the figure
    fig.update_layout(
        xaxis_title="retention time (s)",
        yaxis_title="intensity (counts)",
        plot_bgcolor="rgb(255,255,255)",
        template="plotly_white",
        showlegend=True,
        margin=dict(l=0, r=0, t=0, b=0),
        height=400,
    )
    return fig


@st.cache_resource
def get_feature_intensity_plot(metabolite):
    df = pd.DataFrame(
        {
            "sample": [i for i in metabolite.index if i.endswith(".mzML")]
        }
    )
    df["intensity"] = df["sample"].apply(lambda x: metabolite[x])
    df["sample"] = df["sample"].apply(lambda x: Path(x).stem)
    df = df.sort_values("sample")
    df = add_color_column(df)

    # Create a mapping from sample to color
    color_map = dict(zip(df["sample"], df["color"]))

    # Plot bar chart
    fig = px.bar(
        df, x="sample", y="intensity", color="sample", color_discrete_map=color_map
    )

    fig.update_layout(
        xaxis_title="",
        yaxis_title="intensity (AUC)",
        plot_bgcolor="rgb(255,255,255)",
        template="plotly_white",
        showlegend=True,
        margin=dict(l=0, r=0, t=0, b=0),
        height=400,
    )
    return fig

def spectralmatching_summary(s):
    s.index = [
        i.replace("SpectralMatch_", "").replace("SpectralMatch", "name")
        for i in s.index
    ]
    # Split each value by " ## " and create a DataFrame
    split_data = s.apply(lambda x: x.split(' ## '))

    # Convert to a DataFrame with single entries
    df = pd.DataFrame(split_data.tolist(), index=s.index).T

    # Rename columns to match the original Series index
    df.columns = s.index
    for i, row in df.iterrows():
        row.name = f"Spectra Match #{i+1}"
        c1, c2 = st.columns([0.7, 0.3])
        with c1:
            st.dataframe(row, use_container_width=True)
        with c2:
            if "smiles" in row.index:
                if row["smiles"] != "None":
                    molecule = Chem.MolFromSmiles(row["smiles"])
                    img = Draw.MolToImage(molecule)
                    st.image(img, use_container_width=True)

def sirius_summary(s):
    """s containing SIRIUS, CSI:FingerID and CANOPUS results for selected metabolite"""
    samples = list(set(s.split("_")[1] for s in s.index if s.startswith("SIRIUS_")))
    c1, c2 = st.columns([0.6, 0.4])
    if not samples:
        return
    if len(samples) > 1:
        sample = c1.selectbox("select file", samples)
    else:
        sample = samples[0]

    s = s[[i for i in s.index if f"_{sample}_" in i]]

    s.index = [
        i.replace(f"SIRIUS_{sample}_molecularFormula", "formula (SIRIUS)").replace(f"_{sample}_", "").replace("CANOPUS", "").replace("CSI:FingerID", "").replace("_", " ").replace("molecularFormula", "formula (CSI:FingerID)")
        for i in s.index
    ]
    s.name = f"{sample}; {s.name}"
    s = s.dropna()
    with c1:
        s.index = [i.replace(f"_{sample}_", " ") for i in s.index]
        st.dataframe(s, use_container_width=True)
    with c2:
        if "InChI" in s.index:
            molecule = Chem.MolFromInchi(s["InChI"])
            img = Draw.MolToImage(molecule)
            st.image(img, width=500)


def ms2query_summary(s):
    s = s[s != "nan"]
    s.index = [
        i.lstrip("MS2Query_").replace("_", " ").replace("cf ", "").replace("npc ", "")
        for i in s.index
    ]
    c1, c2 = st.columns([0.6, 0.4])
    with c1:
        st.dataframe(s, use_container_width=True)
    with c2:
        if "smiles" in s.index:
            molecule = Chem.MolFromSmiles(s["smiles"])
            img = Draw.MolToImage(molecule)
            st.image(img, use_container_width=True)


def gnps_summary(s):
    """Display GNPS FBMN library match annotations for selected metabolite."""
    s = s[s != ""]
    if s.empty:
        st.info("No GNPS annotations available for this metabolite.")
        return

    s.index = [
        i.replace("GNPS_", "").replace("_", " ")
        for i in s.index
    ]
    s.name = "GNPS Library Match"

    c1, c2 = st.columns([0.6, 0.4])
    with c1:
        st.dataframe(s, use_container_width=True)
    with c2:
        if "SMILES" in s.index and s["SMILES"]:
            try:
                molecule = Chem.MolFromSmiles(s["SMILES"])
                if molecule:
                    img = Draw.MolToImage(molecule)
                    st.image(img, use_container_width=True)
            except:
                pass
        elif "INCHI" in s.index and s["INCHI"]:
            try:
                molecule = Chem.MolFromInchi(s["INCHI"])
                if molecule:
                    img = Draw.MolToImage(molecule)
                    st.image(img, use_container_width=True)
            except:
                pass


@st.cache_resource
def plot_consensus_map(df):
    fig = go.Figure()

    df["mean"] = df.loc[:, [c for c in df.columns if c.endswith(".mzML")]].mean(axis=1)
    df = df.sort_values("mean")

    # define custom data for hovering
    meta_values = [
        df.index,
        df["mz"].round(5),
        df["RT"].round(),
        df["mean"],
        df["charge"],
        df["quality"],
    ]

    hovertemplate = """
<b>name: %{customdata[0]}<br>
mz: %{customdata[1]}<br>
RT: %{customdata[2]}<br>
intensity: %{customdata[3]}<br>
charge: %{customdata[4]}<br>
quality: %{customdata[5]}<br>
"""

    if "adduct" in df.columns:
        meta_values.append(df["adduct"])
        hovertemplate += "adduct: %{customdata[6]}<br>"

    for s in [col for col in df.columns if col.endswith("mzML")]:
        meta_values.append(df[s])
        hovertemplate += (
            s[:-5] + ": %{customdata[" + str(len(meta_values) - 1) + "]}<br>"
        )

    customdata = np.stack(meta_values, axis=-1)

    fig.add_trace(
        go.Scattergl(
            name="feature",
            x=df["RT"],
            y=df["mz"],
            mode="markers",
            marker_color=df["mean"],
            marker_symbol="square",
            marker_size=12,
            customdata=customdata,
            hovertemplate=hovertemplate,
        )
    )
    fig.update_layout(
        xaxis_title="retention time",
        yaxis_title="m/z",
        showlegend=False,
        plot_bgcolor="rgb(255,255,255)",
        margin=dict(l=0, r=0, t=0, b=0),
    )
    fig.update_traces(
        showlegend=False,
        marker_colorscale=COLOR_SCALE,
        selector=dict(type="scattergl"),
    )

    return fig

def download_section(workflow_dir):
    with st.popover("⬇️ Downloads", use_container_width=True):
        if st.button(
            "Prepare files for download", use_container_width=True
        ):
            with open(
                Path(workflow_dir, "results", "results.zip"), "rb"
            ) as fp:
                st.success("Files are ready.")
                st.download_button(
                    label="⬇️ Download Now",
                    data=fp,
                    file_name="UmetaFlow-results.zip",
                    mime="application/zip",
                    use_container_width=True,
                )
        st.markdown("""**Feature Matrix**:
                        
The main result file with consensus feature meta data (*m/z*, retention time, charge, adduct), intensities, annotations and additional information such as original & consensus feature IDs (for mapping with sample feature map data) and re-quantification status. In **tsv** and **parquet** format.

**Meta Value Table Template**:

A template table for entering sample meta data for downstream processing tools such as the FBMN stats guide. In **tsv** file format.

**Sample Feature Maps:**

Dataframes for feature information on sample level. Containing meta data, intensities and chromatogram data for the monoisotopic mass trace. In two directories: **ffm-df** (from the initial feature detection step with FeatureFinderMetabo) and **ffmid-df** (Re-quantified by FeatureFinderMetaboIdent) In **parquet** file format.
                        
**SIRIUS Input Files**
                        
Input files for SIRIUS in the **sirius-export** directory. Useful to run SIRIUS independently for instance in the GUI app. In **.ms** file format.

**GNPS Input Files**

Input files for GNPS FBMN and IIMN in the **gnps-export** directory.                 
""")
        
def help_section():
    st.markdown("""
**Welcome to the UmetaFlow Result Dashboard!**
                
This page provides a detailed overview of the analysis results for your UmetaFlow run. Here's how to navigate and interpret the information.

**1. Metabolite Selection**

Use checkbox on the left side of the table to select a metabolite from the final FeatureMatrix or take a look at a scatter plot to check out consensus feature distribution.

The table contains an intensity column with bar charts for quick identificantion of differential metabolites. Hover over the column header to see the order of samples in the bar chart.

💡 Use the option to **filter** the data.

**2. 📝 Summary:**

Displays key metrics and statistics for the selected metabolite, such as:

- Retention time
- Mass-to-charge ratio (m/z)
- charge
- adduct
- if it has been re-quantified

**3. Annotations:**

Each of the annotation methods will be displayed in a separate panel (if data is available).

**🎯 Spectral Matching**

MS2 spectral matching based on in-house library. Here, multiple matches can be displayed together with the structure from SMILES (if available). Take a look at the score and mass error in ppm.

**⭐ SIRIUS, CSI:FingerID & CANOPUS**

Formula, structure and compound class predictions. If available, the structure will be displayed based on InChI determined from CSI:FingerID.

**🤖 MS2Query**

Predicted compound names structures and compound classes together with structure drawn from SMILES string.

💡 In general, a score of > 0.7 can be considered a good match. Take a look at the precursor *m/z* difference. While a very small difference indicates an exact match, larger differences are valid and refer to chemical analogues.

**🌐 GNPS FBMN**

Library matches from GNPS Feature-Based Molecular Networking. Import your GraphML results via the GNPS FBMN page to annotate features with compound names, structures (SMILES/InChI), adducts, and match quality scores.

**4. Quantification: 📈 Chromatograms & 📊 Intensities**

This section shows quantitative differences in your samples. Check out the chromatograms to determine valid feature detection settings. E.g. chromatograms should be complete and not cut off.
                    
**Tips**

- Expand or collapse sections using the 🔽 icons
- Hover over charts for detailed tooltips
- Hovering over charts opens a menu in the top right corner to save the image or make them full screen
- Tables offer a search and download button in the top right corner
- If new results don't appear, try to reset the filter (can contain old results)
""")
import streamlit as st
from pathlib import Path
import shutil
import tempfile
import zipfile

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from itertools import cycle

import json

from openms_insight import LinePlot, MirrorPlot, StateManager

from src.common.common import show_fig, load_parquet, export_folder_input
from src.common.insight_tables import cache_id_for, cache_ready, drop_stale_caches
from src.feature_explorer import FeatureExplorer
from src.feature_tags import all_tags, has_any_tag, load_tags, set_feature_tags
from src.ms2_viewer import LIBRARY_SIDE, QUERY_SIDE, load_ms2_data, peaks_frame


def add_color_column(df):
    color_cycle = cycle(px.colors.qualitative.Plotly)
    df["color"] = [next(color_cycle) for _ in range(len(df))]
    return df


@st.dialog("🔎 Filter Feature Matrix")
def filter_dialog(df, tags=None):
    tags = tags or {}
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
    filter_annotation_type = c2.selectbox("annotation type", ["all", "Spectral Matcher", "SIRIUS", "MS2Query"], 0)
    c1, c2 = st.columns(2)
    tag_choice = c1.multiselect(
        "tagged with any of",
        all_tags(tags),
        placeholder="all features" if tags else "no tags yet",
        disabled=not tags,
        help="Tag features under the table: select a feature, then add tags such as contaminant or standard.",
    )
    exclude_tags = c2.multiselect(
        "leave out features tagged",
        all_tags(tags),
        placeholder="none" if tags else "no tags yet",
        disabled=not tags,
    )
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
        cols = ["SpectralMatch", "SIRIUS_", "CSI:FingerID", "CANOPUS", "MS2Query"]
        if filter_annotation_type == "all":
            mask = cols
        elif filter_annotation_type == "Spectral Matcher":
            mask = [cols[0]]
        elif filter_annotation_type == "SIRIUS":
            mask = cols[1:4]
        elif filter_annotation_type == "MS2Query":
            mask = [cols[-1]]
        df_annotation = df[
            [c for c in df.columns if any([c.startswith(k) for k in mask])]
        ].replace('', pd.NA).dropna(how="all")
        df = df.loc[df_annotation.index, :]
    if tag_choice:
        filter_text += f" **tags:** {', '.join(tag_choice)};"
        df = df[has_any_tag(df, tags, tag_choice)]
    if exclude_tags:
        filter_text += f" **without tags:** {', '.join(exclude_tags)};"
        df = df[~has_any_tag(df, tags, exclude_tags)]
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
    """Draw the feature table and return its FeatureExplorer (None if there are no features)."""
    df = load_parquet(
        Path(st.session_state.results_dir, "consensus-dfs", "feature-matrix.parquet")
    )

    if df.empty:
        st.error("FeatureMatrix is empty.")
        return None
    tags = load_tags(st.session_state.results_dir, df)

    c1, c2, c3 = st.columns([0.5, 0.25, 0.25])
    if "feature-matrix-filtered" in st.session_state:
        if c2.button("❌ Reset", use_container_width=True):
            del st.session_state["feature-matrix-filtered"]
            del st.session_state["fm-filter-info"]
            st.rerun()
        st.success(st.session_state["fm-filter-info"])
    if c3.button("🔎 Filter", use_container_width=True):
        filter_dialog(df, tags)
    if "feature-matrix-filtered" in st.session_state:
        if not st.session_state["feature-matrix-filtered"].empty:
            df = st.session_state["feature-matrix-filtered"]
    c1.markdown(f"**Feature Matrix** ({df.shape[0]} metabolites)")
    st.caption(
        "Click a feature in the table or in the feature map (retention time vs *m/z*, "
        "colored by the highest intensity across samples) to see its details. "
        "Column headers sort and filter the table; its download button exports the filtered rows."
    )

    explorer = FeatureExplorer(
        df, cache_dir=Path(st.session_state.results_dir, "insight-cache"), tags=tags
    )
    explorer.table()
    return explorer


def feature_tags_editor(metabolite):
    """Tags of the selected feature (MZmine row tags), saved with the results."""
    results_dir = st.session_state.results_dir
    tags = load_tags(results_dir)
    current = tags.get(str(metabolite.name), [])
    chosen = st.multiselect(
        f"🏷️ Tags of feature {metabolite.name}",
        options=sorted(set(all_tags(tags)) | set(current), key=str.casefold),
        default=current,
        accept_new_options=True,
        placeholder="Add a tag, e.g. contaminant, standard, follow up",
        key=f"feature-tags-{metabolite.name}",
        help="Tags show in the table's tags column and can be filtered on with 🔎 Filter. "
        "They are kept with these results and dropped when the workflow runs again.",
    )
    if chosen != current:
        set_feature_tags(results_dir, metabolite, chosen)
        st.rerun()


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

def _theme_type():
    """Active app theme, "dark" or "light" (light before the first rerun)."""
    try:
        return "dark" if st.context.theme.type == "dark" else "light"
    except Exception:
        return "light"


def draw_structure(molecule, theme="light", size=(400, 300)):
    """PNG bytes of a 2D structure on a transparent background.

    Dark theme uses RDKit's dark-mode palette (light bonds, brightened heteroatom
    colours), so the drawing stays legible on Streamlit's dark background instead
    of showing as a white tile. Returns None for a molecule that failed to parse.
    """
    if molecule is None:
        return None
    drawer = rdMolDraw2D.MolDraw2DCairo(*size)
    opts = drawer.drawOptions()
    if theme == "dark":
        rdMolDraw2D.SetDarkMode(opts)
    opts.setBackgroundColour((0, 0, 0, 0))
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, molecule)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def show_structure(molecule, **image_kwargs):
    """Show a structure drawing in the colours of the active theme."""
    png = draw_structure(molecule, _theme_type())
    if png is None:
        st.caption("Structure could not be drawn.")
        return
    st.image(png, **image_kwargs)


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
                    show_structure(
                        Chem.MolFromSmiles(row["smiles"]), use_container_width=True
                    )

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
            show_structure(Chem.MolFromInchi(s["InChI"]), width=500)


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
            show_structure(Chem.MolFromSmiles(s["smiles"]), use_container_width=True)


RESULT_FILES_INFO = """**Feature Matrix**:
                        
The main result file with consensus feature meta data (*m/z*, retention time, charge, adduct), intensities, annotations and additional information such as original & consensus feature IDs (for mapping with sample feature map data) and re-quantification status. In **tsv** and **parquet** format.

**Meta Value Table Template**:

A template table for entering sample meta data for downstream processing tools such as the FBMN stats guide. In **tsv** file format.

**Sample Feature Maps:**

Dataframes for feature information on sample level. Containing meta data, intensities and chromatogram data for the monoisotopic mass trace. In two directories: **ffm-df** (from the initial feature detection step with FeatureFinderMetabo) and **ffmid-df** (Re-quantified by FeatureFinderMetaboIdent) In **parquet** file format.
                        
**SIRIUS Input Files**
                        
Input files for SIRIUS in the **sirius-export** directory. Useful to run SIRIUS independently for instance in the GUI app. In **.ms** file format.

**GNPS Input Files**

Input files for GNPS FBMN and IIMN in the **gnps-export** directory.                 
"""


def download_section(workflow_dir):
    zip_path = Path(workflow_dir, "results", "results.zip")
    # A local install writes the result files straight into a folder.
    if st.session_state.location == "local":
        with st.popover("📤 Export", use_container_width=True):
            export_dir = export_folder_input("umetaflow-results")
            no_folder = not st.session_state["export_dir"].strip()
            if no_folder:
                st.warning("Enter an export folder.")
            if st.button(
                "📤 Export", type="primary", use_container_width=True,
                disabled=no_folder or not zip_path.exists(),
            ):
                target = Path(export_dir, "UmetaFlow-results")
                try:
                    with zipfile.ZipFile(zip_path) as zf:
                        bad = zf.testzip()
                        if bad is not None:
                            raise zipfile.BadZipFile(f"corrupt member {bad}")
                        files = [i for i in zf.infolist() if not i.is_dir()]
                        # Extract next to the target first, so a failure
                        # leaves the export folder as it was.
                        target.parent.mkdir(parents=True, exist_ok=True)
                        with tempfile.TemporaryDirectory(
                            dir=target.parent, prefix=".umetaflow-export-"
                        ) as tmp:
                            zf.extractall(tmp)
                            shutil.copytree(tmp, target, dirs_exist_ok=True)
                    st.success(f"Exported {len(files)} files to **{target}**.")
                except (OSError, zipfile.BadZipFile) as e:
                    st.error(f"Export failed: {e}")
            st.markdown(RESULT_FILES_INFO)
        return
    with st.popover("⬇️ Downloads", use_container_width=True):
        if st.button(
            "Prepare files for download", use_container_width=True
        ):
            with open(zip_path, "rb") as fp:
                st.success("Files are ready.")
                st.download_button(
                    label="⬇️ Download Now",
                    data=fp,
                    file_name="UmetaFlow-results.zip",
                    mime="application/zip",
                    use_container_width=True,
                )
        st.markdown(RESULT_FILES_INFO)

def help_section():
    st.markdown("""
**Welcome to the UmetaFlow Result Dashboard!**
                
This page provides a detailed overview of the analysis results for your UmetaFlow run. Here's how to navigate and interpret the information.

**1. Feature table and feature map**

Click a row in the feature table, or a point in the feature map (retention time vs *m/z*), to select a metabolite. Both views stay in sync. The table has one intensity column per sample, the number of samples a feature was found in, and which tools annotated it. Sort or filter by any column header, and download the filtered table as CSV.

Tag the selected feature under the table (for example *contaminant*, *standard* or *follow up*). Tags appear in the table's **tags** column, and **🔎 Filter** shows or leaves out features by tag. Tags are saved with the results and dropped when the workflow runs again.

The chromatograms of the selected feature in every sample are shown next to the feature map.

💡 Use the option to **filter** the data.

**2. 📝 Summary:**

Displays key metrics and statistics for the selected metabolite, such as:

- Retention time
- Mass-to-charge ratio (m/z)
- charge
- adduct
- if it has been re-quantified

**🔬 MS2 spectrum**

The MS2 spectrum recorded for the feature. If the in-house library matched it, the library spectrum is mirrored below it, with the peaks that pair up highlighted, the match score, the cosine similarity and the number of matched peaks.

**3. Annotations:**
                    
Each of the three annotation methods will be displayed in a separate panel (if data is available).
                    
**🎯 Spectral Matching**

MS2 spectral matching based on in-house library. Here, multiple matches can be displayed together with the structure from SMILES (if available). Take a look at the score and mass error in ppm.

**⭐ SIRIUS, CSI:FingerID & CANOPUS**

Formula, structure and compound class predictions. If available, the structure will be displayed based on InChI determined from CSI:FingerID.
                    
**🤖 MS2Query**

Predicted compound names structures and compound classes together with structure drawn from SMILES string.

💡 In general, a score of > 0.7 can be considered a good match. Take a look at the precursor *m/z* difference. While a very small difference indicates an exact match, larger differences are valid and refer to chemical analogues.

**4. Quantification: 📈 Chromatograms & 📊 Intensities**

This section shows quantitative differences in your samples. Check out the chromatograms to determine valid feature detection settings. E.g. chromatograms should be complete and not cut off.
                    
**Tips**

- Expand or collapse sections using the 🔽 icons
- Hover over charts for detailed tooltips
- Hovering over charts opens a menu in the top right corner to save the image or make them full screen
- Tables offer a search and download button in the top right corner
- If new results don't appear, try to reset the filter (can contain old results)
""")

@st.cache_data(show_spinner=False)
def _load_ms2_data(results_dir: str, library: str, tol_ppm: float, mtimes: tuple):
    # mtimes only key the cache, so a re-run of the workflow invalidates it
    return load_ms2_data(Path(results_dir), Path(library) if library else None, tol_ppm)


def _ms2_inputs(workflow_dir: Path):
    """MS2 library file and fragment tolerance (ppm) the spectral matching used."""
    library_dir = Path(workflow_dir, "input-files", "ms2-library")
    libraries = sorted(p for p in library_dir.glob("*") if p.is_file() and p.name != "external_files.txt") if library_dir.is_dir() else []
    tol = 500.0  # MetaboliteSpectralMatcher default
    for name in ("params-translated.json", "params.json"):
        path = Path(workflow_dir, name)
        if path.is_file():
            try:
                msm = json.loads(path.read_text()).get("MetaboliteSpectralMatcher", {})
            except ValueError:
                continue
            if msm.get("algorithm:mass_error_unit", "ppm") == "ppm":
                tol = float(msm.get("algorithm:frag_mass_error_value", tol))
            break
    return (libraries[0] if libraries else None), tol


def ms2_section(metabolite, workflow_dir: Path) -> None:
    """MS2 spectrum of the selected feature; library hits as a mirror plot."""
    results_dir = Path(workflow_dir, "results")
    library, tol = _ms2_inputs(Path(workflow_dir))
    sources = [results_dir / "gnps-export" / "MS2.mgf", results_dir / "spectral-matcher" / "MS2-matches.mzTab"]
    if library:
        sources.append(library)
    mtimes = tuple(p.stat().st_mtime_ns if p.exists() else 0 for p in sources)
    data = _load_ms2_data(str(results_dir), str(library or ""), tol, mtimes)
    if data is None:
        return
    with st.expander("🔬 **MS2 spectrum**", expanded=True):
        query_ids = data.feature_spectra.get(str(metabolite.name), [])
        matches = data.matches.get(str(metabolite.name), [])
        if not query_ids:
            st.info("No MS2 spectrum was recorded for this feature.")
            return

        cache_dir = Path(results_dir, "insight-cache")
        peaks = peaks_frame(data)
        cache_id = cache_id_for("ms2_peaks", data.peaks)
        state = StateManager("umetaflow_ms2")
        mirrored = [m for m in matches if m.library_id is not None]
        labels = {row.spectrum_id: row.label for row in data.spectra.itertuples()}

        if mirrored:
            options = list(range(len(mirrored)))
            pick = 0
            if len(mirrored) > 1:
                pick = st.selectbox(
                    "library hit", options,
                    format_func=lambda i: f"{mirrored[i].name} (score {mirrored[i].score:g})",
                    key="ms2-library-hit",
                )
            match = mirrored[pick]
            mirror_id = f"mirror_{cache_id}"
            if not cache_ready(cache_dir, mirror_id):
                drop_stale_caches(cache_dir, "mirror_ms2_peaks", mirror_id)
                MirrorPlot(
                    cache_id=mirror_id, data=peaks.lazy(), cache_path=str(cache_dir),
                    filters_top={QUERY_SIDE: "spectrum_id"}, filters_bottom={LIBRARY_SIDE: "spectrum_id"},
                    filter_defaults_top={QUERY_SIDE: -1}, filter_defaults_bottom={LIBRARY_SIDE: -1},
                    interactivity={"ms2_peak": "peak_id"},
                    x_column="mz", y_column="intensity",
                    highlight_column="matched", annotation_column="label",
                    title_top="feature MS2", title_bottom="library spectrum",
                    x_label="m/z", y_label="relative intensity (%)",
                )
            state.set_selection(QUERY_SIDE, int(match.query_id))
            state.set_selection(LIBRARY_SIDE, int(match.library_id))
            c1, c2, c3, c4 = st.columns(4)
            c1.metric("library hit", match.name)
            c2.metric("match score", f"{match.score:g}" if match.score is not None else "–")
            c3.metric("cosine", f"{match.cosine:.2f}" if match.cosine is not None else "–")
            c4.metric("matched peaks", match.matched_peaks)
            st.caption(
                f"Top: the feature's MS2 spectrum ({labels.get(match.query_id, '')}). Bottom, mirrored: the library "
                f"spectrum of {match.name}. Peaks that pair up within {tol:g} ppm are highlighted (labeled from 5 % intensity). "
                "The cosine uses square-root intensities of the paired peaks."
            )
            MirrorPlot(cache_id=mirror_id, cache_path=str(cache_dir))(
                key="umetaflow_ms2_mirror", state_manager=state, height=420
            )
        else:
            pick = query_ids[0]
            if len(query_ids) > 1:
                pick = st.selectbox("spectrum", query_ids, format_func=lambda i: labels.get(i, str(i)),
                                    key="ms2-query-spectrum")
            line_id = f"ms2line_{cache_id}"
            if not cache_ready(cache_dir, line_id):
                drop_stale_caches(cache_dir, "ms2line_ms2_peaks", line_id)
                LinePlot(
                    cache_id=line_id, data=peaks.lazy(), cache_path=str(cache_dir),
                    filters={QUERY_SIDE: "spectrum_id"}, filter_defaults={QUERY_SIDE: -1},
                    interactivity={"ms2_peak": "peak_id"},
                    x_column="mz", y_column="intensity",
                    title="feature MS2", x_label="m/z", y_label="relative intensity (%)",
                )
            state.set_selection(QUERY_SIDE, int(pick))
            if matches:
                st.caption(f"Library hit {matches[0].name}: its spectrum was not found in the uploaded library.")
            else:
                st.caption("No spectral library hit for this feature.")
            LinePlot(cache_id=line_id, cache_path=str(cache_dir))(
                key="umetaflow_ms2_line", state_manager=state, height=380
            )

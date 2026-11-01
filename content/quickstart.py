import streamlit as st
from src.common.common import page_setup
from pathlib import Path
import json

params = page_setup(page="main")

st.title("A universal metabolomics tool")
cols = st.columns([0.8, 0.1])
cols[0].markdown(
    """
This app offers the powerful UmetaFlow **[1]** pipeline for untargeted metabolomics in an accessible user interface. Raw data pre-processing converts raw data to a feature quantification table by feature detection, alignment, grouping, adduct annotation and optional re-quantification of missing values. Features can be annotated by in-house libraries based on MS1 m/z and retention time matching as well as MS2 fragment spectrum similarity as well as with formula, structure and compound classes by SIRIUS **[2]**, CSI:FingerID **[3]** & CANOPUS **[4]** and chemical analogues by MS2Query **[5]**. Furthermore, required input files for GNPS Feature Based Molecular Networking **[6]** and Ion Identity Molecular Networking **[7]** can be generated. Besides the untargeted pipeline, this app offers some additional features, such as an interface to explore raw data and metabolite identification and quantification via extracted ion chromatograms based on exact m/z values generated conveniently by an included m/z calculator. For downstream processing statistical analysis can be performed within the app or in the popular [FBmn STATS GUIde](https://github.com/axelwalter/streamlit-metabolomics-statistics) for statistical analyis of metabolomics data **[8]**.

UmetaFlow is further implemented as a [snakemake pipeline](https://github.com/NBChub/snakemake-UmetaFlow) and as a Python version in [Jupyter notebooks](https://github.com/eeko-kon/pyOpenMS_UmetaFlow) based on [pyOpenMS](https://pyopenms.readthedocs.io/en/latest/index.html).
    """)

cols[1].image("assets/umetaflow-logo.png")
cols[1].image("assets/pyopenms-logo.png")


if Path("UmetaFlow-Setup.exe").exists():
    st.markdown(
            """
## Installation

### Download for Windows

Download and run the installer. It installs UmetaFlow for your user without admin rights; start it from the desktop icon or the start menu.
        """)
    with open("UmetaFlow-Setup.exe", "rb") as file:
        st.download_button(
                    label="Download for Windows",
                    data=file,
                    file_name="UmetaFlow-Setup.exe",
                    mime="application/vnd.microsoft.portable-executable",
                    type="primary"
                )
# st.image("assets/umetaflow-app-overview.png", width=800)

st.markdown(
    """
## Workspaces
On the left side of this page you can define a workspace where all your data including uploaded `mzML` files will be stored. Entering a workspace will switch to an existing one or create a new one if it does not exist yet. In the web app, you can share your results via the unique workspace ID. Be careful with sensitive data, anyone with access to this ID can view your data.

## 📁 File Handling
Upload `mzML` files via the **File Upload** tab. The data will be stored in your workspace. With the web app you can upload only one file at a time.
Locally there is no limit in files. However, it is recommended to upload large number of files by specifying the path to a directory containing the files.

Your uploaded files will be shown in the sidebar of all tabs dealing with the files, e.g. the **Metabolomics** tab. Checked file names will be used for analysis.

Result files are available via specified download buttons or, if run locally, within the workspace directory.


## UmetaFlow: Untargeted Quantification & Identification

1. **Pre-Processing**
This initial step transforms raw data into a  feature quantification matrix (FeatureMatrix) containing metabolite intensities across all samples. Includes steps such as RT alignment, adduct detection and mapping of MS2 spectra top features for subsequent annotation.

2. **Re-Quantification**
Unique to UmetaFlow, this optional step is designed to detect low-abundance features that may be missed during pre-processing (e.g. features below the noise-threshold intensity). By re-evaluating the data for subtle signals, this step enhances sensitivity, improving the quantification of minor metabolites.

3. **Annotation**
UmetaFlow offers complementary annotation methods in order to maximize meaningful biological interpretation. Best annotations are from spectral matching with in-house libraries, generated from authentic standards, followed by querying publicly available databases for metabolite identification and structure prediction by CSI:FingerID and GNPS. If no specific compound was found, analog search by MS2Query and compound class prediction by CANOPUS can still provide valuable insights.

## Targeted Quantification: Extracted Ion Chromatograms

#### 📟 m/z Calculator

The m/z calculator facilitates the calculation of mass-to-charge ratios (m/z) for metabolites and includes a method to easily combine metabolites into large molecules.

This table can be used as input for the Extracted Ion Chromatograms workflow.

#### 🔍 Extracted Ion Chromatograms

Simple workflow for the extraction of chromatograms by `m/z` (and optionally `RT` range) value. Produces a **Feature Matrix** file with area under the curve intensities as well as a **Meta Data** template and the chromatogram data for each file.

Area intensities of different variants (e.g. adducts or neutral losses) of a metabolite can be combined. Put a `#` with the name first and variant second (e.g. `glucose` and `glucose#[M+Na]+`).  

### Downstream Processing

#### 📈 Statistics

We recommend the [FBMN stats guide](https://fbmn-statsguide.gnps2.org/) for statistical analysis of metabolomics data.

**[1]** Kontou, Eftychia E., et al. "UmetaFlow: an untargeted metabolomics workflow for high-throughput data processing and analysis." Journal of Cheminformatics 15.1 (2023): 52**.

**[2]**	Dührkop K, Fleischauer M, Ludwig M, Aksenov AA, Melnik AV, Meusel M, et al. SIRIUS 4: a rapid tool for turning tandem mass spectra into metabolite structure information. Nat Methods 2019;16:299–302. https://doi.org/10.1038/s41592-019-0344-8.

**[3]**	Dührkop K, Shen H, Meusel M, Rousu J, Böcker S. Searching molecular structure databases with tandem mass spectra using CSI:FingerID. Proc Natl Acad Sci 2015;112:12580–5. https://doi.org/10.1073/pnas.1509788112.

**[4]**	Dührkop K, Nothias L-F, Fleischauer M, Reher R, Ludwig M, Hoffmann MA, et al. Systematic classification of unknown metabolites using high-resolution fragmentation mass spectra. Nat Biotechnol 2021;39:462–71. https://doi.org/10.1038/s41587-020-0740-8.

**[5]** de Jonge NF, Louwen JJR, Chekmeneva E, Camuzeaux S, Vermeir FJ, Jansen RS, et al. MS2Query: reliable and scalable MS2 mass spectra-based analogue search. Nat Commun 2023;14:1752. https://doi.org/10.1038/s41467-023-37446-4.

**[6]**	Nothias L-F, Petras D, Schmid R, Dührkop K, Rainer J, Sarvepalli A, et al. Feature-based molecular networking in the GNPS analysis environment. Nat Methods 2020;17:905–8. https://doi.org/10.1038/s41592-020-0933-6.

**[7]**	Schmid R, Petras D, Nothias L-F, Wang M, Aron AT, Jagels A, et al. Ion identity molecular networking for mass spectrometry-based metabolomics in the GNPS environment. Nat Commun 2021;12:3832. https://doi.org/10.1038/s41467-021-23953-9.

**[8]** Shah, Abzer K. Pakkir, et al. "The Hitchhiker’s Guide to Statistical Analysis of Feature-based Molecular Networks from Non-Targeted Metabolomics Data." (2023).
    """)

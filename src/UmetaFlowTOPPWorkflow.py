import streamlit as st
from pathlib import Path
from .workflow.WorkflowManager import WorkflowManager

import pandas as pd

import shutil
import json
import sys


from src.metabolomicsresults import *


class Workflow(WorkflowManager):
    # Setup pages for upload, parameter, execution and results.
    # For layout use any streamlit components such as tabs (as shown in example), columns, or even expanders.
    def __init__(self, workspace) -> None:
        flag_file = Path(workspace, "umetaflow-expert-flag.txt")

        name = "UmetaFlow"
        self.expert_mode = False
        if flag_file.exists():
            self.expert_mode = True
            name = "UmetaFlow-Expert"

        # Initialize the parent class with the workflow name.
        super().__init__(name, workspace)

    def upload(self) -> None:
        return

    def add_sirius_path_to_session_state(self):
        if "sirius-path" not in st.session_state:
            possible_paths = [  # potential SIRIUS locations in increasing priority
                str(Path("sirius")),  # anywhere
                str(Path(sys.prefix, "bin", "sirius")),  # in current conda environment
                str(
                    Path(".", "sirius", "sirius.exe")
                ),  # in case of Windows executables
            ]
            st.session_state["sirius-path"] = ""
            for path in possible_paths:
                if shutil.which(path) is not None:
                    st.session_state["sirius-path"] = path

    def configure_simple(self) -> None:
        cols = st.columns(4)
        with cols[0]:
            self.ui.input_widget(
                "ion_mode",
                "positive",
                "**ion mode**",
                options=["positive", "negative"],
                help="MS instrument ion mode.",
            )
        with cols[1]:
            self.ui.input_widget(
                "mz_tolerance",
                10.0,
                "***m/z* tolerance**",
                min_value=0.1,
                max_value=50,
                step_size=5,
                help="MS instrument mass accuracy in ppm.",
            )
        with cols[2]:
            self.ui.input_widget(
                "RT_tolerance",
                30.0,
                "**RT tolerance**",
                step_size=5,
                min_value=5,
                max_value=120,
                help="RT tolerance for feature linking and re-quantification.",
            )
        with cols[3]:
            self.ui.input_widget(
                "num_threads",
                1,
                "number of threads",
                step_size=1,
                min_value=1,
                max_value=20,
                help="Number of threads to use for computation. On a standard laptop it is recommended to use the default setting (1 thread). Can be increased on powerful machines to reduce analysis time significantly.",
            )
        tabs = st.tabs(
            ["⚙️ **Pre-Processing**", "🔎 **Re-Quantification**", "🏷️ **Annotation**"]
        )
        with tabs[0]:
            st.markdown(
                "**Feature Detection**",
                help="Untargeted extraction of metabolic features from raw peak maps.",
            )
            cols = st.columns(4)

            with cols[0]:
                self.ui.input_widget(
                    "ffm:algorithm:common:noise_threshold_int",
                    1000.0,
                    "noise_threshold_intensity",
                    min_value=10,
                    max_value=100000,
                    step_size=100,
                    help="Intensity threshold below which peaks are regarded as noise.",
                )
            with cols[1]:
                self.ui.input_widget(
                    "ffm:algorithm:common:chrom_peak_snr",
                    3.0,
                    "chrom_peak_snr",
                    min_value=1,
                    max_value=10,
                    step_size=1,
                    help="Minimum signal-to-noise a mass trace should have.",
                )
            with cols[2]:
                self.ui.input_widget(
                    "ffm:algorithm:common:chrom_fwhm",
                    5.0,
                    "chrom_fwhm",
                    min_value=1,
                    max_value=60,
                    step_size=5,
                    help="Expected chromatographic peak width (in seconds).",
                )
            with cols[3]:
                self.ui.input_widget(
                    "ffm:algorithm:ffm:remove_single_traces",
                    "true",
                    "remove_single_traces",
                    options=["true", "false"],
                    help="Remove unassembled traces (single traces).",
                )

            self.ui.input_widget(
                "adduct-detection",
                False,
                "**Adduct Detection (optional)**",
                help="Detect adduct configuration in features.",
            )
            cols = st.columns(2)
            with cols[0]:
                self.ui.input_widget(
                    "adducts_pos",
                    "H:+:0.6 Na:+:0.1 NH4:+:0.1 H-1O-1:+:0.1 H-3O-2:+:0.1",
                    "adducts_pos",
                    help="Potential adducts in positive mode. Format: adduct:charge:probability.",
                )
            with cols[1]:
                self.ui.input_widget(
                    "adducts_neg",
                    "H-1:-:1 H-2O-1:0:0.05 CH2O2:0:0.5",
                    "adducts_neg",
                    help="Potential adducts in negative mode. Format: adduct:charge:probability.",
                )
            with st.columns(2)[0].container(border=True):
                st.image(str(Path("assets", "metabolomics-preprocessing.png")))
        with tabs[1]:
            self.ui.input_widget(
                "requantify",
                False,
                "**Re-quantify** features with missing values",
                help="Re-quantify consensus features in the FeatureMatrix with at least one missing value.",
            )
            with st.columns(2)[0].container(border=True):
                st.image(str(Path("assets", "requant.png")))
        with tabs[2]:
            self.ui.input_widget(
                "annotate-ms2",
                False,
                "MS2 spectral library matching with **in-house library**",
                help="Based on MS2 spectrum similarity.",
            )
            cols = st.columns([0.75, 0.25])
            with cols[0]:
                self.ui.simple_file_uploader(
                    "ms2-library", "mgf", "MS2 library in mgf format"
                )

            with cols[1]:
                st.image(str(Path("assets", "OpenMS.png")), width=200)
            st.divider()
            cols = st.columns([0.75, 0.25])
            with cols[0]:
                self.ui.input_widget(
                    "export-sirius",
                    True,
                    "export input files for **SIRIUS**",
                    help="Only export input files for SIRIUS, without executing SIRIUS.",
                )
                self.ui.input_widget(
                    "run-sirius",
                    False,
                    "predict compound **chemical formulas** with SIRIUS",
                    help="Run SIRIUS for chemical formula annoation.",
                )
                self.ui.input_widget(
                    "run-fingerid",
                    False,
                    "predict compound **structures and classes** with CSI:FingerID & CANOPUS",
                    help="SIRIUS, CSI:FingerID, and CANOPUS are powerful tools for molecular annotation from high-resolution tandem mass spectrometry (MS/MS) data. SIRIUS focuses on molecular formula elucidation. It uses isotope patterns and fragmentation spectra to construct fragmentation trees, which explain how precursor ions break down into smaller fragments. These trees are scored against theoretical formulas in order to determine the best matching formula. CSI:FingerID predicts molecular structures based on molecular fingerprints derived from the fragmentation trees, which encode chemical properties such as functional groups. These fingerprints are matched against large chemical databases (e.g. PubChem) for compound identification. CANOPUS complements these tools by predicting compound classes directly from MS2 spectra and does not require a database. CANOPUS classifies molecules into hierarchical categories (e.g. lipids, flavonoids), which provides valuable biological context, even when features remain unidentified by CSI:FingerID.",
                )
            with cols[1]:
                st.image(str(Path("assets", "sirius.png")), width=200)
            # SiriusExport
            self.add_sirius_path_to_session_state()

            if st.session_state["sirius-path"]:
                cols = st.columns([0.25, 0.25, 0.5])
                with cols[0]:
                    self.ui.input_widget(
                        "sirius-user-email",
                        "",
                        "SIRIUS user E-mail",
                        help="Email address from a valid SIRIUS account.",
                        widget_type="text",
                    )
                with cols[1]:
                    self.ui.input_widget(
                        "sirius-user-password",
                        "",
                        "SIRIUS user password",
                        help="Password from a valid SIRIUS account. **Not encrypted**, will be stored **unencrypted in plain text** in parameters and show up in log files.",
                        widget_type="password",
                    )
                cols = st.columns(4)
                with cols[0]:
                    self.ui.input_widget(
                        "sirius-profile",
                        name="profile",
                        options=["default", "qtof", "orbitrap", "fticr"],
                        default="default",
                        help="Name of the configuration profile",
                    )
                with cols[1]:
                    self.ui.input_widget(
                        "sirius-maxmz",
                        500,
                        "max *m/z*",
                        min_value=50,
                        max_value=1000,
                        step_size=50,
                        help="Only considers compounds with a precursor m/z lower or equal for **structure predictions**. All other compounds in the input will be skipped. Recommended to be below 600, otherwise very long execution times are expected.",
                    )
                db_options = [
                    "none",
                    "ALL",
                    "ALL_BUT_INSILICO",
                    "BIO",
                    "PUBCHEM",
                    "MESH",
                    "HMDB",
                    "KNAPSACK",
                    "CHEBI",
                    "PUBMED",
                    "KEGG",
                    "HSDB",
                    "MACONDA",
                    "METACYC",
                    "GNPS",
                    "ZINCBIO",
                    "UNDP",
                    "YMDB",
                    "PLANTCYC",
                    "NORMAN",
                    "ADDITIONAL",
                    "PUBCHEMANNOTATIONBIO",
                    "PUBCHEMANNOTATIONDRUG",
                    "PUBCHEMANNOTATIONSAFETYANDTOXIC",
                    "PUBCHEMANNOTATIONFOOD",
                    "KEGGMINE",
                    "ECOCYCMINE",
                    "YMDBMINE",
                ]
                with cols[2]:
                    self.ui.input_widget(
                        "sirius-db",
                        "none",
                        "formula database",
                        options=db_options,
                        help="Search formulas in the given database. If no database is given all possible molecular formulas will be respected (no database is used).",
                    )
                with cols[3]:
                    self.ui.input_widget(
                        "sirius-structure-db",
                        "BIO",
                        "structure database",
                        options=db_options,
                        help="Search structure in the given database.",
                    )
            else:
                st.columns(2)[0].info(
                    "💡 To run SIRIUS within UmetaFlow, install the command line tool."
                )
            st.divider()
            cols = st.columns([0.75, 0.25])
            with cols[0]:
                self.ui.input_widget(
                    "export-gnps",
                    True,
                    "export input files for **GNPS FBMN & IIMN**",
                    help="GNPS (Global Natural Products Social Molecular Networking) is an open-access platform designed for the analysis, sharing, and annotation of MS2 spectra, particularly in natural products research. By leveraging community-contributed spectral libraries, GNPS enables the identification of metabolites and the exploration of related molecules in networks. Feature-based molecular networking (FBMN) builds networks where nodes represent molecular features and edges indicate spectral similarity. This approach provides a visual and interactive way to explore the relationships between metabolites, aiding in the identification of unknown compounds which have identified neighbours. If adduct detection was enabled during pre-processing, Ion identity molecular networking (IIMN) can reduce the complexity in FBMN by collapsing different ion species of the same metabolite, which would have otherwise appeared as separate nodes. UmetaFlow exports all necessary files for GNPS FBMN and IIMN, but these tools must be executed externally. GNPS result files, including a table of library hits and the network graph, can then be used to annotate the FeatureMatrix with library hits and enrich the FBMN network graph with SIRIUS results.",
                )
            with cols[1]:
                st.image(str(Path("assets", "GNPS_logo.png")), width=200)
            st.divider()
            cols = st.columns([0.75, 0.25])
            with cols[0]:
                self.ui.input_widget(
                    "run-ms2query",
                    False,
                    "predict **chemical analogues and compound classes** with MS2Query",
                    help="Unlike traditional library search methods that focus exclusively on exact matches, MS2Query identifies related spectra of analogs with high chemical similarity. It employs machine learning-based chemical similarity predictors trained on existing spectral libraries. UmetaFlow automatically downloads the default models, which are trained on the GNPS library. The FeatureMatrix is annotated with MS2Query hits, including compound names, structures, and classes.",
                )
            with cols[1]:
                st.image(str(Path("assets", "ms2query_logo.png")), width=200)
            with st.columns(2)[0].container(border=True):
                st.image(str(Path("assets", "annotations.png")))

    def configure_expert(self) -> None:
        tabs = st.tabs(
            ["⚙️ **Pre-Processing**", "🔎 **Re-Quantification**", "🏷️ **Annotation**"]
        )
        with tabs[0]:
            t = st.tabs(
                [
                    "Precursor Mass Correction",
                    "Feature Detection",
                    "Adduct Detection",
                    "Map Alignment",
                    "Feature Linking",
                ]
            )
            with t[0]:
                self.ui.input_widget(
                    "correct-precursor",
                    True,
                    "correct precursor mass to highest intensity peak",
                    help="Correct precursor mass to highest intensity MS peak peak.",
                )
                self.ui.input_TOPP(
                    "HighResPrecursorMassCorrector",
                    include_parameters=[
                        "highest_intensity_peak:mz_tolerance",
                    ],
                    custom_defaults={"highest_intensity_peak:mz_tolerance": 100.0},
                )
            with t[1]:
                # Parameters for FeatureFinderMetabo TOPP tool.
                self.ui.input_TOPP(
                    "FeatureFinderMetabo",
                    exclude_parameters=["report_convex_hulls", "quant_method"],
                    custom_defaults={
                        "algorithm:common:noise_threshold_int": 1000.0,
                        "algorithm:mtd:mass_error_ppm": 10.0,
                        "algorithm:ffm:remove_single_traces": "true",
                        "algorithm:ffm:report_convex_hulls": "true",
                    },
                )
            with t[2]:
                # A single checkbox widget for workflow logic.
                self.ui.input_widget(
                    "adduct-detection",
                    False,
                    "enable **adduct detection**",
                    help="Detect feature adducts using the OpenMS TOPP tool *MetaboliteAdductDecharger*.",
                )
                # Paramters for MetaboliteAdductDecharger TOPP tool.
                self.ui.input_TOPP("MetaboliteAdductDecharger")
            with t[3]:
                self.ui.input_widget(
                    "map-alignement",
                    True,
                    "enable **map alignement**",
                    help="Align features to a reference map using the OpenMS TOPP tool *MapAlignerPoseClustering*.",
                )
                self.ui.input_TOPP(
                    "MapAlignerPoseClustering", exclude_parameters=["index"]
                )
            with t[4]:
                self.ui.input_TOPP(
                    "FeatureLinkerUnlabeledKD",
                )
        with tabs[1]:
            self.ui.input_widget(
                "requantify",
                False,
                "**re-quantify** features with missing values",
                help="Re-quantify missing values in consensus features using the OpenMS TOPP tool *FeatureFinderMetaboIdent*.",
            )
            self.ui.input_TOPP("FeatureFinderMetaboIdent")
        with tabs[2]:
            t = st.tabs(
                [
                    "In-house MS2 library",
                    "SIRIUS, CSI:FingerID & CANOPUS",
                    "GNPS FBMN & IIMN File Export",
                    "M2Query",
                ]
            )
            with t[0]:
                self.ui.input_widget(
                    "annotate-ms2",
                    False,
                    "annotate consensus features",
                    help="Based on MS2 spectrum similarity.",
                )
                self.ui.simple_file_uploader(
                    "ms2-library", "mgf", "MS2 library in mgf format"
                )
                self.ui.input_TOPP("MetaboliteSpectralMatcher")
            with t[1]:
                if "SiriusExport-path" not in st.session_state:
                    possible_paths = (
                        [  # potential SIRIUS locations in increasing priority
                            str(Path("SiriusExport")),  # anywhere
                            str(
                                Path(sys.prefix, "bin", "SiriusExport")
                            ),  # in current conda environment
                        ]
                    )
                    st.session_state["SiriusExport-path"] = ""
                    for path in possible_paths:
                        if shutil.which(path) is not None:
                            st.session_state["SiriusExport-path"] = path
                if st.session_state["SiriusExport-path"]:
                    st.markdown("**Pre-processing and file export**")
                    self.ui.input_widget(
                        "export-sirius",
                        False,
                        "export files for SIRIUS",
                        help="Generate input files for SIRIUS from raw data and feature information using the OpenMS TOPP tool *SiriusExport*.",
                    )
                    self.ui.input_TOPP("SiriusExport")
                self.add_sirius_path_to_session_state()

                if st.session_state["sirius-path"]:
                    st.markdown("**SIRIUS user login**")
                    cols = st.columns([0.25, 0.25, 0.5])
                    with cols[0]:
                        self.ui.input_widget(
                            "sirius-user-email",
                            "",
                            "Email",
                            help="Email address from a valid SIRIUS account.",
                            widget_type="text",
                        )
                    with cols[1]:
                        self.ui.input_widget(
                            "sirius-user-password",
                            "",
                            "password **NOT ENCRYPTED!**",
                            help="Password from a valid SIRIUS account. **Not encrypted**, will be stored in **plain text** in parameters and show up in log files.",
                            widget_type="password",
                        )
                    st.markdown("**SIRIUS**")
                    self.ui.input_widget(
                        "run-sirius",
                        False,
                        "predict **sum formulas**",
                        help="Generate input files for SIRIUS from raw data and feature information using the OpenMS TOPP tool *SiriusExport*.",
                    )
                    cols = st.columns(4)
                    with cols[0]:
                        self.ui.input_widget(
                            "sirius-profile",
                            name="profile",
                            options=["default", "qtof", "orbitrap", "fticr"],
                            default="default",
                            help="Name of the configuration profile",
                        )
                    with cols[1]:
                        self.ui.input_widget(
                            "sirius-maxmz",
                            300,
                            "max m/z",
                            min_value=50,
                            max_value=1000,
                            step_size=50,
                            help="Only considers compounds with a precursor m/z lower or equal. All other compounds in the input will be skipped. Recommended to be below 600, otherwise very long execution times are expected.",
                        )
                    with cols[2]:
                        self.ui.input_widget(
                            "sirius-db",
                            "none",
                            "database formula prediction",
                            options=[
                                "none",
                                "ALL",
                                "ALL_BUT_INSILICO",
                                "BIO",
                                "PUBCHEM",
                                "MESH",
                                "HMDB",
                                "KNAPSACK",
                                "CHEBI",
                                "PUBMED",
                                "KEGG",
                                "HSDB",
                                "MACONDA",
                                "METACYC",
                                "GNPS",
                                "ZINCBIO",
                                "UNDP",
                                "YMDB",
                                "PLANTCYC",
                                "NORMAN",
                                "ADDITIONAL",
                                "PUBCHEMANNOTATIONBIO",
                                "PUBCHEMANNOTATIONDRUG",
                                "PUBCHEMANNOTATIONSAFETYANDTOXIC",
                                "PUBCHEMANNOTATIONFOOD",
                                "KEGGMINE",
                                "ECOCYCMINE",
                                "YMDBMINE",
                            ],
                            help="Search formulas in the given database. If no database is given all possible molecular formulas will be respected (no database is used).",
                        )
                    cols = st.columns(4)
                    with cols[0]:
                        self.ui.input_widget(
                            "sirius-elements-considered",
                            "SBrClBSe",
                            "elements considered",
                            help="Set the allowed elements for rare element detection. Example: `SBrClBSe` to allow the elements S,Br,Cl,B and Se.",
                        )
                    with cols[1]:
                        self.ui.input_widget(
                            "sirius-elements-enforced",
                            "CHNOP",
                            "elements enforced",
                            help="Example: CHNOPSCl to allow the elements C, H, N, O, P, S and Cl. Add numbers in brackets to restrict the minimal and maximal allowed occurrence of these elements: CHNOP[5]S[8]Cl[1-2]. When one number is given then it is interpreted as upper bound. Default: C,H,N,O,P",
                        )
                    with cols[2]:
                        self.ui.input_widget(
                            "sirius-ions-considered",
                            "[M+H]+,[M+K]+,[M+Na]+,[M+H-H2O]+,[M+H-H4O2]+,[M+NH4]+,[M-H]-,[M+Cl]-,[M-H2O-H]-,[M+Br]-",
                            "ions considered",
                            help="The ion type/adduct of the MS/MS data in a comma separated list. Default: '[M+H]+,[M+K]+,[M+Na]+,[M+H-H2O]+,[M+H-H4O2]+,[M+NH4]+,[M-H]-,[M+Cl]-,[M-H2O-H]-,[M+Br]-'",
                        )
                    with cols[0]:
                        self.ui.input_widget(
                            "sirius-ppm-max",
                            10.0,
                            "ppm max",
                            help="Maximum allowed mass deviation in ppm for decomposing masses. Default: 10.0 ppm",
                        )
                    with cols[1]:
                        self.ui.input_widget(
                            "sirius-ppm-max-ms2",
                            10.0,
                            "ppm max MS2",
                            help="Maximum allowed mass deviation in ppm for decomposing masses in MS2. Default: 10.0 ppm",
                        )
                    self.ui.input_widget(
                        "run-fingerid",
                        False,
                        "predict **molecular structures**",
                        help="This subtool is dedicated to predicting molecular structures based on tandem mass spectrometry (MS/MS) data. It utilizes a fragmentation tree approach for the annotation of fragment spectra.",
                    )
                    cols = st.columns(4)
                    with cols[0]:
                        self.ui.input_widget(
                            "sirius-structure-db",
                            "BIO",
                            "structure database",
                            options=[
                                "ALL",
                                "ALL_BUT_INSILICO",
                                "BIO",
                                "PUBCHEM",
                                "MESH",
                                "HMDB",
                                "KNAPSACK",
                                "CHEBI",
                                "PUBMED",
                                "KEGG",
                                "HSDB",
                                "MACONDA",
                                "METACYC",
                                "GNPS",
                                "ZINCBIO",
                                "UNDP",
                                "YMDB",
                                "PLANTCYC",
                                "NORMAN",
                                "ADDITIONAL",
                                "PUBCHEMANNOTATIONBIO",
                                "PUBCHEMANNOTATIONDRUG",
                                "PUBCHEMANNOTATIONSAFETYANDTOXIC",
                                "PUBCHEMANNOTATIONFOOD",
                                "KEGGMINE",
                                "ECOCYCMINE",
                                "YMDBMINE",
                            ],
                            help="Search structure in the given database.",
                        )
                    self.ui.input_widget(
                        "run-canopus",
                        False,
                        "predict **compound classes**",
                        help="Predict compound categories for each compound individually based on its predicted molecular fingerprint (CSI:FingerID) using CANOPUS.",
                    )
                else:
                    st.info(
                        "💡 Install SIRIUS as command line tool to annotate features with formula, structural and compound class predictions."
                    )
            with t[2]:
                self.ui.input_widget(
                    "export-gnps",
                    True,
                    "export files for GNPS FBMN and IIMN",
                    help="Generate input files for GNPS feature based molecular networking (FBMN) and ion identity molecular networking (IIMN) from raw data and feature information using the OpenMS TOPP tool *GNPSExport*.",
                )
                self.ui.input_TOPP("GNPSExport")
            with t[3]:
                self.ui.input_widget(
                    "run-ms2query",
                    False,
                    "predict **chemical analogues and compound classes** with MS2Query",
                    help="Unlike traditional library search methods that focus exclusively on exact matches, MS2Query identifies related spectra of analogs with high chemical similarity. It employs machine learning-based chemical similarity predictors trained on existing spectral libraries. UmetaFlow automatically downloads the default models, which are trained on the GNPS library. The FeatureMatrix is annotated with MS2Query hits, including compound names, structures, and classes.",
                )
                self.ui.input_widget(
                    "ion_mode",
                    "positive",
                    name="ion mode for MS2Query",
                    options=["positive", "negative"],
                    help="Ion mode for MS2Query analog search.",
                )

    def configure(self) -> None:
        if not self.expert_mode:
            self.configure_simple()
        # Else, configure expert mode
        else:
            self.configure_expert()

    def format_simple_params(self) -> dict:
        """Takes the paramter file from simple configuration page and translate into one useable in this workflow."""
        with open(
            Path(self.workflow_dir, "..", "umetaflow", "params.json"), "r"
        ) as f:
            simple = json.load(f)

        expert_params_path = Path(
            self.workflow_dir, "..", "umetaflow-expert", "params.json"
        )
        if expert_params_path.exists():
            with open(expert_params_path, "r") as f:
                expert = json.load(f)
        else:
            with open("default-parameters.json", "r") as f:
                expert = json.load(f)["umetaflow-expert"]

        new = expert.copy()

        new["ion_mode"] = "positive"

        for key in ["correct-precursor", "map-alignement"]:
            new[key] = True

        for k, v in simple.items():
            if k in new.keys():
                new[k] = v
            if "ffm:" in k:
                new["FeatureFinderMetabo"][k[4:]] = v

        # mz and RT tolerances
        new["HighResPrecursorMassCorrector"][
            "highest_intensity_peak:mz_tolerance"
        ] = 100.0

        new["FeatureFinderMetabo"]["algorithm:mtd:mass_error_ppm"] = simple[
            "mz_tolerance"
        ]
        new["FeatureFinderMetabo"]["algorithm:ffm:isotope_filtering_model"] = "none"
        new["FeatureFinderMetabo"]["algorithm:ffm:report_convex_hulls"] = "true"

        new["MapAlignerPoseClustering"][
            "algorithm:pairfinder:distance_MZ:max_difference"
        ] = 10.0
        new["MapAlignerPoseClustering"]["algorithm:pairfinder:distance_MZ:unit"] = "ppm"

        new["FeatureLinkerUnlabeledKD"]["algorithm:warp:enabled"] = "false"
        new["FeatureLinkerUnlabeledKD"]["algorithm:link:rt_tol"] = simple[
            "RT_tolerance"
        ]
        new["FeatureLinkerUnlabeledKD"]["algorithm:link:mz_tol"] = simple[
            "mz_tolerance"
        ]

        new["sirius-ppm-max"] = simple["mz_tolerance"]
        new["sirius-ppm-max-ms2"] = simple["mz_tolerance"]

        new["MetaboliteSpectralMatcher"]["algorithm:merge_spectra"] = "false"

        new["FeatureFinderMetaboIdent"]["extract:mz_window"] = simple["mz_tolerance"]
        new["FeatureFinderMetaboIdent"]["extract:rt_window"] = simple["RT_tolerance"]

        # Adduct detection
        new["MetaboliteAdductDecharger"][
            "algorithm:MetaboliteFeatureDeconvolution:potential_adducts"
        ] = (
            simple["adducts_pos"].replace(" ", "\n")
            if simple["ion_mode"] == "positive"
            else simple["adducts_neg"].replace(" ", "\n")
        )
        if simple["ion_mode"] == "negative":
            new["MetaboliteAdductDecharger"][
                "algorithm:MetaboliteFeatureDeconvolution:negative_mode"
            ] = ""

        new["MetaboliteAdductDecharger"][
            "algorithm:MetaboliteFeatureDeconvolution:max_neutrals"
        ] = 1
        new["MetaboliteAdductDecharger"][
            "algorithm:MetaboliteFeatureDeconvolution:charge_span_max"
        ] = 3
        if simple["ion_mode"] == "positive":
            new["MetaboliteAdductDecharger"][
                "algorithm:MetaboliteFeatureDeconvolution:charge_max"
            ] = 3
        else:
            new["MetaboliteAdductDecharger"][
                "algorithm:MetaboliteFeatureDeconvolution:charge_min"
            ] = -3
            new["MetaboliteAdductDecharger"][
                "algorithm:MetaboliteFeatureDeconvolution:max_neutrals"
            ] = 1
            new["MetaboliteAdductDecharger"][
                "algorithm:MetaboliteFeatureDeconvolution:charge_max"
            ] = -1

        # SIRIUS logic
        if not simple["run-sirius"] and simple["run-fingerid"]:
            new["run-sirius"] = True

        if simple["run-fingerid"]:
            new["run-canopus"] = True

        # threads for each TOPP tool
        for k, v in new.items():
            if isinstance(v, dict):
                new[k]["threads"] = simple["num_threads"]

        self.parameter_manager.params_file = Path(
            Path(self.parameter_manager.params_file).parent, "params-translated.json"
        )
        with open(self.parameter_manager.params_file, "w") as f:
            json.dump(new, f, indent=4)
        return new

    def execution(self) -> None:
        # Check if run in expert mode, if not paramters need to be formatted to be compatible with this framework.
        if self.expert_mode:
            self.logger.log("Running UmetaFlow with expert mode configuration.")
        else:
            self.params = self.format_simple_params()

        # Get mzML files
        df_path = Path(self.workflow_dir, '..', "mzML-files.tsv")

        if not df_path.exists():
            mzML = []
        else:
            df = pd.read_csv(df_path, sep="\t")

            # Filter the DataFrame for files where "use in workflow" is True
            selected_files = df[df["use in workflows"] == True]["file name"].tolist()

            # Construct full file paths
            mzML = [
                str(Path(self.workflow_dir, '..', "mzML-files", file_name))
                for file_name in selected_files
            ]

        if len(mzML) == 0:
            self.logger.log(
                "ERROR: Select at leat two mzML files to run this workflow."
            )
            return

        # # Get mzML input files from self.params.
        # mzML = self.file_manager.get_files(self.params["mzML-files"])

        # # Log any messages.
        self.logger.log(f"Number of input mzML files: {len(mzML)}")
        self.logger.log(f"mzML files: {[Path(p).name for p in mzML]}")

        # Precursor m/z correction to highest intensity MS1 peak
        if self.params["correct-precursor"]:
            self.logger.log("Correcting precursor m/z to highest intensity MS1 peak.")
            mzML_pmc = self.file_manager.get_files(mzML, "mzML", "mzML-pmc")
            self.executor.run_topp(
                "HighResPrecursorMassCorrector",
                {"in": mzML, "out": mzML_pmc},
            )
            mzML = mzML_pmc

        # Feature Detection
        self.logger.log("Detecting features.")
        ffm = self.file_manager.get_files(mzML, "featureXML", "ffm-featureXML")
        self.executor.run_topp(
            "FeatureFinderMetabo",
            input_output={
                "in": mzML,
                "out": ffm,
                "out_chrom": self.file_manager.get_files(
                    mzML, set_results_dir="ffm-chroms"
                ),
            },
        )

        # Adduct Detection
        if self.params["adduct-detection"]:
            self.logger.log("Detecting adducts.")
            # Run MetaboliteAdductDecharger for adduct detection, with disabled logs.
            self.executor.run_topp(
                "MetaboliteAdductDecharger",
                {"in": ffm, "out_fm": ffm},
            )

        # Map Alignement
        if self.params["map-alignement"]:
            self.logger.log("Aligning feature maps.")
            trafos = self.file_manager.get_files(
                ffm, "trafoXML", "trafos", collect=True
            )
            # Run MapAlignerPoseClustering for map alignement, with disabled logs.
            self.executor.run_topp(
                "MapAlignerPoseClustering",
                {
                    "in": self.file_manager.get_files(ffm, collect=True),
                    "out": self.file_manager.get_files(ffm, collect=True),
                    "trafo_out": trafos,
                },
            )
            self.logger.log("Transforming mzML files based on map alignement.")
            # Transform mzML files
            self.executor.run_topp(
                "MapRTTransformer",
                {
                    "in": mzML,
                    "out": mzML,
                    "trafo_in": self.file_manager.get_files(trafos),
                },
            )

        # Export FFM feature maps to dataframes (including chromatograms)
        self.executor.run_python("export_ffm_df", {"in": ffm})

        # Feature Linking and Export to pd.DataFrame
        self.logger.log("Linking features.")
        consensusXML = self.file_manager.get_files(
            "feature-matrix-ffm", "consensusXML", "feature-linker"
        )
        self.executor.run_topp(
            "FeatureLinkerUnlabeledKD",
            {"in": self.file_manager.get_files(ffm, collect=True), "out": consensusXML},
        )

        # Export to DataFrame
        consensus_df = self.file_manager.get_files(
            "feature-matrix", "parquet", "consensus-dfs"
        )

        self.executor.run_python(
            "export_consensus_df", {"in": consensusXML, "out": consensus_df}
        )

        # Requantify features with missing values
        if self.params["requantify"]:
            self.logger.log("Re-quantifying features with missing values.")
            # Prepare library
            ffmid_library = self.file_manager.get_files(
                "library", "tsv", "ffmid-library"
            )
            consensus_df_ffm_complete = self.file_manager.get_files(
                "consensus-df-ffm-complete", "parquet", "consensus-dfs"
            )
            self.executor.run_python(
                "generate_FFMID_library",
                {
                    "in": consensus_df,
                    "out": ffmid_library,
                    "out_ffm": consensus_df_ffm_complete,
                },
            )

            # Run FeatureFinderMetaboIdent
            ffmid = self.file_manager.get_files(mzML, "featureXML", "ffmid-featureXML")
            self.executor.run_topp(
                "FeatureFinderMetaboIdent",
                {"in": mzML, "out": ffmid, "id": ffmid_library},
            )

            # Perform Adduct detection on re-quantified features
            if self.params["adduct-detection"]:
                self.logger.log("Detecting adducts for re-quantified features.")
                # Run MetaboliteAdductDecharger for adduct detection.
                self.executor.run_topp(
                    "MetaboliteAdductDecharger",
                    {"in": ffmid, "out_fm": ffmid},
                )

            # Export re-quantified feature maps to dataframes (including chromatograms)
            self.executor.run_python("export_ffmid_df", {"in": ffmid})

            # Link re-quantified features
            consensusXML_ffmid = self.file_manager.get_files(
                "feature-matrix-ffmid", "consensusXML", "feature-linker"
            )
            self.executor.run_topp(
                "FeatureLinkerUnlabeledKD",
                {
                    "in": self.file_manager.get_files(ffmid, collect=True),
                    "out": consensusXML_ffmid,
                },
            )

            # Export to DataFrame
            consensus_df_ffmid = self.file_manager.get_files(
                "feature-matrix-ffmid", "parquet", "consensus-dfs"
            )
            self.executor.run_python(
                "export_consensus_df",
                {"in": consensusXML_ffmid, "out": consensus_df_ffmid},
            )

            # Merge consensus_df and consensus_df_ffmid
            self.executor.run_python(
                "merge_consensus_df",
                {
                    "in": [consensus_df_ffm_complete, consensus_df_ffmid],
                    "out": consensus_df,
                },
            )

            # Merge feature maps from FFM and FFMID from merged consensus table
            self.executor.run_python(
                "merge_ffm_ffmid_df",
                {
                    "in": consensus_df,
                },
            )

            # Re-create feature maps from consensus df
            self.executor.run_python(
                "recreate_feature_maps",
                {"in": str(Path(self.file_manager.workflow_dir, "results"))},
            )

            # Ensure mzML and featureXML file paths are ordered the same for SiriusExport and GNPSExport
            ffm = sorted(
                [
                    str(p)
                    for p in Path(
                        Path(consensus_df[0]).parent.parent,
                        "feature-maps-recreated",
                    ).glob("*.featureXML")
                ]
            )
            mzML = sorted(mzML)

        self.add_sirius_path_to_session_state()
        if st.session_state["sirius-path"]:
            if (
                self.params["run-sirius"]
                or self.params["run-fingerid"]
                or self.params["run-canopus"]
            ):
                if not (
                    self.params["sirius-user-email"]
                    and self.params["sirius-user-password"]
                ):
                    self.logger.log(
                        "WARNING: SIRIUS account info incomplete. SIRIUS will not be executed and features not annotated."
                    )
                    st.session_state["sirius-path"] = ""
            else:
                st.session_state["sirius-path"] = ""

        if self.params["export-sirius"] or st.session_state["sirius-path"]:
            self.logger.log("Exporting input files for SIRIUS.")
            sirius_ms_files = self.file_manager.get_files(mzML, "ms", "sirius-export")
            self.executor.run_topp(
                "SiriusExport",
                {
                    "in": mzML,
                    "in_featureinfo": ffm,
                    "out": sirius_ms_files,
                },
            )
            if st.session_state["sirius-path"]:
                self.logger.log("Logging in to SIRIUS...")
                self.executor.run_command(
                    [
                        st.session_state["sirius-path"],
                        "login",
                        f"--email={self.params['sirius-user-email']}",
                        f"--password={self.params['sirius-user-password']}",
                    ]
                )
                sirius_projects = [
                    Path(
                        self.workflow_dir, "results", "sirius-projects", Path(file).stem
                    )
                    for file in sirius_ms_files
                ]
                commands = []
                for ms, project in zip(sirius_ms_files, sirius_projects):
                    if Path(ms).stat().st_size > 0:
                        project.mkdir(parents=True)
                        command = [
                            st.session_state["sirius-path"],
                            "--input",
                            ms,
                            "--project",
                            str(project),
                            "--no-compression",
                            "--maxmz",
                            self.params["sirius-maxmz"],
                            "formula",
                            "--db",
                            self.params["sirius-db"],
                            "--ions-considered",
                            self.params["sirius-ions-considered"],
                            "--elements-considered",
                            self.params["sirius-elements-considered"],
                            "--elements-enforced",
                            self.params["sirius-elements-enforced"],
                            "--ppm-max",
                            self.params["sirius-ppm-max"],
                            "--ppm-max-ms2",
                            self.params["sirius-ppm-max-ms2"],
                            "--profile",
                            self.params["sirius-profile"],
                            "--candidates",
                            "1",
                        ]
                        if self.params["run-fingerid"] or self.params["run-canopus"]:
                            command.append("fingerprint")
                        if self.params["run-fingerid"]:
                            command += [
                                "structure",
                                "--db",
                                self.params["sirius-structure-db"],
                            ]
                        if self.params["run-canopus"]:
                            command.append("canopus")
                        command.append("write-summaries")
                        commands.append(command)
                if commands:
                    self.logger.log("Running SIRIUS... (might take a VERY long time)")
                    for command in commands:
                        self.logger.log(f"Running file {Path(command[2]).stem}...")
                        self.executor.run_command(command)
                else:
                    self.logger.log("No MS2 data for SIRIUS to process.")

        if (
            self.params["export-gnps"]
            or self.params["annotate-ms2"]
            or self.params["run-ms2query"]
        ):
            self.logger.log("Exporting input files for GNPS.")
            # Map MS2 specs to features
            self.executor.run_topp(
                "IDMapper",
                {
                    "in": ffm,
                    "spectra:in": mzML,
                    "out": ffm,
                    "id": self.file_manager.get_files(
                        str(Path("assets", "empty.idXML"))
                    ),
                },
            )
            # Link features with MS2 info
            gnps_consensus = self.file_manager.get_files(
                "feature-matrix-gnps", "consensusXML", "feature-linker"
            )
            self.executor.run_topp(
                "FeatureLinkerUnlabeledKD",
                {
                    "in": self.file_manager.get_files(ffm, collect=True),
                    "out": gnps_consensus,
                },
            )

            # Filter consensus features which have missing values
            self.executor.run_topp(
                "FileFilter",
                {"in": gnps_consensus, "out": gnps_consensus},
                custom_params={"id:remove_unannotated_features": ""},
            )

            # Export to dataframe
            self.executor.run_python(
                "export_consensus_df",
                {
                    "in": gnps_consensus,
                    "out": self.file_manager.get_files(
                        "feature-matrix-gnps", "parquet", "consensus-dfs"
                    ),
                },
            )

            # Run GNPSExport
            self.executor.run_topp(
                "GNPSExport",
                {
                    "in_cm": gnps_consensus,
                    "in_mzml": self.file_manager.get_files(mzML, collect=True),
                    "out": self.file_manager.get_files("MS2", "mgf", "gnps-export"),
                    "out_quantification": self.file_manager.get_files(
                        "feature-quantification", "txt", "gnps-export"
                    ),
                    "out_pairs": self.file_manager.get_files(
                        "pairs", "csv", "gnps-export"
                    ),
                    "out_meta_values": self.file_manager.get_files(
                        "meta-values", "tsv", "gnps-export"
                    ),
                },
            )

        if self.params["annotate-ms2"]:
            dir_path = Path(self.workflow_dir, "input-files", "ms2-library")
            if dir_path.exists():
                dir_path = Path(self.workflow_dir, "input-files", "ms2-library")
                if dir_path.exists():
                    files = [p for p in dir_path.iterdir()]
                    if files:
                        self.logger.log("Annotating consensus features on MS2 level.")
                        self.executor.run_topp(
                            "FileConverter",
                            {
                                "in": self.file_manager.get_files(
                                    "MS2", "mgf", "gnps-export"
                                ),
                                "out": self.file_manager.get_files(
                                    "MS2", "mzML", "spectral-matcher"
                                ),
                            },
                        )
                        self.executor.run_topp(
                            "MetaboliteSpectralMatcher",
                            {
                                "in": self.file_manager.get_files(
                                    "MS2", "mzML", "spectral-matcher"
                                ),
                                "database": self.file_manager.get_files(str(files[0])),
                                "out": self.file_manager.get_files(
                                    "MS2-matches", "mzTab", "spectral-matcher"
                                ),
                            },
                            custom_params={"algorithm:merge_spectra": "false"},
                        )
                        self.executor.run_python(
                            "annotate-ms2",
                            {
                                "in_mzTab": self.file_manager.get_files(
                                    "MS2-matches", "mzTab", "spectral-matcher"
                                ),
                                "in_mzML": self.file_manager.get_files(
                                    "MS2", "mzML", "spectral-matcher"
                                ),
                                "in_mgf": self.file_manager.get_files(
                                    "MS2", "mgf", "gnps-export"
                                ),
                                "in_gnps_consensus": self.file_manager.get_files(
                                    "feature-matrix-gnps", "parquet", "consensus-dfs"
                                ),
                                "out": consensus_df,
                            },
                        )

        if st.session_state["sirius-path"]:
            self.executor.run_python("annotate-sirius", {"in": consensus_df})

        if self.params["run-ms2query"]:
            self.logger.log(
                "Detecting chemical analogues and compound classes with MS2Query. This will take a while..."
            )
            self.executor.run_python(
                "run_ms2query",
                {
                    "in": consensus_df,
                    "in_mgf": self.file_manager.get_files("MS2", "mgf", "gnps-export"),
                    "out_ms2query_csv": self.file_manager.get_files(
                        "MS2", "csv", "ms2query"
                    ),
                    "ion_mode": self.params["ion_mode"],
                },
            )

        # ZIP all relevant files for Download
        self.executor.run_python("zip-result-files", {"in": consensus_df})

    def results(self) -> None:
        # Set current results directory
        st.session_state.results_dir = Path(self.workflow_dir, "results")

        # Check if results exist at all
        if not st.session_state.results_dir.exists():
            st.info("No results yet.")
            return

        c1, c2 = st.columns(2)
        # Help
        with st.sidebar:
            help_section()
        # Downloads
        with c2:
            download_section(self.workflow_dir)

        # Select a metabolite from the final FeatureMatrix
        metabolite = metabolite_selection()

        if metabolite is None:
            return

        # Metrics
        with st.expander("📝 **Summary**", expanded=True):
            metabolite_metrics(metabolite)

        # Annotations
        spectralmatching = metabolite[
            [i for i in metabolite.index if i.startswith("SpectralMatch")]
        ].replace('', pd.NA).dropna()
        if not spectralmatching.empty:
            with st.expander(f"🎯 **MS2 Spectral Matching**", expanded=True):
                spectralmatching_summary(spectralmatching)
        sirius = metabolite[
            [
                i
                for i in metabolite.index
                if any(i.startswith(k) for k in ("SIRIUS_", "CSI", "CANOPUS"))
            ]
        ].dropna()
        if not sirius.empty:
            with st.expander(f"⭐ **SIRIUS, CSI:FingerID & CANOPUS**", expanded=True):
                sirius_summary(sirius)

        ms2query = metabolite[
            [i for i in metabolite.index if i.startswith("MS2Query")]
        ].dropna()
        if not ms2query.empty:
            with st.expander(f"🤖 **MS2Query**", expanded=True):
                ms2query_summary(ms2query)

        gnps = metabolite[
            [i for i in metabolite.index if i.startswith("GNPS_")]
        ].replace('', pd.NA).dropna()
        if not gnps.empty:
            with st.expander(f"🌐 **GNPS FBMN**", expanded=True):
                gnps_summary(gnps)

        # Chromatograms and Intensities
        chrom_data = get_chroms_for_each_sample(metabolite)
        chrom_fig = get_feature_chromatogram_plot(chrom_data)
        auc_fig = get_feature_intensity_plot(metabolite)

        c1, c2 = st.columns(2)
        with c1:
            with st.expander(f"**📈 Chromatograms**", expanded=True):
                show_fig(chrom_fig, f"chromatograms_{metabolite.name}")

        with c2:
            with st.expander(f"**📊 Intensities**", expanded=True):
                show_fig(auc_fig, f"AUC_{metabolite.name}")



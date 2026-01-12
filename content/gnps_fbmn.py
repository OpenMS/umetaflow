"""
GNPS Feature-Based Molecular Networking (FBMN) Integration Page

Provides 4 tabs:
1. Export - Shows exported GNPS files and link to GNPS upload
2. Import - Import GraphML results (manual upload or Task-ID)
3. Library Matches - View library matches with structure display
4. Network - Interactive network visualization
"""

import streamlit as st
import pandas as pd
from pathlib import Path

from src.common.common import page_setup, show_table, show_fig
from src.gnps import (
    GNPS_FBMN_URL,
    check_gnps_status,
    download_gnps_graphml,
    parse_gnps_graphml,
    get_library_matches,
    map_gnps_to_feature_matrix,
    plot_gnps_network,
    get_network_statistics,
)


params = page_setup()

st.title("GNPS Feature-Based Molecular Networking")

# Check if workflow results exist
workflow_dir = Path(st.session_state.workspace, "umetaflow")
gnps_export_dir = Path(workflow_dir, "results", "gnps-export")
results_dir = Path(workflow_dir, "results")

# Session state for GNPS data
if "gnps_graph" not in st.session_state:
    st.session_state.gnps_graph = None
if "gnps_annotations" not in st.session_state:
    st.session_state.gnps_annotations = None

tabs = st.tabs(["Export", "Import Results", "Library Matches", "Network Visualization"])

# =============================================================================
# Tab 1: Export
# =============================================================================
with tabs[0]:
    st.header("GNPS Export Files")

    st.markdown("""
    UmetaFlow automatically generates the required files for GNPS Feature-Based
    Molecular Networking (FBMN) when the workflow is run with the GNPS export option enabled.
    """)

    if gnps_export_dir.exists():
        # List exported files
        exported_files = list(gnps_export_dir.glob("*"))

        if exported_files:
            st.success(f"Found {len(exported_files)} GNPS export files")

            file_info = []
            for f in exported_files:
                size_kb = f.stat().st_size / 1024
                file_info.append({
                    "File": f.name,
                    "Size": f"{size_kb:.1f} KB" if size_kb < 1024 else f"{size_kb/1024:.1f} MB",
                    "Type": f.suffix.upper()[1:] if f.suffix else "Unknown"
                })

            st.dataframe(pd.DataFrame(file_info), use_container_width=True, hide_index=True)

            # Download buttons
            st.subheader("Download Files")
            cols = st.columns(min(len(exported_files), 4))
            for i, f in enumerate(exported_files):
                with cols[i % 4]:
                    with open(f, "rb") as file:
                        st.download_button(
                            label=f.name,
                            data=file.read(),
                            file_name=f.name,
                            key=f"download_{f.name}"
                        )

            st.divider()

            # GNPS Upload Instructions
            st.subheader("Upload to GNPS")
            st.markdown(f"""
            **Required files for GNPS FBMN:**
            1. **MS2.mgf** - MS/MS spectra file
            2. **feature-quantification.txt** - Feature quantification table

            **Optional files:**
            - **pairs.csv** - Pre-computed spectral pairs (for IIMN)
            - **meta-values.tsv** - Additional metadata

            **Steps to submit:**
            1. Go to [GNPS FBMN Workflow]({GNPS_FBMN_URL})
            2. Upload the MS2.mgf file as "Spectrum Files"
            3. Upload the feature-quantification.txt as "Quantification Table"
            4. Configure parameters and submit
            5. Once complete, download the GraphML file and import it in the next tab
            """)

            st.link_button("Open GNPS FBMN", GNPS_FBMN_URL, type="primary")

        else:
            st.warning("GNPS export directory exists but is empty.")
    else:
        st.info("""
        No GNPS export files found.

        To generate GNPS files:
        1. Go to **Configure** page
        2. Enable "export files for GNPS FBMN and IIMN" in the Annotation section
        3. Run the workflow
        """)

# =============================================================================
# Tab 2: Import Results
# =============================================================================
with tabs[1]:
    st.header("Import GNPS Results")

    st.markdown("""
    Import your GNPS FBMN results to annotate the feature matrix with library matches.
    You can either:
    - **Enter a Task ID** to automatically download results (requires completed GNPS job)
    - **Upload a GraphML file** manually
    """)

    import_method = st.radio(
        "Import method",
        ["Manual Upload", "GNPS Task ID"],
        horizontal=True
    )

    if import_method == "GNPS Task ID":
        st.subheader("Automatic Import via Task ID")

        task_id = st.text_input(
            "GNPS Task ID",
            placeholder="e.g., c9219a98c83f4f14b54f6f7e5e5f1234",
            help="Enter the Task ID from your GNPS FBMN job"
        )

        if task_id:
            col1, col2 = st.columns([1, 3])

            with col1:
                check_status = st.button("Check Status", type="primary")

            if check_status:
                with st.spinner("Checking GNPS job status..."):
                    status = check_gnps_status(task_id)

                if status.get("status") == "DONE":
                    st.success("Job completed! Ready to download results.")

                    if st.button("Download GraphML"):
                        with st.spinner("Downloading GraphML..."):
                            gnps_results_dir = Path(st.session_state.workspace, "gnps-results")
                            graphml_path = download_gnps_graphml(task_id, gnps_results_dir)

                            if graphml_path and graphml_path.exists():
                                st.success(f"Downloaded: {graphml_path.name}")

                                # Parse and store in session state
                                annotations_df, G = parse_gnps_graphml(graphml_path)
                                st.session_state.gnps_graph = G
                                st.session_state.gnps_annotations = annotations_df
                                st.session_state.gnps_graphml_path = graphml_path

                                stats = get_network_statistics(G)
                                st.info(f"""
                                **Network Statistics:**
                                - Total nodes: {stats['total_nodes']}
                                - Total edges: {stats['total_edges']}
                                - Components: {stats['num_components']}
                                - Library matches: {stats['library_matches']}
                                """)
                            else:
                                st.error("Failed to download GraphML file")

                elif status.get("status") == "RUNNING":
                    st.warning("Job is still running. Please wait for completion.")
                elif status.get("status") == "ERROR":
                    st.error(f"Error checking status: {status.get('error', 'Unknown error')}")
                else:
                    st.info(f"Job status: {status.get('status', 'Unknown')}")

    else:  # Manual Upload
        st.subheader("Manual GraphML Upload")

        uploaded_file = st.file_uploader(
            "Upload GraphML file from GNPS",
            type=["graphml"],
            help="Upload the molecular network GraphML file from your GNPS FBMN results"
        )

        if uploaded_file is not None:
            # Save uploaded file
            gnps_results_dir = Path(st.session_state.workspace, "gnps-results")
            gnps_results_dir.mkdir(parents=True, exist_ok=True)
            graphml_path = gnps_results_dir / uploaded_file.name

            with open(graphml_path, "wb") as f:
                f.write(uploaded_file.getbuffer())

            st.success(f"Uploaded: {uploaded_file.name}")

            # Parse GraphML
            with st.spinner("Parsing GraphML..."):
                try:
                    annotations_df, G = parse_gnps_graphml(graphml_path)
                    st.session_state.gnps_graph = G
                    st.session_state.gnps_annotations = annotations_df
                    st.session_state.gnps_graphml_path = graphml_path

                    stats = get_network_statistics(G)
                    st.success("GraphML parsed successfully!")

                    col1, col2, col3, col4 = st.columns(4)
                    col1.metric("Nodes", stats['total_nodes'])
                    col2.metric("Edges", stats['total_edges'])
                    col3.metric("Components", stats['num_components'])
                    col4.metric("Library Matches", stats['library_matches'])

                except Exception as e:
                    st.error(f"Error parsing GraphML: {str(e)}")

    # Annotation section (appears after successful import)
    st.divider()
    st.subheader("Annotate Feature Matrix")

    if st.session_state.gnps_annotations is not None:
        feature_matrix_path = Path(results_dir, "consensus-dfs", "feature-matrix.parquet")
        gnps_matrix_path = Path(results_dir, "consensus-dfs", "feature-matrix-gnps.parquet")

        if feature_matrix_path.exists():
            st.success("Feature matrix found. Ready to annotate.")

            if st.button("Annotate Feature Matrix with GNPS Results", type="primary"):
                with st.spinner("Annotating feature matrix..."):
                    try:
                        feature_matrix = pd.read_parquet(feature_matrix_path)
                        gnps_matrix = pd.read_parquet(gnps_matrix_path) if gnps_matrix_path.exists() else pd.DataFrame()

                        annotated_df = map_gnps_to_feature_matrix(
                            st.session_state.gnps_annotations,
                            feature_matrix,
                            gnps_matrix
                        )

                        # Count annotations
                        gnps_count = (annotated_df['GNPS_Compound_Name'] != '').sum()

                        # Save annotated matrix
                        annotated_df.to_parquet(feature_matrix_path)
                        annotated_df.to_csv(str(feature_matrix_path).replace('.parquet', '.tsv'), sep='\t')

                        st.success(f"Successfully annotated {gnps_count} features with GNPS library matches!")
                        st.info("View updated results in the **Results** page.")

                    except Exception as e:
                        st.error(f"Error during annotation: {str(e)}")
        else:
            st.warning("Feature matrix not found. Please run the UmetaFlow workflow first.")
    else:
        st.info("Import a GraphML file first to enable annotation.")

# =============================================================================
# Tab 3: Library Matches
# =============================================================================
with tabs[2]:
    st.header("GNPS Library Matches")

    if st.session_state.gnps_annotations is not None:
        library_matches = get_library_matches(st.session_state.gnps_annotations)

        if len(library_matches) > 0:
            st.success(f"Found {len(library_matches)} library matches")

            # Filters
            with st.expander("Filters", expanded=True):
                col1, col2 = st.columns(2)

                with col1:
                    # Filter by compound name
                    search_term = st.text_input("Search compound name", "")
                    if search_term:
                        library_matches = library_matches[
                            library_matches['Compound_Name'].str.contains(search_term, case=False, na=False)
                        ]

                with col2:
                    # Filter by MQScore
                    if 'MQScore' in library_matches.columns:
                        try:
                            library_matches['MQScore_float'] = pd.to_numeric(library_matches['MQScore'], errors='coerce')
                            min_score = st.slider(
                                "Minimum MQScore",
                                0.0, 1.0, 0.0,
                                help="Filter by minimum match quality score"
                            )
                            library_matches = library_matches[
                                (library_matches['MQScore_float'] >= min_score) |
                                (library_matches['MQScore_float'].isna())
                            ]
                        except:
                            pass

            # Display columns
            display_cols = [
                'cluster_index', 'Compound_Name', 'SMILES', 'Adduct',
                'MQScore', 'MZErrorPPM', 'precursor_mass', 'RTConsensus'
            ]
            display_cols = [c for c in display_cols if c in library_matches.columns]

            show_table(library_matches[display_cols], "gnps-library-matches")

            # Structure viewer
            st.subheader("Structure Viewer")

            compounds_with_smiles = library_matches[
                (library_matches['SMILES'].notna()) &
                (library_matches['SMILES'] != '')
            ]

            if len(compounds_with_smiles) > 0:
                selected_compound = st.selectbox(
                    "Select compound to view structure",
                    compounds_with_smiles['Compound_Name'].unique()
                )

                if selected_compound:
                    smiles = compounds_with_smiles[
                        compounds_with_smiles['Compound_Name'] == selected_compound
                    ]['SMILES'].iloc[0]

                    try:
                        from rdkit import Chem
                        from rdkit.Chem import Draw
                        import io

                        mol = Chem.MolFromSmiles(smiles)
                        if mol:
                            img = Draw.MolToImage(mol, size=(400, 300))
                            st.image(img, caption=f"{selected_compound}\nSMILES: {smiles}")
                        else:
                            st.warning(f"Could not parse SMILES: {smiles}")
                    except ImportError:
                        st.info(f"SMILES: {smiles}")
                        st.caption("Install RDKit to view molecular structures")
            else:
                st.info("No compounds with SMILES found.")
        else:
            st.info("No library matches found in the imported GraphML.")
    else:
        st.info("Import a GraphML file in the **Import Results** tab to view library matches.")

# =============================================================================
# Tab 4: Network Visualization
# =============================================================================
with tabs[3]:
    st.header("Molecular Network Visualization")

    if st.session_state.gnps_graph is not None:
        G = st.session_state.gnps_graph
        stats = get_network_statistics(G)

        # Network statistics
        col1, col2, col3, col4, col5 = st.columns(5)
        col1.metric("Nodes", stats['total_nodes'])
        col2.metric("Edges", stats['total_edges'])
        col3.metric("Components", stats['num_components'])
        col4.metric("Library Matches", stats['library_matches'])
        col5.metric("Largest Component", stats['largest_component_size'])

        # Component selector
        st.subheader("View Network")

        # Get unique components
        components = set()
        for node, attrs in G.nodes(data=True):
            comp = attrs.get('ComponentIndex', attrs.get('component', ''))
            if comp:
                try:
                    components.add(int(float(comp)))
                except (ValueError, TypeError):
                    pass

        components = sorted(list(components))

        view_option = st.radio(
            "View option",
            ["All components", "Single component"],
            horizontal=True
        )

        selected_component = None
        if view_option == "Single component" and components:
            selected_component = st.selectbox(
                "Select component",
                components,
                help="Components are clusters of related features in the molecular network"
            )

        # Check network size
        max_nodes = 500
        if stats['total_nodes'] > max_nodes and view_option == "All components":
            st.warning(f"""
            Network has {stats['total_nodes']} nodes.
            For better performance, consider viewing individual components.
            Displaying first {max_nodes} nodes only.
            """)

        # Generate network plot
        with st.spinner("Generating network visualization..."):
            try:
                fig = plot_gnps_network(G, selected_component)
                show_fig(fig, f"gnps-network-component-{selected_component or 'all'}")
            except Exception as e:
                st.error(f"Error generating network: {str(e)}")

        # Node details on selection
        st.subheader("Node Details")
        st.info("Hover over nodes in the network to see details. Click to select a node.")

    else:
        st.info("Import a GraphML file in the **Import Results** tab to view the molecular network.")

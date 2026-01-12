"""
GNPS FBMN Integration Module

Helper functions for:
- GNPS API status checking and result download
- GraphML parsing and annotation extraction
- Network visualization
"""

import requests
import pandas as pd
import networkx as nx
from pathlib import Path
from typing import Optional, Tuple
import streamlit as st


GNPS_BASE_URL = "https://gnps.ucsd.edu/ProteoSAFe"
GNPS_FBMN_URL = "https://gnps.ucsd.edu/ProteoSAFe/index.jsp?params=%7B%22workflow%22:%22FEATURE-BASED-MOLECULAR-NETWORKING%22%7D"


def check_gnps_status(task_id: str) -> dict:
    """
    Check GNPS job status via API.

    Args:
        task_id: GNPS task identifier

    Returns:
        dict with status information including 'status' key
    """
    url = f"{GNPS_BASE_URL}/status_json.jsp?task={task_id}"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        return response.json()
    except requests.RequestException as e:
        return {"status": "ERROR", "error": str(e)}


def get_gnps_result_files(task_id: str) -> dict:
    """
    Get available result file URLs for a GNPS task.

    Args:
        task_id: GNPS task identifier

    Returns:
        dict mapping file type to download URL
    """
    base = f"{GNPS_BASE_URL}/DownloadResultFile"
    return {
        "graphml": f"{base}?task={task_id}&file=gnps_molecular_network_graphml/",
        "clustersummary": f"{base}?task={task_id}&file=clusterinfosummarygroup_attributes_withIDs_withcomponentID/",
        "libraryhits": f"{base}?task={task_id}&file=result_specnets_DB/",
    }


def download_gnps_graphml(task_id: str, output_dir: Path) -> Optional[Path]:
    """
    Download GraphML network file from GNPS.

    Args:
        task_id: GNPS task identifier
        output_dir: Directory to save the file

    Returns:
        Path to downloaded file or None if failed
    """
    urls = get_gnps_result_files(task_id)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        response = requests.get(urls["graphml"], timeout=60)
        response.raise_for_status()

        output_path = output_dir / f"{task_id}.graphml"
        with open(output_path, "wb") as f:
            f.write(response.content)
        return output_path
    except requests.RequestException:
        return None


def parse_gnps_graphml(graphml_path: Path) -> Tuple[pd.DataFrame, nx.Graph]:
    """
    Parse GNPS FBMN GraphML file and extract annotations.

    Args:
        graphml_path: Path to GraphML file

    Returns:
        Tuple of (annotations DataFrame, NetworkX Graph)
    """
    G = nx.read_graphml(graphml_path)

    # Extract node attributes for annotation
    annotations = []
    for node_id, attrs in G.nodes(data=True):
        # Get cluster index (scan number for mapping)
        cluster_index = attrs.get('cluster index', attrs.get('shared name', node_id))
        try:
            cluster_index = int(float(cluster_index))
        except (ValueError, TypeError):
            cluster_index = node_id

        annotation = {
            'cluster_index': cluster_index,
            'Compound_Name': attrs.get('Compound_Name', attrs.get('LibraryID', '')),
            'SMILES': attrs.get('SMILES', attrs.get('Smiles', '')),
            'INCHI': attrs.get('INCHI', attrs.get('InChI', attrs.get('inchi', ''))),
            'Adduct': attrs.get('Adduct', attrs.get('Ion_Source', '')),
            'ComponentIndex': attrs.get('ComponentIndex', attrs.get('component', '')),
            'MQScore': attrs.get('MQScore', attrs.get('MassDiff', '')),
            'MZErrorPPM': attrs.get('MZErrorPPM', ''),
            'SharedPeaks': attrs.get('SharedPeaks', ''),
            'LibrarySpectrumID': attrs.get('LibrarySpectrumID', attrs.get('SpectrumID', '')),
            'precursor_mass': attrs.get('parent mass', attrs.get('precursor mass', '')),
            'RTConsensus': attrs.get('RTConsensus', attrs.get('RTMean', '')),
        }
        annotations.append(annotation)

    df = pd.DataFrame(annotations)
    return df, G


def get_library_matches(annotations_df: pd.DataFrame) -> pd.DataFrame:
    """
    Filter annotations to only include library matches.

    Args:
        annotations_df: DataFrame from parse_gnps_graphml

    Returns:
        DataFrame with only features that have library matches
    """
    # Filter for rows with actual compound names (not empty or N/A)
    matches = annotations_df[
        (annotations_df['Compound_Name'].notna()) &
        (annotations_df['Compound_Name'] != '') &
        (annotations_df['Compound_Name'] != 'N/A')
    ].copy()
    return matches


def map_gnps_to_feature_matrix(
    gnps_annotations: pd.DataFrame,
    feature_matrix: pd.DataFrame,
    gnps_feature_matrix: pd.DataFrame
) -> pd.DataFrame:
    """
    Map GNPS annotations to the consensus feature matrix.

    The mapping works via:
    1. cluster_index (from GraphML) matches scan number in GNPS MGF export
    2. gnps_feature_matrix contains consensus_feature_id to scan mapping
    3. feature_matrix is indexed by metabolite name

    Args:
        gnps_annotations: DataFrame from parse_gnps_graphml
        feature_matrix: Main consensus feature matrix
        gnps_feature_matrix: GNPS-specific feature matrix with scan info

    Returns:
        Updated feature matrix with GNPS annotation columns
    """
    # Initialize GNPS columns
    gnps_columns = [
        'GNPS_Compound_Name', 'GNPS_SMILES', 'GNPS_INCHI', 'GNPS_Adduct',
        'GNPS_ComponentIndex', 'GNPS_MQScore', 'GNPS_MZErrorPPM',
        'GNPS_LibrarySpectrumID'
    ]
    for col in gnps_columns:
        feature_matrix[col] = ''

    # Create mapping from cluster_index to annotations
    gnps_annotations = gnps_annotations.set_index('cluster_index')

    # Map via consensus_feature_id
    if 'consensus_feature_id' in gnps_feature_matrix.columns:
        for idx, row in gnps_feature_matrix.iterrows():
            consensus_id = row.get('consensus_feature_id')
            if consensus_id is not None and consensus_id in gnps_annotations.index:
                ann = gnps_annotations.loc[consensus_id]
                if idx in feature_matrix.index:
                    feature_matrix.loc[idx, 'GNPS_Compound_Name'] = ann.get('Compound_Name', '')
                    feature_matrix.loc[idx, 'GNPS_SMILES'] = ann.get('SMILES', '')
                    feature_matrix.loc[idx, 'GNPS_INCHI'] = ann.get('INCHI', '')
                    feature_matrix.loc[idx, 'GNPS_Adduct'] = ann.get('Adduct', '')
                    feature_matrix.loc[idx, 'GNPS_ComponentIndex'] = str(ann.get('ComponentIndex', ''))
                    feature_matrix.loc[idx, 'GNPS_MQScore'] = str(ann.get('MQScore', ''))
                    feature_matrix.loc[idx, 'GNPS_MZErrorPPM'] = str(ann.get('MZErrorPPM', ''))
                    feature_matrix.loc[idx, 'GNPS_LibrarySpectrumID'] = ann.get('LibrarySpectrumID', '')

    return feature_matrix


def plot_gnps_network(G: nx.Graph, selected_component: Optional[int] = None) -> "plotly.graph_objects.Figure":
    """
    Create interactive network visualization using Plotly.

    Args:
        G: NetworkX Graph from GraphML
        selected_component: Optional component index to highlight

    Returns:
        Plotly Figure object
    """
    import plotly.graph_objects as go

    # Filter to selected component if specified
    if selected_component is not None:
        nodes_in_component = [
            n for n, d in G.nodes(data=True)
            if str(d.get('ComponentIndex', '')) == str(selected_component)
        ]
        G = G.subgraph(nodes_in_component)

    if len(G.nodes()) == 0:
        return go.Figure().add_annotation(text="No nodes to display", showarrow=False)

    # Calculate layout
    pos = nx.spring_layout(G, k=1/len(G.nodes())**0.5 if len(G.nodes()) > 1 else 1, iterations=50)

    # Create edges trace
    edge_x, edge_y = [], []
    for edge in G.edges():
        if edge[0] in pos and edge[1] in pos:
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines'
    )

    # Create nodes trace
    node_x, node_y, node_text, node_color = [], [], [], []
    for node in G.nodes():
        if node in pos:
            x, y = pos[node]
            node_x.append(x)
            node_y.append(y)

            attrs = G.nodes[node]
            compound = attrs.get('Compound_Name', 'Unknown')
            mz = attrs.get('parent mass', attrs.get('precursor mass', 'N/A'))
            rt = attrs.get('RTConsensus', 'N/A')
            component = attrs.get('ComponentIndex', 'N/A')

            text = f"ID: {node}<br>Compound: {compound}<br>m/z: {mz}<br>RT: {rt}<br>Component: {component}"
            node_text.append(text)

            # Color by component index
            try:
                node_color.append(int(float(attrs.get('ComponentIndex', 0))))
            except (ValueError, TypeError):
                node_color.append(0)

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        text=node_text,
        marker=dict(
            showscale=True,
            colorscale='Viridis',
            color=node_color,
            size=10,
            colorbar=dict(
                thickness=15,
                title='Component',
                xanchor='left',
                titleside='right'
            ),
            line_width=2
        )
    )

    fig = go.Figure(
        data=[edge_trace, node_trace],
        layout=go.Layout(
            title='GNPS Molecular Network',
            showlegend=False,
            hovermode='closest',
            margin=dict(b=20, l=5, r=5, t=40),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor='white'
        )
    )

    return fig


def get_network_statistics(G: nx.Graph) -> dict:
    """
    Calculate basic network statistics.

    Args:
        G: NetworkX Graph

    Returns:
        dict with network statistics
    """
    # Get component info
    components = {}
    for node, attrs in G.nodes(data=True):
        comp = attrs.get('ComponentIndex', 'Unknown')
        if comp not in components:
            components[comp] = 0
        components[comp] += 1

    # Count library matches
    library_matches = sum(
        1 for _, attrs in G.nodes(data=True)
        if attrs.get('Compound_Name') and attrs.get('Compound_Name') not in ['', 'N/A']
    )

    return {
        'total_nodes': G.number_of_nodes(),
        'total_edges': G.number_of_edges(),
        'num_components': len(components),
        'library_matches': library_matches,
        'largest_component_size': max(components.values()) if components else 0,
    }

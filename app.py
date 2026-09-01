import streamlit as st
from pathlib import Path

if __name__ == '__main__':
    pages = {
        "" : [
            st.Page(Path("content", "quickstart.py"), title="Quickstart", icon="👋"),
            st.Page(Path("content", "file_upload.py"), title="File Upload", icon="📂"),
            st.Page(Path("content", "raw_data_viewer.py"), title="View MS data", icon="👀"),
        ],
        "UmetaFlow: Untargeted Quantification & Identification" : [
            st.Page(Path("content", "umetaflow_configure.py"), title="Configure", icon="⚙️"),
            st.Page(Path("content", "umetaflow_run.py"), title="Run", icon="🚀"),
            st.Page(Path("content", "umetaflow_results.py"), title="Results", icon="📊"),
        ],
        "Targeted Quantification: Extracted Ion Chromatograms": [
            st.Page(Path("content", "mz_calculator.py"), title="m/z Calculator", icon="📟"),
            st.Page(Path("content", "extracted_ion_chromatograms.py"), title="Extracted Ion Chromatograms", icon="🔍"),
        ],
        "Downstream": [
            st.Page(Path("content", "statistics.py"), title="Statistics", icon="📈"),
        ],
        "Downstream (OpenMS-Insight)": [
            st.Page(Path("content", "downstream_filtering.py"), title="Filtering", icon="🧹"),
            st.Page(Path("content", "downstream_imputation.py"), title="Imputation", icon="🧩"),
            st.Page(Path("content", "downstream_normalization.py"), title="Normalization", icon="⚖️"),
            st.Page(Path("content", "downstream_statistics.py"), title="Statistical Inference", icon="🔬"),
            st.Page(Path("content", "downstream_results.py"), title="Volcano / PCA / Heatmap", icon="📊"),
            st.Page(Path("content", "downstream_kegg_pathway.py"), title="KEGG Pathway", icon="🧬"),
        ]
    }

    pg = st.navigation(pages)
    pg.run()
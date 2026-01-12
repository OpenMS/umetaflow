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
            st.Page(Path("content", "gnps_fbmn.py"), title="GNPS FBMN", icon="🌐"),
        ]
    }

    pg = st.navigation(pages)
    pg.run()
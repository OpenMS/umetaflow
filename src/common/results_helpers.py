"""Adapter layer that exposes umetaflow's feature-matrix as openms_insight-compatible tables.

Mirrors the role of quantms-web's `src/common/results_helpers.py`, but reads the
metabolomics feature-matrix (rows = features, columns = mzML samples) instead of a
proteomics quant_results CSV.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st

ID_COLUMN = "metabolite"
GROUPS_FILE = "sample-groups.json"


def get_umetaflow_dir(workspace: Path) -> Path:
    return Path(workspace, "umetaflow")


def get_feature_matrix_path(workspace: Path) -> Path | None:
    path = Path(
        get_umetaflow_dir(workspace), "results", "consensus-dfs", "feature-matrix.parquet"
    )
    return path if path.exists() else None


@st.cache_data
def _load_abundance_data(workspace_path: str, mtime: float) -> tuple:
    path = Path(
        workspace_path, "umetaflow", "results", "consensus-dfs", "feature-matrix.parquet"
    )
    df = pd.read_parquet(path)

    sample_cols = [c for c in df.columns if c.endswith(".mzML")]

    pivot_df = df[sample_cols].copy()
    pivot_df.insert(0, ID_COLUMN, df.index)
    pivot_df = pivot_df.reset_index(drop=True)

    expr_df = pivot_df.set_index(ID_COLUMN)[sample_cols].replace(0, np.nan)
    expr_df = np.log2(expr_df + 1).dropna()

    return pivot_df, expr_df, sample_cols


def get_abundance_data(workspace: Path) -> tuple | None:
    """Load (pivot_df, expr_df, sample_cols) from the workspace's feature-matrix.parquet.

    Returns None if no untargeted UmetaFlow run has produced a feature-matrix yet.
    """
    path = get_feature_matrix_path(workspace)
    if path is None:
        return None
    return _load_abundance_data(str(workspace), path.stat().st_mtime)


def get_id_column(*_args, **_kwargs) -> str:
    """Row identifier column. Kept as a function (not a constant) to mirror quantms-web's API."""
    return ID_COLUMN


def get_sample_group_map(workspace: Path, sample_cols: list[str]) -> dict[str, str]:
    """Load {sample_name: group_name} assignments saved via the group-assignment widget.

    Samples with no saved assignment map to "" (unassigned).
    """
    path = Path(get_umetaflow_dir(workspace), GROUPS_FILE)
    saved: dict[str, str] = {}
    if path.exists():
        with open(path, "r", encoding="utf-8") as f:
            saved = json.load(f)
    return {s: saved.get(s, "") for s in sample_cols}


def save_sample_group_map(workspace: Path, group_map: dict[str, str]) -> None:
    path = Path(get_umetaflow_dir(workspace), GROUPS_FILE)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(group_map, f, indent=2, ensure_ascii=False)


def render_group_assignment(workspace: Path, sample_cols: list[str]) -> dict[str, str]:
    """Editable table for assigning each sample to a biological group.

    Renders a data editor pre-filled with any previously saved assignments, lets the
    user edit and save them, and returns the (possibly just-saved) group map.
    """
    st.subheader("Sample Groups")
    st.caption(
        "Assign each sample to a biological group. This is required for filtering, "
        "imputation, normalization and statistical testing below."
    )

    current = get_sample_group_map(workspace, sample_cols)
    edit_df = pd.DataFrame(
        {"sample": list(current.keys()), "group": list(current.values())}
    )
    edited = st.data_editor(
        edit_df,
        column_config={
            "sample": st.column_config.TextColumn("Sample", disabled=True),
            "group": st.column_config.TextColumn("Group"),
        },
        hide_index=True,
        use_container_width=True,
        key="sample_group_editor",
    )

    # Reflect the table as currently edited (even if not yet saved to disk) so that
    # filtering/imputation/normalization/statistics further down THIS SAME page run
    # use what's actually in the table, not a stale on-disk copy.
    live_map = dict(zip(edited["sample"], edited["group"].fillna("")))

    just_saved = False
    if st.button("Save group assignment"):
        save_sample_group_map(workspace, live_map)
        just_saved = True
        st.success("Saved. This assignment will now be used on every Downstream page.")

    if not just_saved and live_map != current:
        st.caption(
            "⚠️ Unsaved changes — used below on this page, but click "
            "**Save group assignment** so other pages (Imputation, Normalization, "
            "Statistics...) see them too."
        )

    return live_map

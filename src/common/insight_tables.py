"""openms-insight helpers: streamed tables and content-keyed plot caches.

``st.dataframe`` ships the whole table to the browser, which stalls on
full-scale feature matrices. An openms-insight ``Table`` caches the data on
disk once and sends one page at a time; its download button exports every row
that passes the table's filters. Same approach as quantms-web's helper.
"""
import hashlib
import shutil
import tempfile
from pathlib import Path

import pandas as pd
import polars as pl
import streamlit as st

ROW_FIELD = "row"
# Rows sampled when fingerprinting a dataframe
_FINGERPRINT_SAMPLE_ROWS = 2000


def fingerprint_dataframe(df: pd.DataFrame) -> str:
    """Cheap content hash of a dataframe (shape, schema, column sums, sampled rows)."""
    digest = hashlib.sha256()
    digest.update(repr((df.shape, [str(c) for c in df.columns], [str(t) for t in df.dtypes])).encode())
    digest.update(df.select_dtypes(include="number").sum().to_numpy(dtype="float64").tobytes())
    if len(df):
        step = max(len(df) // _FINGERPRINT_SAMPLE_ROWS, 1)
        rows = df.iloc[::step].astype(str)
        digest.update(pd.util.hash_pandas_object(rows, index=False).values.tobytes())
        digest.update(df.iloc[[-1]].astype(str).to_csv(index=False).encode())
    return digest.hexdigest()


def default_cache_dir() -> Path:
    """``<workspace>/insight_tables``, or a temp dir outside a workspace."""
    workspace = st.session_state.get("workspace")
    base = Path(workspace) if workspace else Path(tempfile.gettempdir(), "umetaflow")
    return base / "insight_tables"


def cache_id_for(name: str, df: pd.DataFrame) -> str:
    """Cache id that changes with the table's content, so a new run never shows old rows.

    Rebuilding on every rerun (``regenerate_cache=True``) instead hands the
    browser new data each time and the component never leaves "Loading...".
    """
    return f"{name}_{fingerprint_dataframe(df)[:16]}"


def drop_stale_caches(cache_dir: Path, name: str, keep: str) -> None:
    """Remove earlier caches of the component ``name`` (other content fingerprints)."""
    for stale in cache_dir.glob(f"{name}_" + "?" * 16):
        if stale.name != keep:
            shutil.rmtree(stale, ignore_errors=True)


def to_polars(df: pd.DataFrame, column_overrides: dict | None = None) -> tuple[pl.DataFrame, list[dict]]:
    """Convert to polars with Tabulator-safe field names and column definitions.

    Tabulator reads a "." in a field name as a nested lookup (every sample
    column ends in ".mzML"), so such columns get a safe field name and keep
    their original name as the column title. A ``row`` index is added for
    linking components.
    """
    df = df.reset_index(drop=True)
    renames, column_definitions, used = {}, [], {ROW_FIELD}
    for column in df.columns:
        field = str(column).replace(".", "_").replace(" ", "_")
        while field in used:
            field = f"{field}_"
        used.add(field)
        renames[column] = field
        col_def = {"field": field, "title": str(column), "headerTooltip": True}
        if pd.api.types.is_bool_dtype(df[column]):
            col_def["sorter"] = "boolean"
            col_def["formatter"] = "tickCross"
            col_def["formatterParams"] = {"crossElement": False}
        elif pd.api.types.is_numeric_dtype(df[column]):
            col_def["sorter"] = "number"
            col_def["hozAlign"] = "right"
        col_def.update((column_overrides or {}).get(column, {}))
        column_definitions.append(col_def)
    df = df.rename(columns=renames)
    for column in df.columns:
        if df[column].dtype == object:
            df[column] = df[column].map(lambda v: None if v is None or v != v else str(v))
    return pl.from_pandas(df).with_row_index(ROW_FIELD), column_definitions

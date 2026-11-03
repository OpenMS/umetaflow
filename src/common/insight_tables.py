"""openms-insight helpers: streamed tables and content-keyed plot caches.

``st.dataframe`` ships the whole table to the browser, which stalls on
full-scale feature matrices. An openms-insight ``Table`` caches the data on
disk once and sends one page at a time; its download button exports every row
that passes the table's filters. Same approach as quantms-web's helper.
"""
import hashlib
import json
import shutil
import tempfile
from pathlib import Path

import pandas as pd
import polars as pl
import streamlit as st

ROW_FIELD = "row"


def fingerprint_dataframe(df: pd.DataFrame) -> str:
    """Content hash of a dataframe (shape, schema and every cell)."""
    digest = hashlib.sha256()
    digest.update(repr((df.shape, [str(c) for c in df.columns], [str(t) for t in df.dtypes])).encode())
    # Every row: feature tables are small enough, and a sampled hash can miss
    # a change in an unsampled row and reuse a stale cache
    digest.update(pd.util.hash_pandas_object(df.astype(str), index=False).values.tobytes())
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


def cache_ready(cache_dir: Path, cache_id: str) -> bool:
    """True once openms-insight finished a cache: it writes manifest.json last.

    A manifest that does not parse (e.g. truncated by a crash mid-write) means
    the cache is unusable; it is removed so the caller rebuilds it.
    """
    cache = Path(cache_dir) / cache_id
    try:
        with open(cache / "manifest.json", encoding="utf-8") as f:
            return isinstance(json.load(f), dict)
    except FileNotFoundError:
        return False
    except (OSError, ValueError):
        shutil.rmtree(cache, ignore_errors=True)
        return False


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
            df[column] = df[column].map(lambda v: None if pd.isna(v) else str(v))
    return pl.from_pandas(df).with_row_index(ROW_FIELD), column_definitions

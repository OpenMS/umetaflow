"""Linked feature table and feature map for the UmetaFlow results page.

The MZmine / MS-DIAL core screen: a feature list with one area column per
sample, a 2D feature map (RT vs m/z), and the selected feature's shapes and
intensities below. Table and map are openms-insight components linked through
one ``feature`` selection, so clicking either selects the feature in both.
"""
from pathlib import Path

import pandas as pd
import streamlit as st
from openms_insight import Heatmap, StateManager, Table

from src.feature_tags import tags_column
from src.common.insight_tables import (
    ROW_FIELD,
    cache_id_for,
    cache_ready,
    default_cache_dir,
    drop_stale_caches,
    to_polars,
)

SELECTION = "feature"

# Annotation sources and the column prefixes their tools write into the matrix
ANNOTATION_SOURCES = {
    "spectral library": ("SpectralMatch",),
    "SIRIUS": ("SIRIUS_",),
    "CSI:FingerID": ("CSI:FingerID",),
    "CANOPUS": ("CANOPUS",),
    "MS2Query": ("MS2Query",),
}


def sample_columns(df: pd.DataFrame) -> list[str]:
    """Per-sample intensity columns of a UmetaFlow feature matrix, sorted by name."""
    return sorted(c for c in df.columns if str(c).endswith(".mzML"))


def _has_value(series: pd.Series) -> pd.Series:
    return series.notna() & (series.astype(str).str.strip() != "")


def annotation_sources(df: pd.DataFrame) -> pd.Series:
    """Comma-separated list of the tools that annotated each feature ('' if none)."""
    found = pd.DataFrame(index=df.index)
    for label, prefixes in ANNOTATION_SOURCES.items():
        cols = [c for c in df.columns if any(str(c).startswith(p) for p in prefixes)]
        found[label] = df[cols].apply(_has_value).any(axis=1) if cols else False
    return found.apply(lambda row: ", ".join(k for k, v in row.items() if v), axis=1)


def build_feature_table(df: pd.DataFrame, tags: dict[str, list[str]] | None = None) -> pd.DataFrame:
    """Display table for the feature explorer, one row per consensus feature.

    Columns: feature id, m/z, RT, charge, adduct (if detected), number of
    samples the feature was found in, maximum intensity, annotation sources,
    user tags (``tags``: feature id -> list of tags), re-quantified flag and
    one intensity column per sample (named by file stem).
    Row order follows ``df``, so row ``i`` of the result is ``df.iloc[i]``.
    """
    samples = sample_columns(df)
    intensities = df[samples].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    table = pd.DataFrame(
        {
            "feature": df.index.astype(str),
            "m/z": df["mz"].astype(float).round(4).to_numpy(),
            "RT (s)": df["RT"].astype(float).round(1).to_numpy(),
            "charge": df["charge"].to_numpy(),
        }
    )
    if "adduct" in df.columns:
        table["adduct"] = df["adduct"].fillna("").astype(str).to_numpy()
    table["found in"] = (intensities > 0).sum(axis=1).to_numpy()
    table["max intensity"] = intensities.max(axis=1).round(0).to_numpy()
    table["annotated by"] = annotation_sources(df).to_numpy()
    table["tags"] = tags_column(df, tags or {}).to_numpy()
    if "re-quantified" in df.columns:
        table["re-quantified"] = df["re-quantified"].fillna(False).astype(bool).to_numpy()
    # Sample columns last: when the table is too wide, Tabulator folds the
    # trailing columns into a sub-row, and those should be the per-sample values
    for col in samples:
        name = Path(col).stem
        while name in table.columns:  # a sample named e.g. "charge" must not replace that column
            name += " (sample)"
        table[name] = intensities[col].round(0).to_numpy()
    return table


def _decimals(n: int) -> dict:
    return {"formatter": "money", "formatterParams": {"thousand": False, "precision": n}}


def _ensure_caches(table: pd.DataFrame, cache_dir: Path) -> tuple[str, str]:
    """Build the table and feature-map caches for ``table`` if they don't exist yet."""
    table_id = cache_id_for("feature_table", table)
    map_id = cache_id_for("feature_map", table)
    if cache_ready(cache_dir, table_id) and cache_ready(cache_dir, map_id):
        return table_id, map_id
    cache_dir.mkdir(parents=True, exist_ok=True)
    drop_stale_caches(cache_dir, "feature_table", table_id)
    drop_stale_caches(cache_dir, "feature_map", map_id)

    overrides = {"m/z": _decimals(4), "RT (s)": _decimals(1), "max intensity": _decimals(0)}
    overrides.update({c: _decimals(0) for c in table.columns[table.columns.get_loc("max intensity") + 1:]
                      if pd.api.types.is_float_dtype(table[c]) and not pd.api.types.is_bool_dtype(table[c])})
    pl_df, column_definitions = to_polars(table, overrides)
    Table(
        cache_id=table_id,
        data=pl_df.lazy(),
        cache_path=str(cache_dir),
        column_definitions=column_definitions,
        index_field=ROW_FIELD,
        interactivity={SELECTION: ROW_FIELD},
        initial_sort=[{"column": "max_intensity", "dir": "desc"}],
        title="feature-matrix",
    )
    Heatmap(
        cache_id=map_id,
        data=pl_df.select([ROW_FIELD, "RT_(s)", "m/z", "max_intensity"]).lazy(),
        cache_path=str(cache_dir),
        x_column="RT_(s)",
        y_column="m/z",
        intensity_column="max_intensity",
        interactivity={SELECTION: ROW_FIELD},
        title="Feature map",
        x_label="retention time (s)",
        y_label="m/z",
        intensity_label="max intensity",
    )
    return table_id, map_id


class FeatureExplorer:
    """Linked feature table and feature map; draw each where the page wants it.

    Both components share one ``StateManager``, so clicking a row or a point
    selects that feature in both. ``selected()`` returns the selected row of
    ``df`` (the table preselects its first row, the most intense feature).
    """

    def __init__(
        self,
        df: pd.DataFrame,
        cache_dir: Path | None = None,
        key: str = "umetaflow",
        tags: dict[str, list[str]] | None = None,
    ):
        self.df = df
        self.key = key
        self.cache_dir = Path(cache_dir) if cache_dir else default_cache_dir()
        table = build_feature_table(df, tags)
        self.table_id, self.map_id = _ensure_caches(table, self.cache_dir)
        # One selection state per set of rows: a selected row index of the
        # unfiltered matrix must not select a different feature after filtering.
        # Tags are left out, so tagging the selected feature keeps it selected.
        rows_id = cache_id_for("feature_rows", table.drop(columns="tags"))
        self.state = StateManager(f"{key}_feature_explorer_{rows_id}")

    def table(self, height: int = 420) -> None:
        Table(cache_id=self.table_id, cache_path=str(self.cache_dir))(
            key=f"{self.key}_feature_table", state_manager=self.state, height=height
        )

    def feature_map(self, height: int = 420) -> None:
        Heatmap(cache_id=self.map_id, cache_path=str(self.cache_dir))(
            key=f"{self.key}_feature_map", state_manager=self.state, height=height
        )

    def selected(self) -> pd.Series | None:
        row = self.state.get_selection(SELECTION)
        if row is None or not 0 <= int(row) < len(self.df):
            return None
        return self.df.iloc[int(row)]

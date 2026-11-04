"""User tags on consensus features (MZmine's feature list row tags).

Tags are free-text labels such as "contaminant", "standard" or "follow up"
that a user puts on features while browsing the results, then filters the
feature table by. They live in ``<results_dir>/feature-tags.json`` next to the
feature matrix, so a new run of the workflow (which replaces ``results/``)
starts without tags. Each entry also keeps the feature's m/z and RT: a tag
whose feature no longer has that m/z and RT is ignored instead of landing on a
different feature.

Streamlit-free.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

TAGS_FILE = "feature-tags.json"
MZ_TOLERANCE = 1e-4
RT_TOLERANCE = 1e-2


def tags_path(results_dir: str | Path) -> Path:
    return Path(results_dir, TAGS_FILE)


def clean_tags(tags) -> list[str]:
    """Strip, drop empty and duplicate tags (case-insensitive), keep the first spelling, sort."""
    seen: dict[str, str] = {}
    for tag in tags or []:
        if tag is None:
            continue
        text = " ".join(str(tag).split())
        if text and text.casefold() not in seen:
            seen[text.casefold()] = text
    return sorted(seen.values(), key=str.casefold)


def load_tags(results_dir: str | Path, df: pd.DataFrame | None = None) -> dict[str, list[str]]:
    """Tags per feature id (``str(df.index)``).

    With ``df`` (the feature matrix), entries whose feature is missing or has
    a different m/z or RT are left out. An unreadable file counts as no tags.
    """
    path = tags_path(results_dir)
    try:
        raw = json.loads(path.read_text())
    except (OSError, ValueError):
        return {}
    if not isinstance(raw, dict):
        return {}
    tags: dict[str, list[str]] = {}
    for feature, entry in raw.items():
        if not isinstance(entry, dict):
            continue
        cleaned = clean_tags(entry.get("tags"))
        if not cleaned:
            continue
        if df is not None and not _matches(df, str(feature), entry):
            continue
        tags[str(feature)] = cleaned
    return tags


def _matches(df: pd.DataFrame, feature: str, entry: dict) -> bool:
    index = df.index.astype(str)
    if feature not in set(index):
        return False
    row = df.loc[df.index[index == feature][0]]
    try:
        return (
            abs(float(row["mz"]) - float(entry.get("mz"))) <= MZ_TOLERANCE
            and abs(float(row["RT"]) - float(entry.get("RT"))) <= RT_TOLERANCE
        )
    except (TypeError, ValueError, KeyError):
        return False


def set_feature_tags(results_dir: str | Path, feature: pd.Series, tags) -> dict[str, list[str]]:
    """Replace the tags of one feature (a row of the feature matrix) and save the file.

    An empty ``tags`` removes the feature's entry. Returns the tags of all
    features as ``load_tags`` (without a matrix) would.
    """
    path = tags_path(results_dir)
    try:
        raw = json.loads(path.read_text())
        if not isinstance(raw, dict):
            raw = {}
    except (OSError, ValueError):
        raw = {}
    key = str(feature.name)
    cleaned = clean_tags(tags)
    if cleaned:
        raw[key] = {"tags": cleaned, "mz": float(feature["mz"]), "RT": float(feature["RT"])}
    else:
        raw.pop(key, None)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(raw, indent=1, sort_keys=True))
    tmp.replace(path)
    return load_tags(results_dir)


def all_tags(tags: dict[str, list[str]]) -> list[str]:
    """Every tag in use, sorted case-insensitively."""
    return clean_tags(t for feature_tags in tags.values() for t in feature_tags)


def tags_column(df: pd.DataFrame, tags: dict[str, list[str]]) -> pd.Series:
    """Comma-separated tags per row of ``df`` ('' for untagged features)."""
    return pd.Series(
        [", ".join(tags.get(str(i), [])) for i in df.index], index=df.index, dtype=object
    )


def has_any_tag(df: pd.DataFrame, tags: dict[str, list[str]], wanted) -> pd.Series:
    """Boolean mask: rows tagged with at least one of ``wanted`` (case-insensitive)."""
    wanted_keys = {t.casefold() for t in clean_tags(wanted)}
    return pd.Series(
        [bool(wanted_keys & {t.casefold() for t in tags.get(str(i), [])}) for i in df.index],
        index=df.index,
    )

"""Tests for user tags on features (src/feature_tags.py)."""
import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.feature_explorer import build_feature_table  # noqa: E402
from src.feature_tags import (  # noqa: E402
    TAGS_FILE,
    all_tags,
    clean_tags,
    has_any_tag,
    load_tags,
    set_feature_tags,
    tags_column,
)


def _matrix():
    return pd.DataFrame(
        {
            "charge": [1, 1, 2],
            "RT": [301.228, 120.04, 540.12],
            "mz": [235.118873, 180.06, 755.39309],
            "A.mzML": [5.9e6, 0.0, 1.05e5],
        },
        index=pd.Index(["f1", "f2", "f3"], name="metabolite"),
    )


def test_clean_tags_strips_dedupes_case_insensitively_and_sorts():
    assert clean_tags([" follow  up ", "Contaminant", "contaminant", "", None and "x", "blank"]) == [
        "blank", "Contaminant", "follow up",
    ]


def test_set_and_load_round_trip(tmp_path):
    df = _matrix()
    set_feature_tags(tmp_path, df.loc["f1"], ["contaminant"])
    tags = set_feature_tags(tmp_path, df.loc["f3"], ["standard", "follow up"])
    assert tags == {"f1": ["contaminant"], "f3": ["follow up", "standard"]}
    assert load_tags(tmp_path, df) == tags
    assert all_tags(tags) == ["contaminant", "follow up", "standard"]


def test_empty_tags_remove_the_entry(tmp_path):
    df = _matrix()
    set_feature_tags(tmp_path, df.loc["f1"], ["contaminant"])
    assert set_feature_tags(tmp_path, df.loc["f1"], []) == {}
    assert json.loads((tmp_path / TAGS_FILE).read_text()) == {}


def test_tags_of_a_changed_or_missing_feature_are_ignored(tmp_path):
    df = _matrix()
    set_feature_tags(tmp_path, df.loc["f1"], ["contaminant"])
    set_feature_tags(tmp_path, df.loc["f2"], ["standard"])
    moved = df.copy()
    moved.loc["f1", "RT"] = 400.0  # same id, different feature after a new run
    assert load_tags(tmp_path, moved) == {"f2": ["standard"]}
    assert load_tags(tmp_path, df.drop(index="f2")) == {"f1": ["contaminant"]}


def test_unreadable_file_means_no_tags(tmp_path):
    (tmp_path / TAGS_FILE).write_text("{not json")
    assert load_tags(tmp_path) == {}
    # and saving replaces it
    assert set_feature_tags(tmp_path, _matrix().loc["f2"], ["x"]) == {"f2": ["x"]}


def test_tags_column_and_mask():
    df = _matrix()
    tags = {"f1": ["contaminant"], "f3": ["follow up", "standard"]}
    assert tags_column(df, tags).tolist() == ["contaminant", "", "follow up, standard"]
    assert has_any_tag(df, tags, ["Standard", "nope"]).tolist() == [False, False, True]
    assert has_any_tag(df, tags, []).tolist() == [False, False, False]


def test_feature_table_has_a_tags_column_before_the_samples():
    df = _matrix()
    table = build_feature_table(df, {"f2": ["blank"]})
    assert table["tags"].tolist() == ["", "blank", ""]
    cols = list(table.columns)
    assert cols.index("tags") < cols.index("A")
    assert build_feature_table(df)["tags"].tolist() == ["", "", ""]

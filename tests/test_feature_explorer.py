"""Tests for the feature explorer's table building (src/feature_explorer.py)."""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.feature_explorer import annotation_sources, build_feature_table, sample_columns  # noqa: E402
from src.common.insight_tables import ROW_FIELD, cache_id_for, to_polars  # noqa: E402


def _matrix():
    df = pd.DataFrame(
        {
            "charge": [1, 1, 2],
            "RT": [301.228, 120.04, 540.12],
            "mz": [235.118873, 180.06, 755.39309],
            "adduct": ["[M+H]+", None, "[M+Na]+"],
            "Treatment.mzML": [5.9e6, 0.0, 1.05e5],
            "Blank.mzML": [5.3e6, 2.0e3, None],
            "Blank.mzML_IDs": ["1", "2", None],
            "re-quantified": [False, True, False],
            "SpectralMatch": ["caffeine", "", None],
            "SIRIUS_Blank_molecularFormula": [None, "C6H12O6", None],
        },
        index=pd.Index(["235.1189@301.23", "180.06@120.04", "755.3931@540.12"], name="metabolite"),
    )
    return df


def test_sample_columns_sorted_and_without_id_columns():
    assert sample_columns(_matrix()) == ["Blank.mzML", "Treatment.mzML"]


def test_annotation_sources_ignores_empty_strings():
    sources = annotation_sources(_matrix()).tolist()
    assert sources == ["spectral library", "SIRIUS", ""]


def test_build_feature_table_rows_follow_matrix():
    table = build_feature_table(_matrix())
    assert table["feature"].tolist() == ["235.1189@301.23", "180.06@120.04", "755.3931@540.12"]
    assert table.columns.tolist() == [
        "feature", "m/z", "RT (s)", "charge", "adduct", "found in", "max intensity",
        "annotated by", "tags", "re-quantified", "Blank", "Treatment",
    ]
    # Missing and zero intensities both count as not found
    assert table["found in"].tolist() == [2, 1, 1]
    assert table["max intensity"].tolist() == [5.9e6, 2.0e3, 1.05e5]
    assert table["adduct"].tolist() == ["[M+H]+", "", "[M+Na]+"]


def test_to_polars_makes_tabulator_safe_fields():
    table = build_feature_table(_matrix())
    pl_df, defs = to_polars(table)
    fields = {d["title"]: d["field"] for d in defs}
    assert fields["m/z"] == "m/z"
    assert fields["RT (s)"] == "RT_(s)"
    assert all("." not in f for f in fields.values())
    assert pl_df[ROW_FIELD].to_list() == [0, 1, 2]
    assert pl_df.columns[0] == ROW_FIELD


def test_cache_id_changes_with_content():
    table = build_feature_table(_matrix())
    other = table.copy()
    other.loc[0, "max intensity"] = 1.0
    assert cache_id_for("feature_table", table) == cache_id_for("feature_table", table.copy())
    assert cache_id_for("feature_table", table) != cache_id_for("feature_table", other)


def test_sample_named_like_a_display_column_gets_its_own_column():
    df = _matrix().rename(columns={"Treatment.mzML": "charge.mzML"})
    table = build_feature_table(df)
    assert table["charge"].tolist() == [1, 1, 2]
    assert table["charge (sample)"].tolist() == [5.9e6, 0.0, 1.05e5]


def test_to_polars_handles_pd_na_in_object_columns():
    df = pd.DataFrame({"adduct": pd.array(["[M+H]+", pd.NA], dtype=object)})
    pl_df, _ = to_polars(df)
    assert pl_df["adduct"].to_list() == ["[M+H]+", None]


def test_cache_ready_needs_manifest(tmp_path):
    from src.common.insight_tables import cache_ready
    (tmp_path / "t").mkdir()
    assert not cache_ready(tmp_path, "t")
    (tmp_path / "t" / "manifest.json").write_text("{}")
    assert cache_ready(tmp_path, "t")


def test_cache_ready_rejects_truncated_manifest(tmp_path):
    from src.common.insight_tables import cache_ready

    (tmp_path / "t").mkdir()
    (tmp_path / "t" / "manifest.json").write_text('{"version": ')
    assert not cache_ready(tmp_path, "t")
    assert not (tmp_path / "t").exists()

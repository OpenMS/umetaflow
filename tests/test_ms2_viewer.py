"""Tests for the MS2 / library mirror data loading (src/ms2_viewer.py)."""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.ms2_viewer import load_ms2_data, match_peaks, peaks_frame  # noqa: E402

QUERY_MGF = """BEGIN IONS
SCANS=1
FEATURE_ID=e_111
PEPMASS=508.003
100.0 50
136.0615 1000
200.0 10
END IONS
BEGIN IONS
SCANS=2
FEATURE_ID=e_222
PEPMASS=191.0767
191.0767 400
END IONS
"""

LIBRARY_MGF = """BEGIN IONS
NAME=ADP
PEPMASS=428.0367
136.0615 100
END IONS
BEGIN IONS
NAME=ATP
PEPMASS=508.003
136.0616 100
300.0 20
END IONS
"""

MZTAB = "\n".join([
    "MTD\tmzTab-version\t1.0.0",
    "",
    "SMH\tidentifier\tsmiles\tdescription\texp_mass_to_charge\topt_ppm_error\topt_match_score\topt_spec_native_id",
    "SML\tnull\tNC1=NC=NC2=C1N=CN2\tATP\t508.003\t-0.01\t48.2\tindex=0",
    "",
])


def _results(tmp_path, with_matches=True):
    res = tmp_path / "results"
    (res / "gnps-export").mkdir(parents=True)
    (res / "consensus-dfs").mkdir()
    (res / "gnps-export" / "MS2.mgf").write_text(QUERY_MGF)
    pd.DataFrame(
        {"consensus_feature_id": ["111", "222", "333"]},
        index=pd.Index(["508.003@264.84", "191.0767@123.83", "99.0@10.0"], name="metabolite"),
    ).to_parquet(res / "consensus-dfs" / "feature-matrix-gnps.parquet")
    if with_matches:
        (res / "spectral-matcher").mkdir()
        (res / "spectral-matcher" / "MS2-matches.mzTab").write_text(MZTAB)
    library = tmp_path / "library.mgf"
    library.write_text(LIBRARY_MGF)
    return res, library


def test_match_peaks_pairs_within_tolerance_once():
    mq, ml, cosine = match_peaks(
        np.array([100.0, 136.0615]), np.array([5.0, 100.0]),
        np.array([136.0616, 136.07, 300.0]), np.array([100.0, 90.0, 20.0]), tol_ppm=10,
    )
    assert mq.tolist() == [False, True]
    assert ml.tolist() == [True, False, False]
    assert 0 < cosine < 1


def test_match_peaks_identical_spectra_cosine_one():
    mz, inty = np.array([10.0, 20.0]), np.array([1.0, 4.0])
    _, _, cosine = match_peaks(mz, inty, mz, inty, tol_ppm=1)
    assert abs(cosine - 1.0) < 1e-9


def test_load_ms2_data_links_features_spectra_and_library(tmp_path):
    res, library = _results(tmp_path)
    data = load_ms2_data(res, library, frag_tol_ppm=10)
    assert data.feature_spectra == {"508.003@264.84": [0], "191.0767@123.83": [1], "99.0@10.0": []}
    [match] = data.matches["508.003@264.84"]
    assert (match.name, match.query_id, match.score) == ("ATP", 0, 48.2)
    # The library spectrum named ATP, not the first spectrum with that peak
    lib = data.spectra.set_index("spectrum_id").loc[match.library_id]
    assert (lib["kind"], lib["label"], lib["precursor_mz"]) == ("library", "ATP", 508.003)
    assert match.matched_peaks == 1
    # Intensities are relative to the base peak
    q = data.peaks[data.peaks["spectrum_id"] == 0]
    assert q["intensity"].max() == 100.0
    assert q.loc[q["matched"], "mz"].tolist() == [136.0615]


def test_load_ms2_data_without_matches_or_export(tmp_path):
    res, library = _results(tmp_path, with_matches=False)
    data = load_ms2_data(res, library)
    assert data.matches == {}
    assert not data.peaks["matched"].any()
    assert load_ms2_data(tmp_path / "missing", library) is None


def test_peaks_frame_labels_only_strong_matched_peaks(tmp_path):
    res, library = _results(tmp_path)
    frame = peaks_frame(load_ms2_data(res, library, frag_tol_ppm=10))
    labelled = frame.filter(frame["label"] != "")
    assert sorted(labelled["mz"].to_list()) == [136.0615, 136.0616]
    assert frame["peak_id"].to_list() == list(range(frame.height))

"""MS2 spectra of features and library-match mirror plots for the results page.

MZmine, MS-DIAL and GNPS show a spectral library hit as a mirror plot: the
feature's MS2 spectrum on top, the library spectrum below, matched peaks
highlighted. UmetaFlow keeps everything needed in its results folder:

- ``gnps-export/MS2.mgf``: one MS2 spectrum per feature and file (GNPSExport,
  ``FEATURE_ID=e_<consensus id>`` of ``feature-matrix-gnps``),
- ``consensus-dfs/feature-matrix-gnps.parquet``: those consensus ids per
  feature (same ``m/z@RT`` feature names as the final matrix),
- ``spectral-matcher/MS2-matches.mzTab``: MetaboliteSpectralMatcher hits; the
  query is spectrum ``index=N`` of ``MS2.mzML``, which is the N-th spectrum of
  ``MS2.mgf`` (FileConverter keeps the order); ``description`` is the library
  spectrum's ``NAME``,
- the MS2 library (mgf) the user uploaded.

So existing results get the viewer without re-running the workflow.
"""
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
from pyteomics import mgf

QUERY_SIDE = "ms2_query"
LIBRARY_SIDE = "ms2_library"


@dataclass
class Match:
    query_id: int
    library_id: int | None
    name: str
    score: float | None
    ppm_error: float | None
    smiles: str = ""
    cosine: float | None = None
    matched_peaks: int = 0


@dataclass
class MS2Data:
    peaks: pd.DataFrame  # spectrum_id, mz, intensity (% of base peak), matched
    spectra: pd.DataFrame  # spectrum_id, kind, label, precursor_mz
    feature_spectra: dict = field(default_factory=dict)  # feature name -> [query spectrum ids]
    matches: dict = field(default_factory=dict)  # feature name -> [Match]

    def has_feature(self, feature: str) -> bool:
        return bool(self.feature_spectra.get(feature))


def _read_mgf(path: Path) -> list[dict]:
    with mgf.MGF(str(path), use_header=False, convert_arrays=1, read_charges=False) as reader:
        return list(reader)


def _float(value) -> float | None:
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    return None if np.isnan(f) else f


def _precursor(params: dict) -> float | None:
    pepmass = params.get("pepmass")
    if isinstance(pepmass, (tuple, list)):
        pepmass = pepmass[0]
    return _float(pepmass)


def match_peaks(mz_a: np.ndarray, int_a: np.ndarray, mz_b: np.ndarray, int_b: np.ndarray, tol_ppm: float):
    """Greedy one-to-one peak matching by intensity within ``tol_ppm``.

    Returns (matched mask of a, matched mask of b, cosine of the matched
    intensities with square-root weighting, as in GNPS-style comparisons).
    """
    matched_a = np.zeros(len(mz_a), dtype=bool)
    matched_b = np.zeros(len(mz_b), dtype=bool)
    if len(mz_a) == 0 or len(mz_b) == 0:
        return matched_a, matched_b, 0.0
    pairs = []
    for i in np.argsort(-int_a):
        tol = mz_a[i] * tol_ppm * 1e-6
        candidates = np.where(~matched_b & (np.abs(mz_b - mz_a[i]) <= tol))[0]
        if len(candidates):
            j = candidates[np.argmax(int_b[candidates])]
            matched_a[i] = matched_b[j] = True
            pairs.append((i, j))
    wa, wb = np.sqrt(int_a), np.sqrt(int_b)
    norm = np.linalg.norm(wa) * np.linalg.norm(wb)
    cosine = float(sum(wa[i] * wb[j] for i, j in pairs) / norm) if norm else 0.0
    return matched_a, matched_b, cosine


def _read_matches(mztab: Path) -> pd.DataFrame:
    """Small-molecule section of a MetaboliteSpectralMatcher mzTab."""
    rows, header = [], None
    with open(mztab, encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "SMH":
                header = parts[1:]
            elif parts[0] == "SML" and header:
                rows.append(dict(zip(header, parts[1:])))
    return pd.DataFrame(rows, columns=header or [])


def load_ms2_data(results_dir: Path, library: Path | None, frag_tol_ppm: float = 500.0) -> MS2Data | None:
    """Collect query/library spectra, feature -> spectra and feature -> matches.

    Returns None when the run did not export MS2 spectra (GNPS export,
    spectral matching and MS2Query all off).
    """
    results_dir = Path(results_dir)
    mgf_path = results_dir / "gnps-export" / "MS2.mgf"
    gnps_matrix = results_dir / "consensus-dfs" / "feature-matrix-gnps.parquet"
    if not mgf_path.is_file() or not gnps_matrix.is_file():
        return None

    queries = _read_mgf(mgf_path)
    peaks, spectra = [], []
    by_consensus: dict[str, list[int]] = {}
    for idx, spec in enumerate(queries):
        params = spec["params"]
        mz, inty = spec["m/z array"], spec["intensity array"]
        base = inty.max() if len(inty) else 1.0
        peaks.append(pd.DataFrame({"spectrum_id": idx, "mz": mz, "intensity": 100.0 * inty / base, "matched": False}))
        feature_id = str(params.get("feature_id", "")).removeprefix("e_")
        by_consensus.setdefault(feature_id, []).append(idx)
        spectra.append({"spectrum_id": idx, "kind": "query", "label": f"scan {params.get('scans', idx)}",
                        "precursor_mz": _precursor(params)})

    gnps = pd.read_parquet(gnps_matrix, columns=["consensus_feature_id"])
    feature_spectra = {
        str(name): by_consensus.get(str(cid), [])
        for name, cid in zip(gnps.index, gnps["consensus_feature_id"].astype(str))
    }

    matches: dict[str, list[Match]] = {}
    mztab = results_dir / "spectral-matcher" / "MS2-matches.mzTab"
    if mztab.is_file():
        library_spectra = _read_mgf(library) if library and Path(library).is_file() else []
        query_to_feature = {q: f for f, qs in feature_spectra.items() for q in qs}
        next_id = len(queries)
        library_ids: dict[int, int] = {}  # library index -> spectrum_id
        for _, hit in _read_matches(mztab).iterrows():
            native = str(hit.get("opt_spec_native_id", ""))
            if not native.startswith("index="):
                continue
            query_id = int(native.removeprefix("index="))
            feature = query_to_feature.get(query_id)
            if feature is None:
                continue
            name = str(hit.get("description", ""))
            precursor = _float(hit.get("exp_mass_to_charge"))
            # Library spectrum: same NAME, closest precursor m/z
            candidates = [i for i, s in enumerate(library_spectra) if str(s["params"].get("name", "")) == name]
            lib_idx = None
            if candidates:
                lib_idx = min(candidates, key=lambda i: abs((_precursor(library_spectra[i]["params"]) or 0) - (precursor or 0)))
            match = Match(
                query_id=query_id, library_id=None, name=name,
                score=_float(hit.get("opt_match_score")), ppm_error=_float(hit.get("opt_ppm_error")),
                smiles="" if str(hit.get("smiles", "null")) == "null" else str(hit.get("smiles")),
            )
            if lib_idx is not None:
                if lib_idx not in library_ids:
                    spec = library_spectra[lib_idx]
                    inty = spec["intensity array"]
                    base = inty.max() if len(inty) else 1.0
                    library_ids[lib_idx] = next_id
                    peaks.append(pd.DataFrame({"spectrum_id": next_id, "mz": spec["m/z array"],
                                               "intensity": 100.0 * inty / base, "matched": False}))
                    spectra.append({"spectrum_id": next_id, "kind": "library", "label": name,
                                    "precursor_mz": _precursor(spec["params"])})
                    next_id += 1
                match.library_id = library_ids[lib_idx]
            matches.setdefault(feature, []).append(match)

    peaks_df = pd.concat(peaks, ignore_index=True) if peaks else pd.DataFrame(
        columns=["spectrum_id", "mz", "intensity", "matched"])
    peaks_df["spectrum_id"] = peaks_df["spectrum_id"].astype(int)
    data = MS2Data(peaks=peaks_df, spectra=pd.DataFrame(spectra), feature_spectra=feature_spectra, matches=matches)
    _mark_matched_peaks(data, frag_tol_ppm)
    return data


def _mark_matched_peaks(data: MS2Data, tol_ppm: float) -> None:
    """Set ``matched`` on peaks that pair up in a match, and each match's cosine."""
    groups = {sid: idx for sid, idx in data.peaks.groupby("spectrum_id").groups.items()}
    matched = data.peaks["matched"].to_numpy().copy()
    for feature_matches in data.matches.values():
        for m in feature_matches:
            if m.library_id is None or m.query_id not in groups or m.library_id not in groups:
                continue
            qa, la = groups[m.query_id], groups[m.library_id]
            q, lib = data.peaks.loc[qa], data.peaks.loc[la]
            mq, ml, m.cosine = match_peaks(q["mz"].to_numpy(), q["intensity"].to_numpy(),
                                            lib["mz"].to_numpy(), lib["intensity"].to_numpy(), tol_ppm)
            m.matched_peaks = int(mq.sum())
            matched[data.peaks.index.get_indexer(qa[mq])] = True
            matched[data.peaks.index.get_indexer(la[ml])] = True
    data.peaks["matched"] = matched


def peaks_frame(data: MS2Data) -> pl.DataFrame:
    """Peaks for the insight caches, with an m/z label on matched peaks."""
    df = data.peaks.copy()
    df["peak_id"] = np.arange(len(df))
    # Label matched peaks of at least 5 % relative intensity; more labels overlap
    df["label"] = np.where(df["matched"] & (df["intensity"] >= 5.0), df["mz"].round(4).astype(str), "")
    return pl.from_pandas(df[["peak_id", "spectrum_id", "mz", "intensity", "matched", "label"]])

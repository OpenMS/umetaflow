import json
import sys
from pathlib import Path
import pandas as pd
import shutil

from ms2query.run_ms2query import (
    run_complete_folder,
    download_zenodo_files,
)
from ms2query.ms2library import create_library_object_from_one_dir
from ms2query.utils import SettingsRunMS2Query

DEFAULTS = [
    {"key": "in", "value": [], "help": "Feature Matrix tsv file", "hide": True},
    {"key": "in_mgf", "value": [], "help": "GNPS mgf file", "hide": True},
    {
        "key": "out_m2query_csv",
        "value": [],
        "help": "MS2Query output file.",
        "hide": True,
    },
    {
        "key": "ion_mode",
        "value": "positive",
        "help": "Ion mode for MS2Query.",
        "hide": True,
    },
]


def allow_ms2deepscore_model_types():
    """Let torch's weights-only load read the MS2Deepscore model file.

    From torch 2.6 on, ``torch.load`` only unpickles tensors and plain Python
    types unless told otherwise. The MS2Deepscore model that MS2Query loads
    also stores its settings, which hold numpy arrays. Allowing exactly those
    numpy types keeps the check on for everything else, instead of loading
    the file with ``weights_only=False``.
    """
    try:
        import numpy as np
        import torch
    except ImportError:
        return
    add_safe_globals = getattr(torch.serialization, "add_safe_globals", None)
    if add_safe_globals is None:  # torch before 2.4 loads without the check
        return
    multiarray = np.core.multiarray
    add_safe_globals(
        [multiarray._reconstruct, multiarray.scalar, np.ndarray, np.dtype]
        + [
            type(np.dtype(t))
            for t in ("float64", "float32", "float16", "int64", "int32", "int16", "int8", "uint8", "bool")
        ]
    )


def get_params():
    if len(sys.argv) > 1:
        with open(sys.argv[1], "r") as f:
            return json.load(f)
    else:
        return {}


def download_ms2query_libraries(ms2query_library_files_directory, ion_mode, flag_file):
    print("Downloading MS2Query Models...")
    if ms2query_library_files_directory.exists():
        shutil.rmtree(ms2query_library_files_directory)

    ms2query_library_files_directory.mkdir(exist_ok=True, parents=True)

    # Downloads pretrained models and files for MS2Query (>2GB download)
    download_zenodo_files(ion_mode, ms2query_library_files_directory)

    Path(flag_file).touch()


def ms2query_annotations(feature_matrix, ms2query_csv):
    df_gnps = pd.read_parquet(
        Path(Path(feature_matrix).parent, "feature-matrix-gnps.parquet"))

    df_ms2query = pd.read_csv(ms2query_csv)
    df_ms2query["feature_id"] = df_ms2query["feature_id"].apply(lambda x: int(x[2:]))
    df_ms2query = df_ms2query.set_index("feature_id")

    ms2query_columns = [
        "ms2query_model_prediction",
        "precursor_mz_difference",
        "inchikey",
        "analog_compound_name",
        "smiles",
        "cf_kingdom",
        "cf_superclass",
        "cf_class",
        "cf_subclass",
        "cf_direct_parent",
        "npc_class_results",
        "npc_superclass_results",
        "npc_pathway_results",
    ]

    df = pd.read_parquet(Path(feature_matrix).with_suffix(".parquet"))

    for i in df_gnps["consensus_feature_id"]:
        if i in df_ms2query.index:
            for col in ms2query_columns:
                df.loc[
                    df_gnps.index[df_gnps["consensus_feature_id"] == i].tolist()[0],
                    f"MS2Query_{col}",
                ] = str(df_ms2query.loc[i, col])

    df.to_parquet(feature_matrix)
    df.to_csv(Path(feature_matrix).with_suffix(".tsv"), sep="\t")


if __name__ == "__main__":
    params = get_params()

    consensus_file = params["in"][0]
    mgf_spectra = params["in_mgf"][0]
    results_file = params["out_ms2query_csv"][0]

    # Set the location where downloaded library and model files are stored (directly in workspace)
    ms2query_library_files_directory = Path(
        Path(Path(results_file).parent.parent.parent.parent, "ms2query-models"),
        params["ion_mode"],
    )
    flag = Path(ms2query_library_files_directory, "download-complete")
    if not flag.exists():
        download_ms2query_libraries(
            ms2query_library_files_directory, params["ion_mode"], flag
        )

    allow_ms2deepscore_model_types()

    # Create a MS2Library object
    ms2library = create_library_object_from_one_dir(ms2query_library_files_directory)

    if Path(results_file).exists():
        Path(results_file).unlink()

    run_complete_folder(
        ms2library=ms2library,
        folder_with_spectra=Path(mgf_spectra).parent,
        results_folder=Path(results_file).parent,
        settings=SettingsRunMS2Query(additional_metadata_columns=("FEATURE_ID",)),
    )

    ms2query_annotations(consensus_file, results_file)
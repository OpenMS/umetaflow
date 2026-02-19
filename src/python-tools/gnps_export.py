import json
import sys
import pyopenms as poms

############################
# default paramter values #
###########################
#
# Mandatory keys for each parameter
# key: a unique identifier
# value: the default value
#
# Optional keys for each parameter
# name: the name of the parameter
# hide: don't show the parameter in the parameter section (e.g. for input/output files)
# options: a list of valid options for the parameter
# min: the minimum value for the parameter (int and float)
# max: the maximum value for the parameter (int and float)
# step_size: the step size for the parameter (int and float)
# help: a description of the parameter
# widget_type: the type of widget to use for the parameter (default: auto)
# advanced: whether or not the parameter is advanced (default: False)

DEFAULTS = [
    {"key": "in_cm", "value": [], "help": "Input consensusXML file.", "hide": True},
    {"key": "in_mzml", "value": [], "help": "Original mzML files containing the MS2 spectra.", "hide": True},
    {"key": "out", "value": [], "help": "Output MGF file.", "hide": True},
    {"key": "out_quantification", "value": [], "help": "Output feature quantification table.", "hide": True},
    {"key": "out_pairs", "value": [], "help": "Output supplementary pairs table for IIMN.", "hide": True},
    {"key": "out_meta_values", "value": [], "help": "Output meta value file.", "hide": True},
    {
        "key": "output_type",
        "name": "output type",
        "value": "most_intense",
        "options": ["most_intense", "merged_spectra"],
        "help": "Specificity of MGF output information. 'most_intense' returns the most intense MS2 scan per consensus feature. 'merged_spectra' merges all MS2 scans per consensus feature into one representative spectrum.",
    },
    {
        "key": "peptide_cutoff",
        "name": "peptide cutoff",
        "value": 5,
        "min": -1,
        "step_size": 1,
        "help": "Number of most intense peptides to consider per consensus element. Use -1 to consider all identifications.",
    },
    {
        "key": "ms2_bin_size",
        "name": "MS2 bin size (Da)",
        "value": 0.02,
        "min": 0.0,
        "step_size": 0.01,
        "help": "Bin size in Da for fragment ions when merging MS2 scans.",
    },
    {
        "key": "cos_similarity",
        "name": "cosine similarity",
        "value": 0.9,
        "min": 0.0,
        "max": 1.0,
        "step_size": 0.1,
        "help": "Cosine similarity threshold for merged_spectra output.",
        "advanced": True,
    },
]


def get_params():
    if len(sys.argv) > 1:
        with open(sys.argv[1], "r") as f:
            return json.load(f)
    else:
        return {}


if __name__ == "__main__":
    params = get_params()

    consensus_file = params["in_cm"][0]
    mzml_files = params["in_mzml"]
    # Handle nested list from collect=True in FileManager
    if mzml_files and isinstance(mzml_files[0], list):
        mzml_files = mzml_files[0]
    output_mgf = params["out"][0]
    output_quant = params["out_quantification"][0]
    output_pairs = params["out_pairs"][0]
    output_meta = params["out_meta_values"][0]

    # Load consensus map
    consensus_map = poms.ConsensusMap()
    poms.ConsensusXMLFile().load(consensus_file, consensus_map)

    # Conditionally annotate for IIMN if linked groups metadata exists
    for feature in consensus_map:
        if feature.metaValueExists("IIMN_linked_groups"):
            poms.IonIdentityMolecularNetworking.annotateConsensusMap(consensus_map)
            break

    # Set up GNPSMGFFile with algorithm parameters
    gnps_mgf = poms.GNPSMGFFile()
    p = gnps_mgf.getParameters()
    p.setValue(b"output_type", params.get("output_type", "most_intense"))
    p.setValue(b"peptide_cutoff", int(params.get("peptide_cutoff", 5)))
    p.setValue(b"ms2_bin_size", float(params.get("ms2_bin_size", 0.02)))
    p.setValue(
        b"merged_spectra:cos_similarity", float(params.get("cos_similarity", 0.9))
    )
    gnps_mgf.setParameters(p)

    # Export MGF file
    gnps_mgf.store(
        poms.String(consensus_file),
        [f.encode() for f in mzml_files],
        poms.String(output_mgf),
    )

    # Export supplementary pairs table for IIMN
    poms.IonIdentityMolecularNetworking.writeSupplementaryPairTable(
        consensus_map, output_pairs
    )

    # Export feature quantification table
    poms.GNPSQuantificationFile().store(consensus_map, output_quant)

    # Export meta values
    poms.GNPSMetaValueFile().store(consensus_map, output_meta)

    print("GNPS export completed successfully.")

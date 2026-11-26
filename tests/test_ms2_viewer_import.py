"""The workflow imports pyteomics only when it reads an mgf file.

Every UmetaFlow workflow page imports src.UmetaFlowTOPPWorkflow, which imports
src.metabolomicsresults and with it src.ms2_viewer. Importing pyteomics.mgf at
the top of src.ms2_viewer cost about 0.2 s of the first of those pages after a
restart, although only the MS2 view of a results page reads an mgf file.
"""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_workflow_modules_do_not_import_pyteomics_mgf():
    code = (
        "import sys\n"
        "import src.ms2_viewer, src.metabolomicsresults, src.UmetaFlowTOPPWorkflow\n"
        "print('pyteomics.mgf' in sys.modules)\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=ROOT,
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().splitlines()[-1] == "False"

"""The workflow imports plotly express only where a feature's intensities are drawn.

Every UmetaFlow workflow page imports src.UmetaFlowTOPPWorkflow, which imports
src.metabolomicsresults. Importing plotly.express at its top cost about 0.3 s of
the first of those pages after a restart, although only the bar chart of a
feature's intensities on the results page uses it.
"""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_workflow_modules_do_not_import_plotly_express():
    code = (
        "import sys\n"
        "import src.metabolomicsresults, src.UmetaFlowTOPPWorkflow\n"
        "print('plotly.express' in sys.modules)\n"
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

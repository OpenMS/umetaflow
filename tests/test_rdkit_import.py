"""rdkit is imported only where a structure is drawn, not by every page."""
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize(
    "module", ["src.metabolomicsresults", "src.UmetaFlowTOPPWorkflow"]
)
def test_module_does_not_import_rdkit(module):
    code = (
        "import importlib, sys\n"
        f"importlib.import_module({module!r})\n"
        "print('rdkit' in sys.modules)\n"
    )
    out = subprocess.run(
        [sys.executable, "-c", code], cwd=ROOT, capture_output=True, text=True
    )
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip().splitlines()[-1] == "False"

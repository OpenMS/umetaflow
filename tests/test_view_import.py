"""The raw data viewer's module does not load plotly express, which it never uses."""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_view_does_not_import_plotly_express():
    code = (
        "import sys\n"
        "import src.view\n"
        "print('plotly.express' in sys.modules)\n"
    )
    out = subprocess.run(
        [sys.executable, "-c", code], cwd=ROOT, capture_output=True, text=True
    )
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip().splitlines()[-1] == "False"

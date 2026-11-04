"""Smoke tests for UmetaFlow's pages (streamlit.testing AppTest).

App-owned (sync.app_owned in app.yaml): the template's test_gui.py launches
the template's example pages, which UmetaFlow does not ship.
"""
import json
import shutil
from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest

# Imported at collection, before tests/ runs: several unit tests there swap
# `streamlit` for a mock while importing src.common, and the pages must find
# the module already loaded with the real one (as in the template's test_gui.py).
from src import fileupload  # noqa: F401


@pytest.fixture
def launch(request):
    test = AppTest.from_file(request.param)

    ## Initialize session state ##
    with open("settings.json", "r") as f:
        test.session_state.settings = json.load(f)
    test.session_state.settings["test"] = True
    test.secrets["workspace"] = "test"
    return test


# Test launching of all pages
@pytest.mark.parametrize(
    "launch",
    (
        # "content/quickstart.py", # NOTE: this page does not work due to streamlit.errors.StreamlitPageNotFoundError error
        "content/file_upload.py",
        "content/raw_data_viewer.py",
        "content/umetaflow_configure.py",
        "content/umetaflow_run.py",
        "content/umetaflow_results.py",
        "content/mz_calculator.py",
        "content/extracted_ion_chromatograms.py",
        "content/statistics.py",
        "content/workspace_sync.py",
    ),
    indirect=True,
)
def test_launch(launch):
    """Test if all pages can be launched without errors."""
    launch.run(timeout=60)
    assert not launch.exception


@pytest.mark.parametrize("launch", ["content/file_upload.py"], indirect=True)
def test_file_upload_load_example(launch):
    launch.run(timeout=60)
    for i in launch.tabs:
        if i.label == "Example Data":
            i.button[0].click().run(timeout=60)
            assert not launch.exception


@pytest.mark.parametrize(
    "launch,example",
    [
        ("content/raw_data_viewer.py", "Blank.mzML"),
        ("content/raw_data_viewer.py", "Treatment.mzML"),
    ],
    indirect=["launch"],
)
def test_view_raw_ms_data(launch, example):
    launch.run(timeout=60)

    mzML_dir = Path(launch.session_state.workspace, "mzML-files")
    for f in Path("example-data", "mzML").glob("*.mzML"):
        try:
            shutil.copy(f, mzML_dir)
        except shutil.SameFileError:
            pass  # already a symlink to the same source (on Linux)
    launch.run(timeout=60)

    launch.selectbox[0].select(example).run(timeout=60)
    assert not launch.exception


@pytest.mark.parametrize("launch", ["content/topp_workflow_parameter.py"], indirect=True)
def test_parameter_page_has_no_advanced_toggle(launch):
    """Advanced parameters are reached by search (TOPP) or an expander (python
    tools), not by a page-wide toggle."""
    launch.run(timeout=30)
    assert not launch.exception
    assert not [t for t in launch.toggle if t.label == "Show advanced parameters"]
    # src/python-tools/example.py has one advanced parameter
    # (an expander with an icon is listed as a status block by AppTest)
    blocks = list(launch.expander) + list(launch.status)
    assert [e for e in blocks if e.label == "1 advanced parameter"]

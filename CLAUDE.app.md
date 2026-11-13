# UmetaFlow

UmetaFlow's own guidance, imported by the template-owned `CLAUDE.md` (which describes the framework).

## What This Is

**UmetaFlow** is an untargeted metabolomics app built on the OpenMS streamlit-template: feature detection, alignment, linking and re-quantification with TOPP tools, export for SIRIUS and GNPS, and MS2 annotation (spectral library matching, SIRIUS/CSI:FingerID/CANOPUS, MS2Query). Everything below `src/workflow/`, `src/common/`, the Dockerfiles, `k8s/base`, `k8s/components` and the CI workflows is the template's and arrives through the weekly template sync; `sync.app_owned` in `app.yaml` lists what UmetaFlow owns.

## UmetaFlow's own code

- `src/UmetaFlowTOPPWorkflow.py`: the workflow (slug `umetaflow`), with a simple and an expert parameter mode (`expert_mode` in its `params.json`). It finds SIRIUS on `PATH` (`_get_sirius_path`) and caches the path in `<workspace>/runtime-config.json`.
- Pages: `content/umetaflow_*.py` (configure, run, results), `content/extracted_ion_chromatograms.py`, `content/statistics.py`, `content/mz_calculator.py`, plus the template's file upload and raw data viewer.
- Results and plots: `src/metabolomicsresults.py`, `src/dataframes.py`, `src/eic.py`, `src/stats.py`, `src/feature_explorer.py`, `src/ms2_viewer.py`, `src/plotting/`, `src/view.py`, `src/fileupload.py`, `src/masscalculator.py`.
- Python tools: `src/python-tools/` (export, merge, annotation and the MS2Query run).

## Dependencies

- `requirements.txt` lists direct dependencies only. `ms2query==1.5.3` pins `matchms<=0.26.4`, which needs `numpy<1.25`, so `numpy==1.24.4` and `pyopenms==3.2.0` (the last pyOpenMS that allows it) stay. PyPI has no linux/aarch64 wheel of pyopenms 3.2.0, so `deploy.images` in `app.yaml` leaves out the simple ARM image; the full ARM image builds pyOpenMS from source. The full Docker image builds pyOpenMS 3.6 from the OpenMS source; its wheel asks for numpy>=2, but `requirements.txt` installs numpy 1.24.4 over it (pip warns), and pyOpenMS 3.6 works with it (nanobind does not use NumPy's C API). The image's import check covers this.
- `torch` comes in through `ms2deepscore` (MS2Query). From torch 2.6 on, `torch.load` refuses the numpy arrays in the MS2Deepscore model file; `src/python-tools/run_ms2query.py` allows exactly those numpy types (`torch.serialization.add_safe_globals`) instead of turning the check off.
- `conda-requirements.txt` holds SIRIUS (`sirius-ms`, conda-forge), marked `# [x86_64]` because SIRIUS 5 has no ARM build: the ARM images ship without it. The template's Dockerfiles install it into `streamlit-env`; the full image's `mono` runs ThermoRawFileParser.

## Deployment

`k8s/overlays/prod` deploys to `umetaflow.webapps.openms.de` / `.org` (name prefix `umetaflow-`, image `ghcr.io/openms/umetaflow`, release tags `vX.Y.Z-full`). Its Argo CD Application and ImageUpdater are cluster configuration, kept outside this repository. App values (hosts, Redis URL, share size, resources) are patches in the overlay, never edits in `k8s/base` or `k8s/components`.

## Conventions & gotchas

- Tests: `python -m pytest test_gui.py tests/` as in the template; `test_gui.py` launches UmetaFlow's pages. `test.py` (run by `workflow-tests.yml`) is a placeholder.
- The workflow's presets live in `presets.json` under `umetaflow`; renaming the workflow orphans them.
- Template syncs: merge those PRs with a merge commit, never squash.

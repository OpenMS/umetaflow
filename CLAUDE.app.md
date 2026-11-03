# UmetaFlow

UmetaFlow (untargeted metabolomics) is built on the OpenMS streamlit-template. Its own code is `src/UmetaFlowTOPPWorkflow.py` (the workflow, slug `umetaflow`), the `content/umetaflow_*.py` pages, `src/python-tools/`, `src/plotting/`, `src/dataframes.py`, `src/eic.py`, `src/stats.py` and `src/metabolomicsresults.py`. UmetaFlow does not deploy to Kubernetes (`deploy.k8s.enabled: false` in `app.yaml`), so the template's `k8s/` is not part of this repository. The weekly `template-sync.yml` merges the template; merge those PRs with a merge commit, never squash.

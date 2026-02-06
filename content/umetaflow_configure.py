import streamlit as st
from src.common.common import page_setup
from src.UmetaFlowTOPPWorkflow import Workflow

# The rest of the page can, but does not have to be changed
params = page_setup()

wf = Workflow(st.session_state["workspace"])

# Expert mode from params.json
expert_mode = wf.params.get("expert_mode", False)

def update_expert_mode():
    wf.params["expert_mode"] = st.session_state["umetaflow-expert-mode"]
    wf.parameter_manager.save_parameters_to_json()

st.toggle(
    "**Expert Mode**",
    expert_mode,
    key="umetaflow-expert-mode",
    on_change=update_expert_mode
)

wf.show_parameter_section()

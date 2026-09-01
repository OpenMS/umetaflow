"""Metabolite -> KEGG compound -> KEGG pathway lookups via the KEGG REST API.

Unlike openms_insight's `calculate_kegg_enrichment` (which maps UniProt/gene IDs to
KEGG pathways via MyGene.info), this module resolves *compound names* to KEGG
compound IDs directly against KEGG's own REST API (https://rest.kegg.jp), which is
the right lookup path for metabolomics feature annotations (from SIRIUS/CSI:FingerID,
GNPS or in-house library matches) rather than genes/proteins.
"""
from collections import defaultdict

import pandas as pd
import requests
import streamlit as st

KEGG_BASE = "https://rest.kegg.jp"


@st.cache_data(show_spinner=False, ttl=60 * 60 * 24)
def find_kegg_compound(name: str) -> str | None:
    """Best-effort match of a compound name to a single KEGG compound ID (e.g. 'C00031')."""
    name = (name or "").strip()
    if not name:
        return None
    try:
        resp = requests.get(f"{KEGG_BASE}/find/compound/{name}", timeout=10)
        resp.raise_for_status()
    except requests.RequestException:
        return None

    lines = [line for line in resp.text.splitlines() if line.strip()]
    if not lines:
        return None

    # Prefer an exact (case-insensitive) name match among the hits; otherwise take the first hit.
    name_lower = name.lower()
    best = None
    for line in lines:
        kegg_id, _, names_field = line.partition("\t")
        candidate_names = [n.strip().lower() for n in names_field.split(";")]
        if name_lower in candidate_names:
            best = kegg_id
            break
    if best is None:
        best = lines[0].split("\t")[0]

    return best.replace("cpd:", "")


@st.cache_data(show_spinner=False, ttl=60 * 60 * 24)
def get_pathways_for_compound(compound_id: str) -> list[tuple[str, str]]:
    """Return [(pathway_id, pathway_name)] for a KEGG compound ID."""
    try:
        resp = requests.get(f"{KEGG_BASE}/link/pathway/cpd:{compound_id}", timeout=10)
        resp.raise_for_status()
    except requests.RequestException:
        return []

    pathway_ids = []
    for line in resp.text.splitlines():
        if not line.strip():
            continue
        _, _, pathway_ref = line.partition("\t")
        pathway_ids.append(pathway_ref.strip().replace("path:", ""))

    if not pathway_ids:
        return []

    names = _pathway_names(tuple(sorted(set(pathway_ids))))
    return [(pid, names.get(pid, pid)) for pid in pathway_ids]


@st.cache_data(show_spinner=False, ttl=60 * 60 * 24)
def _pathway_names(pathway_ids: tuple[str, ...]) -> dict[str, str]:
    try:
        resp = requests.get(f"{KEGG_BASE}/list/pathway", timeout=15)
        resp.raise_for_status()
    except requests.RequestException:
        return {}

    all_names = {}
    for line in resp.text.splitlines():
        if not line.strip():
            continue
        pid, _, name = line.partition("\t")
        all_names[pid.replace("path:", "")] = name.strip()
    return {pid: all_names.get(pid, pid) for pid in pathway_ids}


def map_names_to_pathways(names: list[str], progress_label: str = "") -> pd.DataFrame:
    """Resolve a list of compound names to KEGG compound IDs + pathway memberships.

    Returns a long-format DataFrame with columns: name, kegg_compound_id, pathway_id,
    pathway_name. Names that don't resolve to a KEGG compound are dropped.
    """
    rows = []
    progress = st.progress(0.0, text=progress_label) if progress_label else None
    unique_names = sorted({n for n in names if n and str(n).strip()})

    for i, name in enumerate(unique_names):
        if progress:
            progress.progress((i + 1) / max(len(unique_names), 1), text=f"{progress_label}: {name}")
        compound_id = find_kegg_compound(name)
        if not compound_id:
            continue
        for pathway_id, pathway_name in get_pathways_for_compound(compound_id):
            rows.append(
                {
                    "name": name,
                    "kegg_compound_id": compound_id,
                    "pathway_id": pathway_id,
                    "pathway_name": pathway_name,
                }
            )

    if progress:
        progress.empty()

    return pd.DataFrame(rows, columns=["name", "kegg_compound_id", "pathway_id", "pathway_name"])


def run_pathway_enrichment(
    mapping_df: pd.DataFrame, foreground_names: set[str], background_names: set[str]
) -> pd.DataFrame:
    """Fisher's exact test per KEGG pathway: is it over-represented in the foreground set?

    Args:
        mapping_df: Output of `map_names_to_pathways`.
        foreground_names: Significant compound names (e.g. passing p-adj/log2FC cutoffs).
        background_names: All annotated compound names considered (superset of foreground).

    Returns:
        DataFrame sorted by p-value ascending, columns:
        pathway_id, pathway_name, fg_count, bg_count, p_value.
    """
    from scipy.stats import fisher_exact

    pathway_to_fg = defaultdict(set)
    pathway_to_bg = defaultdict(set)

    for row in mapping_df.itertuples(index=False):
        if row.name in background_names:
            pathway_to_bg[row.pathway_id].add(row.name)
        if row.name in foreground_names:
            pathway_to_fg[row.pathway_id].add(row.name)

    total_fg = len(foreground_names)
    total_bg = len(background_names)
    pathway_names = mapping_df.drop_duplicates("pathway_id").set_index("pathway_id")["pathway_name"]

    results = []
    for pathway_id, fg_set in pathway_to_fg.items():
        a = len(fg_set)
        b = total_fg - a
        c = len(pathway_to_bg[pathway_id]) - a
        d = total_bg - total_fg - c
        if a == 0:
            continue
        _, p_value = fisher_exact([[a, b], [max(c, 0), max(d, 0)]], alternative="greater")
        results.append(
            {
                "pathway_id": pathway_id,
                "pathway_name": pathway_names.get(pathway_id, pathway_id),
                "fg_count": a,
                "bg_count": len(pathway_to_bg[pathway_id]),
                "p_value": p_value,
            }
        )

    return pd.DataFrame(results).sort_values("p_value") if results else pd.DataFrame(
        columns=["pathway_id", "pathway_name", "fg_count", "bg_count", "p_value"]
    )

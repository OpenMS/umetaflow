"""Metabolite -> KEGG compound -> KEGG pathway lookups.

Unlike openms_insight's `calculate_kegg_enrichment` (which maps UniProt/gene IDs to
KEGG pathways via MyGene.info), this module resolves metabolomics feature annotations
(from SIRIUS/CSI:FingerID, GNPS, MS2Query or in-house library matches) to KEGG compound
IDs, preferring a structure-based match (SMILES -> PubChem CID -> KEGG cross-reference)
over a name-based one (KEGG's own REST API, https://rest.kegg.jp), since compound names
from spectral-library search tools are often systematic/IUPAC names that don't exist
verbatim in KEGG, while the underlying structure usually does.
"""
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
import pubchempy as pcp
import requests
import streamlit as st

KEGG_BASE = "https://rest.kegg.jp"
PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


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
def find_kegg_compound_by_smiles(smiles: str) -> str | None:
    """Structure-based match: SMILES -> PubChem CID (exact structure) -> KEGG compound ID.

    More reliable than name search for compounds with systematic/IUPAC names (common in
    MS2Query/GNPS analog-search output), since it doesn't depend on that exact string
    existing as a KEGG synonym - only on the structure being registered in both PubChem
    and KEGG.
    """
    smiles = (smiles or "").strip()
    if not smiles:
        return None

    try:
        compounds = pcp.get_compounds(smiles, "smiles")
    except Exception:
        return None
    if not compounds:
        return None
    cid = compounds[0].cid
    if not cid:
        return None

    # PubChem tracks external-database registry IDs (incl. KEGG compound IDs) as xrefs.
    # Retry a couple of times: PubChem's PUG-REST occasionally returns a transient 5xx.
    for attempt in range(3):
        try:
            resp = requests.get(
                f"{PUBCHEM_BASE}/compound/cid/{cid}/xrefs/RegistryID/JSON", timeout=15
            )
            resp.raise_for_status()
            break
        except requests.RequestException:
            if attempt == 2:
                return None
            time.sleep(1.5 * (attempt + 1))
    else:
        return None

    try:
        registry_ids = resp.json()["InformationList"]["Information"][0].get("RegistryID", [])
    except (KeyError, IndexError, ValueError):
        return None

    # KEGG compound IDs look like "C" + 5 digits (e.g. "C06672"); other registries use
    # different formats, so this pattern is enough to pick them out of the mixed xref list.
    kegg_ids = [r for r in registry_ids if len(r) == 6 and r[0] == "C" and r[1:].isdigit()]
    return kegg_ids[0] if kegg_ids else None


def resolve_compound_id(name: str, smiles: str | None) -> tuple[str | None, str]:
    """Resolve a KEGG compound ID, preferring structure (SMILES) over name.

    Returns (compound_id, method) where method is "smiles", "name", or "none" (no match
    by either route) - surfaced in the mapping table so match quality stays visible.
    """
    if smiles:
        compound_id = find_kegg_compound_by_smiles(smiles)
        if compound_id:
            return compound_id, "smiles"
    compound_id = find_kegg_compound(name)
    if compound_id:
        return compound_id, "name"
    return None, "none"


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


def map_compounds_to_pathways(
    entries: list[tuple[str, str | None]], progress_label: str = "", max_workers: int = 4
) -> pd.DataFrame:
    """Resolve (name, smiles) pairs to KEGG compound IDs + pathway memberships.

    For each entry, a structure-based match (SMILES -> PubChem -> KEGG) is tried first;
    if that fails (no SMILES, or no match), falls back to a name-based KEGG search.

    The compound-resolution step (one PubChem + KEGG round trip per entry, with retries
    on failures) is the slow part and is I/O-bound, so it's run concurrently across
    `max_workers` threads - a large background set resolved one entry at a time in the
    Streamlit script thread can run for hours with no UI feedback, which looks to the
    browser like the connection dropped. `max_workers=4` keeps this under PubChem's
    request-rate guidance (<=5 req/s) even with retries in flight.

    Returns a long-format DataFrame with columns: name, kegg_compound_id, match_method
    ("smiles" or "name"), pathway_id, pathway_name. Entries that don't resolve by either
    route are dropped.
    """
    # De-duplicate by name, keeping the first non-empty SMILES seen for that name.
    unique_entries: dict[str, str | None] = {}
    for name, smiles in entries:
        if not name or not str(name).strip():
            continue
        if name not in unique_entries or (not unique_entries[name] and smiles):
            unique_entries[name] = smiles
    sorted_names = sorted(unique_entries)

    progress = st.progress(0.0, text=progress_label) if progress_label else None
    resolved: dict[str, tuple[str | None, str]] = {}
    completed = 0
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        future_to_name = {
            pool.submit(resolve_compound_id, name, unique_entries[name]): name
            for name in sorted_names
        }
        for future in as_completed(future_to_name):
            name = future_to_name[future]
            try:
                resolved[name] = future.result()
            except Exception:
                resolved[name] = (None, "none")
            completed += 1
            if progress:
                progress.progress(
                    completed / max(len(sorted_names), 1), text=f"{progress_label}: {name}"
                )
    if progress:
        progress.empty()

    # Pathway lookups are cheap/cached (only ~19k KEGG pathways total) - run sequentially
    # in name order so the output stays deterministic.
    rows = []
    for name in sorted_names:
        compound_id, method = resolved[name]
        if not compound_id:
            continue
        for pathway_id, pathway_name in get_pathways_for_compound(compound_id):
            rows.append(
                {
                    "name": name,
                    "kegg_compound_id": compound_id,
                    "match_method": method,
                    "pathway_id": pathway_id,
                    "pathway_name": pathway_name,
                }
            )

    return pd.DataFrame(
        rows, columns=["name", "kegg_compound_id", "match_method", "pathway_id", "pathway_name"]
    )


def run_pathway_enrichment(
    mapping_df: pd.DataFrame, foreground_names: set[str], background_names: set[str]
) -> pd.DataFrame:
    """Fisher's exact test per KEGG pathway: is it over-represented in the foreground set?

    Args:
        mapping_df: Output of `map_compounds_to_pathways`.
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

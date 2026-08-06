"""
reagents.py - Functional methods for analyzing TCR antibody reagents.
Loads data dynamically from the external JSON registry tracked in config.py.
"""

import json
from typing import Dict, List, NamedTuple, Optional
from .config import REFERENCE_PATHS  # Import configuration constants

class AntibodyReagent(NamedTuple):
    clone: str
    reference: str
    protein_name: str
    primary_target: str
    all_reactive_genes: List[str]
    non_reactive_genes: List[str]


def _load_registry() -> Dict[str, List[AntibodyReagent]]:
    """Helper to lazily load and parse the JSON registry into structured NamedTuples."""
    json_path = REFERENCE_PATHS.get('antibody_registry')
    
    if not json_path or not json_path.exists():
        # Fallback safeguard or informative error
        raise FileNotFoundError(f"Antibody database missing at configured path: {json_path}")
        
    with open(json_path, 'r') as f:
        raw_data = json.load(f)
        
    # Reconstruct dictionary into strong Python types
    registry = {}
    for gene, mabs in raw_data.items():
        registry[gene] = [AntibodyReagent(**mab) for mab in mabs]
    return registry


def find_antibodies_for_gene(gene: str) -> List[AntibodyReagent]:
    """Returns all antibody clones targeting a specific alpha or beta chain gene."""
    registry = _load_registry()
    # Normalize input string (e.g., 'trbv13-2*01' -> 'TRBV13-2')
    clean_gene = gene.split('*')[0].strip().upper()
    return registry.get(clean_gene, [])


def evaluate_ab_pair(trav_gene: str, trbv_gene: str) -> dict:
    """
    Evaluates whether a given alpha/beta pair can be targeted by any known reagents.
    """
    alpha_mabs = find_antibodies_for_gene(trav_gene)
    beta_mabs = find_antibodies_for_gene(trbv_gene)
    
    if alpha_mabs and beta_mabs:
        status = "Full Detection (Clones available for both Alpha and Beta chains)"
    elif alpha_mabs:
        status = "Partial Detection (Alpha chain covered only)"
    elif beta_mabs:
        status = "Partial Detection (Beta chain covered only)"
    else:
        status = "No Direct Reagents Found in Panel"
        
    return {
        "pair": f"{trav_gene} / {trbv_gene}",
        "status": status,
        "alpha_antibodies": [mab.clone for mab in alpha_mabs],
        "beta_antibodies": [mab.clone for mab in beta_mabs]
    }
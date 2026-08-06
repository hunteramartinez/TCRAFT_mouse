"""
config.py - Configuration and constants for the package
"""

import os
from pathlib import Path

# Get package root directory
PACKAGE_ROOT = Path(__file__).parent  # TCRAFT/scripts/
REFERENCES_ROOT = PACKAGE_ROOT.parent / 'references'  # TCRAFT/references/
MOUSE_REF = REFERENCES_ROOT / 'mouse'

# Reference file paths
REFERENCE_PATHS = {
    'trav_genes': [
        MOUSE_REF / 'imgt_genes_aa' / 'IMGT_TRAV_aa.fa',
        MOUSE_REF / 'imgt_genes_aa' / 'IMGT_TRAV_extra_aa.fa',
    ],
    'trbv_genes': [
        MOUSE_REF / 'imgt_genes_aa' / 'IMGT_TRBV_aa.fa',
        MOUSE_REF / 'imgt_genes_aa' / 'IMGT_TRBV_extra_aa.fa',
    ],
    'traj_genes': [MOUSE_REF / 'imgt_genes_aa' / 'IMGT_TRAJ_aa.fa'],
    'trbj_genes': [MOUSE_REF / 'imgt_genes_aa' / 'IMGT_TRBJ_aa.fa'],
    'trbc_constant': MOUSE_REF / 'imgt_genes_aa' / 'IMGT_TRBC_LVL-Cys_aa.fa',
    'trac_constant': MOUSE_REF / 'imgt_genes_aa' / 'IMGT_TRAC_LVL-Cys_aa.fa',
    'p2a_linker': REFERENCES_ROOT / 'linkers_P2A_T2A.txt',
    'codon_table': MOUSE_REF / 'mouse_codon_usage.json',
    'antibody_registry': MOUSE_REF / 'tcr_antibodies' / 'tcr_mab_registry.json',
    'hdr_sequences': MOUSE_REF / 'mm_trac_HDR_sequences.fa',
    'neb_bbsi': REFERENCES_ROOT / 'neb_tables' / 'BbsI_ligation_efficiency.csv',
    'neb_bsai': REFERENCES_ROOT / 'neb_tables' / 'BsaI_ligation_efficiency.csv',
    'neb_bsmbi': REFERENCES_ROOT / 'neb_tables' / 'BsmBI_ligation_efficiency.csv',
}

# Chain configuration - Updated to reflect AIRR-style assembly
CHAIN_CONFIG = {
    'alpha': {
        'locus': 'TRA',
        # V_CALL provides fwr1-fwr3; JUNCTION_AA provides the CDR3; J_CALL provides fwr4
        'components': {
            'variable': 'v_call',
            'junction': 'junction_aa',
            'joining': 'j_call',
            'constant': 'c_call'
        },
        'constant_region': 'trac_constant',
    },
    'beta': {
        'locus': 'TRB',
        'components': {
            'variable': 'v_call',
            'junction': 'junction_aa',
            'joining': 'j_call',
            'constant': 'c_call'
        },
        'constant_region': 'trbc_constant',
    }
}

# Linker configuration
LINKER_CONFIG = {
    'p2a': {
        'name': 'P2A',
        'type': 'nucleotide',
        'reference_file': 'p2a_linker',
        'description': 'Porcine Teschovirus-1 2A self-cleaving peptide'
    }
}

# Codon optimization
CODON_OPTIMIZATION = {
    'organism': 'mouse',
    'reference_file': 'codon_table'
}

# Default columns in AIRR file (Standardized to match Cell Ranger/MiXCR output)
AIRR_COLUMNS = {
    'clone_id': 'clone_id',
    'locus': 'locus',
    'v_call': 'v_call',
    'd_call': 'd_call',
    'j_call': 'j_call',
    'c_call': 'c_call',
    'junction_aa': 'junction_aa', # Amino acid sequence of the junction
    'sequence': 'sequence',
}
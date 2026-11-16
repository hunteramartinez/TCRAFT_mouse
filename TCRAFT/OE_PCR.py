"""
OE-PCR gBlock splitter.

Splits a fully assembled HDR template into two gBlocks that share a
central overlapping region.  During overlap-extension PCR, HiFi polymerase
anneals the complementary overlapping strands and extends outward to
reconstitute the full-length template, which is then amplified with the
existing terminal primers.

Public API
----------
design_gblocks(seq, ...)  -> dict
write_gblocks_fasta(result, filepath, clone_id='')
"""

import math
from pathlib import Path

# ── SantaLucia 1998 nearest-neighbor parameters ──────────────────────────────
# Table 2 values; complement-symmetric pairs share the same entries.

_NN_DH = {   # kcal/mol
    'AA': -7.9,  'TT': -7.9,
    'AT': -7.2,
    'TA': -7.2,
    'CA': -8.5,  'TG': -8.5,
    'GT': -8.4,  'AC': -8.4,
    'CT': -7.8,  'AG': -7.8,
    'GA': -8.2,  'TC': -8.2,
    'CG': -10.6,
    'GC': -9.8,
    'GG': -8.0,  'CC': -8.0,
}

_NN_DS = {   # cal/(mol·K)
    'AA': -22.2, 'TT': -22.2,
    'AT': -20.4,
    'TA': -21.3,
    'CA': -22.7, 'TG': -22.7,
    'GT': -22.4, 'AC': -22.4,
    'CT': -21.0, 'AG': -21.0,
    'GA': -22.2, 'TC': -22.2,
    'CG': -27.2,
    'GC': -24.4,
    'GG': -19.9, 'CC': -19.9,
}

# Initiation corrections (applied once per terminal base pair)
_INIT_DH = {'GC': 0.1, 'AT': 2.3}   # kcal/mol
_INIT_DS = {'GC': -2.8, 'AT': 4.1}  # cal/(mol·K)

_R = 1.987  # cal/(mol·K), gas constant


def _tm_nn(seq, na_conc=0.05, oligo_conc=250e-9):
    """
    Nearest-neighbor Tm (SantaLucia 1998) with linear NaCl correction.

    Parameters
    ----------
    seq       : uppercase DNA string
    na_conc   : NaCl concentration in M (default 50 mM)
    oligo_conc: total strand concentration in M (default 250 nM)

    Returns
    -------
    Tm in °C, rounded to 1 decimal place.
    """
    seq = seq.upper()
    if len(seq) < 2:
        raise ValueError('Sequence must be at least 2 nt for Tm calculation')

    dH = 0.0  # kcal/mol
    dS = 0.0  # cal/(mol·K)

    for i in range(len(seq) - 1):
        pair = seq[i:i + 2]
        dH += _NN_DH.get(pair, 0.0)
        dS += _NN_DS.get(pair, 0.0)

    # Initiation correction at both terminal base pairs
    for end_base in (seq[0], seq[-1]):
        key = 'GC' if end_base in 'GC' else 'AT'
        dH += _INIT_DH[key]
        dS += _INIT_DS[key]

    dH_cal = dH * 1000  # kcal → cal
    tm_1m = dH_cal / (dS + _R * math.log(oligo_conc / 4.0)) - 273.15
    return round(tm_1m + 16.6 * math.log10(na_conc), 1)


def _has_gc_clamp(seq):
    """True if both the first and last base of seq are G or C."""
    return seq[0] in 'GC' and seq[-1] in 'GC'


def find_oe_split(seq, target_tm=65.0, search_window=50,
                  min_overlap=15, max_overlap=35):
    """
    Find the optimal overlap position and length near the sequence midpoint.

    Scans a ±search_window bp window around the midpoint, trying overlap
    lengths from min_overlap to max_overlap.  Requires a GC clamp at both
    ends of the overlap (since HiFi extension initiates from both ends of
    the junction).  Returns the candidate whose Tm is closest to target_tm.

    Parameters
    ----------
    seq           : uppercase DNA string (fully assembled HDR template)
    target_tm     : desired overlap Tm in °C (default 65.0)
    search_window : bp to search either side of midpoint (default 50)
    min_overlap   : minimum overlap length in nt (default 15)
    max_overlap   : maximum overlap length in nt (default 35)

    Returns
    -------
    dict with keys: overlap_start, overlap_end, overlap_seq, overlap_tm,
                    overlap_len
    """
    seq = seq.upper()
    midpoint = len(seq) // 2

    best = None
    best_score = float('inf')

    for center in range(midpoint - search_window, midpoint + search_window + 1):
        for length in range(min_overlap, max_overlap + 1):
            start = center - length // 2
            end   = start + length
            if start < 0 or end > len(seq):
                continue
            overlap = seq[start:end]
            if not _has_gc_clamp(overlap):
                continue
            tm = _tm_nn(overlap)
            score = abs(tm - target_tm)
            if score < best_score:
                best_score = score
                best = {
                    'overlap_start': start,
                    'overlap_end':   end,
                    'overlap_seq':   overlap,
                    'overlap_tm':    tm,
                    'overlap_len':   length,
                }

    if best is None:
        raise ValueError(
            f'No GC-clamped overlap found within ±{search_window} bp of midpoint '
            f'with length {min_overlap}–{max_overlap} nt. '
            'Try widening search_window or relaxing min/max_overlap.'
        )

    return best


def design_gblocks(seq, target_tm=65.0, search_window=50,
                   min_overlap=15, max_overlap=35):
    """
    Split a fully assembled HDR template into two overlapping gBlocks.

    gBlock 1 covers seq[0 : overlap_end]   (ends with the overlap)
    gBlock 2 covers seq[overlap_start : ]  (starts with the overlap)

    The two gBlocks are ordered from IDT and combined in an OE-PCR reaction:
    HiFi polymerase anneals the shared overlap and extends outward to give the
    full-length template, which is then amplified with existing terminal primers.

    Parameters
    ----------
    seq           : fully assembled HDR sequence (any case; uppercased internally)
    target_tm     : target overlap Tm in °C (default 65.0)
    search_window : bp window either side of midpoint to search (default 50)
    min_overlap   : minimum overlap length in nt (default 15)
    max_overlap   : maximum overlap length in nt (default 35)

    Returns
    -------
    dict with keys:
        gblock1, gblock2     : uppercase DNA strings
        overlap_seq          : shared overlap sequence
        overlap_tm           : nearest-neighbor Tm of the overlap (°C)
        overlap_len          : overlap length in nt
        overlap_start        : 0-based start of overlap in full sequence
        overlap_end          : 0-based end (exclusive) of overlap in full sequence
        total_len            : length of the full assembled sequence
    """
    seq = seq.upper()
    split = find_oe_split(seq, target_tm=target_tm, search_window=search_window,
                          min_overlap=min_overlap, max_overlap=max_overlap)

    overlap_start = split['overlap_start']
    overlap_end   = split['overlap_end']

    return {
        'gblock1':       seq[:overlap_end],
        'gblock2':       seq[overlap_start:],
        'overlap_seq':   split['overlap_seq'],
        'overlap_tm':    split['overlap_tm'],
        'overlap_len':   split['overlap_len'],
        'overlap_start': overlap_start,
        'overlap_end':   overlap_end,
        'total_len':     len(seq),
    }


def write_gblocks_tsv(result, filepath, clone_id=''):
    """
    Write a two-row TSV (Name, Sequence) for the two gBlocks.

    Mirrors the column structure of constant_gene_blocks.csv / variable_gene_blocks.csv
    so the sequences can be copied directly into an IDT order form.

    Parameters
    ----------
    result   : dict returned by design_gblocks()
    filepath : path to write (created or overwritten)
    clone_id : optional label used in the Name column
    """
    suffix = f'_{clone_id}' if clone_id else ''
    filepath = Path(filepath)
    with open(filepath, 'w') as fh:
        fh.write(f'gBlock1{suffix}\t{result["gblock1"]}\n')
        fh.write(f'gBlock2{suffix}\t{result["gblock2"]}\n')


def write_gblocks_fasta(result, filepath, clone_id=''):
    """
    Write a two-entry FASTA file for the two gBlocks.

    FASTA headers include gBlock length, overlap length, and overlap Tm so
    the ordering sheet is self-documenting.

    Parameters
    ----------
    result   : dict returned by design_gblocks()
    filepath : path to write (created or overwritten)
    clone_id : optional label appended to FASTA entry names
    """
    suffix = f'_{clone_id}' if clone_id else ''
    g1  = result['gblock1']
    g2  = result['gblock2']
    ol  = result['overlap_len']
    tm  = result['overlap_tm']

    filepath = Path(filepath)
    with open(filepath, 'w') as fh:
        fh.write(f'>gBlock1{suffix}  len={len(g1)}  overlap_nt={ol}  overlap_tm={tm:.1f}C\n')
        fh.write(g1 + '\n')
        fh.write(f'>gBlock2{suffix}  len={len(g2)}  overlap_nt={ol}  overlap_tm={tm:.1f}C\n')
        fh.write(g2 + '\n')

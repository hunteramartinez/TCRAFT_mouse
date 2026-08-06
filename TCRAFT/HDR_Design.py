#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
HDR_Design.py — Mouse TCR HDR template designer.

Designs all DNA fragments needed for Cas9-RNP mediated knock-in of a paired
TCR alpha/beta into the mouse TRAC locus using two-step golden gate assembly.

Assembly architecture (5' → 3' in the final HDR template):
  [LH] - [TRBV] - [CDR3β] - [TRBC-LVL-Cys] - [T2A] - [TRAV] - [CDR3α] - [TRAC-LVL-Cys] - [RH]

Fragment types:
  Constant (ordered once):
    • 24 LH-TRBV gene blocks  (one per functional TRBV gene)
    • 1  TRAC-RH gene block
    • 1  TRBC-T2A gene block

  Variable (ordered per TCR):
    • 1  CDR3β-TRAV-CDR3α gene block

Golden gate scheme:
  Step 1 — BbsI digest + ligation:
    [LH-TRBV]──1a──[CDR3β·TRAV·CDR3α]──1b──[TRAC-RH]
    • 1a = TRBV-specific 4nt overhang (unique per gene, enables TRBV-pooled assembly)
    • 1b = constant 4nt from the start of TRAC-LVL-Cys

  Step 2 — BsmBI digest + ligation:
    assembled-step1 ──2a──[TRBC-T2A]──2b──
    • 2a = constant 4nt from the start of TRBC-LVL-Cys
    • 2b = constant 4nt from the end of T2A

Usage:
  CLI:   TCRAFT-hdr input_airr.tsv --trbv TRBV5 --output_dir ./hdr_output
  API:   from TCRAFT.HDR_Design import generate_hdr_design
"""

import argparse
import json
import os
import random
import sys
from pathlib import Path

import importlib.resources as pkg_resources
import numpy as np
import pandas as pd
from Bio.Restriction import BbsI, BsaI, BsmBI, SapI, XbaI, EcoRI, XhoI, EcoRV, NotI, Analysis
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from tqdm import tqdm


# ── Reference paths ─────────────────────────────────────────────────────────

_REF_ROOT = Path(__file__).parent / 'references'
_MOUSE    = _REF_ROOT / 'mouse'

_PATHS = {
    'trbv_aa':   [_MOUSE / 'imgt_genes_aa' / 'IMGT_TRBV_aa.fa',
                  _MOUSE / 'imgt_genes_aa' / 'IMGT_TRBV_extra_aa.fa'],
    'trav_aa':   [_MOUSE / 'imgt_genes_aa' / 'IMGT_TRAV_aa.fa',
                  _MOUSE / 'imgt_genes_aa' / 'IMGT_TRAV_extra_aa.fa'],
    'trbj_aa':   [_MOUSE / 'imgt_genes_aa' / 'IMGT_TRBJ_aa.fa'],
    'traj_aa':   [_MOUSE / 'imgt_genes_aa' / 'IMGT_TRAJ_aa.fa'],
    'trac_lvlcys': _MOUSE / 'imgt_genes_aa' / 'IMGT_TRAC_LVL-Cys_aa.fa',
    'trbc_lvlcys': _MOUSE / 'imgt_genes_aa' / 'IMGT_TRBC_LVL-Cys_aa.fa',
    'hdr_seqs':  _MOUSE / 'mm_trac_HDR_sequences.fa',
    'codon_table': _MOUSE / 'mouse_codon_usage.json',
    'linkers':   _REF_ROOT / 'linkers_P2A_T2A.txt',
    'neb_bbsi':  _REF_ROOT / 'neb_tables' / 'BbsI_ligation_efficiency.csv',
    'neb_bsmbi': _REF_ROOT / 'neb_tables' / 'BsmBI_ligation_efficiency.csv',
}

# Restriction sites to screen out of all designed sequences
_FORBIDDEN_RE_SITES = [BbsI.site, BsmBI.site, BsaI.site, SapI.site,
                        XbaI.site, EcoRI.site, XhoI.site, EcoRV.site, NotI.site]

# BbsI golden gate primers: GAAGAC + 2nt spacer | 2nt spacer + GTCTTC
_BBSI_FWD = 'GAAGACnn'   # appended to 3' end of fragment; nn = spacer nucleotides
_BBSI_REV = 'nnGTCTTC'   # prepended to 5' end of fragment (cuts to expose 5' OH)

# BsmBI golden gate: CGTCTC + 1nt spacer | 1nt spacer + GAGACG
_BSMBI_FWD = 'CGTCTCn'
_BSMBI_REV = 'nGAGACG'


# ── FASTA / file parsers ─────────────────────────────────────────────────────

def _parse_fasta(path):
    """Return list of (header, sequence) from a FASTA file."""
    records, header, parts = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(parts)))
                header, parts = line[1:], []
            else:
                parts.append(line.replace(' ', '').replace('-', ''))
    if header is not None:
        records.append((header, ''.join(parts)))
    return records


def _load_fasta_dict(paths, filter_functional=True, preferred_allele='*01'):
    """
    Load FASTA files into {gene_name: aa_seq} dict.
    When filter_functional=True, only '|F|' entries are kept.
    For genes with multiple alleles, *01 is preferred; otherwise the first entry is used.
    """
    gene_dict = {}
    for path in (paths if isinstance(paths, list) else [paths]):
        for header, seq in _parse_fasta(path):
            if filter_functional and '|F|' not in header:
                continue
            parts = header.split('|')
            if len(parts) < 2:
                continue
            allele = parts[1].strip()          # e.g. TRBV5*01
            gene = allele.split('*')[0]         # e.g. TRBV5
            # Prefer *01 allele; only overwrite if this is *01 or gene not yet seen
            if gene not in gene_dict or allele.endswith(preferred_allele):
                gene_dict[gene] = seq.upper()
    return gene_dict


def _load_hdr_sequences(path):
    """Return {entry_name: sequence} from mm_trac_HDR_sequences.fa."""
    seqs = {}
    for header, seq in _parse_fasta(path):
        seqs[header.strip()] = seq
    return seqs


def _load_linkers(path):
    """Return {linker_name: nt_sequence} from linkers_P2A_T2A.txt."""
    linkers = {}
    for header, seq in _parse_fasta(path):
        linkers[header.strip()] = seq.upper()
    return linkers


def _load_codon_table(path):
    """Load mouse codon usage JSON: {AA: [[codon, freq], ...]}."""
    with open(path) as fh:
        return json.load(fh)


def _load_neb_table(path):
    """
    Load NEB ligation efficiency table into a DataFrame indexed by overhang.

    Table convention (Pryor et al. 2020):
      row = overhang X on one fragment's sticky end
      col = overhang Y on the adjacent fragment's sticky end
      value = observed ligation count

    Self-ligation (correct):  df.loc[X, RC(X)]
    Cross-ligation (wrong):   df.loc[X, RC(Y)] where Y is another OH in the same reaction
    Physical misligation risk: Hamming(X, Y) — same as Hamming(RC(X), RC(Y)) — predicts
      cross-ligation.  Hamming = 1 → ~7-13% misligation.  Hamming ≥ 2 → ~0%.

    encoding='utf-8-sig' strips the BOM that is present in the NEB CSV files.
    """
    df = pd.read_csv(path, index_col=0, encoding='utf-8-sig')
    df.index   = df.index.str.strip().str.upper()
    df.columns = df.columns.str.strip().str.upper()
    return df


def _oh_self_score(oh, neb_df):
    """Self-ligation efficiency for a 4nt overhang: df[oh, RC(oh)]."""
    oh = oh.upper()
    partner = rev_comp(oh).upper()
    try:
        return int(neb_df.loc[oh, partner])
    except (KeyError, ValueError):
        return 0


# ── Reference data (module-level, loaded once) ──────────────────────────────

_codon_table  = _load_codon_table(_PATHS['codon_table'])
_hdr_seqs     = _load_hdr_sequences(_PATHS['hdr_seqs'])
_linkers      = _load_linkers(_PATHS['linkers'])

# P2A leader placed between the left homology arm and the TRBV gene in every
# LH-TRBV gene block.  It appears to contain an in-frame stop (TGA) when read
# in isolation, but the reading frame is established by the upstream endogenous
# TRAC ORF carried across the left homology arm — the LH arm positions the P2A
# so that TGA falls out-of-frame relative to the genomic ORF.  Do NOT alter
# or codon-optimise this sequence; it is a fixed, frame-verified design element.
_p2a_leader_nt = _hdr_seqs.get('mm_trac_p2a_leader_sequence', '').upper()
_trbv_aa_dict = _load_fasta_dict(_PATHS['trbv_aa'])
_trav_aa_dict = _load_fasta_dict(_PATHS['trav_aa'])
_trbj_aa_dict = _load_fasta_dict(_PATHS['trbj_aa'])
_traj_aa_dict = _load_fasta_dict(_PATHS['traj_aa'])

# LVL-Cys constant regions (strip leading X = J-C junction residue)
_trac_lvlcys_aa = _parse_fasta(_PATHS['trac_lvlcys'])[0][1].lstrip('X').upper()
_trbc_lvlcys_aa = _parse_fasta(_PATHS['trbc_lvlcys'])[0][1].lstrip('X').upper()

# Homology arms (use uppercase portion as the "interface" sequence)
_lh_nt = _hdr_seqs.get('mm_trac_left_homology_arm', '')
_rh_nt = _hdr_seqs.get('mm_trac_right_homology_arm', '')

# T2A linker NT sequence
_t2a_nt = _linkers.get('furin-GSG-T2A', _linkers.get('T2A', ''))
if not _t2a_nt:
    # Fall back to P2A if T2A not found under that name
    _t2a_nt = next((v for k, v in _linkers.items() if 'T2A' in k.upper()), '')

_neb_bbsi  = _load_neb_table(_PATHS['neb_bbsi'])
_neb_bsmbi = _load_neb_table(_PATHS['neb_bsmbi'])


# ── Sequence utilities ───────────────────────────────────────────────────────

_COMPLEMENT = str.maketrans('ACGTacgt', 'TGCAtgca')

def rev_comp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def translate(nt_seq):
    codon_table = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }
    aa = []
    for i in range(0, len(nt_seq) - len(nt_seq) % 3, 3):
        aa.append(codon_table.get(nt_seq[i:i+3].upper(), 'X'))
    return ''.join(aa)


# ── Codon optimization ───────────────────────────────────────────────────────

def _best_codon(aa):
    """Most-frequent codon for an amino acid (deterministic)."""
    options = _codon_table.get(aa, [])
    if not options:
        raise ValueError(f'No codon entry for amino acid: {aa!r}')
    return max(options, key=lambda x: x[1])[0]


def back_translate_deterministic(aa_seq):
    """Back-translate using the single most-frequent codon per amino acid."""
    return ''.join(_best_codon(aa) for aa in aa_seq)


def _sample_codon(aa, freq_cutoff=0.15):
    options = [(c, f) for c, f in _codon_table.get(aa, []) if f > freq_cutoff]
    if not options:
        options = _codon_table.get(aa, [])
    codons, freqs = zip(*options)
    freqs = np.array(freqs, dtype=float)
    freqs /= freqs.sum()
    return random.choices(codons, weights=freqs)[0]


def back_translate_random(aa_seq, freq_cutoff=0.15):
    """Back-translate with weighted-random codon sampling."""
    seq = ''.join(_sample_codon(aa, freq_cutoff) for aa in aa_seq)
    assert Seq(seq).translate() == aa_seq
    return seq


# ── Restriction site / homopolymer screening ────────────────────────────────

def check_re_sites(seq, sites=None):
    """Return (position, site) of first forbidden RE site found, or (-1, None)."""
    if sites is None:
        sites = _FORBIDDEN_RE_SITES
    seq_u = seq.upper()
    for site in sites:
        s = site.upper()
        rc_s = rev_comp(s)
        for needle in (s, rc_s):
            idx = seq_u.find(needle)
            if idx != -1:
                return idx, site
    return -1, None


def check_homopolymer(seq, min_len=5):
    """Return (position, run) of first homopolymer >= min_len, or (-1, None)."""
    for nt in 'ACGT':
        run = nt * min_len
        idx = seq.upper().find(run)
        if idx != -1:
            return idx, run
    return -1, None


def _remove_re_sites(seq, aa, sites=None):
    if sites is None:
        sites = _FORBIDDEN_RE_SITES
    iters = 0
    while check_re_sites(seq, sites)[0] != -1:
        pos, site = check_re_sites(seq, sites)
        start = max(0, pos // 3 - 2)
        end   = min((pos + len(site)) // 3 + 2, len(aa) - 1)
        replace = None
        for opt in _dfs_seq_builder(aa[start:end+1]):
            if check_re_sites(opt, sites)[0] == -1:
                replace = opt
                break
        if replace is None or iters >= 20:
            raise ValueError(f'Could not remove RE site in: {seq}')
        seq = seq[:start*3] + replace + seq[(end+1)*3:]
        iters += 1
    return seq


def _remove_homopolymers(seq, aa, sites=None):
    if sites is None:
        sites = _FORBIDDEN_RE_SITES
    iters = 0
    while check_homopolymer(seq)[0] != -1:
        pos, run = check_homopolymer(seq)
        start = max(0, pos // 3 - 2)
        end   = min((pos + len(run)) // 3 + 2, len(aa) - 1)
        replace = None
        for opt in _dfs_seq_builder(aa[start:end+1]):
            if check_re_sites(opt, sites)[0] == -1 and check_homopolymer(opt)[0] == -1:
                replace = opt
                break
        if replace is None or iters >= 20:
            raise ValueError(f'Could not remove homopolymer in: {seq}')
        seq = seq[:start*3] + replace + seq[(end+1)*3:]
        iters += 1
    return seq


def _dfs_seq_builder(aa_seq):
    """Yield all NT sequences that translate to aa_seq (uses all codons, no freq cutoff)."""
    def _rec(pos, current):
        if pos >= len(aa_seq):
            yield current
            return
        for codon, _ in _codon_table.get(aa_seq[pos], []):
            yield from _rec(pos + 1, current + codon)
    yield from _rec(0, '')


def clean_nt_sequence(nt_seq, aa_seq, label='sequence'):
    """Remove RE sites and homopolymers from a codon-optimized sequence."""
    try:
        nt_seq = _remove_re_sites(nt_seq, aa_seq)
        nt_seq = _remove_homopolymers(nt_seq, aa_seq)
        assert Seq(nt_seq).translate() == aa_seq, 'Translation mismatch after cleaning'
    except Exception as e:
        raise ValueError(f'Failed to clean {label}: {e}')
    return nt_seq


def generate_clean_nt(aa_seq, label='', max_tries=30, freq_cutoff=0.15):
    """
    Generate a codon-optimized, restriction-site-free, homopolymer-free NT sequence.
    Uses random sampling; retries up to max_tries times.
    """
    for attempt in range(max_tries):
        try:
            seq = back_translate_random(aa_seq, freq_cutoff)
            seq = clean_nt_sequence(seq, aa_seq, label)
            assert check_re_sites(seq)[0] == -1
            assert check_homopolymer(seq)[0] == -1
            return seq
        except Exception:
            if attempt == max_tries - 1:
                raise ValueError(f'generate_clean_nt failed for {label!r} after {max_tries} tries')
    raise ValueError(f'generate_clean_nt failed for {label!r}')


# ── V-gene trimming at CDR3 boundary ────────────────────────────────────────

def trim_vgene_to_cys104(aa_seq, gene_name=''):
    """
    Trim a full L+V AA sequence to the conserved Cys104 (the last C of FWR3).
    Returns the trimmed AA sequence including the terminal C.

    IMGT V-gene sequences sometimes extend into the CDR3 (e.g. TRBV5 ends in CASSQ).
    We truncate at the last C that represents Cys104 — identified as the C that is
    followed by residues with CDR3-like composition before the sequence ends.

    Heuristic: find the LAST uppercase C in the sequence that is at least 10 residues
    from the N-terminus and at least 2 residues from the end (some V genes end with
    a few CDR3-templated residues after C104).
    """
    # Strip leading * or X that appear in some IMGT entries
    aa = aa_seq.lstrip('X*').rstrip('*')
    # Walk backwards to find the last C that looks like Cys104
    # Cys104 is always followed by something (CDR3 or end-of-sequence)
    # We accept any C that occurs in the latter half of the sequence and is
    # followed by 0-5 residues (V-templated CDR3 start)
    last_c_idx = None
    for i in range(len(aa) - 1, max(10, len(aa) // 2), -1):
        if aa[i] == 'C':
            last_c_idx = i
            break
    if last_c_idx is None:
        raise ValueError(f'Could not find conserved Cys104 in {gene_name}: {aa}')
    return aa[:last_c_idx + 1]   # include the C


# ── NEB overhang utilities ───────────────────────────────────────────────────

def oh_ligation_score(oh, neb_df):
    """Diagonal (self-ligation) score for a 4nt overhang from a NEB table."""
    oh = oh.upper()
    if oh in neb_df.index and oh in neb_df.columns:
        return neb_df.loc[oh, oh]
    return 0


def oh_cross_score(oh1, oh2, neb_df):
    """Off-diagonal ligation score between two overhangs (higher = more crosstalk)."""
    o1, o2 = oh1.upper(), oh2.upper()
    rc2 = rev_comp(o2).upper()
    score = 0
    for a, b in [(o1, o2), (o1, rc2), (rev_comp(o1).upper(), o2), (rev_comp(o1).upper(), rc2)]:
        if a in neb_df.index and b in neb_df.columns:
            score = max(score, neb_df.loc[a, b])
    return score


def choose_best_oh(candidates, reserved, neb_df, min_self=100, max_cross=30):
    """
    From a list of candidate 4nt OHs, return the one with highest self-ligation
    and acceptable cross-ligation relative to reserved OHs.
    """
    ranked = sorted(
        candidates,
        key=lambda oh: oh_ligation_score(oh, neb_df),
        reverse=True,
    )
    for oh in ranked:
        if oh_ligation_score(oh, neb_df) < min_self:
            continue
        if rev_comp(oh).upper() == oh.upper():
            continue  # palindrome — would self-anneal
        if any(oh_cross_score(oh, r, neb_df) > max_cross for r in reserved):
            continue
        return oh
    # Relax constraints if nothing qualifies
    return ranked[0] if ranked else candidates[0]


# ── Golden gate site builders ────────────────────────────────────────────────

_GG_TAIL = 'ACAT'  # stub outside RE sites so BioPython can detect the cut position

# Empirically verified spacers (see design notes):
#   BbsI  (GAAGAC, 2/6): 3nt spacer GGA for 3' site; 1nt spacer A for 5' site
#   BsmBI (CGTCTC, 1/5): 2nt spacer CA for both 3' and 5' sites (both use CGTCTC on top)
#
# In every case, after enzyme digest the RIGHT downstream piece starts with oh_4nt
# as a 5' top-strand overhang.  For the LEFT upstream piece, the oh appears as a
# complementary 5' bottom-strand overhang (which anneals with the partner right piece).


def _bbsi_3prime_site(oh_4nt):
    """
    Append to a fragment's 3' end.  After BbsI digest, the downstream waste piece
    starts with oh_4nt (5' top-strand overhang); the functional upstream piece has
    the complementary 5' bottom-strand overhang at its right end.

    Layout: [...body...] GAAGAC GGA [oh_4nt] ACAT
    """
    return f'GAAGACGGA{oh_4nt}{_GG_TAIL}'


def _bbsi_5prime_site(oh_4nt):
    """
    Prepend to a fragment's 5' end.  After BbsI digest, the downstream functional
    piece starts with oh_4nt as a 5' top-strand overhang.

    Layout: ACAT [oh_4nt] A GTCTTC [...body...]
    """
    return f'{_GG_TAIL}{oh_4nt}AGTCTTC'


def _bsmbi_3prime_site(oh_4nt):
    """
    Append to a fragment's 3' end (step 2 BsmBI reaction).
    After BsmBI digest, the downstream waste piece starts with oh_4nt; the
    functional upstream piece has the complementary 5' bottom-strand overhang.

    Layout: [...body...] CGTCTC CA [oh_4nt] ACAT
    """
    return f'CGTCTCCA{oh_4nt}{_GG_TAIL}'


def _bsmbi_5prime_site(oh_4nt):
    """
    Prepend to a fragment's 5' end (step 2 BsmBI reaction).
    After BsmBI digest, the downstream functional piece starts with oh_4nt.

    Layout: ACAT CGTCTC CA [oh_4nt] [...body...]
    """
    return f'{_GG_TAIL}CGTCTCCA{oh_4nt}'


# ── Overhang design ──────────────────────────────────────────────────────────

def _last_4nt(nt_seq):
    return nt_seq[-4:].upper()


def _first_4nt(nt_seq):
    return nt_seq[:4].upper()


def _hamming(a, b):
    return sum(x != y for x, y in zip(a.upper(), b.upper()))


def _oh_conflicts(oh, occupied, min_dist=2):
    """Return True if oh or RC(oh) is within min_dist of any OH in occupied."""
    rc = rev_comp(oh).upper()
    oh = oh.upper()
    for existing in occupied:
        e = existing.upper()
        erc = rev_comp(e).upper()
        for a, b in [(oh, e), (oh, erc), (rc, e), (rc, erc)]:
            if _hamming(a, b) < min_dist:
                return True
    return False


def _best_neb_oh_from_sequence(nt_seq, neb_df, scan_window=16,
                                reserved=(), min_hamming=2, min_score=50):
    """
    Scan the first `scan_window` positions of nt_seq for a 4nt window that:
      - is not palindromic
      - has NEB self-ligation score >= min_score
      - satisfies Hamming >= min_hamming from all OHs in `reserved`
    Returns (oh_4nt, position) of the highest-scoring qualifying window.
    Falls back to best available if min_score cannot be met.
    """
    candidates = []
    for k in range(min(scan_window, len(nt_seq) - 4)):
        oh = nt_seq[k:k+4].upper()
        if oh == rev_comp(oh):             # palindrome
            continue
        if any(oh == c * 4 for c in 'ACGT'):  # homopolymer
            continue
        if _oh_conflicts(oh, reserved, min_hamming):
            continue
        score = _oh_self_score(oh, neb_df)
        candidates.append((score, k, oh))

    if not candidates:
        return _first_4nt(nt_seq), 0     # last-resort fallback

    # Prefer high score; among ties, prefer position closer to 0 (biology-natural)
    best = max(candidates, key=lambda t: t[0])
    if best[0] < min_score:
        # Nothing meets the threshold — return best available with a note
        pass
    return best[2], best[1]


def _best_neb_oh_from_pool(neb_df, reserved=(), min_hamming=2, min_score=100):
    """
    Return the highest-NEB-self-ligation-scoring 4nt OH that satisfies the
    Hamming and palindrome constraints.  Falls back to lower scores if needed.
    """
    import itertools
    all4 = [''.join(x) for x in itertools.product('ACGT', repeat=4)]
    scored = sorted(
        [(_oh_self_score(oh, neb_df), oh) for oh in all4],
        reverse=True,
    )
    for score, oh in scored:
        if oh == rev_comp(oh):
            continue
        if any(oh == c * 4 for c in 'ACGT'):
            continue
        if _oh_conflicts(oh, reserved, min_hamming):
            continue
        return oh
    raise RuntimeError('_best_neb_oh_from_pool: exhausted all 4-mers')


def design_constant_overhangs(trac_nt, trbc_nt):
    """
    Select 1b, 2a, 2b overhangs using NEB ligation efficiency tables.

    Reaction constraints:
      Step 1 (BbsI): 1a OHs + 1b — only within this reaction need Hamming >= 2.
      Step 2 (BsmBI): 2a + 2b — only these two need Hamming >= 2 from each other.

    1b: scan first 16nt of TRAC-LVL-Cys, pick the window with the highest BbsI
        self-ligation score.  Position k is stored so build_trac_rh_fragment knows
        to put TRAC_nt[k+4:] in the gene-block body (TRAC_nt[:k] goes into the
        variable fragment CDR3α bridge).

    2a: same logic using TRBC-LVL-Cys and the BsmBI table.
        TRBC_nt[k+4:] goes in the TRBC-T2A gene block; TRBC_nt[:k] goes in the
        variable fragment CDR3β bridge.

    2b: chosen from the full 4-mer pool ranked by BsmBI self-ligation score,
        with Hamming >= 2 from 2a.

    Returns (oh_1b, split_1b, oh_2a, split_2a, oh_2b).
    split_X = the byte offset into the sequence where the OH window starts.
    """
    oh_1b, split_1b = _best_neb_oh_from_sequence(trac_nt, _neb_bbsi,  scan_window=16)
    oh_2a, split_2a = _best_neb_oh_from_sequence(trbc_nt, _neb_bsmbi, scan_window=16)

    # oh_2b MUST be derived from within the T2A sequence — not a synthetic linker.
    #
    # Reading-frame rationale: every other OH in this design comes from the gene
    # sequence it flanks, so after ligation the full sequence is reconstituted with
    # no extra nucleotides.  A synthetic 4-nt oh_2b would add 4 nt between T2A and
    # TRAV; since 4 mod 3 = 1, this creates a +1 frameshift into TRAV.
    #
    # By taking oh_2b = T2A[split_2b : split_2b+4], the assembled product is:
    #   TRBC + T2A[:split_2b] + T2A[split_2b:split_2b+4] + T2A[split_2b+4:] + TRAV
    #   = TRBC + T2A + TRAV  (T2A fully reconstituted, TRAV in correct frame)
    #
    # T2A[:split_2b] lives in the TRBC-T2A gene block.
    # T2A[split_2b+4:] (the T2A tail) lives in the variable gene block, immediately
    # before TRAV.
    t2a = _t2a_nt.upper()
    t2a_candidates = []
    for n in range(len(t2a) - 4):
        oh = t2a[n:n + 4]
        if (oh == rev_comp(oh)
                or any(oh == c * 4 for c in 'ACGT')
                or _oh_conflicts(oh, [oh_2a])):
            continue
        t2a_candidates.append((_oh_self_score(oh, _neb_bsmbi), n, oh))

    if not t2a_candidates:
        raise RuntimeError('No valid oh_2b found in T2A — adjust T2A or oh_2a')
    _, split_2b, oh_2b = max(t2a_candidates, key=lambda t: t[0])

    return oh_1b, split_1b, oh_2a, split_2a, oh_2b, split_2b


def design_1a_overhangs(trbv_nt_dict, reserved_ohs=(), min_hamming=2):
    """
    For each TRBV gene, select a unique 4nt 1a overhang with Hamming distance
    >= min_hamming from every other assigned OH (and its RC) in the same BbsI
    reaction.

    Returns a dict:
      {gene_name: {'oh': oh_4nt, 'split_from_end': int}}

    where split_from_end is the number of nucleotides from the END of the TRBV
    coding sequence where the 4nt window starts.  The LH-TRBV gene block body
    uses trbv_nt[:-split_from_end], and the variable fragment bridge contains
    trbv_nt[-split_from_end+4:] (the tail after the OH, before CDR3β).

    Strategy:
      Phase 1 — scan the last 30nt of the TRBV NT sequence for a qualifying window.
      Phase 2 — exhaustive balanced-GC search across all 256 4-mers if phase 1 fails.
                 This means the OH is a designed junction sequence (not gene-derived)
                 and no bridge is needed (split_from_end=4).
      Genes that cannot satisfy the constraint are flagged as separate-pool candidates.
    """
    import itertools

    assigned = list(oh.upper() for oh in reserved_ohs)

    oh_map = {}
    needs_separate_pool = []

    def _ok(oh):
        oh = oh.upper()
        return (len(oh) == 4
                and oh != rev_comp(oh)
                and not any(oh == c * 4 for c in 'ACGT')
                and not _oh_conflicts(oh, assigned, min_hamming))

    # Pre-rank all 256 4-mers by BbsI self-ligation efficiency (highest first).
    all4_by_score = sorted(
        [''.join(x) for x in itertools.product('ACGT', repeat=4)],
        key=lambda s: _oh_self_score(s, _neb_bbsi),
        reverse=True,
    )

    def _phase1_candidates(nt):
        """All (score, start, oh, split) from the last 30nt of nt satisfying _ok."""
        n = len(nt)
        out = []
        for start in range(n - 4, max(n - 31, 3), -1):
            oh = nt[start:start + 4].upper()
            if _ok(oh):
                out.append((_oh_self_score(oh, _neb_bbsi), start, oh, n - start))
        return out

    # Process genes hardest-first: genes with fewest sequence-derived candidates
    # (estimated against only the reserved OHs) are processed first so they claim
    # good OHs before easier genes exhaust the pool.
    initial_candidate_counts = {}
    for gene, nt in trbv_nt_dict.items():
        n = len(nt)
        count = 0
        for start in range(n - 4, max(n - 31, 3), -1):
            oh = nt[start:start + 4].upper()
            if (oh != rev_comp(oh)
                    and not any(oh == c * 4 for c in 'ACGT')
                    and not _oh_conflicts(oh, list(assigned), min_hamming)):
                count += 1
        initial_candidate_counts[gene] = count

    gene_order = sorted(trbv_nt_dict.keys(),
                        key=lambda g: (initial_candidate_counts[g], g))

    for gene in gene_order:
        nt = trbv_nt_dict[gene]
        found = False

        # Phase 1: collect all valid windows, pick highest NEB score.
        candidates = _phase1_candidates(nt)
        if candidates:
            best = max(candidates, key=lambda t: (t[0], -t[3]))
            _, _, oh, split = best
            oh_map[gene] = {'oh': oh, 'split_from_end': split,
                            'neb_score': _oh_self_score(oh, _neb_bbsi)}
            assigned.append(oh.upper())
            found = True

        if not found:
            # Phase 2: synthetic OH ranked by NEB score.
            for oh in all4_by_score:
                if _ok(oh):
                    oh_map[gene] = {'oh': oh, 'split_from_end': 4,
                                    'neb_score': _oh_self_score(oh, _neb_bbsi)}
                    assigned.append(oh.upper())
                    found = True
                    break

        if not found:
            needs_separate_pool.append(gene)
            # Force-assign the highest-scoring OH that doesn't literally duplicate or
            # reverse-complement an already-assigned OH (best effort).
            assigned_set = set(assigned)
            for oh in all4_by_score:
                if oh == rev_comp(oh) or any(oh == c*4 for c in 'ACGT'):
                    continue
                if oh in assigned_set or rev_comp(oh) in assigned_set:
                    continue
                oh_map[gene] = {'oh': oh, 'split_from_end': 4,
                                'neb_score': _oh_self_score(oh, _neb_bbsi)}
                assigned.append(oh)
                break

    if needs_separate_pool:
        print(f'  [WARN] {len(needs_separate_pool)} gene(s) could not satisfy Hamming '
              f'>= {min_hamming} — flag for a separate BbsI reaction pool: '
              f'{needs_separate_pool}')

    return oh_map


# ── Gene sequence builders ───────────────────────────────────────────────────

def build_trbv_gene_nt(trbv_name):
    """
    Return the deterministically codon-optimised NT sequence for the TRBV
    gene up to and including Cys104.
    """
    aa_full = _trbv_aa_dict.get(trbv_name)
    if aa_full is None:
        raise ValueError(f'TRBV gene not found in references: {trbv_name}')
    aa_trimmed = trim_vgene_to_cys104(aa_full, trbv_name)
    nt = back_translate_deterministic(aa_trimmed)
    nt = clean_nt_sequence(nt, aa_trimmed, label=trbv_name)
    return nt


def build_trav_gene_nt(trav_name):
    """
    Return the deterministically codon-optimised full TRAV NT sequence (with leader)
    up to and including Cys104.
    """
    aa_full = _trav_aa_dict.get(trav_name)
    if aa_full is None:
        raise ValueError(f'TRAV gene not found in references: {trav_name}')
    aa_trimmed = trim_vgene_to_cys104(aa_full, trav_name)
    nt = back_translate_deterministic(aa_trimmed)
    nt = clean_nt_sequence(nt, aa_trimmed, label=trav_name)
    return nt


def build_trac_nt():
    """Return cleaned NT sequence for TRAC-LVL-Cys (deterministic codon optimisation)."""
    nt = back_translate_deterministic(_trac_lvlcys_aa)
    return clean_nt_sequence(nt, _trac_lvlcys_aa, label='TRAC-LVL-Cys')


def build_trbc_nt():
    """Return cleaned NT sequence for TRBC-LVL-Cys (deterministic codon optimisation)."""
    nt = back_translate_deterministic(_trbc_lvlcys_aa)
    return clean_nt_sequence(nt, _trbc_lvlcys_aa, label='TRBC-LVL-Cys')


# ── Fragment assemblers ──────────────────────────────────────────────────────

def build_lh_trbv_fragment(trbv_name, oh_1a, split_from_end):
    """
    Build the LH-TRBV constant gene block.

    The 4nt oh_1a is the golden gate junction overhang.  It lies within the TRBV
    coding sequence at position [-split_from_end : -split_from_end+4].  The gene
    block body contains TRBV_nt up to (but not including) the oh_1a window; the oh
    appears only as the sticky end in the BbsI site so that after ligation the full
    TRBV sequence is reconstituted exactly once.

    Gene-block layout:
      [LH] [TRBV_nt[:-split_from_end]] GAAGACGGA [oh_1a] ACAT

    After BbsI digest and ligation with the variable fragment:
      [LH] [P2A_leader] [TRBV_nt[:-split_from_end]] + oh_1a + [bridge + CDR3β + ...]
      = [LH] [P2A_leader] [TRBV_nt_full] [CDR3β ...]
    """
    lh   = _lh_nt.upper()
    p2a  = _p2a_leader_nt          # frame-verified by upstream TRAC ORF via LH arm
    trbv_nt = build_trbv_gene_nt(trbv_name)
    body = trbv_nt[:-split_from_end] if split_from_end > 0 else trbv_nt
    gene_block = lh + p2a + body + _bbsi_3prime_site(oh_1a)
    _validate_gene_block_re(gene_block, expected_bbsi=1, expected_bsmbi=0, label=trbv_name)
    return gene_block


def build_trac_rh_fragment(oh_1b, split_1b):
    """
    Build the TRAC-LVL-Cys + right homology arm constant gene block.

    oh_1b = trac_nt[split_1b : split_1b+4] — the BbsI 1b junction overhang.
    The gene-block body starts at trac_nt[split_1b+4:]; the first split_1b nt of
    TRAC go into the variable fragment as a bridge (may be empty if split_1b == 0).

    Gene-block layout:
      ACAT [oh_1b] A GTCTTC [TRAC_nt[split_1b+4:]] [RH]

    After BbsI digest, the right functional piece = oh_1b + TRAC_nt[split_1b+4:] + RH.
    """
    trac_nt = build_trac_nt()
    rh = _rh_nt.upper()
    gene_block = _bbsi_5prime_site(oh_1b) + trac_nt[split_1b + 4:] + rh
    _validate_gene_block_re(gene_block, expected_bbsi=1, expected_bsmbi=0, label='TRAC-RH')
    return gene_block


def build_trbc_t2a_fragment(oh_2a, split_2a, oh_2b, split_2b):
    """
    Build the TRBC-LVL-Cys + partial-T2A constant gene block.

    oh_2a = trbc_nt[split_2a : split_2a+4]  — BsmBI 2a junction.
    oh_2b = t2a_nt[split_2b : split_2b+4]   — BsmBI 2b junction, derived from
            WITHIN T2A so the assembled product reconstitutes full T2A with no
            extra nucleotides (avoids the +1 frameshift that a synthetic oh_2b
            would introduce between T2A and TRAV).

    Gene-block body = TRBC_nt[split_2a+4:] + T2A_nt[:split_2b]
    The T2A tail (T2A_nt[split_2b+4:]) lives in the variable gene block.

    Gene-block layout:
      ACAT CGTCTCCA [oh_2a] [TRBC[split_2a+4:]] [T2A[:split_2b]] CGTCTCCA [oh_2b] ACAT

    After BsmBI digest + ligation, the assembled sequence is:
      TRBC + T2A[:split_2b] + oh_2b(=T2A[split_2b:split_2b+4]) + T2A[split_2b+4:] + TRAV
      = TRBC + T2A + TRAV  (correct reading frame throughout)
    """
    trbc_nt = build_trbc_nt()
    t2a     = _t2a_nt.upper()
    body    = trbc_nt[split_2a + 4:] + t2a[:split_2b]
    gene_block = _bsmbi_5prime_site(oh_2a) + body + _bsmbi_3prime_site(oh_2b)
    _validate_gene_block_re(gene_block, expected_bbsi=0, expected_bsmbi=2, label='TRBC-T2A')
    return gene_block


def _generate_variable_parts(trbv_name, trav_name, cdr3b_aa, cdr3a_aa,
                               oh_1a, trbv_bridge,
                               oh_1b, trac_bridge,
                               oh_2a, trbc_bridge,
                               oh_2b, t2a_tail,
                               max_tries=30):
    """
    Internal: generate and return all nucleotide components of the variable fragment.
    Returns a dict with keys:
      cdr3b_nt, cdr3a_nt, trav_nt, gene_block
    Raises ValueError if restriction sites or homopolymers cannot be resolved.
    """
    trav_nt = build_trav_gene_nt(trav_name)
    cdr3b_body = cdr3b_aa[1:]
    cdr3a_body = cdr3a_aa[1:]

    for attempt in range(max_tries):
        try:
            cdr3b_nt = generate_clean_nt(cdr3b_body, label=f'CDR3b({trbv_name})')
            cdr3a_nt = generate_clean_nt(cdr3a_body, label=f'CDR3a({trav_name})')
            # t2a_tail: T2A[split_2b+4:] — the remaining T2A after the 2b OH window.
            #   Placed before TRAV so the assembled product reconstitutes full T2A
            #   (TRBC-T2A fragment holds T2A[:split_2b]; variable holds T2A[split_2b+4:]).
            gene_block = (
                _bbsi_5prime_site(oh_1a)
                + trbv_bridge
                + cdr3b_nt
                + trbc_bridge
                + _bsmbi_3prime_site(oh_2a)
                + _bsmbi_5prime_site(oh_2b)
                + t2a_tail
                + trav_nt
                + cdr3a_nt
                + trac_bridge
                + _bbsi_3prime_site(oh_1b)
            )
            _validate_gene_block_re(gene_block, expected_bbsi=2, expected_bsmbi=2,
                                    label=f'variable({trbv_name}/{trav_name})')
            return {
                'gene_block': gene_block,
                'cdr3b_nt':   cdr3b_nt,
                'cdr3a_nt':   cdr3a_nt,
                'trav_nt':    trav_nt,
            }
        except ValueError:
            if attempt == max_tries - 1:
                raise
    raise ValueError('_generate_variable_parts: exceeded max_tries')


def build_variable_fragment(trbv_name, trav_name, cdr3b_aa, cdr3a_aa,
                              oh_1a, trbv_bridge,
                              oh_1b, trac_bridge,
                              oh_2a, trbc_bridge,
                              oh_2b, t2a_tail,
                              max_tries=30):
    """
    Build the per-TCR variable gene block.

    Gene-block layout (5'→3'):
      ACAT [oh_1a] A GTCTTC           ← BbsI 5' site exposing oh_1a
      [trbv_bridge]                   ← TRBV tail after oh_1a window (may be empty)
      [CDR3β_nt]                      ← codon-optimised CDR3β (minus leading C)
      [trbc_bridge]                   ← TRBC prefix before oh_2a (may be empty)
      CGTCTCCA [oh_2a] ACAT           ← BsmBI 3' site: CDR3β/TRBC junction
      ACAT CGTCTCCA [oh_2b]           ← BsmBI 5' site: T2A/TRAV junction
      [TRAV_nt]                       ← full TRAV with leader, codon-optimised
      [CDR3α_nt]                      ← codon-optimised CDR3α (minus leading C)
      [trac_bridge]                   ← TRAC prefix before oh_1b (may be empty)
      GAAGACGGA [oh_1b] ACAT          ← BbsI 3' site exposing oh_1b

    trbv_bridge: the TRBV nucleotides AFTER the oh_1a window and before CDR3β
                 (may be empty if the oh window ends exactly at Cys104).

    After step-1 BbsI ligation:
      full_TRBV_nt + CDR3β_nt + [inner BsmBI sites] + TRAV_nt + CDR3α_nt + TRAC

    After step-2 BsmBI ligation (TRBC-T2A inserted):
      full_TRBV_nt + CDR3β_nt + TRBC_nt + T2A_nt + TRAV_nt + CDR3α_nt + TRAC
    """
    parts = _generate_variable_parts(trbv_name, trav_name, cdr3b_aa, cdr3a_aa,
                                     oh_1a, trbv_bridge,
                                     oh_1b, trac_bridge,
                                     oh_2a, trbc_bridge,
                                     oh_2b, t2a_tail, max_tries)
    return parts['gene_block']


def _validate_gene_block_re(seq, expected_bbsi, expected_bsmbi, label=''):
    """Assert the gene block has exactly the expected number of BbsI and BsmBI sites."""
    from Bio.Restriction import Analysis as RA
    rb = BbsI + BsmBI
    result = RA(rb, Seq(seq.upper())).full()
    n_bbsi  = len(result.get(BbsI,  []))
    n_bsmbi = len(result.get(BsmBI, []))
    if n_bbsi != expected_bbsi:
        raise ValueError(
            f'{label}: expected {expected_bbsi} BbsI site(s), found {n_bbsi} in:\n{seq[:80]}...')
    if n_bsmbi != expected_bsmbi:
        raise ValueError(
            f'{label}: expected {expected_bsmbi} BsmBI site(s), found {n_bsmbi} in:\n{seq[:80]}...')


# ── GenBank output ───────────────────────────────────────────────────────────

def _feat(start, end, feat_type, strand=1, **qualifiers):
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    return SeqFeature(FeatureLocation(start, end, strand=strand),
                      type=feat_type, qualifiers=qualifiers)


def _write_gb(seq_nt, record_id, description, features, filepath):
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
    record = SeqRecord(
        Seq(seq_nt.upper()),
        id=record_id[:16],
        name=record_id[:16],
        description=description,
        features=features,
    )
    record.annotations['molecule_type'] = 'DNA'
    SeqIO.write(record, filepath, 'genbank')


def _annotate_lh_trbv_gb(lh_nt, trbv_body_nt, oh_1a, trbv_name, filepath):
    """Write annotated GenBank for an LH-TRBV constant gene block."""
    p2a  = _p2a_leader_nt
    bbsi = _bbsi_3prime_site(oh_1a)          # GAAGACGGA + oh + ACAT
    seq  = lh_nt.upper() + p2a + trbv_body_nt.upper() + bbsi

    lh_end   = len(lh_nt)
    p2a_end  = lh_end + len(p2a)
    trbv_end = p2a_end + len(trbv_body_nt)
    oh_start = trbv_end + 9   # after GAAGACGGA (9 nt)

    features = [
        _feat(0,        lh_end,     'misc_feature', label=['Left Homology Arm'],
              note=['Homology arm for Cas9-RNP HDR at mTRAC locus']),
        _feat(lh_end,   p2a_end,    'misc_feature', label=['P2A leader'],
              note=['P2A self-cleaving peptide. Reading frame is established by the '
                    'upstream endogenous TRAC ORF across the left homology arm — '
                    'the apparent in-frame TGA is out-of-frame in the genomic context. '
                    'Do NOT codon-optimise or alter this sequence.']),
        _feat(p2a_end,  trbv_end,   'gene',         label=[f'TRBV ({trbv_name})'],
              gene=[trbv_name],
              note=['Codon-optimised, includes signal peptide; ends at BbsI 1a junction']),
        _feat(trbv_end, len(seq)-4, 'protein_bind', label=['BbsI recognition (1a)'],
              bound_moiety=['BbsI'],
              note=[f'Step-1 golden gate site; 1a overhang = {oh_1a}']),
        _feat(oh_start, oh_start+4, 'misc_feature', label=[f'1a OH ({oh_1a})'],
              note=['TRBV-specific BbsI overhang for step-1 pooled assembly']),
    ]
    _write_gb(seq, f'LH-{trbv_name}', f'LH-{trbv_name} gene block | 1a={oh_1a}',
              features, filepath)


def _annotate_trac_rh_gb(trac_nt, rh_nt, oh_1b, filepath):
    """Write annotated GenBank for the TRAC-RH constant gene block."""
    bbsi = _bbsi_5prime_site(oh_1b)          # ACAT + oh + A + GTCTTC
    seq  = bbsi + trac_nt[4:].upper() + rh_nt.upper()
    trac_start = len(bbsi)
    rh_start   = trac_start + len(trac_nt) - 4

    features = [
        _feat(0,          len(bbsi)-6,    'protein_bind', label=['BbsI recognition (1b)'],
              bound_moiety=['BbsI'],
              note=[f'Step-1 golden gate site; 1b overhang = {oh_1b}']),
        _feat(4,          8,              'misc_feature', label=[f'1b OH ({oh_1b})'],
              note=['Constant BbsI overhang at CDR3α/TRAC junction']),
        _feat(trac_start, rh_start,       'CDS',          label=['TRAC (LVL-Cys)'],
              gene=['Trac'], note=['LVL-Cys engineered constant region (improved TCR pairing)']),
        _feat(rh_start,   len(seq),       'misc_feature', label=['Right Homology Arm'],
              note=['Homology arm for Cas9-RNP HDR at mTRAC locus']),
    ]
    _write_gb(seq, 'TRAC-RH', f'TRAC-RH gene block | 1b={oh_1b}', features, filepath)


def _annotate_trbc_t2a_gb(trbc_nt, t2a_nt, oh_2a, oh_2b, filepath):
    """Write annotated GenBank for the TRBC-T2A constant gene block."""
    left  = _bsmbi_5prime_site(oh_2a)        # ACAT + CGTCTCCA + oh_2a
    right = _bsmbi_3prime_site(oh_2b)        # CGTCTCCA + oh_2b + ACAT
    seq   = left + trbc_nt[4:].upper() + t2a_nt.upper() + right

    trbc_start = len(left)
    t2a_start  = trbc_start + len(trbc_nt) - 4
    right_start = t2a_start + len(t2a_nt)

    features = [
        _feat(0,           len(left)-4,   'protein_bind', label=['BsmBI recognition (2a)'],
              bound_moiety=['BsmBI'],
              note=[f'Step-2 golden gate site; 2a overhang = {oh_2a}']),
        _feat(4,           8,             'misc_feature', label=[f'2a OH ({oh_2a})'],
              note=['Constant BsmBI overhang at CDR3β/TRBC junction']),
        _feat(trbc_start,  t2a_start,     'CDS',          label=['TRBC (LVL-Cys)'],
              gene=['Trbc'], note=['LVL-Cys engineered constant region']),
        _feat(t2a_start,   right_start,   'misc_feature', label=['furin-GSG-T2A'],
              note=['Self-cleaving T2A peptide separating beta and alpha chains']),
        _feat(right_start, right_start+8, 'protein_bind', label=['BsmBI recognition (2b)'],
              bound_moiety=['BsmBI'],
              note=[f'Step-2 golden gate site; 2b overhang = {oh_2b}']),
        _feat(right_start+8, right_start+12, 'misc_feature', label=[f'2b OH ({oh_2b})'],
              note=['Constant BsmBI overhang at T2A/TRAV junction']),
    ]
    _write_gb(seq, 'TRBC-T2A', f'TRBC-T2A gene block | 2a={oh_2a} 2b={oh_2b}',
              features, filepath)


def _annotate_variable_gb(trbv_name, trav_name, cdr3b_aa, cdr3a_aa,
                           gene_block, trbv_bridge, cdr3b_nt, trav_nt, cdr3a_nt,
                           oh_1a, oh_1b, oh_2a, oh_2b, filepath):
    """Write annotated GenBank for a per-TCR variable gene block."""
    # Site lengths
    bbsi5  = len(_bbsi_5prime_site(oh_1a))   # = 4+4+1+6 = 15
    bsmbi3 = len(_bsmbi_3prime_site(oh_2a))  # = 6+2+4+4 = 16
    bsmbi5 = len(_bsmbi_5prime_site(oh_2b))  # = 4+6+2+4 = 16
    bbsi3  = len(_bbsi_3prime_site(oh_1b))   # = 6+3+4+4 = 17

    p = bbsi5
    bridge_end  = p + len(trbv_bridge)
    cdr3b_end   = bridge_end + len(cdr3b_nt)
    bsmbi3_end  = cdr3b_end + bsmbi3
    bsmbi5_end  = bsmbi3_end + bsmbi5
    trav_end    = bsmbi5_end + len(trav_nt)
    cdr3a_end   = trav_end + len(cdr3a_nt)

    features = [
        _feat(0,         bbsi5,       'protein_bind', label=['BbsI recognition (1a)'],
              bound_moiety=['BbsI'],
              note=[f'Step-1 GG site; 1a OH = {oh_1a} (TRBV-specific)']),
        _feat(4,         8,           'misc_feature', label=[f'1a OH ({oh_1a})'],
              note=[f'Pairs with LH-{trbv_name} gene block']),
    ]
    if trbv_bridge:
        features.append(
            _feat(p, bridge_end, 'misc_feature', label=[f'TRBV bridge ({trbv_name})'],
                  note=['Germline TRBV tail after the 1a overhang, before CDR3β'])
        )
    features += [
        _feat(bridge_end,  cdr3b_end,  'misc_feature', label=[f'CDR3β ({cdr3b_aa})'],
              note=[f'Codon-optimised CDR3β (aa: {cdr3b_aa[1:]}); excludes leading C in TRBV']),
        _feat(cdr3b_end,   bsmbi3_end, 'protein_bind', label=['BsmBI recognition (2a)'],
              bound_moiety=['BsmBI'],
              note=[f'Step-2 GG site; 2a OH = {oh_2a}']),
        _feat(cdr3b_end+8, cdr3b_end+12, 'misc_feature', label=[f'2a OH ({oh_2a})'],
              note=['Junction to TRBC-T2A constant fragment']),
        _feat(bsmbi3_end,  bsmbi5_end, 'protein_bind', label=['BsmBI recognition (2b)'],
              bound_moiety=['BsmBI'],
              note=[f'Step-2 GG site; 2b OH = {oh_2b}']),
        _feat(bsmbi3_end+4, bsmbi3_end+8, 'misc_feature', label=[f'2b OH ({oh_2b})'],
              note=['Junction from TRBC-T2A constant fragment']),
        _feat(bsmbi5_end,  trav_end,   'gene',         label=[f'TRAV ({trav_name})'],
              gene=[trav_name],
              note=['Full codon-optimised TRAV, includes signal peptide, ends at Cys104']),
        _feat(trav_end,    cdr3a_end,  'misc_feature', label=[f'CDR3α ({cdr3a_aa})'],
              note=[f'Codon-optimised CDR3α (aa: {cdr3a_aa[1:]}); excludes leading C in TRAV']),
        _feat(cdr3a_end,   cdr3a_end+bbsi3, 'protein_bind', label=['BbsI recognition (1b)'],
              bound_moiety=['BbsI'],
              note=[f'Step-1 GG site; 1b OH = {oh_1b}']),
        _feat(cdr3a_end+9, cdr3a_end+13, 'misc_feature', label=[f'1b OH ({oh_1b})'],
              note=['Junction to TRAC-RH constant fragment']),
    ]

    desc = (f'Variable fragment | {trbv_name}/{trav_name} | '
            f'CDR3b={cdr3b_aa} CDR3a={cdr3a_aa}')
    _write_gb(gene_block, 'variable', desc, features, filepath)


def _annotate_assembled_gb(trbv_name, trav_name, cdr3b_aa, cdr3a_aa,
                            trbv_nt, cdr3b_nt, trav_nt, cdr3a_nt,
                            trac_nt, trbc_nt, t2a_nt, lh_nt, rh_nt,
                            oh_2b, clone_id, filepath):
    """
    Write the fully assembled HDR template as an annotated GenBank record.
    This represents the DNA that integrates into the genome after HDR.
    No golden gate recognition sequences — only the functional insert.
    """
    # Build assembled sequence (no GG sites; T2A reconstituted fully, no oh_2b linker)
    p2a = _p2a_leader_nt
    parts = [lh_nt, p2a, trbv_nt, cdr3b_nt, trbc_nt, t2a_nt, trav_nt, cdr3a_nt, trac_nt, rh_nt]
    seq = ''.join(p.upper() for p in parts)

    pos = {}
    cursor = 0
    for label, s in zip(
        ['lh','p2a','trbv','cdr3b','trbc','t2a','trav','cdr3a','trac','rh'], parts
    ):
        pos[label] = (cursor, cursor + len(s))
        cursor += len(s)

    features = [
        _feat(*pos['lh'],   'misc_feature', label=['Left Homology Arm'],
              note=['~300bp homology for Cas9-RNP HDR at mTRAC locus']),
        _feat(*pos['p2a'],  'misc_feature', label=['P2A leader'],
              note=['P2A self-cleaving peptide. Reading frame established by upstream '
                    'TRAC ORF across the LH arm: last nt of LH + first 2 nt of P2A '
                    'complete a shared TCT=Ser codon.']),
        _feat(*pos['trbv'], 'gene',         label=[f'TRBV ({trbv_name})'],
              gene=[trbv_name],
              note=[f'Codon-optimised {trbv_name}; includes signal peptide']),
        _feat(*pos['cdr3b'],'misc_feature', label=[f'CDR3β'],
              note=[f'aa: {cdr3b_aa} | codon-optimised | {len(cdr3b_nt)//3} codons']),
        _feat(*pos['trbc'], 'CDS',          label=['TRBC (LVL-Cys)'],
              gene=['Trbc'],
              note=['LVL-Cys engineered constant region for improved TCRαβ pairing']),
        _feat(*pos['t2a'],  'misc_feature', label=['furin-GSG-T2A'],
              note=['Self-cleaving T2A. oh_2b is derived from within this sequence '
                    'so T2A is fully reconstituted with no extra nt at the T2A/TRAV '
                    'junction — preserving the reading frame into TRAV.']),
        _feat(*pos['trav'], 'gene',         label=[f'TRAV ({trav_name})'],
              gene=[trav_name],
              note=[f'Codon-optimised {trav_name}; includes signal peptide']),
        _feat(*pos['cdr3a'],'misc_feature', label=[f'CDR3α'],
              note=[f'aa: {cdr3a_aa} | codon-optimised | {len(cdr3a_nt)//3} codons']),
        _feat(*pos['trac'], 'CDS',          label=['TRAC (LVL-Cys)'],
              gene=['Trac'],
              note=['LVL-Cys engineered constant region for improved TCRαβ pairing']),
        _feat(*pos['rh'],   'misc_feature', label=['Right Homology Arm'],
              note=['~300bp homology for Cas9-RNP HDR at mTRAC locus']),
    ]

    desc = (f'HDR template | {clone_id} | {trbv_name}/{trav_name} | '
            f'CDR3b={cdr3b_aa} CDR3a={cdr3a_aa} | {len(seq)} bp')
    _write_gb(seq, clone_id[:16], desc, features, filepath)


def write_gb_outputs(output_dir, constant_df, variable_df, oh_1a_map,
                     trbv_nt_cache, oh_2b, trac_nt, trbc_nt, t2a_nt,
                     lh_nt, rh_nt, verbose=True):
    """
    Write GenBank files for all constant and variable fragments.

    Output layout:
      output_dir/
        gb/
          constants/
            LH-TRBV1.gb … LH-TRBV31.gb
            TRAC-RH.gb
            TRBC-T2A.gb
          clonotypes/
            <clone_id>/
              variable.gb    ← gene block with GG site annotations
              assembled.gb   ← final HDR template, no GG sites
    """
    gb_root = Path(output_dir) / 'gb'
    const_dir = gb_root / 'constants'
    clone_dir = gb_root / 'clonotypes'
    const_dir.mkdir(parents=True, exist_ok=True)
    clone_dir.mkdir(parents=True, exist_ok=True)

    if verbose:
        print('\nWriting GenBank files...')

    # ── Constant fragments ────────────────────────────────────────────────────
    # TRBC-T2A
    row = constant_df[constant_df['Fragment'] == 'TRBC-T2A']
    if not row.empty:
        oh_2a = row.iloc[0]['OH_left']
        _annotate_trbc_t2a_gb(trbc_nt, t2a_nt, oh_2a, oh_2b,
                               const_dir / 'TRBC-T2A.gb')

    # TRAC-RH
    row = constant_df[constant_df['Fragment'] == 'TRAC-RH']
    if not row.empty:
        oh_1b = row.iloc[0]['OH_left']
        _annotate_trac_rh_gb(trac_nt, rh_nt, oh_1b, const_dir / 'TRAC-RH.gb')

    # LH-TRBV per gene
    for _, row in constant_df[constant_df['Fragment'].str.startswith('LH-')].iterrows():
        gene = row['Gene']
        oh_1a = row['OH_right']
        if gene not in trbv_nt_cache or gene not in oh_1a_map:
            continue
        trbv_nt = trbv_nt_cache[gene]
        split   = oh_1a_map[gene]['split_from_end']
        body    = trbv_nt[:-split] if split > 0 else trbv_nt
        _annotate_lh_trbv_gb(lh_nt, body, oh_1a, gene,
                              const_dir / f'LH-{gene}.gb')

    if verbose:
        print(f'  Wrote {len(list(const_dir.glob("*.gb")))} constant .gb files')

    # ── Variable + assembled per clonotype ───────────────────────────────────
    oh_1b_seq = trac_nt[:4]   # used for assembled reconstruction

    written = 0
    for _, row in variable_df.iterrows():
        clone_id  = str(row['Clone'])
        trbv_name = row['TRBV']
        trav_name = row['TRAV']
        cdr3b_aa  = row['CDR3b']
        cdr3a_aa  = row['CDR3a']
        gene_block = row['Sequence']

        if trbv_name not in oh_1a_map or trbv_name not in trbv_nt_cache:
            continue

        oh_1a_entry = oh_1a_map[trbv_name]
        oh_1a       = oh_1a_entry['oh']
        split       = oh_1a_entry['split_from_end']
        trbv_nt_full = trbv_nt_cache[trbv_name]
        trbv_bridge  = trbv_nt_full[-split+4:] if split > 4 else ''

        # Reconstruct cdr3b_nt, cdr3a_nt, trav_nt from the gene block
        # (parse back out of the stored gene block string)
        bbsi5_len  = len(_bbsi_5prime_site(oh_1a))
        bsmbi3_len = len(_bsmbi_3prime_site(oh_1a_entry['oh']))  # same formula
        bsmbi5_len = len(_bsmbi_5prime_site(oh_2b))
        bbsi3_len  = len(_bbsi_3prime_site(oh_1b_seq))
        try:
            trav_nt_len = len(build_trav_gene_nt(trav_name))
        except ValueError:
            continue
        # CDR3 lengths from gene block
        cdr3b_start = bbsi5_len + len(trbv_bridge)
        cdr3b_end   = len(gene_block) - bsmbi3_len - bsmbi5_len - trav_nt_len - len(cdr3a_aa[1:])*3 - bbsi3_len
        cdr3b_nt    = gene_block[cdr3b_start:cdr3b_end]
        trav_start  = cdr3b_end + bsmbi3_len + bsmbi5_len
        trav_end    = trav_start + trav_nt_len
        cdr3a_nt    = gene_block[trav_end:len(gene_block)-bbsi3_len]
        trav_nt     = gene_block[trav_start:trav_end]

        # Write clonotype folder
        cdir = clone_dir / clone_id
        cdir.mkdir(exist_ok=True)

        _annotate_variable_gb(
            trbv_name, trav_name, cdr3b_aa, cdr3a_aa,
            gene_block, trbv_bridge, cdr3b_nt, trav_nt, cdr3a_nt,
            oh_1a, oh_1b_seq, trbc_nt[:4], oh_2b,
            cdir / 'variable.gb',
        )
        _annotate_assembled_gb(
            trbv_name, trav_name, cdr3b_aa, cdr3a_aa,
            trbv_nt_full, cdr3b_nt, trav_nt, cdr3a_nt,
            trac_nt, trbc_nt, t2a_nt, lh_nt, rh_nt,
            oh_2b, clone_id,
            cdir / 'assembled.gb',
        )
        written += 1

    if verbose:
        print(f'  Wrote {written} clonotype folders ({written*2} .gb files)')
    return gb_root


# ── AIRR parsing ─────────────────────────────────────────────────────────────

def _normalise_gene(name):
    """Strip allele suffix and normalise to uppercase: 'TRBV5*01' → 'TRBV5'."""
    if not name or str(name).lower() in ('nan', 'none', ''):
        return ''
    return str(name).split('*')[0].strip().upper()


def load_airr(path):
    """Load an AIRR rearrangement TSV. Return a DataFrame."""
    df = pd.read_csv(path, sep='\t', low_memory=False)
    return df


def pair_airr_records(df):
    """
    Group AIRR records into alpha/beta pairs by clone_id.
    Returns list of dicts with keys: clone_id, tra_row, trb_row.
    Only fully paired (one TRA + one TRB, both productive) clones are returned.
    """
    pairs = []
    required_cols = {'v_call', 'j_call', 'junction_aa', 'productive'}
    if not required_cols.issubset(df.columns):
        missing = required_cols - set(df.columns)
        raise ValueError(f'AIRR file missing required columns: {missing}')

    group_col = 'clone_id' if 'clone_id' in df.columns else 'cell_id'
    for clone_id, grp in df.groupby(group_col, dropna=False):
        prod = grp[grp['productive'].astype(str).str.upper() == 'T']
        tra = prod[prod['v_call'].str.contains('TRAV', na=False)]
        trb = prod[prod['v_call'].str.contains('TRBV', na=False)]
        if len(tra) == 1 and len(trb) == 1:
            pairs.append({
                'clone_id': clone_id,
                'tra': tra.iloc[0],
                'trb': trb.iloc[0],
            })
    return pairs


# ── Top-level design function ────────────────────────────────────────────────

def generate_hdr_design(airr_path_or_df, trbv_filter=None, output_dir=None, verbose=True):
    """
    Main entry point: design all HDR fragments for TCRs in an AIRR file.

    Parameters
    ----------
    airr_path_or_df : str | Path | DataFrame
        AIRR rearrangement file (TSV) or pre-loaded DataFrame.
    trbv_filter : str | None
        If given, only process clones using this TRBV gene.
        Also controls which LH-TRBV constant fragment(s) to build.
    output_dir : str | Path | None
        Where to write output files. Defaults to ./TCRAFT-hdr_<date>.
    verbose : bool

    Returns
    -------
    dict with keys:
        'constant_fragments'  : DataFrame (LH-TRBV, TRAC-RH, TRBC-T2A sequences)
        'variable_fragments'  : DataFrame (per-TCR variable gene block sequences)
        'overhangs'           : dict of OH sequences (1a per TRBV, 1b, 2a, 2b)
        'invalid'             : DataFrame (skipped clones with error messages)
    """
    if output_dir is None:
        date_str = pd.Timestamp.now().strftime('%Y%m%d')
        output_dir = Path(f'TCRAFT-hdr_{date_str}')
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if isinstance(airr_path_or_df, (str, Path)):
        df = load_airr(airr_path_or_df)
    else:
        df = airr_path_or_df

    pairs = pair_airr_records(df)
    if verbose:
        print(f'Found {len(pairs)} productive paired clonotypes in AIRR file.')
    if trbv_filter:
        pairs = [p for p in pairs if _normalise_gene(p['trb']['v_call']) == trbv_filter.upper()]
        if verbose:
            print(f'  Filtered to {len(pairs)} clones with TRBV = {trbv_filter}.')

    # ── Build constant NT sequences (shared across all TCRs) ──────────────
    if verbose:
        print('\nBuilding constant reference sequences...')
    trac_nt = build_trac_nt()
    trbc_nt = build_trbc_nt()
    t2a_nt  = _t2a_nt.upper()

    oh_1b, split_1b, oh_2a, split_2a, oh_2b, split_2b = design_constant_overhangs(trac_nt, trbc_nt)
    # Bridges / tails: sequence parts that travel in the variable fragment
    trac_bridge = trac_nt[:split_1b].upper()          # TRAC prefix before oh_1b
    trbc_bridge = trbc_nt[:split_2a].upper()          # TRBC prefix before oh_2a
    t2a_tail    = t2a_nt[split_2b + 4:].upper()      # T2A suffix after oh_2b window

    # ── Design 1a overhangs for all functional TRBV genes ─────────────────
    if verbose:
        print('Designing TRBV-specific 1a overhangs...')
    trbv_genes_to_build = (
        [trbv_filter.upper()] if trbv_filter
        else sorted(_trbv_aa_dict.keys())
    )
    trbv_nt_cache = {}
    for gene in trbv_genes_to_build:
        try:
            trbv_nt_cache[gene] = build_trbv_gene_nt(gene)
        except ValueError as e:
            if verbose:
                print(f'  [SKIP] {gene}: {e}')

    # 1a overhangs share the BbsI reaction with 1b only.
    oh_1a_map = design_1a_overhangs(
        trbv_nt_cache,
        reserved_ohs=(oh_1b,),
    )
    # oh_1a_map = {gene: {'oh': str, 'split_from_end': int}}

    if verbose:
        _report_oh_conflicts({g: v['oh'] for g, v in oh_1a_map.items()},
                              oh_1b, oh_2a, oh_2b)  # reported OHs only

    # ── Build constant gene blocks ─────────────────────────────────────────
    if verbose:
        print('\nBuilding constant gene blocks...')
    constant_rows = []

    # TRBC-T2A (one block)
    try:
        trbc_t2a_gb = build_trbc_t2a_fragment(oh_2a, split_2a, oh_2b, split_2b)
        constant_rows.append({
            'Fragment': 'TRBC-T2A',
            'Gene': 'TRBC-LVL-Cys + T2A',
            'OH_left': oh_2a,
            'OH_right': oh_2b,
            'Length': len(trbc_t2a_gb),
            'GC': round(gc_fraction(trbc_t2a_gb), 3),
            'Sequence': trbc_t2a_gb,
        })
        if verbose:
            print(f'  [OK] TRBC-T2A: {len(trbc_t2a_gb)} nt')
    except ValueError as e:
        if verbose:
            print(f'  [FAIL] TRBC-T2A: {e}')

    # TRAC-RH (one block)
    try:
        trac_rh_gb = build_trac_rh_fragment(oh_1b, split_1b)
        constant_rows.append({
            'Fragment': 'TRAC-RH',
            'Gene': 'TRAC-LVL-Cys + RH',
            'OH_left': oh_1b,
            'OH_right': '',
            'Length': len(trac_rh_gb),
            'GC': round(gc_fraction(trac_rh_gb), 3),
            'Sequence': trac_rh_gb,
        })
        if verbose:
            print(f'  [OK] TRAC-RH: {len(trac_rh_gb)} nt')
    except ValueError as e:
        if verbose:
            print(f'  [FAIL] TRAC-RH: {e}')

    # LH-TRBV (one per gene)
    for gene in tqdm(trbv_genes_to_build, desc='LH-TRBV gene blocks', disable=not verbose):
        if gene not in oh_1a_map:
            continue
        oh_1a    = oh_1a_map[gene]['oh']
        split    = oh_1a_map[gene]['split_from_end']
        try:
            gb = build_lh_trbv_fragment(gene, oh_1a, split)
            constant_rows.append({
                'Fragment': f'LH-{gene}',
                'Gene': gene,
                'OH_left': '',
                'OH_right': oh_1a,
                'Length': len(gb),
                'GC': round(gc_fraction(gb), 3),
                'Sequence': gb,
            })
        except ValueError as e:
            if verbose:
                print(f'  [FAIL] LH-{gene}: {e}')

    constant_df = pd.DataFrame(constant_rows)

    # ── Build per-TCR variable gene blocks ────────────────────────────────
    if verbose:
        print(f'\nBuilding variable gene blocks for {len(pairs)} TCRs...')
    var_rows, invalid_rows = [], []

    for pair in tqdm(pairs, desc='Variable fragments', disable=not verbose):
        trb = pair['trb']
        tra = pair['tra']
        clone_id = pair['clone_id']

        trbv = _normalise_gene(trb.get('v_call', ''))
        trav = _normalise_gene(tra.get('v_call', ''))
        cdr3b = str(trb.get('junction_aa', '')).upper().strip()
        cdr3a = str(tra.get('junction_aa', '')).upper().strip()

        # Basic validation
        issues = []
        if not trbv or trbv not in oh_1a_map:
            issues.append(f'TRBV {trbv!r} not in 1a OH map')
        if not trav or trav not in _trav_aa_dict:
            issues.append(f'TRAV {trav!r} not in reference')
        if not cdr3b or not cdr3b.startswith('C'):
            issues.append(f'CDR3b {cdr3b!r} does not start with C')
        if not cdr3a or not cdr3a.startswith('C'):
            issues.append(f'CDR3a {cdr3a!r} does not start with C')
        if issues:
            invalid_rows.append({'clone_id': clone_id, 'TRBV': trbv, 'TRAV': trav,
                                  'CDR3b': cdr3b, 'CDR3a': cdr3a, 'Issue': '; '.join(issues)})
            continue

        oh_1a_entry = oh_1a_map[trbv]
        oh_1a        = oh_1a_entry['oh']
        split        = oh_1a_entry['split_from_end']
        trbv_nt_full = trbv_nt_cache[trbv]
        # Bridge: TRBV sequence AFTER the oh_1a window and before CDR3β
        trbv_bridge  = trbv_nt_full[-split+4:] if split > 4 else ''
        try:
            gb = build_variable_fragment(trbv, trav, cdr3b, cdr3a,
                                         oh_1a, trbv_bridge,
                                         oh_1b, trac_bridge,
                                         oh_2a, trbc_bridge,
                                         oh_2b, t2a_tail)
            var_rows.append({
                'Clone': clone_id,
                'TRBV': trbv,
                'TRAV': trav,
                'CDR3b': cdr3b,
                'CDR3a': cdr3a,
                'OH_1a': oh_1a,
                'OH_1b': oh_1b,
                'OH_2a': oh_2a,
                'OH_2b': oh_2b,
                'Length': len(gb),
                'GC': round(gc_fraction(gb), 3),
                'Sequence': gb,
            })
        except ValueError as e:
            invalid_rows.append({'clone_id': clone_id, 'TRBV': trbv, 'TRAV': trav,
                                  'CDR3b': cdr3b, 'CDR3a': cdr3a, 'Issue': str(e)})

    var_df     = pd.DataFrame(var_rows)
    invalid_df = pd.DataFrame(invalid_rows)

    # ── Save outputs ──────────────────────────────────────────────────────
    oh_summary = {
        'oh_1a': {g: v['oh'] for g, v in oh_1a_map.items()},
        'oh_1b': oh_1b,
        'oh_2a': oh_2a,
        'oh_2b': oh_2b,
    }

    constant_df.to_csv(output_dir / 'constant_gene_blocks.csv', index=False)
    if not var_df.empty:
        var_df.to_csv(output_dir / 'variable_gene_blocks.csv', index=False)
    if not invalid_df.empty:
        invalid_df.to_csv(output_dir / 'invalid_clones.csv', index=False)
    with open(output_dir / 'overhangs.json', 'w') as fh:
        json.dump(oh_summary, fh, indent=2)

    _write_assembly_report(output_dir, oh_summary, constant_df, var_df, invalid_df)

    # GenBank output
    if not var_df.empty:
        write_gb_outputs(
            output_dir, constant_df, var_df,
            oh_1a_map, trbv_nt_cache, oh_2b,
            trac_nt, trbc_nt, t2a_nt,
            _lh_nt.upper(), _rh_nt.upper(),
            verbose=verbose,
        )

    if verbose:
        print(f'\n✓  Designed {len(var_df)} variable fragments.')
        if not invalid_df.empty:
            print(f'   Skipped {len(invalid_df)} clones — see invalid_clones.csv.')
        print(f'   Output written to: {output_dir}')

    return {
        'constant_fragments': constant_df,
        'variable_fragments': var_df,
        'overhangs': oh_summary,
        'invalid': invalid_df,
    }


# ── Reporting helpers ─────────────────────────────────────────────────────────

def _report_oh_conflicts(oh_1a_map, oh_1b, oh_2a, oh_2b):
    s1b = _oh_self_score(oh_1b, _neb_bbsi)
    s2a = _oh_self_score(oh_2a, _neb_bsmbi)
    s2b = _oh_self_score(oh_2b, _neb_bsmbi)
    print(f'  OH summary (score = NEB self-ligation efficiency):')
    print(f'    1b={oh_1b} (BbsI score={s1b})   '
          f'2a={oh_2a} (BsmBI score={s2a})   2b={oh_2b} (BsmBI score={s2b})')
    print(f'  1a OHs ({len(oh_1a_map)} unique):')
    for g, oh in sorted(oh_1a_map.items()):
        score = _oh_self_score(oh, _neb_bbsi)
        print(f'    {g:<28} {oh}  (BbsI score={score})')


def _write_assembly_report(output_dir, ohs, const_df, var_df, invalid_df):
    report_path = output_dir / 'assembly_report.txt'
    with open(report_path, 'w') as fh:
        fh.write('TCRAFT Mouse TCR HDR Design Report\n')
        fh.write('=' * 60 + '\n\n')
        fh.write('Golden Gate Overhang Summary\n')
        fh.write('-' * 40 + '\n')
        fh.write(f"  Step 1 enzyme: BbsI\n")
        fh.write(f"  Step 2 enzyme: BsmBI\n")
        fh.write(f"  1b overhang (CDR3α/TRAC junction): {ohs['oh_1b']}\n")
        fh.write(f"  2a overhang (CDR3β/TRBC junction): {ohs['oh_2a']}\n")
        fh.write(f"  2b overhang (T2A/TRAV junction):   {ohs['oh_2b']}\n\n")
        fh.write('Per-TRBV 1a Overhangs\n')
        fh.write('-' * 40 + '\n')
        for gene, oh in sorted(ohs['oh_1a'].items()):
            fh.write(f"  {gene:<20} {oh}\n")
        fh.write('\nConstant Fragment Summary\n')
        fh.write('-' * 40 + '\n')
        for _, row in const_df.iterrows():
            fh.write(f"  {row['Fragment']:<20} {row['Length']} nt  GC={row['GC']:.0%}\n")
        fh.write(f'\nVariable Fragments: {len(var_df)}\n')
        if not var_df.empty:
            fh.write(f"  Lengths: {var_df['Length'].min()}–{var_df['Length'].max()} nt\n")
            fh.write(f"  GC:      {var_df['GC'].min():.0%}–{var_df['GC'].max():.0%}\n")
        if not invalid_df.empty:
            fh.write(f'\nSkipped Clones: {len(invalid_df)}\n')
            for _, row in invalid_df.iterrows():
                fh.write(f"  {row.get('clone_id','?')}: {row.get('Issue','')}\n")


# ── CLI ───────────────────────────────────────────────────────────────────────

def run_CLI():
    parser = argparse.ArgumentParser(
        description='TCRAFT-hdr: Design mouse TCR HDR template fragments.')
    parser.add_argument('airr', type=str,
                        help='AIRR rearrangement TSV file.')
    parser.add_argument('--trbv', type=str, default=None,
                        help='Restrict to clones using this TRBV gene (e.g. TRBV5).')
    parser.add_argument('--output_dir', type=str, default=None,
                        help='Output directory (default: ./TCRAFT-hdr_<date>).')
    args = parser.parse_args()
    generate_hdr_design(args.airr, trbv_filter=args.trbv, output_dir=args.output_dir)


if __name__ == '__main__':
    run_CLI()

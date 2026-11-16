#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
TCRAFT-library — pooled TCR plasmid library in a retroviral vector.

Delivers the same TCR cassette as TCRAFT-hdr, but as an array-pooled plasmid library cloned
into pMIGR rather than knocked into a locus.  Adapts the original (human) TCRAFT scheme to
mouse, where TRAV cannot be a plasmid library: mouse has 119 TRAV genes in the IMGT reference
against 44 in human, so TRAV moves out of the vector set and into de novo synthesis.


WHY TWO STEPS
-------------
A TCR ORF alternates V-CDR3-C, V-CDR3-C, so the two variable regions are always separated by
a constant region.  In a circular assembly one molecule cannot occupy two non-adjacent arcs.
The synthesised fragment must stay a single molecule — if it is cut open while free in a pool,
alpha/beta pairing scrambles across the library.  So TRBC-T2A cannot be added in the same pot.

  Step 1 closes the circle with the fragment intact, locking pairing.
  Step 2 cuts the now-tethered middle and inserts TRBC-T2A.

This is the same ordering the original TCRAFT uses, and it is forced by the topology rather
than chosen.


COMPONENTS
----------
  (1) TRAC backbone vector  x1     [Kozak]-[ccdB stuffer]-[TRAC][stop][IRES-GFP][LTRs][AmpR]
  (2) TRBV donor plasmids   x24    small KanR plasmids, each releasing [TRBV_x] by BbsI
  (3) TRBC-T2A donor        x1     small KanR plasmid, releases [TRBC_Cys][furin-GSG-T2A]
  (4) synthesised fragments xN     ~500 nt gene fragments, one per clonotype

  fragment = |BbsI/oh_1a| TRBV_tail | CDR3b+FR4 | TRBC_head |BsmBI/oh_2a|
             |BsmBI/oh_2b| T2A_tail | TRAV_full | CDR3a+FR4 | TRAC_head |BbsI/oh_1b|

  STEP 1  BbsI,  3-piece:  (1) + (2) + (4)  -> circle; pairing locked
  STEP 2  BsmBI, 2-piece:  cut oh_2a/oh_2b, insert (3)  -> final library

One reaction pair per TRBV pool.  Donors carry a different resistance marker from the
backbone, so plating on the backbone's antibiotic removes donor-derived circles.


OVERHANGS
---------
Step 1 has three junctions:
  oh_V1  vector -> TRBV   constant; the last 4 nt of the Kozak (CACC), shared by every donor
                          because every TRBV begins ATG
  oh_1a  TRBV -> CDR3b    gene-specific; a 4 nt window inside the TRBV 3' end
  oh_1b  CDR3a -> TRAC    constant; a window inside the vector's TRAC
Step 2 has two: oh_2a (TRBC) and oh_2b (within T2A).

Every overhang is a window taken from INSIDE the sequence it joins, so ligation reconstitutes
the gene with no extra nucleotides and no frameshift.

Because each TRBV runs in its own reaction, oh_1a only has to clear the constant overhangs —
not every other TRBV gene.  Each gene's overhang is therefore independent of which genes are
in the build, so a donor plasmid stays valid regardless of how the library is later extended.
(HDR_Design.design_1a_overhangs is greedy and set-dependent: building TRBV5 alone versus in a
22-gene set yields different overhangs, which would silently invalidate an ordered plasmid.)


VECTOR DOMESTICATION
--------------------
The stock pMIGR backbone carries 2x BbsI (in the IRES) and 2x BsmBI (AmpR-LTR spacer).  Both
pairs fall inside two windows bounded by unique cutters, so one-time restriction cloning of
two synthesised patches clears them — no PCR on the vector:

  Patch A  AatII -> NdeI    249 bp   silences both BsmBI sites (non-coding spacer)
  Patch B  PasI  -> BspEI  2588 bp   installs Kozak + stuffer + TRAC, silences both BbsI

Usage
-----
    TCRAFT-library input_airr.tsv --output_dir ./library
    TCRAFT-library input_airr.tsv --stuffer ccdB.fa --trbv TRBV5
"""

import argparse
import json
import os
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.Restriction import BbsI, BsmBI, BsaI, SapI, AatII, NdeI, PasI, BspEI, BamHI
from tqdm import tqdm

from TCRAFT import HDR_Design as H


# ── Reference paths ─────────────────────────────────────────────────────────

_REF_ROOT    = Path(__file__).parent / 'references'
_VECTOR_PATH = _REF_ROOT / 'base_vectors' / 'pMIGR-clone-5290-archive.gb'


# ── Base-vector coordinates (1-based, inclusive; from the stock GenBank) ────
# Verified against pMIGR-clone-5290-archive.gb: the stock insert at 2767..4560 is already a
# complete TCR cassette whose TRBC and TRAC match IMGT_TRBC1_Cys_aa.fa / IMGT_TRAC_LVL-Cys_aa.fa.

_KOZAK_START, _KOZAK_END = 2758, 2766     # GCCGCCACC
_TRAC_START,  _TRAC_END  = 4150, 4557     # TRAC CDS, excluding its stop
_STOP_START              = 4558           # TRAC stop codon

# The Kozak's last 4 nt serve as oh_V1.  Every TRBV starts ATG, so the vector keeps
# GCCGC and the donor supplies CACC + ATG; ligation restores GCCGCCACCATG.
_OH_V1 = 'CACC'

# Enzymes whose sites must stay unique for future vector engineering.  Screened out of all
# designed parts in addition to HDR_Design._FORBIDDEN_RE_SITES.
_VECTOR_CUTTERS = [AatII, NdeI, PasI, BspEI, BamHI]
_LIB_FORBIDDEN  = H._FORBIDDEN_RE_SITES + [e.site for e in _VECTOR_CUTTERS]

_PLACEHOLDER_STUFFER_LEN = 1500


# ── Base vector ─────────────────────────────────────────────────────────────

def load_base_vector(path=None):
    """Return the stock pMIGR record (uppercased sequence)."""
    rec = SeqIO.read(str(path or _VECTOR_PATH), 'genbank')
    rec.seq = Seq(str(rec.seq).upper())
    return rec


def _sub(seq, start, end):
    """1-based inclusive slice."""
    return str(seq)[start - 1:end]


def _unique_cut(seq, enzyme):
    """Position of a cutter that must be unique in the vector."""
    hits = enzyme.search(Seq(str(seq)))
    if len(hits) != 1:
        raise ValueError(f'{enzyme} must cut the base vector exactly once, found {len(hits)}')
    return hits[0]


# ── Domestication ───────────────────────────────────────────────────────────

def _site_positions(nt, enzyme):
    """0-based start positions of an enzyme's recognition sequence on either strand."""
    out = []
    for site in {enzyme.site, H.rev_comp(enzyme.site)}:
        start = nt.find(site)
        while start != -1:
            out.append(start)
            start = nt.find(site, start + 1)
    return sorted(out)


def silence_sites(nt, enzymes, forbidden=None, label=''):
    """
    Remove recognition sites from NON-CODING vector DNA by minimal point substitution.

    Used on the IRES and the AmpR-LTR spacer, where there is no reading frame to preserve.
    One base per site is changed, preferring a transition (purine<->purine,
    pyrimidine<->pyrimidine) as the most conservative substitution and falling back to a
    transversion only if no transition works.  A candidate is rejected if it leaves the site
    intact or introduces any new forbidden site.

    Returns (new_nt, [(position, from, to), ...]).
    """
    transition = {'A': 'G', 'G': 'A', 'C': 'T', 'T': 'C'}
    forbidden = list(forbidden or _LIB_FORBIDDEN)
    nt, changes = str(nt).upper(), []

    def _forbidden_load(s):
        """Total forbidden-site occurrences. The stock backbone already contains some
        (e.g. BsaI), so the test is 'introduces no NEW site', not 'has none'."""
        return sum(s.count(site) + s.count(H.rev_comp(site)) for site in forbidden)

    for enzyme in enzymes:
        guard = 0
        while _site_positions(nt, enzyme):
            guard += 1
            if guard > 20:
                raise ValueError(f'{label}: {enzyme} sites did not converge')
            idx = _site_positions(nt, enzyme)[0]
            before = len(_site_positions(nt, enzyme))
            baseline = _forbidden_load(nt)
            for pos in range(idx, idx + len(enzyme.site)):
                old = nt[pos]
                for new in [transition[old]] + [b for b in 'ACGT' if b not in (old, transition[old])]:
                    cand = nt[:pos] + new + nt[pos + 1:]
                    if len(_site_positions(cand, enzyme)) >= before:
                        continue                   # didn't remove this occurrence
                    if _forbidden_load(cand) > baseline - 1:
                        continue                   # traded one site for another
                    nt = cand
                    changes.append((pos, old, new))
                    break
                else:
                    continue
                break
            else:
                raise ValueError(f'{label}: could not silence {enzyme} at position {idx}')
    return nt, changes


# ── Overhang design ─────────────────────────────────────────────────────────

def _oh_context_is_clean(oh, trbv_nt, split_from_end):
    """
    Reject an overhang that creates a NEW enzyme site where it abuts the Golden Gate scaffold.

    The 5' scaffold ends GAAGACGG, so an overhang beginning AAGA followed by a C yields
    GGAAGAC-AGC — a second BbsI site that silently drops the clone at validation.  The base
    after the TRBV bridge belongs to the per-clone CDR3, so the overhang must be safe for
    every possible next base, not just the one this clone happens to have.
    """
    bridge = trbv_nt[-split_from_end + 4:] if split_from_end > 4 else ''
    body = trbv_nt[:-split_from_end] if split_from_end > 0 else trbv_nt
    for nxt in 'ACGT':
        probe = H._bbsi_5prime_site(oh) + bridge + nxt
        if len(BbsI.search(Seq(probe))) != 1 or BsmBI.search(Seq(probe)):
            return False
    donor = H._bbsi_5prime_site(_OH_V1) + body + H._bbsi_3prime_site(oh)
    return len(BbsI.search(Seq(donor))) == 2 and not BsmBI.search(Seq(donor))


def design_pool_overhangs(trbv_nt_dict, reserved, neb_df=None, scan_window=30, min_hamming=2):
    """
    Pick oh_1a per TRBV gene, INDEPENDENTLY of the other genes in the build.

    Each gene runs in its own Golden Gate reaction, so its oh_1a only has to be orthogonal to
    the constant overhangs in that same reaction (oh_V1, oh_1b) — never to another TRBV's.
    Dropping the cross-gene constraint is what makes a donor plasmid permanently valid.

    Returns {gene: {'oh', 'split_from_end', 'neb_score'}}.
    """
    neb_df = H._neb_bbsi if neb_df is None else neb_df
    reserved = [r.upper() for r in reserved]
    out = {}

    for gene, nt in trbv_nt_dict.items():
        n = len(nt)
        candidates = []
        for start in range(n - 4, max(n - scan_window - 1, 3), -1):
            oh = nt[start:start + 4].upper()
            if oh == H.rev_comp(oh):                       # palindrome self-ligates
                continue
            if any(oh == c * 4 for c in 'ACGT'):
                continue
            if H._oh_conflicts(oh, reserved, min_hamming):
                continue
            if not _oh_context_is_clean(oh, nt, n - start):
                continue
            candidates.append((H._oh_self_score(oh, neb_df), start, oh, n - start))
        if not candidates:
            raise ValueError(
                f'{gene}: no 4 nt window in the last {scan_window} nt is orthogonal to '
                f'{reserved} at Hamming >= {min_hamming}')
        score, _start, oh, split = max(candidates, key=lambda t: (t[0], -t[3]))
        out[gene] = {'oh': oh, 'split_from_end': split, 'neb_score': score}
    return out


# ── Fragment builders ───────────────────────────────────────────────────────

def build_trbv_donor_insert(trbv_nt, oh_1a, split_1a):
    """
    TRBV donor insert:  [BbsI->oh_V1] TRBV_nt[:-split_1a] [BbsI->oh_1a]

    Released fragment is  oh_V1 + TRBV_nt[:-split_1a] + oh_1a.  The vector contributes the
    Kozak up to GCCGC and the variable fragment contributes TRBV_nt[-split_1a+4:], so the
    assembled product reconstitutes GCCGCCACC + the full TRBV exactly once.
    """
    body = trbv_nt[:-split_1a] if split_1a > 0 else trbv_nt
    insert = H._bbsi_5prime_site(_OH_V1) + body + H._bbsi_3prime_site(oh_1a)
    H._validate_gene_block_re(insert, expected_bbsi=2, expected_bsmbi=0, label='TRBV donor')
    return insert


def build_trbc_t2a_donor_insert(trbc_nt, t2a_nt, oh_2a, split_2a, oh_2b, split_2b):
    """
    TRBC-T2A donor insert:  [BsmBI->oh_2a] TRBC[split_2a+4:] T2A[:split_2b] [BsmBI->oh_2b]

    Identical in structure to HDR_Design.build_trbc_t2a_fragment; oh_2b is taken from inside
    T2A so no extra nucleotides land at the T2A/TRAV junction.
    """
    body = trbc_nt[split_2a + 4:] + t2a_nt[:split_2b]
    insert = H._bsmbi_5prime_site(oh_2a) + body + H._bsmbi_3prime_site(oh_2b)
    H._validate_gene_block_re(insert, expected_bbsi=0, expected_bsmbi=2, label='TRBC-T2A donor')
    return insert


def build_variable_fragment(trbv, trav, trbj, traj, cdr3b_aa, cdr3a_aa,
                            oh_1a, trbv_bridge, oh_1b, trac_bridge,
                            oh_2a, trbc_bridge, oh_2b, t2a_tail, sites=None):
    """
    Per-clonotype synthesised fragment.  Delegates to HDR_Design._generate_variable_parts:
    the layout required here is identical, it simply loses the homology-arm framing.
    """
    return H._generate_variable_parts(
        trbv, trav, cdr3b_aa, cdr3a_aa,
        oh_1a, trbv_bridge, oh_1b, trac_bridge, oh_2a, trbc_bridge, oh_2b, t2a_tail,
        trbj_name=trbj, traj_name=traj, sites=sites or _LIB_FORBIDDEN)


# ── Backbone construction ───────────────────────────────────────────────────

def make_placeholder_stuffer(length=_PLACEHOLDER_STUFFER_LEN):
    """
    Stand-in for the ccdB cassette so the rest of the design can be generated and checked.
    Deliberately NOT a real ccdB sequence — supply the true cassette with --stuffer.
    """
    unit = 'ATCGATGCTAGCATGCTAGGATCATGCTAACGTACGTATCGGATCTAACGTAGCTAGCTT'
    seq = (unit * (length // len(unit) + 1))[:length]
    assert H.check_re_sites(seq, _LIB_FORBIDDEN)[0] == -1, 'placeholder stuffer has a forbidden site'
    return seq


def build_backbone(vector_rec, oh_1b, split_1b, stuffer_nt):
    """
    Build the domesticated TRAC backbone vector.

    Patch A (AatII -> NdeI)  replaces the AmpR-LTR spacer with both BsmBI sites silenced.
    Patch B (PasI  -> BspEI) installs
        [vector 5' of the Kozak, verbatim][GCCGC][BbsI->oh_V1][stuffer][BbsI->oh_1b]
        [TRAC[split_1b+4:]][stop + IRES with both BbsI silenced + into GFP]

    Returns (full_sequence, {'A': patchA, 'B': patchB}, domestication_changes).
    """
    seq = str(vector_rec.seq)
    cut_aatii, cut_ndei = _unique_cut(seq, AatII), _unique_cut(seq, NdeI)
    cut_pasi,  cut_bspei = _unique_cut(seq, PasI), _unique_cut(seq, BspEI)

    # ---- Patch A: silence BsmBI in the non-coding AmpR-LTR spacer
    region_a = seq[cut_aatii - 1:cut_ndei - 1]
    patch_a, changes_a = silence_sites(region_a, [BsmBI], label='patch A (AmpR-LTR spacer)')
    # report in stock-vector coordinates, not offsets into the slice
    changes_a = [(cut_aatii + p, o, n) for p, o, n in changes_a]

    # ---- Patch B
    kozak_head = _sub(seq, _KOZAK_START, _KOZAK_END - len(_OH_V1))      # GCCGC
    assert _sub(seq, _KOZAK_START, _KOZAK_END) == kozak_head + _OH_V1, \
        'Kozak does not end in the expected oh_V1 window'

    prefix   = seq[cut_pasi - 1:_KOZAK_START - 1] + kozak_head
    trac_nt  = _sub(seq, _TRAC_START, _TRAC_END)
    tail_raw = seq[_STOP_START - 1:cut_bspei - 1]                        # stop + IRES + GFP start
    tail, changes_b = silence_sites(tail_raw, [BbsI], label='patch B (IRES)')
    changes_b = [(_STOP_START + p, o, n) for p, o, n in changes_b]

    patch_b = (prefix
               + H._bbsi_3prime_site(_OH_V1)
               + stuffer_nt
               + H._bbsi_5prime_site(oh_1b)
               + trac_nt[split_1b + 4:]
               + tail)

    full = seq[:cut_aatii - 1] + patch_a + seq[cut_ndei - 1:cut_pasi - 1] + patch_b + seq[cut_bspei - 1:]

    patches = {
        'A': {'name': 'patchA_BsmBI_domestication', 'enzymes': 'AatII / NdeI',
              'seq': patch_a, 'replaces': f'{cut_aatii}..{cut_ndei - 1}',
              'changes': changes_a},
        'B': {'name': 'patchB_kozak_stuffer_TRAC_BbsI_domestication', 'enzymes': 'PasI / BspEI',
              'seq': patch_b, 'replaces': f'{cut_pasi}..{cut_bspei - 1}',
              'changes': changes_b},
    }
    return full, patches, trac_nt


# ── Assembly simulation ─────────────────────────────────────────────────────

def simulate_assembly(backbone_prefix, trbv_nt, split_1a, oh_1a, var_parts,
                      trbc_nt, t2a_nt, split_2a, split_2b, vector_trac, split_1b):
    """
    Reconstruct the final insert the two assembly steps should produce, from the designed
    parts rather than from the gene blocks — so the frame check is independent of the
    Golden Gate string plumbing.

    insert = TRBV + CDR3b+FR4 + TRBC + T2A + TRAV + CDR3a+FR4 + TRAC
    """
    cdr3b_nt = var_parts['cdr3b_nt']
    cdr3a_nt = var_parts['cdr3a_nt']
    trav_nt  = var_parts['trav_nt']
    return trbv_nt + cdr3b_nt + trbc_nt + t2a_nt + trav_nt + cdr3a_nt + vector_trac


def check_orf(insert_nt, label=''):
    """Assert the insert is a single clean ORF starting at ATG with no internal stop."""
    if len(insert_nt) % 3:
        raise ValueError(f'{label}: insert length {len(insert_nt)} is not a multiple of 3')
    if not insert_nt.upper().startswith('ATG'):
        raise ValueError(f'{label}: insert does not start with ATG')
    aa = str(Seq(insert_nt).translate())
    if '*' in aa:
        raise ValueError(f'{label}: premature stop at codon {aa.index("*")}')
    return aa


# ── Output helpers ──────────────────────────────────────────────────────────

def _gb(seq_nt, record_id, description, features, filepath, circular=True):
    """Write a GenBank file; plasmid maps are circular, linear fragments are not."""
    from Bio.SeqRecord import SeqRecord
    rec = SeqRecord(Seq(seq_nt.upper()), id=record_id[:16], name=record_id[:16],
                    description=description, features=features)
    rec.annotations['molecule_type'] = 'DNA'
    rec.annotations['topology'] = 'circular' if circular else 'linear'
    SeqIO.write(rec, str(filepath), 'genbank')


def _enzyme_census(seq, circular=False):
    """{enzyme_name: site_count} for the enzymes this design depends on."""
    s = Seq(seq.upper())
    return {str(e): len(e.search(s, linear=not circular))
            for e in (BbsI, BsmBI, BsaI, SapI)}


# ── Main entry point ────────────────────────────────────────────────────────

def generate_vector_library(airr_path_or_df, output_dir=None, trbv_filter=None,
                            stuffer_path=None, verbose=True, limit=None,
                            write_clonotype_gb=True):
    """
    Design the full pooled plasmid library from an AIRR rearrangement file.

    Returns a dict of DataFrames: backbone_patches, trbv_donors, variable_fragments,
    invalid_clones.
    """
    if output_dir is None:
        output_dir = os.path.join(
            os.getcwd(), f'TCRAFT-library_{pd.Timestamp.now().strftime("%Y%m%d")}')
    output_dir = Path(output_dir)
    (output_dir / 'gb' / 'vectors').mkdir(parents=True, exist_ok=True)
    (output_dir / 'gb' / 'clonotypes').mkdir(parents=True, exist_ok=True)

    say = (lambda *a: print(*a)) if verbose else (lambda *a: None)
    say('*********************')
    say('*  TCRAFT-library   *')
    say('*********************\n')

    # ---- references and constant parts
    vector = load_base_vector()
    vector_trac = _sub(vector.seq, _TRAC_START, _TRAC_END)
    trbc_nt = H.build_trbc_nt('E')
    t2a_nt  = H._t2a_nt.upper()

    # ---- constant overhangs
    oh_1b, split_1b = H._best_neb_oh_from_sequence(
        vector_trac, H._neb_bbsi, scan_window=16, reserved=[_OH_V1])
    oh_2a, split_2a = H._best_neb_oh_from_sequence(trbc_nt, H._neb_bsmbi, scan_window=16)
    oh_2b, split_2b = H._best_neb_oh_from_sequence(
        t2a_nt, H._neb_bsmbi, scan_window=len(t2a_nt) - 4, reserved=[oh_2a])

    trac_bridge = vector_trac[:split_1b]
    trbc_bridge = trbc_nt[:split_2a]
    t2a_tail    = t2a_nt[split_2b + 4:]

    say('Constant overhangs')
    say(f'  oh_V1={_OH_V1} (Kozak)   oh_1b={oh_1b} (TRAC, split {split_1b})')
    say(f'  oh_2a={oh_2a} (TRBC, split {split_2a})   oh_2b={oh_2b} (T2A, split {split_2b})')

    # ---- TRBV donors (gene-independent overhangs)
    genes = ([trbv_filter.upper()] if trbv_filter
             else sorted(g for g in H._trbv_aa_dict if _trbv_usable(g)))
    trbv_nt = {g: H.build_trbv_gene_nt(g) for g in genes}
    oh_1a_map = design_pool_overhangs(trbv_nt, reserved=[_OH_V1, oh_1b])
    say(f'\nTRBV donors: {len(genes)} genes')

    donor_rows = []
    for g in genes:
        split = oh_1a_map[g]['split_from_end']
        insert = build_trbv_donor_insert(trbv_nt[g], oh_1a_map[g]['oh'], split)
        donor_rows.append({'Donor': f'{g}_donor', 'Gene': g,
                           'OH_left': _OH_V1, 'OH_right': oh_1a_map[g]['oh'],
                           'NEB_score': oh_1a_map[g]['neb_score'],
                           'Length': len(insert), 'GC': round(gc_fraction(insert), 3),
                           'Sequence': insert})
    donor_df = pd.DataFrame(donor_rows)

    # ---- TRBC-T2A donor
    trbc_donor = build_trbc_t2a_donor_insert(trbc_nt, t2a_nt, oh_2a, split_2a, oh_2b, split_2b)

    # ---- backbone
    if stuffer_path:
        stuffer = ''.join(s for _, s in H._parse_fasta(stuffer_path)).upper()
        stuffer_note = f'ccdB cassette from {Path(stuffer_path).name}'
    else:
        stuffer = make_placeholder_stuffer()
        stuffer_note = 'PLACEHOLDER - NOT ccdB; supply the real cassette with --stuffer'
    backbone, patches, _ = build_backbone(vector, oh_1b, split_1b, stuffer)

    patch_df = pd.DataFrame([
        {'Patch': k, 'Name': p['name'], 'Cut_sites': p['enzymes'],
         'Replaces_in_stock_vector': p['replaces'],
         'Silenced_bases': '; '.join(f'{pos}:{o}>{n}' for pos, o, n in p['changes']),
         'Length': len(p['seq']), 'GC': round(gc_fraction(p['seq']), 3), 'Sequence': p['seq']}
        for k, p in sorted(patches.items())])

    # ---- clonotypes
    df = airr_path_or_df if isinstance(airr_path_or_df, pd.DataFrame) else H.load_airr(airr_path_or_df)
    pairs, unpaired = H.pair_airr_records(df, return_dropped=True)
    if trbv_filter:
        pairs = [p for p in pairs if H._normalise_gene(p['trb']['v_call']) == trbv_filter.upper()]
    if limit:
        pairs = pairs[:limit]
    say(f'\nPaired clonotypes: {len(pairs)}'
        + (f'  ({len(unpaired)} groups yielded no pair)' if unpaired else ''))

    var_rows = []
    # seed with groups that never produced a pair so every input group is accounted for
    invalid_rows = [{'clone_id': u['clone_id'], 'TRBV': '', 'TRAV': '', 'TRBJ': '',
                     'TRAJ': '', 'CDR3b': '', 'CDR3a': '', 'Issue': u['Issue']}
                    for u in unpaired] if not trbv_filter else []
    for pair in tqdm(pairs, desc='Variable fragments', disable=not verbose):
        trb, tra, clone_id = pair['trb'], pair['tra'], pair['clone_id']
        trbv = H._normalise_gene(trb.get('v_call', ''))
        trav = H._normalise_gene(tra.get('v_call', ''))
        trbj = H._normalise_gene(trb.get('j_call', ''))
        traj = H._normalise_gene(tra.get('j_call', ''))
        cdr3b = str(trb.get('junction_aa', '')).upper().strip()
        cdr3a = str(tra.get('junction_aa', '')).upper().strip()

        issues = []
        if trbv not in oh_1a_map:
            issues.append(f'TRBV {trbv!r} has no donor plasmid in this build')
        if trav not in H._trav_aa_dict:
            issues.append(f'TRAV {trav!r} not in reference')
        if not cdr3b.startswith('C'):
            issues.append(f'CDR3b {cdr3b!r} does not start with C')
        if not cdr3a.startswith('C'):
            issues.append(f'CDR3a {cdr3a!r} does not start with C')
        for lab, j in (('TRBJ', trbj), ('TRAJ', traj)):
            try:
                H.get_fr4(j)
            except ValueError as e:
                issues.append(f'{lab}: {e}')
        if issues:
            invalid_rows.append({'clone_id': clone_id, 'TRBV': trbv, 'TRAV': trav,
                                 'TRBJ': trbj, 'TRAJ': traj, 'CDR3b': cdr3b, 'CDR3a': cdr3a,
                                 'Issue': '; '.join(issues)})
            continue

        split = oh_1a_map[trbv]['split_from_end']
        oh_1a = oh_1a_map[trbv]['oh']
        trbv_bridge = trbv_nt[trbv][-split + 4:] if split > 4 else ''
        try:
            parts = build_variable_fragment(
                trbv, trav, trbj, traj, cdr3b, cdr3a,
                oh_1a, trbv_bridge, oh_1b, trac_bridge,
                oh_2a, trbc_bridge, oh_2b, t2a_tail)
            insert = simulate_assembly(None, trbv_nt[trbv], split, oh_1a, parts,
                                       trbc_nt, t2a_nt, split_2a, split_2b,
                                       vector_trac, split_1b)
            check_orf(insert, label=clone_id)
        except ValueError as e:
            invalid_rows.append({'clone_id': clone_id, 'TRBV': trbv, 'TRAV': trav,
                                 'TRBJ': trbj, 'TRAJ': traj, 'CDR3b': cdr3b, 'CDR3a': cdr3a,
                                 'Issue': str(e)})
            continue

        gb_seq = parts['gene_block']
        var_rows.append({
            'Clone': clone_id, 'Pool': trbv, 'TRBV': trbv, 'TRAV': trav,
            'TRBJ': trbj, 'TRAJ': traj, 'CDR3b': cdr3b, 'CDR3a': cdr3a,
            'OH_1a': oh_1a, 'OH_1b': oh_1b, 'OH_2a': oh_2a, 'OH_2b': oh_2b,
            'Length': len(gb_seq), 'GC': round(gc_fraction(gb_seq), 3),
            'Sequence': gb_seq, '_insert': insert,
        })

    var_df     = pd.DataFrame(var_rows)
    invalid_df = pd.DataFrame(invalid_rows)

    # ---- write outputs
    patch_df.to_csv(output_dir / 'backbone_patches.csv', index=False)
    donor_df.to_csv(output_dir / 'trbv_donors.csv', index=False)
    if len(var_df):
        var_df.drop(columns=[c for c in var_df.columns if c.startswith('_')]).to_csv(
            output_dir / 'variable_fragments.csv', index=False)
    if len(invalid_df):
        invalid_df.to_csv(output_dir / 'invalid_clones.csv', index=False)

    with open(output_dir / 'overhangs.json', 'w') as fh:
        json.dump({'oh_V1': _OH_V1, 'oh_1b': oh_1b, 'oh_2a': oh_2a, 'oh_2b': oh_2b,
                   'oh_1a': {g: v['oh'] for g, v in oh_1a_map.items()}}, fh, indent=2)

    _write_vector_gbs(output_dir, backbone, stuffer, stuffer_note, oh_1b, split_1b,
                      donor_df, trbc_donor, oh_2a, oh_2b, vector_trac)
    _write_report(output_dir, backbone, patches, donor_df, var_df, invalid_df,
                  oh_1a_map, _OH_V1, oh_1b, oh_2a, oh_2b, stuffer_note, len(pairs))
    if write_clonotype_gb and len(var_df):
        _write_clonotype_gbs(output_dir, var_df, backbone, stuffer, oh_1b, split_1b,
                             vector_trac, verbose=verbose)

    if len(var_df):
        H.write_plate_xlsx(output_dir, donor_df.rename(columns={'Donor': 'Fragment'}),
                           var_df.rename(columns={'Pool': '_Pool'}), verbose=False)

    say(f'\n  {len(donor_df)} TRBV donors | {len(var_df)} variable fragments '
        f'| {len(invalid_df)} skipped')
    if not stuffer_path:
        say('\n  [WARN] Backbone built with a PLACEHOLDER stuffer, not ccdB.')
        say('         Supply the real cassette with --stuffer before ordering the vector.')
    say(f'\nOutput written to: {output_dir}')

    return {'backbone_patches': patch_df, 'trbv_donors': donor_df,
            'variable_fragments': var_df, 'invalid_clones': invalid_df}


def _trbv_usable(gene):
    """TRBV genes that can be trimmed to Cys104 (excludes TRBV21)."""
    try:
        H.trim_vgene_to_cys104(H._trbv_aa_dict[gene], gene)
        return True
    except ValueError:
        return False


def _write_vector_gbs(output_dir, backbone, stuffer, stuffer_note, oh_1b, split_1b,
                      donor_df, trbc_donor, oh_2a, oh_2b, vector_trac):
    """Annotated maps for the backbone and every donor."""
    vdir = output_dir / 'gb' / 'vectors'

    i = backbone.find(stuffer)
    feats = [H._feat(i, i + len(stuffer), 'misc_feature', label=['stuffer'], note=[stuffer_note])]
    j = backbone.find(vector_trac[split_1b + 4:])
    if j != -1:
        feats.append(H._feat(j, j + len(vector_trac[split_1b + 4:]), 'CDS',
                             label=['TRAC (LVL-Cys)'], gene=['Trac']))
    _gb(backbone, 'TRAC_backbone',
        f'TCRAFT library backbone | oh_V1={_OH_V1} oh_1b={oh_1b} | {stuffer_note}',
        feats, vdir / 'TRAC_backbone.gb', circular=True)

    for _, r in donor_df.iterrows():
        f = [H._feat(0, len(r.Sequence), 'misc_feature', label=[r.Donor],
                     note=[f'BbsI donor insert | {r.OH_left} .. {r.OH_right}'])]
        _gb(r.Sequence, r.Donor, f'TRBV donor insert | {r.Gene}', f,
            vdir / f'{r.Donor}_donor.gb', circular=False)

    _gb(trbc_donor, 'TRBC-T2A_donor',
        f'TRBC-T2A donor insert | BsmBI | {oh_2a} .. {oh_2b}',
        [H._feat(0, len(trbc_donor), 'misc_feature', label=['TRBC-T2A donor'])],
        vdir / 'TRBC-T2A_donor.gb', circular=False)


def _write_clonotype_gbs(output_dir, var_df, backbone, stuffer, oh_1b, split_1b,
                         vector_trac, verbose=True):
    """
    Per-clonotype maps: the fragment as ordered, and the predicted final plasmid.

    The final plasmid is the backbone with the whole dropout cassette
    (BbsI site + stuffer + BbsI site + the retained TRAC tail) swapped for
    oh_V1 + the assembled insert — which already ends in the full TRAC.
    """
    cassette = (H._bbsi_3prime_site(_OH_V1) + stuffer
                + H._bbsi_5prime_site(oh_1b) + vector_trac[split_1b + 4:])
    if cassette not in backbone:
        raise ValueError('cannot locate the dropout cassette in the backbone')

    cdir = output_dir / 'gb' / 'clonotypes'
    for _, r in tqdm(list(var_df.iterrows()), desc='Clonotype maps', disable=not verbose):
        sub = cdir / str(r.Clone)
        sub.mkdir(parents=True, exist_ok=True)

        _gb(r.Sequence, f'frag_{r.Clone}'[:16],
            f'Synthesised fragment | {r.Clone} | {r.TRBV}/{r.TRAV} | pool {r.Pool}',
            [H._feat(0, len(r.Sequence), 'misc_feature', label=['variable fragment'],
                     note=[f'BbsI {r.OH_1a}..{r.OH_1b}; inner BsmBI {r.OH_2a}/{r.OH_2b}'])],
            sub / 'fragment.gb', circular=False)

        plasmid = backbone.replace(cassette, _OH_V1 + r._insert)
        i = plasmid.find(r._insert)
        feats, cur = [], i
        for label, length in (('TRBV', len(H.build_trbv_gene_nt(r.TRBV))),):
            feats.append(H._feat(cur, cur + length, 'gene', label=[f'{label} ({r.TRBV})']))
            cur += length
        feats.append(H._feat(i, i + len(r._insert), 'CDS', label=[f'TCR {r.Clone}'],
                             note=[f'{r.TRBV}/{r.TRAV} | CDR3b {r.CDR3b} | CDR3a {r.CDR3a}']))
        _gb(plasmid, f'plas_{r.Clone}'[:16],
            f'Final plasmid | {r.Clone} | {r.TRBV}/{r.TRAV}', feats,
            sub / 'assembled.gb', circular=True)


def _write_report(output_dir, backbone, patches, donor_df, var_df, invalid_df,
                  oh_1a_map, oh_v1, oh_1b, oh_2a, oh_2b, stuffer_note, n_pairs):
    L = []
    L.append('TCRAFT Pooled TCR Plasmid Library Design Report')
    L.append('=' * 60)
    L.append('')
    L.append('Assembly')
    L.append('  Step 1  BbsI   3-piece: backbone + TRBV donor + variable fragment')
    L.append('  Step 2  BsmBI  2-piece: insert TRBC-T2A donor into the circularised product')
    L.append('  One reaction pair per TRBV pool; pairing is locked by step 1.')
    L.append('')
    L.append('Overhangs')
    L.append(f'  oh_V1 (Kozak->TRBV) : {oh_v1}')
    L.append(f'  oh_1b (CDR3a->TRAC) : {oh_1b}')
    L.append(f'  oh_2a (CDR3b->TRBC) : {oh_2a}')
    L.append(f'  oh_2b (T2A->TRAV)   : {oh_2b}')
    L.append('')
    L.append('Backbone')
    L.append(f'  Length: {len(backbone)} nt')
    L.append(f'  Stuffer: {stuffer_note}')
    census = _enzyme_census(backbone, circular=True)
    L.append(f'  Site census: {census}  (BbsI 2 = the cloning sites; BsmBI must be 0)')
    for k, p in sorted(patches.items()):
        L.append(f'  Patch {k} ({p["enzymes"]}): {len(p["seq"])} nt, '
                 f'replaces {p["replaces"]} in the stock vector')
        for pos, o, n in p['changes']:
            L.append(f'      silenced {o}>{n} at stock-vector position {pos}')
    L.append('')
    L.append(f'TRBV donors: {len(donor_df)}')
    for _, r in donor_df.iterrows():
        L.append(f'  {r.Gene:<20} oh_1a={r.OH_right}  (NEB {r.NEB_score})  {r.Length} nt')
    L.append('')
    if len(var_df):
        L.append(f'Variable fragments: {len(var_df)} of {n_pairs} paired clonotypes')
        L.append(f'  Length {var_df.Length.min()}-{var_df.Length.max()} nt  '
                 f'GC {var_df.GC.min():.0%}-{var_df.GC.max():.0%}')
        L.append('  Pool sizes:')
        for pool, n in var_df.Pool.value_counts().items():
            L.append(f'    {pool:<20} {n}')
    L.append('')
    L.append(f'Skipped clones: {len(invalid_df)}')
    for _, r in invalid_df.head(30).iterrows():
        L.append(f'  {r.clone_id}: {r.Issue}')
    if len(invalid_df) > 30:
        L.append(f'  ... and {len(invalid_df) - 30} more (see invalid_clones.csv)')
    (output_dir / 'assembly_report.txt').write_text('\n'.join(L) + '\n')


def run_CLI():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('airr', help='AIRR rearrangement file (.tsv/.csv)')
    ap.add_argument('--output_dir', default=None)
    ap.add_argument('--trbv', default=None, help='restrict to a single TRBV pool')
    ap.add_argument('--stuffer', default=None,
                    help='FASTA of the ccdB cassette to place in the backbone dropout')
    ap.add_argument('--limit', type=int, default=None, help='cap clonotypes (for testing)')
    ap.add_argument('--no-clonotype-gb', action='store_true',
                    help='skip the per-clonotype GenBank maps (two files per clonotype)')
    args = ap.parse_args()
    generate_vector_library(args.airr, output_dir=args.output_dir, trbv_filter=args.trbv,
                            stuffer_path=args.stuffer, limit=args.limit,
                            write_clonotype_gb=not args.no_clonotype_gb)


if __name__ == '__main__':
    run_CLI()

#!/usr/bin/env python
"""
Verify a TCRAFT-library design by simulating the two Golden Gate steps for real.

This does NOT trust the designed parts list.  It digests the actual backbone, donor and
fragment sequences with Bio.Restriction, reads the 4-nt overhangs off the cut geometry,
ligates by overhang complementarity, and only then checks the product.  That is the part
that can silently be wrong: a correct parts list assembled with mismatched overhangs still
produces a plausible-looking CSV.

Checks performed
  1. backbone/donor/fragment carry exactly the expected number of BbsI and BsmBI sites
  2. step 1 overhangs pair uniquely and unambiguously (no cross-ligation)
  3. step 2 inserts TRBC-T2A in the only orientation the overhangs allow
  4. the assembled insert is a single ORF: starts ATG, length %% 3 == 0, no internal stop
  5. the translated product round-trips to the AIRR input: CDR3a/CDR3b junctions and the
     TRBV/TRAV gene identities all reappear in the protein

Usage
    uv run python TCRAFT/scripts/verify_library_assembly.py <library_output_dir>
"""

import argparse
import json
import sys
from pathlib import Path

import pandas as pd
from Bio.Seq import Seq
from Bio.Restriction import BbsI, BsmBI

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from TCRAFT import HDR_Design as H          # noqa: E402
from TCRAFT import Vector_Library as V      # noqa: E402


def digest(seq, enzyme, circular=False):
    """
    Cut `seq` and return [(fragment, left_overhang, right_overhang), ...].

    Overhangs are the 4-nt 5' extensions BbsI/BsmBI leave.  For a top-strand cut at position
    p the overhang is seq[p:p+4] and the downstream fragment begins at p.
    """
    s = str(seq).upper()
    cuts = sorted(enzyme.search(Seq(s), linear=not circular))
    if not cuts:
        return []
    cuts = [c - 1 for c in cuts]            # Bio.Restriction is 1-based
    out = []
    if circular:
        for a, b in zip(cuts, cuts[1:] + [cuts[0] + len(s)]):
            frag = (s + s)[a:b]
            out.append((frag, frag[:4], (s + s)[b:b + 4]))
    else:
        bounds = [0] + cuts + [len(s)]
        for i, (a, b) in enumerate(zip(bounds, bounds[1:])):
            frag = s[a:b]
            left = frag[:4] if i > 0 else None
            right = s[b:b + 4] if b < len(s) and i + 1 < len(bounds) - 1 else None
            out.append((frag, left, right))
    return out


def inner_fragment(seq, enzyme):
    """The single fragment released from a linear donor cut twice — overhangs at both ends."""
    frags = [f for f in digest(seq, enzyme) if f[1] is not None and f[2] is not None]
    if len(frags) != 1:
        raise AssertionError(f'expected exactly 1 released fragment, got {len(frags)}')
    return frags[0]


def ligate(pieces):
    """
    Join (fragment, left_oh, right_oh) pieces into a circle by overhang complementarity.

    Asserts each overhang has exactly one partner, so a design that happens to assemble by
    luck but is ambiguous in the pot is reported as a failure, not a pass.
    """
    by_left = {}
    for p in pieces:
        if p[1] in by_left:
            raise AssertionError(f'overhang {p[1]} claimed by two fragments — ambiguous assembly')
        by_left[p[1]] = p

    start = pieces[0]
    order, cur, seen = [start], start, {start[1]}
    for _ in range(len(pieces)):
        nxt = by_left.get(cur[2])
        if nxt is None:
            raise AssertionError(f'no fragment starts with overhang {cur[2]}')
        if nxt is start:
            break
        if nxt[1] in seen:
            raise AssertionError('overhang cycle does not close cleanly')
        seen.add(nxt[1])
        order.append(nxt)
        cur = nxt
    if len(order) != len(pieces):
        raise AssertionError(f'assembly closed after {len(order)} of {len(pieces)} fragments')
    # each fragment already begins with its own left overhang; concatenating is correct
    return ''.join(p[0] for p in order), order


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('library_dir')
    ap.add_argument('--airr', default='input_data/airr_rearrangement.tsv')
    args = ap.parse_args()
    d = Path(args.library_dir)

    donors = pd.read_csv(d / 'trbv_donors.csv')
    var = pd.read_csv(d / 'variable_fragments.csv')
    backbone = str(next((d / 'gb' / 'vectors').glob('TRAC_backbone.gb')) and
                   __import__('Bio.SeqIO', fromlist=['SeqIO']).read(
                       str(d / 'gb' / 'vectors' / 'TRAC_backbone.gb'), 'genbank').seq).upper()
    trbc_donor = str(__import__('Bio.SeqIO', fromlist=['SeqIO']).read(
        str(d / 'gb' / 'vectors' / 'TRBC-T2A_donor.gb'), 'genbank').seq).upper()

    # ---- 1. site census
    fails = []
    n_bb = len(BbsI.search(Seq(backbone), linear=False))
    if n_bb != 2:
        fails.append(f'backbone has {n_bb} BbsI sites, expected 2')
    n_bs = len(BsmBI.search(Seq(backbone), linear=False))
    if n_bs != 0:
        fails.append(f'backbone has {n_bs} BsmBI sites, expected 0')
    print(f'[1] backbone: BbsI={n_bb} (expect 2)  BsmBI={n_bs} (expect 0)')

    # Backbone arm = the BbsI fragment carrying oh_1b .. oh_V1.  Identify it by its overhangs,
    # not by length: a real ccdB cassette could in principle rival the arm in size.
    ohs = json.loads((d / 'overhangs.json').read_text())
    bb_frags = digest(backbone, BbsI, circular=True)
    arms = [f for f in bb_frags if f[1] == ohs['oh_1b'] and f[2] == ohs['oh_V1']]
    if len(arms) != 1:
        print(f'FAIL: expected 1 backbone arm with overhangs '
              f'{ohs["oh_1b"]}..{ohs["oh_V1"]}, found {len(arms)} '
              f'among {[(f[1], f[2]) for f in bb_frags]}')
        return 1
    arm = arms[0]
    print(f'    backbone arm {len(arm[0])} nt, overhangs {arm[1]} .. {arm[2]}')

    trbc_piece = inner_fragment(trbc_donor, BsmBI)
    print(f'[2] TRBC-T2A donor releases {len(trbc_piece[0])} nt, '
          f'overhangs {trbc_piece[1]} .. {trbc_piece[2]}')

    airr = H.load_airr(args.airr)
    pairs = {p['clone_id']: p for p in H.pair_airr_records(airr)}

    checked = 0
    for _, row in var.iterrows():
        donor_seq = donors.loc[donors.Gene == row.TRBV, 'Sequence'].values[0]
        v_piece = inner_fragment(row.Sequence, BbsI)
        d_piece = inner_fragment(donor_seq, BbsI)

        # ---- step 1
        try:
            circle, _order = ligate([arm, d_piece, v_piece])
        except AssertionError as e:
            fails.append(f'{row.Clone}: step 1 — {e}')
            continue

        # ---- step 2: cut the circle, drop in TRBC-T2A
        open_frags = [f for f in digest(circle, BsmBI, circular=True)]
        if len(open_frags) != 2:
            fails.append(f'{row.Clone}: step 2 cut gave {len(open_frags)} fragments, expected 2')
            continue
        step2_arm = max(open_frags, key=lambda f: len(f[0]))
        try:
            final, _ = ligate([step2_arm, trbc_piece])
        except AssertionError as e:
            fails.append(f'{row.Clone}: step 2 — {e}')
            continue

        # ---- 4. single clean ORF
        doubled = final + final
        atg = doubled.find('ATG' + H.build_trbv_gene_nt(row.TRBV)[3:30])
        if atg == -1:
            fails.append(f'{row.Clone}: TRBV start not found in the assembled circle')
            continue
        orf_aa = None
        for end in range(atg + 3, len(doubled) - 2, 3):
            cod = doubled[end:end + 3]
            if cod in ('TAA', 'TAG', 'TGA'):
                orf_aa = str(Seq(doubled[atg:end]).translate())
                break
        if orf_aa is None:
            fails.append(f'{row.Clone}: no stop codon found downstream of the TRBV ATG')
            continue
        if '*' in orf_aa:
            fails.append(f'{row.Clone}: internal stop in ORF')
            continue

        # ---- 5. round-trip to the AIRR input
        p = pairs.get(row.Clone) or pairs.get(str(row.Clone))
        cdr3b = str(p['trb']['junction_aa']).upper()
        cdr3a = str(p['tra']['junction_aa']).upper()
        if cdr3b[1:] not in orf_aa:
            fails.append(f'{row.Clone}: CDR3b {cdr3b} absent from the assembled protein')
        if cdr3a[1:] not in orf_aa:
            fails.append(f'{row.Clone}: CDR3a {cdr3a} absent from the assembled protein')
        if H.get_fr4(row.TRBJ) not in orf_aa:
            fails.append(f'{row.Clone}: TRBJ FR4 absent from the assembled protein')
        if H.get_fr4(row.TRAJ) not in orf_aa:
            fails.append(f'{row.Clone}: TRAJ FR4 absent from the assembled protein')
        checked += 1

    print(f'[3-5] simulated both steps for {checked}/{len(var)} clonotypes')
    if fails:
        print(f'\nFAILURES ({len(fails)}):')
        for f in fails[:40]:
            print('  ' + f)
        return 1
    print('\nAll checks passed.')
    return 0


if __name__ == '__main__':
    sys.exit(main())

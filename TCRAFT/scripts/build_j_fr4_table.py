#!/usr/bin/env python
"""
Build the mouse J-gene FR4 reference table.

WHY THIS EXISTS
---------------
IMGT's `junction_aa` runs from the conserved Cys104 to the conserved F/W at position 118
*inclusive*.  The J segment continues past 118 with FR4 — the final beta-strand of the V
domain, which physically connects it to the constant region.  Any design that joins
`junction_aa` straight to the constant region silently drops FR4 (9 aa on beta, 10-11 on
alpha) and will not fold correctly.

To append FR4 we need, for each J gene, the index of the position-118 residue within that
gene's amino-acid sequence.  A motif regex is not reliable: the canonical FGXG motif is
degenerate in the ORF genes (TRAJ7 `LGKG`, TRAJ44 `LGAG`, TRBJ1-6 `FAAG`), which is exactly
why IMGT classifies them as ORF rather than F.

METHOD
------
Primary (data-driven): for every productive rearrangement, find the longest suffix of
`junction_aa` that is a prefix of the called J gene.  Recombination trims the J's 5' end and
adds N nucleotides, so this only succeeds for a minority of clones -- but when it succeeds it
pins the anchor exactly.  Across a real repertoire every clone of a given J gene agrees, so
the modal anchor is taken as that gene's canonical value and unanimity is asserted.

Fallback (motif): for J genes with no resolvable clone in the input, take the last match of
`[FWL]G.G` and validate the resulting FR4 length against the per-locus expectation.

Genes are loaded WITHOUT the `|F|` functional filter, so ORF genes (TRBJ1-6, TRAJ7, TRAJ44 --
~350 rows in the reference AIRR) are included rather than silently dropped.

USAGE
-----
    uv run python TCRAFT/scripts/build_j_fr4_table.py input_data/airr_rearrangement.tsv
"""

import argparse
import collections
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from TCRAFT import HDR_Design as H  # noqa: E402


OUT_PATH = Path(__file__).resolve().parents[1] / 'references' / 'mouse' / 'mouse_J_FR4.csv'

# FR4 lengths observed across the mouse IMGT J sets.
_EXPECTED_FR4_LEN = {'TRBJ': (9,), 'TRAJ': (10, 11)}
# The single length used by the fallback when a gene has no clone support.
_CANONICAL_FR4_LEN = {'TRBJ': 9, 'TRAJ': 10}


def _derive_anchors_from_data(sub, j_dict):
    """{gene: Counter(anchor_index)} from longest-suffix/prefix overlap."""
    acc = collections.defaultdict(collections.Counter)
    for junction, j_call in zip(sub['junction_aa'], sub['j_call']):
        gene = str(j_call).split(',')[0].split('*')[0]
        if gene not in j_dict:
            continue
        junction, j_aa = str(junction), j_dict[gene]
        for k in range(min(len(junction), len(j_aa)), 2, -1):
            if junction[-k:] == j_aa[:k]:
                acc[gene][k - 1] += 1   # index in j_aa of the junction's final residue
                break
    return acc


def _anchor_from_length(j_aa, locus):
    """
    Fallback for J genes with no resolvable clone: FR4 has a constant length per locus,
    so the anchor sits that far from the 3' end.  Validated against the data-derived
    anchors -- it reproduces 51 of 52, the lone exception (TRAJ35, a genuine 11-aa FR4)
    being data-derived and so never reaching this path.
    """
    anchor = len(j_aa) - _CANONICAL_FR4_LEN[locus] - 1
    return anchor if anchor >= 0 else None


def _looks_like_fr4(fr4):
    """
    FR4 begins with the J's conserved G-x-G-[TSA] turn into the constant region.
    Pseudogene and frame-shifted entries fail this, which is how low-confidence rows
    from the length fallback are caught (e.g. TRBJ2-6 -> 'WQPIEQPMR').
    """
    return re.match(r'^[GA].G[TSA]', fr4) is not None


def build_table(airr_path):
    df = H.load_airr(airr_path)
    df = df[df['productive'].astype(str).str.upper().isin(['T', 'TRUE'])]

    rows, problems = [], []
    for locus, paths, prefix in (('TRBJ', H._PATHS['trbj_aa'], 'TRB'),
                                 ('TRAJ', H._PATHS['traj_aa'], 'TRA')):
        # filter_functional=False so ORF genes are kept
        j_dict = H._load_fasta_dict(paths, filter_functional=False)
        sub = df[df['v_call'].astype(str).str.startswith(prefix)].dropna(
            subset=['junction_aa', 'j_call'])
        acc = _derive_anchors_from_data(sub, j_dict)

        for gene in sorted(j_dict):
            j_aa = j_dict[gene]
            if gene in acc:
                counter = acc[gene]
                if len(counter) > 1:
                    problems.append(f'{gene}: ambiguous anchors {dict(counter)}')
                anchor, source, support = counter.most_common(1)[0][0], 'data', sum(counter.values())
            else:
                anchor, source, support = _anchor_from_length(j_aa, locus), 'length', 0
                if anchor is None:
                    problems.append(f'{gene}: J sequence shorter than FR4 -> {j_aa!r}')
                    continue

            fr4 = j_aa[anchor + 1:]
            if len(fr4) not in _EXPECTED_FR4_LEN[locus]:
                problems.append(
                    f'{gene}: FR4 length {len(fr4)} outside expected '
                    f'{_EXPECTED_FR4_LEN[locus]} (source={source}) -> {fr4!r}')
            rows.append({'locus': locus, 'J_gene': gene, 'anchor_index': anchor,
                         'anchor_residue': j_aa[anchor], 'FR4_aa': fr4,
                         'FR4_len': len(fr4), 'J_aa': j_aa,
                         'source': source, 'clone_support': support,
                         'usable': _looks_like_fr4(fr4)})

    return pd.DataFrame(rows), problems


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('airr', help='AIRR rearrangement file used to derive anchors')
    ap.add_argument('--out', default=str(OUT_PATH))
    args = ap.parse_args()

    table, problems = build_table(args.airr)
    table = table.sort_values(['locus', 'J_gene']).reset_index(drop=True)
    table.to_csv(args.out, index=False)

    print(f'Wrote {len(table)} J genes to {args.out}')
    for locus, grp in table.groupby('locus'):
        by_src = grp['source'].value_counts().to_dict()
        print(f'  {locus}: {len(grp)} genes  sources={by_src}  usable={int(grp.usable.sum())}  '
              f'FR4 lengths={sorted(grp[grp.usable].FR4_len.unique())}')
    noncanon = table[table.usable & ~table.anchor_residue.isin(['F', 'W'])]
    if len(noncanon):
        print('\nNon-canonical anchor residues (expected for IMGT ORF genes):')
        for _, r in noncanon.iterrows():
            print(f'  {r.J_gene:10s} anchor={r.anchor_residue} FR4={r.FR4_aa}  ({r.source})')
    unusable = table[~table.usable]
    if len(unusable):
        print('\nMarked unusable (FR4 lacks the conserved G.G[TSA] turn; pseudogene/frame-shifted).')
        print('Design modules must skip these and report affected clones:')
        for _, r in unusable.iterrows():
            print(f'  {r.J_gene:10s} FR4={r.FR4_aa!r}  J={r.J_aa}')
    if problems:
        print('\nPROBLEMS:')
        for p in problems:
            print('  ' + p)
        return 1
    print('\nNo problems.')
    return 0


if __name__ == '__main__':
    sys.exit(main())

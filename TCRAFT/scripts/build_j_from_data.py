#!/usr/bin/env python
"""
Derive the J-segment protein used after each CDR3 directly from a paired clonotype table.

For each clone, the J contribution is the stretch of {chain}_sequence_aa between the end of
junction_aa (the conserved F/W) and the start of the constant region motif:
    TRAJ: ... FR4 + <split-codon residue> | I-Q-N...    (residue depends on the J gene: N/D/Y/H)
    TRBJ: ... FR4                         | E-D-L-R-N...
The modal string per J gene is written to references/mouse/mouse_J_from_data.csv. This mirrors the
human design, where the residue between TRAJ and TRAC is carried by the TRAJ sequence and TRAC starts at IQN.
(Beta clones whose contig continues into the unspliced TRBC1 extra-exon stretch carry no EDLRN marker
and are ignored automatically.)

usage: uv run python TCRAFT/scripts/build_j_from_data.py clonotypes.csv
"""
import sys
from pathlib import Path

import pandas as pd

REF = Path(__file__).resolve().parents[1] / "references" / "mouse"
MARK = {"tra": ("TRAJ", "IQN"), "trb": ("TRBJ", "EDLRN")}
d = pd.read_csv(sys.argv[1], low_memory=False,
                usecols=[f"{c}_{k}" for c in MARK for k in ("j_call", "junction_aa", "sequence_aa")]).dropna()
rows = []
for ch, (locus, marker) in MARK.items():
    jc = d[f"{ch}_j_call"].str.split("*").str[0]
    for g, idx in jc.groupby(jc).groups.items():
        seen = []
        for s, jx in zip(d.loc[idx, f"{ch}_sequence_aa"], d.loc[idx, f"{ch}_junction_aa"]):
            p = s.find(jx) + len(jx)
            q = s.find(marker, p)
            if 0 <= q - p <= 24:
                seen.append(s[p:q])
        if seen:
            top = pd.Series(seen).value_counts()
            rows.append(dict(locus=locus, J_gene=g, seq=top.index[0], n=len(seen), n_total=len(idx),
                             frac=round(top.iloc[0] / len(seen), 4), runner_up=top.index[1] if len(top) > 1 else ""))
out = pd.DataFrame(rows).sort_values(["locus", "J_gene"])
out.to_csv(REF / "mouse_J_from_data.csv", index=False)
print(f"{len(out)} J genes ({out.locus.value_counts().to_dict()}); min support n={out.n.min()}; min consistency={out.frac.min()}")
print(out[out.frac < 0.95].to_string(index=False))
old = pd.read_csv(REF / "mouse_J_FR4.csv").set_index("J_gene")
changed = []
for _, r in out.iterrows():
    if r.J_gene in old.index:
        o = old.loc[r.J_gene]
        exp = r.seq if r.locus == "TRBJ" else r.seq[:-1]  # strip the TRAJ split-codon residue to compare with the old FR4
        if o.FR4_aa != exp:
            changed.append((r.J_gene, o.FR4_aa, exp, o.source, int(o.clone_support)))
    else:
        changed.append((r.J_gene, "(absent)", r.seq, "", 0))
print("FR4 differs from the old mouse_J_FR4.csv for", len(changed), "genes:")
for c in changed: print("  ", c)
print("in old table, no data:", sorted(set(old.index) - set(out.J_gene)))

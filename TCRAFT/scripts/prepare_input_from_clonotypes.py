#!/usr/bin/env python
"""
Convert a paired clonotype table (tra_*/trb_* columns, e.g. HAM-07-007-012) into the CSV that
TCRAFT-generate expects: V_alpha, V_beta, J_alpha, J_beta, CDR3_alpha, CDR3_beta.
Gene names are passed through unchanged (alleles, if any, stripped). Rows are written in
file order, so output row i corresponds to input row i.

usage: uv run python TCRAFT/scripts/prepare_input_from_clonotypes.py IN.csv OUT.csv [--sample N --seed S]
"""
import argparse
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("infile"); ap.add_argument("outfile")
ap.add_argument("--sample", type=int, default=None); ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args()

cols = ["clone_id", "clone_size", "tra_v_call", "tra_j_call", "tra_junction_aa", "trb_v_call", "trb_j_call", "trb_junction_aa"]
d = pd.read_csv(a.infile, usecols=cols, low_memory=False).dropna()
if a.sample:
    d = d.sample(a.sample, random_state=a.seed).sort_index()
strip = lambda s: s.astype(str).str.split("*").str[0]
out = pd.DataFrame({
    "V_alpha": strip(d.tra_v_call), "V_beta": strip(d.trb_v_call),
    "J_alpha": strip(d.tra_j_call), "J_beta": strip(d.trb_j_call),
    "CDR3_alpha": d.tra_junction_aa, "CDR3_beta": d.trb_junction_aa,
    "clone_id": d.clone_id, "clone_size": d.clone_size,
})
out.to_csv(a.outfile, index=False)
print(f"wrote {a.outfile}: {len(out)} TCRs")

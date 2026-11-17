#!/usr/bin/env python
"""
Design the final overhang sets for the mouse TCRAFT library and save them to
references/mouse/overhang_sets.json.

  Step 1 (BbsI):        TRBV vectors  - overhangs (TRBV 3' window, TRAC N-term window)
  Step 2 (BsmBI/BsaI):  TRAV vectors  - overhangs (TRBC1 N-term window, TRAV 3' window)

Windows copy the human ranges: constants = first 10 codons (30 nt); TRBV = last 33 nt and
TRAV = last 21 nt of the leader->Cys104 block (the Cys104 codon is preserved: recoding is synonymous).
All TRBV share one pool (A); TRAV is split into two contiguous natural-order halves
(A = BsmBI, B = BsaI), mirroring the human A/B design.

usage: uv run python TCRAFT/scripts/design_overhangs.py [--seeds 3]
"""
import argparse
import json
import re
import statistics
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from TCRAFT.overhang_tools import Fidelity, optimize, rc, window_candidates  # noqa: E402

REF = ROOT / "TCRAFT" / "references" / "mouse"
V_WIN = {"TRB": 33, "TRA": 21}
C_WIN = 30
# Overhangs already used by the (human-derived) vector flanks / middle sites; never reuse them.
FIXED = {"TAAT", "AAAC", "GTCA", "GAAT_PLACEHOLDER"} - {"GAAT_PLACEHOLDER"}
FIXED_ALL = FIXED | {rc(x) for x in FIXED}


def natural_key(s):
    return [int(t) if t.isdigit() else t for t in re.findall(r"\d+|\D+", s)]


def read_fasta_seq(path, key):
    name, seq = None, {}
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:]
            seq[name] = ""
        elif name:
            seq[name] += line
    return [s for n, s in seq.items() if key in n][0].lower()


def backtranslate(aa, usage):
    best = {a: max(c, key=lambda t: t[1])[0].lower() for a, c in usage.items()}
    return "".join(best[a] for a in aa)


def usage_score(span, freq):
    return sum(freq.get(span[i:i + 3], 0) for i in range(0, len(span), 3))


def pick(cands, oh, freq, prefer):
    """Choose one (start, codon_idx, span) for overhang `oh`. prefer: 'early' or 'late' start."""
    opts = cands[oh]
    key = (lambda o: (o[0], -usage_score(o[2], freq))) if prefer == "early" else (lambda o: (-o[0], -usage_score(o[2], freq)))
    return min(opts, key=key)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=3)
    args = ap.parse_args()

    usage = json.load(open(REF / "mouse_codon_usage.json"))
    freq = {c.lower(): f for lst in usage.values() for c, f in lst}
    cons = pd.read_csv(REF / "mouse_V_consensus.csv")
    trac = read_fasta_seq(REF / "mouse_constants_cys_nt.fa", "TRAC_Cys")
    trbc = read_fasta_seq(REF / "mouse_constants_cys_nt.fa", "TRBC1_Cys")

    def filt(d):
        return {o: v for o, v in d.items() if o not in FIXED_ALL}

    trac_d = filt(window_candidates(trac, 0, C_WIN))
    trbc_d = filt(window_candidates(trbc, 0, C_WIN))
    trac_c, trbc_c = list(trac_d), list(trbc_d)
    print(f"constant windows: TRAC {len(trac_c)} candidates, TRBC1 {len(trbc_c)}")

    def v_block(chain):
        g = cons[cons.chain == chain]
        groups = {}
        for _, r in g.iterrows():
            groups.setdefault(r.vector_group, {"aa": r.v_aa, "genes": []})["genes"].append(r.gene)
        names = sorted(groups, key=natural_key)
        dicts = []
        for n in names:
            nt = backtranslate(groups[n]["aa"], usage)
            dicts.append(filt(window_candidates(nt, max(0, len(nt) - V_WIN[chain]), len(nt))))
        return names, groups, dicts

    trbv_names, trbv_groups, trbv_d = v_block("TRB")
    trav_names, trav_groups, trav_d = v_block("TRA")
    for lab, names, d in (("TRBV", trbv_names, trbv_d), ("TRAV", trav_names, trav_d)):
        k = [len(x) for x in d]
        print(f"{lab}: {len(names)} vector groups | candidates/gene min {min(k)} median {statistics.median(k)} max {max(k)}")

    result = {"TRBV": {}, "TRAV": {}, "members": {}, "scores": {}}

    def best_of(fid, lc, rc_, tag):
        best = None
        for sd in range(args.seeds):
            L, R, sc = optimize(fid, lc, rc_, seed=100 + sd, patience=1000)
            print(f"   {tag} seed {sd}: f1={sc['f1']:.4f} mean_p={sc['p_vec'].mean():.4f}")
            if best is None or sc["f"] > best[2]["f"]:
                best = (L, R, sc)
        return best

    # ---- Step 1: TRBV (single pool A), BbsI ----
    print("\nStep 1 TRBV-TRAC (BbsI), one pool:")
    L, R, sc = best_of(Fidelity("BbsI"), [list(d) for d in trbv_d], [trac_c] * len(trbv_names), "TRBV")
    result["scores"]["TRBV_A"] = dict(enzyme="BbsI", n=len(L), f1=sc["f1"], mean_p=float(sc["p_vec"].mean()))
    for i, n in enumerate(trbv_names):
        ls, lc0, lspan = pick(trbv_d[i], L[i], freq, "late")
        rs, rc0, rspan = pick(trac_d, R[i], freq, "early")
        result["TRBV"][n] = dict(set="A", left=L[i], left_info=[ls, lc0, lspan], right=R[i], right_info=[rs, rc0, rspan])
        result["members"][n] = trbv_groups[n]["genes"]

    # ---- Step 2: TRAV, two contiguous natural-order halves ----
    half = (len(trav_names) + 1) // 2
    pools = {"A": (list(range(0, half)), "BsmBI"), "B": (list(range(half, len(trav_names))), "BsaI")}
    for pool, (idx, enz) in pools.items():
        print(f"\nStep 2 TRBC1-TRAV pool {pool} ({enz}), n={len(idx)}: {trav_names[idx[0]]} .. {trav_names[idx[-1]]}")
        L, R, sc = best_of(Fidelity(enz), [trbc_c] * len(idx), [list(trav_d[i]) for i in idx], f"TRAV-{pool}")
        result["scores"][f"TRAV_{pool}"] = dict(enzyme=enz, n=len(idx), f1=sc["f1"], mean_p=float(sc["p_vec"].mean()))
        for j, i in enumerate(idx):
            n = trav_names[i]
            ls, lc0, lspan = pick(trbc_d, L[j], freq, "early")
            rs, rc0, rspan = pick(trav_d[i], R[j], freq, "late")
            result["TRAV"][n] = dict(set=pool, left=L[j], left_info=[ls, lc0, lspan], right=R[j], right_info=[rs, rc0, rspan])
            result["members"][n] = trav_groups[n]["genes"]

    result["v_aa"] = {n: g["aa"] for grp in (trbv_groups, trav_groups) for n, g in grp.items()}
    out = REF / "overhang_sets.json"
    json.dump(result, open(out, "w"), indent=1)
    print("\nwrote", out)
    print(json.dumps(result["scores"], indent=1))


if __name__ == "__main__":
    main()

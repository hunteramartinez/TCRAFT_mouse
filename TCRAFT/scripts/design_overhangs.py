#!/usr/bin/env python
"""
Design overhang sets for the mouse TCRAFT library.

  Step 1 (BbsI):   TRBV vectors,  overhangs = (TRBV 3' window, TRAC N-term window)
  Step 2 (BsmBI/BsaI): TRAV vectors, overhangs = (TRBC1 N-term window, TRAV 3' window)

Windows copy the human ranges: constants = first 10 codons (30 nt); TRBV = last 33 nt,
TRAV = last 21 nt of the leader->Cys104 block (Cys104 codon is preserved by construction).
usage: uv run python TCRAFT/scripts/design_overhangs.py [--pools 1 2 3 4]
"""
import argparse
import json
import statistics
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from TCRAFT.overhang_tools import Fidelity, optimize, window_candidates  # noqa: E402

REF = ROOT / "TCRAFT" / "references" / "mouse"
V_WIN = {"TRB": 33, "TRA": 21}
C_WIN = 30


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pools", type=int, nargs="+", default=[1, 2, 3, 4])
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    usage = json.load(open(REF / "mouse_codon_usage.json"))
    cons = pd.read_csv(REF / "mouse_V_consensus.csv")
    trac = read_fasta_seq(REF / "mouse_constants_cys_nt.fa", "TRAC_Cys")
    trbc = read_fasta_seq(REF / "mouse_constants_cys_nt.fa", "TRBC1_Cys")
    trac_c = list(window_candidates(trac, 0, C_WIN))
    trbc_c = list(window_candidates(trbc, 0, C_WIN))
    print(f"constant-window candidates (palindromes + site-creating recodings removed): "
          f"TRAC {len(trac_c)}, TRBC1 {len(trbc_c)}")

    def v_cands(chain):
        g = cons[cons.chain == chain].drop_duplicates("vector_group").sort_values("vector_group")
        names, cands = [], []
        for _, r in g.iterrows():
            nt = backtranslate(r.v_aa, usage)
            c = list(window_candidates(nt, max(0, len(nt) - V_WIN[chain]), len(nt)))
            names.append(r.vector_group)
            cands.append(c)
        return names, cands

    trbv_names, trbv_c = v_cands("TRB")
    trav_names, trav_c = v_cands("TRA")
    for lab, names, c in (("TRBV", trbv_names, trbv_c), ("TRAV", trav_names, trav_c)):
        k = [len(x) for x in c]
        print(f"{lab}: {len(names)} vector groups | candidates/gene min {min(k)} median {statistics.median(k)} max {max(k)}"
              f" | genes with <10: {[n for n, x in zip(names, k) if x < 10]}")

    # ---- Step 1: all TRBV in one pool, BbsI ----
    f_bbs = Fidelity("BbsI")
    L, R, sc = optimize(f_bbs, trbv_c, [trac_c] * len(trbv_c), seed=args.seed, patience=600,
                        stage1=(4, 8000, 1e-4), stage2=(4, 4000, 1e-3))
    print(f"\nStep 1 TRBV-TRAC, 1 pool of {len(L)} (BbsI): f1={sc['f1']:.3f} mean_p={sc['p_vec'].mean():.3f} "
          f"distinct L/R={len(set(L))}/{len(set(R))}")

    # ---- Step 2: TRAV split into P pools (alternating BsmBI / BsaI tables) ----
    fids = {"BsmBI": Fidelity("BsmBI"), "BsaI": Fidelity("BsaI")}
    for P in args.pools:
        print(f"\nStep 2 TRBC1-TRAV with {P} pool(s) (balanced, gene order):")
        n = len(trav_names)
        for k in range(P):
            idx = list(range(k, n, P))  # interleave so each pool spans all families
            enz = "BsmBI" if k % 2 == 0 else "BsaI"
            L, R, sc = optimize(fids[enz], [trbc_c] * len(idx), [trav_c[i] for i in idx], seed=args.seed + k,
                                patience=600, stage1=(4, 8000, 1e-4), stage2=(4, 4000, 1e-3))
            print(f"  pool {k+1} ({enz}) n={len(idx):2d}: f1={sc['f1']:.3f} mean_p={sc['p_vec'].mean():.3f} "
                  f"distinct L/R={len(set(L))}/{len(set(R))}")


if __name__ == "__main__":
    main()

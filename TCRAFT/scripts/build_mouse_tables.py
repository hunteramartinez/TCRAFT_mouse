#!/usr/bin/env python
"""
Build the mouse equivalents of every human TCRAFT reference table from overhang_sets.json:

  TRBV_TRAC_OHset.csv     TRBC_TRAV_OHset.csv       overhang pair + pool for each V gene
  TRBV_TRAC_OH_frags.csv  TRBC_TRAV_OH_frags.csv    oligo-side / vector-side fragments
  TRBV_vector_inserts.csv TRAV_vector_inserts.csv   gene blocks for the Step-1 / Step-2 vectors
  TRJ_seqs.csv                                       J-segment FR4 protein (after the CDR3-ending F/W)
  Oligo_Components.json, P2A.txt                     copied unchanged from the human design

Schemas match the human files so Generate.py needs only path / codon-table changes.
Overhangs are UPPERCASE, everything else lowercase (human convention).

The flanks around the inserts (XhoI/NotI ends, BsaI/BsmBI/BbsI sites, TAAT/AAAC, the stuffer cassette)
are copied verbatim from the human design as PLACEHOLDERS: the cloning vectors for the mouse
library are still to be designed.

usage: uv run python TCRAFT/scripts/build_mouse_tables.py
"""
import json
import shutil
import sys
from pathlib import Path

import pandas as pd
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from TCRAFT.overhang_tools import SYN, clean_nt, count_problems  # noqa: E402

REF = ROOT / "TCRAFT" / "references"
MREF = REF / "mouse"
C_WIN = 30

# --- placeholders copied from the human vector inserts (see docstring) ---
TRAV_HEAD = "CTCGAGTAATaGAGACCCGTCTCa"
TRAV_TAIL = "gGAGACGCGGTCTCaAAACCTCGAGGCGGCCGC"
TRBV_CASSETTE = "taGTCTTCGTCAcGAGACGTAATgGAGACCatgcGATATCcgtaGGTCTCgAAACCGTCTCgGGGTGAAGACat"


def fasta_seq(path, key):
    cur, seqs = None, {}
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            cur = line[1:]
            seqs[cur] = ""
        elif cur:
            seqs[cur] += line
    return [s for n, s in seqs.items() if key in n][0].lower()


def apply_span(nt, c0, span):
    return nt[: c0 * 3] + span + nt[c0 * 3 + len(span):]


def with_oh_upper(nt, start):
    """lowercase everything, uppercase the 4-nt overhang at [start, start+4)."""
    nt = nt.lower()
    return nt[:start] + nt[start:start + 4].upper() + nt[start + 4:]


def design_v(aa, oh, info, usage):
    start, c0, span = info
    best = {a: max(c, key=lambda t: t[1])[0].lower() for a, c in usage.items()}
    nt = apply_span("".join(best[a] for a in aa), c0, span)
    nt = clean_nt(nt, usage, protect=range(c0, c0 + len(span) // 3))
    assert str(Seq(nt).translate()) == aa, "V translation changed"
    assert nt[start:start + 4].upper() == oh, "overhang lost during cleaning"
    assert count_problems(nt) == 0
    return nt


def main():
    usage = json.load(open(MREF / "mouse_codon_usage.json"))
    sets = json.load(open(MREF / "overhang_sets.json"))
    trac = fasta_seq(MREF / "mouse_constants_cys_nt.fa", "TRAC_Cys")
    trbc = fasta_seq(MREF / "mouse_constants_cys_nt.fa", "TRBC1_Cys")
    p2a = open(REF / "P2A.txt").read().strip()
    assert str(Seq(p2a).translate()).startswith("GSG"), "P2A must start with the GSG linker"

    # ------------------------------------------------------------------ TRBV / TRAC
    ohset_b, frags_b, ins_b = [], [], []
    n = 0
    for grp, d in sets["TRBV"].items():
        vnt = design_v(sets["v_aa"][grp], d["left"], d["left_info"], usage)
        ls = d["left_info"][0]
        rs, rc0, rspan = d["right_info"]
        trac_v = apply_span(trac, rc0, rspan)
        assert str(Seq(trac_v).translate()) == str(Seq(trac).translate())
        assert trac_v[rs:rs + 4].upper() == d["right"]
        for gene in sets["members"][grp]:
            ohset_b.append(dict(N=n, TRBV_name=gene, TRBV=d["left"], TRAC=d["right"], set=d["set"]))
            n += 1
            left = with_oh_upper(vnt[: ls + 4], ls)
            frags_b.append(dict(
                trbv=gene,
                trbv_left_frag=left,
                trbv_right_frag=with_oh_upper(vnt[ls:], 0),
                trac_left_frag=with_oh_upper(trac_v[: rs + 4], rs),
                trac_right_frag=with_oh_upper(trac_v[rs:C_WIN], 0),
            ))
            ins_b.append(dict(TRBV=gene, Sequence=left + TRBV_CASSETTE + with_oh_upper(trac_v[rs:C_WIN], 0)))

    # ------------------------------------------------------------------ TRBC / TRAV
    ohset_a, frags_a, ins_a = [], [], []
    n = 0
    for grp, d in sets["TRAV"].items():
        vnt = design_v(sets["v_aa"][grp], d["right"], d["right_info"], usage)
        rs = d["right_info"][0]
        ls, lc0, lspan = d["left_info"]
        trbc_v = apply_span(trbc, lc0, lspan)
        assert str(Seq(trbc_v).translate()) == str(Seq(trbc).translate())
        assert trbc_v[ls:ls + 4].upper() == d["left"]
        for gene in sets["members"][grp]:
            ohset_a.append(dict(N=n, TRAV_name=gene, TRBC=d["left"], TRAV=d["right"], set=d["set"]))
            n += 1
            trbc_right = with_oh_upper(trbc_v[ls:], 0)
            trav_left = with_oh_upper(vnt[: rs + 4], rs)
            frags_a.append(dict(
                trav=gene,
                trbc_left_frag=with_oh_upper(trbc_v[: ls + 4], ls),
                trbc_right_frag=trbc_right,
                trav_left_frag=trav_left,
                trav_right_frag=with_oh_upper(vnt[rs:], 0),
            ))
            ins_a.append(dict(TRAV=gene, Sequence=TRAV_HEAD + trbc_right + p2a + trav_left + TRAV_TAIL))

    # ------------------------------------------------------------------ J segments (FR4)
    j = pd.read_csv(MREF / "mouse_J_FR4.csv")
    jj = j[j.usable.astype(str).str.lower() == "true"][["J_gene", "FR4_aa"]].rename(columns={"J_gene": "TRJ", "FR4_aa": "seq"})

    out = {
        "TRBV_TRAC_OHset.csv": pd.DataFrame(ohset_b),
        "TRBC_TRAV_OHset.csv": pd.DataFrame(ohset_a),
        "TRBV_TRAC_OH_frags.csv": pd.DataFrame(frags_b),
        "TRBC_TRAV_OH_frags.csv": pd.DataFrame(frags_a),
        "TRBV_vector_inserts.csv": pd.DataFrame(ins_b),
        "TRAV_vector_inserts.csv": pd.DataFrame(ins_a),
        "TRJ_seqs.csv": jj,
    }
    for name, df in out.items():
        df.to_csv(MREF / name, index=False)
        print(f"{name:28s} {len(df):3d} rows")
    shutil.copy(REF / "Oligo_Components.json", MREF / "Oligo_Components.json")
    shutil.copy(REF / "P2A.txt", MREF / "P2A.txt")

    L = pd.DataFrame(frags_a)
    print("TRAV oligo-side lengths (trbc_left + trav_right):",
          (L.trbc_left_frag.str.len() + L.trav_right_frag.str.len()).describe()[["min", "50%", "max"]].to_dict())
    print("TRAV insert lengths:", pd.DataFrame(ins_a).Sequence.str.len().describe()[["min", "max"]].to_dict())
    print("TRBV insert lengths:", pd.DataFrame(ins_b).Sequence.str.len().describe()[["min", "max"]].to_dict())
    print("J segments:", jj.TRJ.str[:4].value_counts().to_dict())


if __name__ == "__main__":
    main()

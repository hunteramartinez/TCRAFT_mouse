#!/usr/bin/env python
"""
Build nucleotide sequences for the Cys-engineered mouse constant regions used by the
mouse TCRAFT cassette (TRBV-TRBC-P2A-TRAV-TRAC), starting at the first full codon after the
J segment (TRAC: D-I-Q..., TRBC: E-D-L...), and report how many 4-nt overhangs can be written
into the N-terminal window of each constant by synonymous recoding.

Engineering (LVL transmembrane mutation intentionally NOT applied):
  TRAC : T48C   (...DK[T>C]VLDM...)
  TRBC : S57C   (...GV[S>C]TDP...)
Protein output is checked against references/mouse/imgt_genes_aa/IMGT_TRAC_LVL-Cys_aa.fa
(upstream of the TM, where the LVL edits are absent) and IMGT_TRBC1_Cys_aa.fa.

usage: uv run python TCRAFT/scripts/build_mouse_constants.py
"""
from itertools import product
from pathlib import Path

from Bio.Data.CodonTable import standard_dna_table
from Bio.Seq import Seq

REF = Path(__file__).resolve().parents[1] / "references" / "mouse"
FWD = standard_dna_table.forward_table
SYN = {}
for cod, aa in FWD.items():
    SYN.setdefault(aa, []).append(cod.lower())


def fasta(path):
    recs, name, seq = [], None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name:
                recs.append((name, "".join(seq)))
            name, seq = line[1:], []
        else:
            seq.append(line)
    if name:
        recs.append((name, "".join(seq)))
    return recs


def rc(s):
    return str(Seq(s).reverse_complement())


def get_constants():
    # TRAC: IMGT nt starts with a partial first codon 'n' -> G (G-AC = Asp)
    trac = [s for n, s in fasta(REF / "imgt_genes_nt" / "IMGT_TRAC_nt.fa") if "|TRAC*01|" in n][0]
    trac = "g" + trac[1:]
    # TRBC1: begin at the E-D-L codons (GAG GAT CTG), dropping the upstream EX0 residues
    trbc = [s for n, s in fasta(REF / "imgt_genes_nt" / "IMGT_TRBC_nt.fa") if "|TRBC1*01|" in n][0]
    trbc = trbc[trbc.index("gaggatctgagaaat"):]
    # trim both to complete codons, drop stop / TM tail not needed for window analysis
    trac = trac[: len(trac) // 3 * 3]
    trbc = trbc[: len(trbc) // 3 * 3]
    return trac, trbc


def apply_cys(nt, motif_aa, pos_in_motif):
    prot = str(Seq(nt).translate())
    i = prot.index(motif_aa) + pos_in_motif
    assert prot[i] in "TS", (prot[i], motif_aa)
    codon_start = i * 3
    return nt[:codon_start] + "tgc" + nt[codon_start + 3:]


def windows_overhangs(nt, n_codons=10):
    """All 4-mers writable by synonymous recoding within the first n_codons codons.
    Returns {4mer: [offsets]} (offset = nt index of the 4-mer start inside the constant)."""
    prot = str(Seq(nt[: n_codons * 3]).translate())
    out = {}
    for s in range(0, n_codons * 3 - 3):
        c0, c1 = s // 3, (s + 3) // 3
        span = range(c0, c1 + 1)
        for combo in product(*[SYN[prot[c]] for c in span]):
            joined = "".join(combo)
            off = s - c0 * 3
            oh = joined[off: off + 4]
            out.setdefault(oh.upper(), set()).add(s)
    return out


def main():
    trac, trbc = get_constants()
    trac_cys = apply_cys(trac, "DKTVLDM", 2)
    trbc_cys = apply_cys(trbc, "GVSTDP", 2)

    p_trac, p_trbc = str(Seq(trac_cys).translate()), str(Seq(trbc_cys).translate())
    ref_trac = fasta(REF / "imgt_genes_aa" / "IMGT_TRAC_LVL-Cys_aa.fa")[0][1].replace("X", "D", 1)
    ref_trbc = fasta(REF / "imgt_genes_aa" / "IMGT_TRBC1_Cys_aa.fa")[0][1].replace("X", "E", 1)
    n = 100  # compare upstream of the transmembrane region (LVL edits sit in the TM)
    assert p_trac[:n] == ref_trac[:n], "TRAC protein mismatch vs LVL-Cys reference (N-term)"
    assert p_trbc[: len(ref_trbc)] == ref_trbc, "TRBC1 protein mismatch vs TRBC1_Cys reference"
    print(f"TRAC Cys: {len(p_trac)} aa  starts {p_trac[:12]}  | TRBC1 Cys: {len(p_trbc)} aa  starts {p_trbc[:12]}")

    out = REF / "mouse_constants_cys_nt.fa"
    with open(out, "w") as fh:
        fh.write(">TRAC_Cys_T48C_no-stop  (starts at D codon after J; LVL not applied)\n" + trac_cys + "\n")
        fh.write(">TRBC1_Cys_S57C_no-stop (starts at E codon after J; LVL n/a)\n" + trbc_cys + "\n")
    print("wrote", out)

    for name, nt in (("TRAC", trac_cys), ("TRBC1", trbc_cys)):
        ohs = windows_overhangs(nt, n_codons=10)
        nonpal = {o: s for o, s in ohs.items() if o != rc(o).upper()}
        print(f"{name}: {len(ohs)} distinct 4-mers writable in first 10 codons "
              f"({len(nonpal)} non-palindromic) of 256 possible")


if __name__ == "__main__":
    main()

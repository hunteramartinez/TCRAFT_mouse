#!/usr/bin/env python
"""
Build per-gene V-region consensus protein sequences (leader -> Cys104, inclusive) from
a paired clonotype table (tra_*/trb_* columns, as in the HAM-07-007-012 file), and
compare each consensus to the IMGT mouse references.

The V block ends at the conserved Cys104 (the first residue of junction_aa), matching the
human TCRAFT convention where the first C of the CDR3 is supplied by the V gene.
Genes with fewer than --min-clones clones fall back to the IMGT reference (C57BL/6J *01
preferred, else the first functional allele).

usage:
  uv run python TCRAFT/scripts/build_v_consensus.py <clonotypes.csv> [--min-clones 30] [--out PATH]
"""
import argparse
import re
from collections import Counter
from pathlib import Path

import pandas as pd

REF_DIR = Path(__file__).resolve().parents[1] / "references" / "mouse" / "imgt_genes_aa"
FASTAS = {
    "TRA": ["IMGT_TRAV_aa.fa", "IMGT_TRAV_extra_aa.fa"],
    "TRB": ["IMGT_TRBV_aa.fa", "IMGT_TRBV_extra_aa.fa"],
}


def read_fasta(path):
    name, seq = None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield name, "".join(seq)
            name, seq = line[1:], []
        else:
            seq.append(line)
    if name:
        yield name, "".join(seq)


def load_refs(chain):
    """-> {gene: [(allele, strain, functionality, aa)]}"""
    refs = {}
    for fn in FASTAS[chain]:
        for hdr, seq in read_fasta(REF_DIR / fn):
            f = hdr.split("|")
            allele, strain, func = f[1], f[2], f[3]
            gene = allele.split("*")[0]
            refs.setdefault(gene, []).append((allele, strain, func, seq.replace(".", "").replace("*", "")))
    return refs


def cut_at_cys104(aa):
    """IMGT V-exons run a few residues past Cys104; cut at the last C within the final 8 aa."""
    i = aa.rfind("C", max(0, len(aa) - 8))
    return aa[: i + 1] if i != -1 else None


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b)) + abs(len(a) - len(b))


def fallback_ref(entries):
    ranked = sorted(entries, key=lambda e: (e[2] != "F", "C57BL/6" not in e[1], e[0]))
    for allele, strain, func, aa in ranked:
        v = cut_at_cys104(aa)
        if v:
            return allele, v
    return None, None


def consensus(seqs):
    """Majority residue per position among the modal-length sequences."""
    lens = Counter(len(s) for s in seqs)
    mode_len, n_mode = lens.most_common(1)[0]
    used = [s for s in seqs if len(s) == mode_len]
    cons = "".join(Counter(col).most_common(1)[0][0] for col in zip(*used))
    return cons, n_mode, lens


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("clonotypes")
    ap.add_argument("--min-clones", type=int, default=30)
    ap.add_argument("--out", default=str(Path(__file__).resolve().parents[1] / "references" / "mouse" / "mouse_V_consensus.csv"))
    args = ap.parse_args()

    cols = [f"{c}_{k}" for c in ("tra", "trb") for k in ("v_call", "junction_aa", "sequence_aa")]
    d = pd.read_csv(args.clonotypes, usecols=cols, low_memory=False)

    rows = []
    for chain in ("TRA", "TRB"):
        c = chain.lower()
        refs = load_refs(chain)
        sub = d[[f"{c}_v_call", f"{c}_junction_aa", f"{c}_sequence_aa"]].dropna()
        sub.columns = ["v", "jx", "seq"]
        # V region = everything before the junction (junction_aa starts at Cys104)
        pos = [s.find(j) for s, j in zip(sub.seq, sub.jx)]
        sub = sub.assign(vaa=[s[:p + 1] if p > 0 and j.startswith("C") else None for s, j, p in zip(sub.seq, sub.jx, pos)])
        sub = sub.dropna(subset=["vaa"])
        for gene, g in sub.groupby("v"):
            n = len(g)
            cons, n_mode, lens = consensus(g.vaa.tolist())
            frac_exact = float((g.vaa == cons).mean())
            entries = refs.get(gene, [])
            best = min(((hamming(cons, cut_at_cys104(a) or a), al) for al, _, _, a in entries), default=(None, None))
            src, final, allele = "consensus", cons, best[1]
            if n < args.min_clones or not re.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", cons):
                fb_allele, fb = fallback_ref(entries)
                src, final, allele = "imgt_fallback", fb, fb_allele
            rows.append(dict(chain=chain, gene=gene, n_clones=n, n_modal_len=n_mode,
                             n_distinct_lengths=len(lens), frac_exact=round(frac_exact, 4),
                             nearest_imgt_allele=best[1], mismatches_vs_imgt=best[0],
                             in_imgt=bool(entries), source=src, used_allele=allele, v_aa=final))
    out = pd.DataFrame(rows).sort_values(["chain", "gene"])
    # genes with identical leader+V protein share one vector; label each by its first member
    out["vector_group"] = out.groupby(["chain", "v_aa"]).gene.transform("first")
    out.to_csv(args.out, index=False)
    print("unique V blocks:", out.groupby("chain").vector_group.nunique().to_dict())
    print(f"wrote {args.out}  ({len(out)} genes)")
    print(out.groupby(["chain", "source"]).size())
    print("not in IMGT:", out.loc[~out.in_imgt, "gene"].tolist())
    print("mismatches vs nearest IMGT allele (consensus genes):")
    print(out[out.source == "consensus"].groupby("chain").mismatches_vs_imgt.describe()[["count", "mean", "max"]])


if __name__ == "__main__":
    main()

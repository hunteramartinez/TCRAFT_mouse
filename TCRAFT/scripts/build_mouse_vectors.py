#!/usr/bin/env python
"""
Build the complete plasmid sequences for the mouse TCRAFT cloning vectors (to be synthesized directly).

  Step 1 (22 TRBV-TRAC, KanR):  Kan+ori skeleton + [SapI-f | TRBV...OH | BbsI stuffer | OH TRAC ... TGA | SapI-r]
  Step 2 (85 TRBC1-TRAV, CamR): Cam+ori skeleton + [pool flank | OH TRBC1 ... P2A ... TRAV OH | pool flank]
  Step 3 destination (pMSCV, AmpR): existing vector with the TCR ORF replaced by a SapI-flanked GFP dropout
        ...GCCGCCACC ATG a GAAGAGC [promoter-GFP dropout] GCTCTTC a TGA GGAATT...   (SapI overhangs ATG / TGA)
     and the SapI site in the pUC origin removed by a single point change (as Gaglione et al. did).

Skeletons: Kan+ori = nt 3442-end of the platform KanR plasmid; Cam+ori = the cat + pUC-ori block of the platform
CamR plasmid with its single EcoRV site (5 nt after the cat stop) mutated GATATC->GATTTC.
Inputs: references/mouse/{TRBV,TRAV}_vector_inserts.csv, *_OH_frags.csv, TRAC_nt_conserved.txt, vectors/source/*.fa
Outputs: references/mouse/vectors/{Step1_TRBV_plasmids,Step2_TRAV_plasmids}.csv, skeleton_*.fa,
         pMSCV_destination_SapI.{fa,gb}, vector_QC.csv, gb/*.gb

usage: uv run python TCRAFT/scripts/build_mouse_vectors.py
"""
import json
import re
import sys
import textwrap
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from TCRAFT.overhang_tools import rc  # noqa: E402

REF = ROOT / "TCRAFT" / "references" / "mouse"
VEC = REF / "vectors"
SRC = VEC / "source"
(VEC / "gb").mkdir(exist_ok=True)

SITES = {"BbsI": "GAAGAC", "BsaI": "GGTCTC", "BsmBI": "CGTCTC", "SapI": "GCTCTTC", "AarI": "CACCTGC",
         "XhoI": "CTCGAG", "NotI": "GCGGCCGC", "EcoRI": "GAATTC", "EcoRV": "GATATC", "XbaI": "TCTAGA"}
# placeholders inherited from the human design (stuffer that carries the BbsI sites + remnant-destruction sites)
CASSETTE = "taGTCTTCGTCAcGAGACGTAATgGAGACCatgcGATATCcgtaGGTCTCgAAACCGTCTCgGGGTGAAGACat"
SAPI_F, SAPI_R = "GCTCTTCa", "aGAAGAGC"  # step-3 flanks: ...GCTCTTC a ATG ... TGA a GAAGAGC
# E. coli constitutive promoter (BBa_J23119) + RBS (BBa_B0034) for the GFP dropout
J23119, RBS = "TTGACAGCTAGCTCAGTCCTAGGTATAATGCTAGC", "AAAGAGGAGAAA"


def read_fa(p):
    return "".join(l.strip() for l in open(p) if not l.startswith(">")).upper()


def write_fa(p, name, s):
    Path(p).write_text(f">{name} | {len(s)} bp\n" + "\n".join(textwrap.wrap(s, 70)) + "\n")


def count(s, site):
    s = s.upper() + s.upper()[:len(site) - 1]  # circular wrap
    return s.count(site) + (s.count(rc(site)) if rc(site) != site else 0)


def audit(s):
    return {k: count(s, v) for k, v in SITES.items()}


def feat(a, b, label, ftype="misc_feature", strand=None):
    return SeqFeature(FeatureLocation(a, b, strand=strand), type=ftype, qualifiers={"label": [label]})


def record(name, short, desc, seq, feats):
    r = SeqRecord(Seq(seq), id=name, name=short, description=desc,
                  annotations={"molecule_type": "DNA", "topology": "circular"})
    r.features = feats
    return r


def main():
    plat = read_fa(SRC / "platform_KanR_URA3_GFP.fa")
    cam = read_fa(SRC / "platform_CamR_part.fa")
    base = read_fa(SRC / "pMSCV_TCR_base.fa")

    # ------------------------------------------------------------------ skeletons
    skelK = plat[plat.index("CCGAACGGTTATCCACAGAATCAG"):]
    skelC = cam[cam.index("GAGACCGACC") + 10:]
    assert skelC.count("GATATC") == 1
    skelC = skelC.replace("GATATC", "GATTTC")

    def orf_aa(s):
        best = ""
        for st in (s, rc(s)):
            for f in range(3):
                p = str(Seq(st[f:len(st) - ((len(st) - f) % 3)]).translate())
                for m in re.finditer(r"M[^*]{150,}\*", p):
                    best = max(best, m.group(), key=len)
        return best
    assert orf_aa(skelC).startswith("MEKKITGYTTVDISQ") and len(orf_aa(skelC)) == 220, "CamR ORF damaged"
    assert orf_aa(skelK).startswith("MSHIQRETSCSRPRLNS"), "KanR ORF damaged"
    kan_i = skelK.index("ATGAGCCATATTCAACGGGAAACG")
    kan_len = (len(orf_aa(skelK))) * 3
    for nm, sk in (("skeleton_KanR_ori", skelK), ("skeleton_CamR_ori", skelC)):
        a = audit(sk)
        assert not any(a.values()), f"{nm} has sites: {a}"
        write_fa(VEC / f"{nm}.fa", nm, sk)
        print(f"{nm}: {len(sk)} bp, site-free for {list(SITES)}")

    # ------------------------------------------------------------------ Step 1: TRBV-TRAC (KanR)
    ins1 = pd.read_csv(REF / "TRBV_vector_inserts.csv").set_index("TRBV")
    fr1 = pd.read_csv(REF / "TRBV_TRAC_OH_frags.csv").set_index("trbv")
    oh1 = pd.read_csv(REF / "TRBV_TRAC_OHset.csv").set_index("TRBV_name")
    conserved = (REF / "TRAC_nt_conserved.txt").read_text().strip().upper()
    assert conserved.endswith("TGA")
    members = json.load(open(REF / "overhang_sets.json"))["members"]
    rows1, qc, gbn = [], [], 0
    for gene in ins1.index:
        left, tracr = fr1.loc[gene, "trbv_left_frag"], fr1.loc[gene, "trac_right_frag"]
        assert (left + CASSETTE + tracr).upper() == ins1.loc[gene, "Sequence"].upper()
        parts = [("SapI site (step-3 cut, fwd) + spacer", SAPI_F), ("TRBV leader->Cys104-side (ATG..OH)", left),
                 ("BbsI stuffer (+BsmBI/BsaI/EcoRV remnant-kill sites)", CASSETTE),
                 ("TRAC first 30 nt (OH..)", tracr), ("TRAC remainder + stop (TGA)", conserved),
                 ("spacer + SapI site (step-3 cut, rev)", SAPI_R)]
        insert = "".join(p for _, p in parts).upper()
        seq = skelK + insert
        pos, feats = len(skelK), [feat(0, len(skelK), "Kan+ori backbone"), feat(kan_i, kan_i + kan_len, "KanR", "CDS", 1)]
        for lab, p in parts:
            feats.append(feat(pos, pos + len(p), f"{gene}: {lab}"))
            pos += len(p)
        gbn += 1
        name = f"pTCRAFT-S1-{gene}"
        SeqIO.write(record(name, f"pTCRAFT_S1_{gbn:02d}", f"TCRAFT step-1 vector {gene} (set {oh1.loc[gene,'set']}); KanR", seq, feats),
                    VEC / "gb" / f"{name}.gb", "genbank")
        rows1.append(dict(TRBV=gene, plasmid=name, Sequence=seq))
        qc.append(dict(plasmid=name, step=1, length=len(seq), **audit(seq)))
    pd.DataFrame(rows1).to_csv(VEC / "Step1_TRBV_plasmids.csv", index=False)

    # ------------------------------------------------------------------ Step 2: TRBC1-P2A-TRAV (CamR)
    ins2 = pd.read_csv(REF / "TRAV_vector_inserts.csv").set_index("TRAV")
    fr2 = pd.read_csv(REF / "TRBC_TRAV_OH_frags.csv").set_index("trav")
    oh2 = pd.read_csv(REF / "TRBC_TRAV_OHset.csv").set_index("TRAV_name")
    p2a = (REF / "P2A.txt").read_text().strip().upper()
    grp_of = {g: grp for grp, gl in members.items() for g in gl}
    rows2, made = [], {}
    for gene in ins2.index:
        grp = grp_of.get(gene, gene)
        ins = ins2.loc[gene, "Sequence"].upper()
        trbc_r, trav_l = fr2.loc[gene, "trbc_right_frag"].upper(), fr2.loc[gene, "trav_left_frag"].upper()
        h = ins.index(trbc_r)
        assert ins[h + len(trbc_r):].startswith(p2a) and ins[h + len(trbc_r) + len(p2a):].startswith(trav_l)
        head, tail = ins[:h], ins[h + len(trbc_r) + len(p2a) + len(trav_l):]
        name = f"pTCRAFT-S2-{grp}"
        seq = skelC + ins
        if name not in made:
            parts = [("pool flank (cloning sites)", head), ("TRBC1 from OH to end", trbc_r), ("GSG-P2A", p2a),
                     ("TRAV leader->OH", trav_l), ("pool flank (cloning sites)", tail)]
            pos, feats = len(skelC), [feat(0, len(skelC), "Cam+ori backbone")]
            for lab, p in parts:
                feats.append(feat(pos, pos + len(p), f"{grp}: {lab}"))
                pos += len(p)
            members_txt = ",".join(members.get(grp, [gene]))
            SeqIO.write(record(name, f"pTCRAFT_S2_{len(made)+1:03d}", f"TCRAFT step-2 donor {grp} (set {oh2.loc[gene,'set']}; genes: {members_txt}); CamR", seq, feats),
                        VEC / "gb" / f"{name}.gb", "genbank")
            made[name] = seq
            qc.append(dict(plasmid=name, step=2, length=len(seq), **audit(seq)))
        rows2.append(dict(TRAV=gene, plasmid=name, Sequence=seq))
    pd.DataFrame(rows2).to_csv(VEC / "Step2_TRAV_plasmids.csv", index=False)

    # ------------------------------------------------------------------ Step 3: pMSCV destination
    site = "GCTCTTCCGCTTCC"
    assert base.count(site) == 1
    b = base.replace(site, "GCTCTTTCGCTTCC")  # remove the pUC-ori SapI site (single C->T)
    atg = b.index("GCCGCCACCATGGCCCCCAGG") + 9
    prot = str(Seq(b[atg:]).translate(to_stop=True)) if len(b[atg:]) % 3 == 0 else str(Seq(b[atg:len(b) - (len(b) - atg) % 3]).translate(to_stop=True))
    stop_end = atg + 3 * len(prot) + 3
    assert b[stop_end - 3:stop_end] == "TGA" and b[stop_end:stop_end + 6] == "GGAATT", (b[stop_end - 3:stop_end + 6])
    g0 = plat.index("ATGCGTAAAGGCGAAGAACTG")
    gfp = plat[g0:plat.index("GGATCC", g0)]
    assert str(Seq(gfp).translate()).startswith("MRKGEELFTGVV") and str(Seq(gfp).translate()).endswith("ELYK")
    drop = J23119 + "TACTAGAG" + RBS + "TACTAG" + gfp + "TAATAA"
    dest = b[:atg] + "ATG" + "A" + "GAAGAGC" + drop + "GCTCTTC" + "A" + "TGA" + b[stop_end:]
    a = audit(dest)
    assert a["SapI"] == 2, a
    assert len(re.findall("ATG.GAAGAGC", dest)) == 1 and len(re.findall("GCTCTTC.TGA", dest)) == 1
    write_fa(VEC / "pMSCV_destination_SapI.fa", "pMSCV_destination_SapI (SapI-flanked GFP dropout; ori SapI removed)", dest)
    i0 = atg
    dfeats = [feat(i0, i0 + 3, "ATG / SapI overhang"), feat(i0 + 4, i0 + 11, "SapI site (rev)"),
              feat(i0 + 11, i0 + 11 + len(drop), "GFP dropout (J23119-RBS-GFP)"),
              feat(i0 + 11 + len(drop), i0 + 18 + len(drop), "SapI site (fwd)"),
              feat(i0 + 19 + len(drop), i0 + 22 + len(drop), "TGA / SapI overhang")]
    SeqIO.write(record("pMSCV-destination-SapI", "pMSCV_dest_SapI", "pMSCV_TCR base with SapI GFP dropout; AmpR", dest, dfeats),
                VEC / "pMSCV_destination_SapI.gb", "genbank")
    qc.append(dict(plasmid="pMSCV_destination_SapI", step=3, length=len(dest), **a))

    q = pd.DataFrame(qc)
    q.to_csv(VEC / "vector_QC.csv", index=False)
    s1, s2 = q[q.step == 1], q[q.step == 2]
    print(f"\nStep 1: {len(s1)} plasmids, {s1.length.min()}-{s1.length.max()} bp")
    print("  site counts (min/max):", {k: (int(s1[k].min()), int(s1[k].max())) for k in SITES})
    print(f"Step 2: {len(s2)} unique plasmids ({len(rows2)} gene rows), {s2.length.min()}-{s2.length.max()} bp")
    print("  site counts (min/max):", {k: (int(s2[k].min()), int(s2[k].max())) for k in SITES})
    print(f"Destination: {len(dest)} bp; sites:", a)
    assert (s1.SapI == 2).all() and (s1.BbsI == 2).all() and (s2.SapI == 0).all() and (s2.BbsI == 0).all()


if __name__ == "__main__":
    main()

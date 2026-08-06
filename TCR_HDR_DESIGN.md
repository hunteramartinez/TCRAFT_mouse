# Mouse TCR HDR Design — Specification

This document specifies the assembly strategy and design rules for Cas9-RNP mediated knock-in of a paired TCRαβ at the endogenous mouse TRAC locus. The implementation lives in `TCRAFT/HDR_Design.py`.

---

## Assembly architecture

The final HDR template (5'→3') after two-step golden gate assembly:

```
[LH] [P2A] [TRBV] [CDR3β] [TRBC-LVL-Cys] [T2A] [TRAV] [CDR3α] [TRAC-LVL-Cys] [RH]
```

- **Beta chain:** TRBV · CDR3β · TRBC-LVL-Cys  
- **Linker:** furin-GSG-T2A (self-cleaving)  
- **Alpha chain:** TRAV · CDR3α · TRAC-LVL-Cys  
- **LVL-Cys mutations** in both constant regions promote TCRαβ pairing over mispairing with endogenous chains.

---

## Genomic context and reading frame

The sgRNA (`GTTAAAAGATCCTCGGTCTC` + PAM `AGG`) cuts within TRAC exon 1, which encodes the beginning of the TRAC constant region (protein sequence: **IQN**PEPAVYQLKDPRS...).

The left homology arm (300 nt) ends with 1 nt of a shared codon at the cut site. The P2A leader sequence begins with 2 nt (`CT`) that complete this codon:

```
...LH_arm_last_nt(T) + P2A_first_2nt(CT) = TCT = Ser  ← shared codon
```

This is intentional: the reading frame is established by the upstream endogenous TRAC ORF across the homology arm. The apparent in-frame `TGA` stop codon in the P2A sequence is **out-of-frame in the genomic context** — it only appears as a stop when reading the P2A sequence in isolation.

**Reading frame verification** (from the start of TRAC exon 1 IQN coding):

| Region | Length | Cumulative from IQN | In frame? |
|---|---|---|---|
| TRAC coding to cut site | 46 nt | 46 nt | — |
| P2A leader | 68 nt | 114 nt = 38 codons | ✓ TRBV ATG in frame |
| TRBV (full, incl. signal peptide) | 3n nt | 114 + 3n nt | ✓ |
| CDR3β | 3n nt | — | ✓ |
| TRBC-LVL-Cys | 3n nt | — | ✓ |
| T2A (75 nt = 25 codons) | 75 nt | — | ✓ TRAV in frame |
| TRAV (full, incl. signal peptide) | 3n nt | — | ✓ |
| CDR3α | 3n nt | — | ✓ |
| TRAC-LVL-Cys | 3n nt | — | ✓ |

---

## Fragment types

### Constant fragments (order once per lab)

| Fragment | Count | Key sequences |
|---|---|---|
| **LH-TRBV** gene blocks | 24 (one per functional TRBV gene) | Left homology arm + P2A leader + full TRBV (with signal peptide) through Cys104 |
| **TRAC-RH** gene block | 1 | TRAC-LVL-Cys + right homology arm |
| **TRBC-T2A** gene block | 1 | TRBC-LVL-Cys + furin-GSG-T2A (partial, up to 2b split) |

### Variable fragment (order per TCR)

One gene block per clonotype, containing:

```
[BbsI 5' / oh_1a] [TRBV bridge] [CDR3β] [TRBC bridge] [BsmBI 3' / oh_2a]
[BsmBI 5' / oh_2b] [T2A tail] [TRAV full] [CDR3α] [TRAC bridge] [BbsI 3' / oh_1b]
```

The bridges and tail are sequence segments derived from the adjacent constant gene (TRBV, TRBC, T2A, TRAC) that fall between the OH window and the true boundary, ensuring exact sequence reconstitution after ligation.

---

## Golden gate scheme

### Step 1 — BbsI reaction

Assembles: **LH-TRBV** + **variable fragment** + **TRAC-RH**

| Junction | OH | Source |
|---|---|---|
| TRBV / CDR3β | **1a** (TRBV-specific) | Derived from within the TRBV NT sequence; unique per gene |
| CDR3α / TRAC | **1b** | Derived from within the TRAC-LVL-Cys NT sequence |

All 24 `1a` overhangs + `1b` must satisfy **Hamming distance ≥ 2** from every other OH in this reaction (including RC comparisons) with **0% cross-ligation** (verified against NEB BbsI ligation efficiency table). One TRBV gene (currently TRBV14) cannot satisfy this constraint and must be run in a **separate BbsI reaction pool**.

### Step 2 — BsmBI reaction

Inserts **TRBC-T2A** into the step-1 product.

| Junction | OH | Source |
|---|---|---|
| CDR3β / TRBC | **2a** | Derived from within the TRBC-LVL-Cys NT sequence |
| T2A / TRAV | **2b** | **Derived from within the T2A sequence** (see critical note below) |

Only `2a` and `2b` need to satisfy Hamming ≥ 2 from each other. No constraint between step-1 and step-2 OHs (different reactions).

---

## Critical design rules

### 1. No internal type IIS sites
All sequences — LH arm, RH arm, P2A leader, TRBV genes, TRAV genes, TRBC, TRAC, T2A, CDR3 regions — must be free of BbsI (`GAAGAC`), BsmBI (`CGTCTC`), BsaI (`GGTCTC`), SapI, XbaI, EcoRI, XhoI, EcoRV, and NotI recognition sequences. The only permitted occurrences of these sites are the designed 1a, 1b, 2a, 2b golden gate junctions.

Codon-optimised sequences (TRBV, TRAV, TRBC, TRAC) are actively screened and resampled to remove any forbidden sites. Fixed sequences (LH, RH, P2A, T2A) are verified clean at module load.

### 2. Overhangs must come from the flanking gene sequence
Every OH is a 4 nt window taken from within the sequence of the gene it flanks:

- **oh_1a** = `TRBV_nt[split_1a : split_1a+4]` → after ligation, full TRBV reconstituted
- **oh_1b** = `TRAC_nt[split_1b : split_1b+4]` → after ligation, full TRAC reconstituted
- **oh_2a** = `TRBC_nt[split_2a : split_2a+4]` → after ligation, full TRBC reconstituted
- **oh_2b** = `T2A_nt[split_2b : split_2b+4]` → after ligation, full T2A reconstituted

**This is the key rule that prevents frameshifts.** A synthetic 4-nt junction (4 mod 3 = 1) between T2A and TRAV would shift the TRAV reading frame by +1. Because oh_2b is derived from within T2A, the assembled product has no extra nucleotides at the T2A/TRAV boundary.

### 3. Overhang selection uses NEB ligation efficiency tables
OH candidates are ranked by NEB self-ligation efficiency (`df[OH, RC(OH)]`). The algorithm processes TRBV genes in order of fewest sequence-derived candidates first (hardest-first greedy assignment). Minimum Hamming distance = 2 between all OHs in the same reaction.

**Cross-ligation interpretation**: `df[X, RC(Y)]` gives the count of ligation events where overhang X mispaired with the fragment intended for overhang Y. A cross-rate below ~1% is acceptable; the step-1 OH set is designed to achieve 0% cross-ligation for all pairs except TRBV14 (flagged for separate pool).

### 4. Codon optimisation rules
- CDR3 regions: weighted-random sampling from mouse codon table (freq > 0.15), resampled up to 30 times to remove forbidden RE sites and homopolymers (≥5 nt)
- Constant V/J/C gene bodies: deterministic best-codon selection (single most-frequent codon per AA), then cleaned for RE sites and homopolymers
- P2A leader: **never codon-optimised** — used exactly as provided in `mm_trac_HDR_sequences.fa`

### 5. No homopolymers ≥ 5 nt in any synthesised sequence

### 6. GC content 25–65% for CDR3 codon-optimised regions

---

## Reference sequences

All reference sequences are in `TCRAFT/references/mouse/`:

| File | Contents |
|---|---|
| `mm_trac_HDR_sequences.fa` | LH arm, RH arm, sgRNA, PAM, P2A leader, full genomic context |
| `imgt_genes_aa/IMGT_TRBV_aa.fa` + `IMGT_TRBV_extra_aa.fa` | 24 functional mouse TRBV AA sequences |
| `imgt_genes_aa/IMGT_TRAV_aa.fa` + `IMGT_TRAV_extra_aa.fa` | Mouse TRAV AA sequences |
| `imgt_genes_aa/IMGT_TRAC_LVL-Cys_aa.fa` | TRAC with LVL-Cys mutation |
| `imgt_genes_aa/IMGT_TRBC_LVL-Cys_aa.fa` | TRBC with LVL-Cys mutation |
| `linkers_P2A_T2A.txt` | furin-GSG-T2A NT sequence |
| `mouse_codon_usage.json` | Mouse codon frequency table |
| `neb_tables/BbsI_ligation_efficiency.csv` | NEB BbsI overhang ligation counts |
| `neb_tables/BsmBI_ligation_efficiency.csv` | NEB BsmBI overhang ligation counts |

---

## TRBV gene coverage

23 functional TRBV genes from IMGT + 1 additional allele from Cell Ranger GRCm38 VDJ v7.0.0 = **24 LH-TRBV gene blocks**. TRBV21 is excluded (only UTR sequence, no L-REGION+V-REGION in any reference). The extra Cell Ranger genes include TRBV10, TRBV26, TRBV30, and the TRBV12-2+TRBV13-2 alternate splice variant.

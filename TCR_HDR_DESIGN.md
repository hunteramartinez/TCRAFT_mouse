# Mouse TCR HDR Design — Specification

This document specifies the assembly strategy and design rules for Cas9-RNP mediated knock-in of a paired TCRαβ at an endogenous mouse TCR locus. The implementation lives in `TCRAFT/HDR_Design.py`.

Two target loci are supported, selected with the `locus=` API parameter or `--locus` CLI flag:

| Locus | Target cell type | Leader | Key difference |
|---|---|---|---|
| **TRAC** (default) | Mature T cells | P2A self-cleaving peptide | Insert rides in-frame on the upstream TRAC ORF |
| **TRBC1** | DN3 double-negative thymocytes | Kozak context + TRBV native ATG | Independent translation start at the TRBC1 locus |

---

## Assembly architecture

### TRAC locus (mature T cells)

The final HDR template (5'→3') after two-step golden gate assembly:

```
[TRAC-LH] [P2A] [TRBV] [CDR3β] [TRBC-Cys] [T2A] [TRAV] [CDR3α] [TRAC-LVL-Cys] [BGH polyA] [TRAC-RH]
```

- **Beta chain:** TRBV · CDR3β · TRBC-Cys
- **Linker:** furin-GSG-T2A (self-cleaving)
- **Alpha chain:** TRAV · CDR3α · TRAC-LVL-Cys
- **LVL-Cys mutations** in both constant regions promote TCRαβ pairing over mispairing with endogenous chains.

### TRBC1 locus (DN3 thymocytes)

```
[TRBC1-LH] [Kozak] [TRBV] [CDR3β] [TRBC-Cys] [T2A] [TRAV] [CDR3α] [TRAC-LVL-Cys] [BGH polyA] [TRBC1-RH]
```

The insert is identical from TRBV onward. Only the flanking elements differ:
- **Left:** TRBC1 homology arm + Kozak context (`GCCGCCACC`) immediately 5' of the TRBV native ATG
- **Right:** TRBC1 homology arm (replaces TRAC-RH); the TRAC-LVL-Cys coding sequence and BGH polyA remain unchanged

In DN3 cells, the Vβ promoter + Eβ enhancer at the TRBC1 locus drives transcription. The Kozak sequence (`GCCGCCACC`) provides optimal ribosome recognition context immediately upstream of the TRBV's own ATG start codon (Met, first residue of the native signal peptide). There is no upstream ORF dependency, unlike the TRAC P2A mode.

---

## Genomic context and reading frame

### TRAC mode

The sgRNA (`GTTAAAAGATCCTCGGTCTC` + PAM `AGG`) cuts within TRAC exon 1. The left homology arm (300 nt) ends at the complete TCT (Ser) codon immediately upstream of the cut site. The P2A leader begins at a clean codon boundary:

```
...LH_arm_last_3nt(TCT=Ser) | P2A_first_codon(GGA=Gly)  ← clean codon boundary
```

The P2A encodes the standard GSG-P2A self-cleaving sequence (22 codons, 66 nt). Every component from P2A onward is an exact multiple of 3 nt.

**Reading frame verification (TRAC mode):**

| Region | Length | In frame? |
|---|---|---|
| P2A leader | 66 nt | ✓ TRBV ATG in frame |
| TRBV (full, incl. signal peptide) | 3n nt | ✓ |
| CDR3β | 3n nt | ✓ |
| TRBC-Cys | 3n nt | ✓ |
| T2A (75 nt = 25 codons) | 75 nt | ✓ TRAV in frame |
| TRAV (full, incl. signal peptide) | 3n nt | ✓ |
| CDR3α | 3n nt | ✓ |
| TRAC-LVL-Cys | 3n nt | ✓ |

### TRBC1 mode

The sgRNA (`GAGAAATGTGACTCCACCCA`) cuts within the TRBC1 locus. The Kozak context (`GCCGCCACC`, 9 nt) is in the 5'UTR and is not translated. The TRBV's own ATG is the start codon for the entire polyprotein.

**Reading frame verification (TRBC1 mode):**

| Region | Length | In frame? |
|---|---|---|
| Kozak (`GCCGCCACC`) | 9 nt | 5'UTR — not translated |
| TRBV (full, incl. signal peptide, ATG-start) | 3n nt | ✓ ORF starts here |
| CDR3β | 3n nt | ✓ |
| TRBC-Cys | 3n nt | ✓ |
| T2A (75 nt = 25 codons) | 75 nt | ✓ TRAV in frame |
| TRAV (full, incl. signal peptide) | 3n nt | ✓ |
| CDR3α | 3n nt | ✓ |
| TRAC-LVL-Cys | 3n nt | ✓ |

Frame is verified at run time by `_check_insert_frame()`, which is locus-aware.

---

## Fragment types

### Constant fragments (order once per lab)

The golden gate scheme is **identical for both loci**. Only the LH-TRBV blocks and TRAC-RH block differ between TRAC and TRBC1 runs.

| Fragment | Count | TRAC contents | TRBC1 contents |
|---|---|---|---|
| **LH-TRBV** gene blocks | N (one per functional TRBV gene) | TRAC-LH + P2A + TRBV through Cys104 | TRBC1-LH + Kozak + TRBV through Cys104 |
| **TRAC-RH** gene block | 1 | TRAC-LVL-Cys + spacer + BGH polyA + **TRAC-RH** | TRAC-LVL-Cys + spacer + BGH polyA + **TRBC1-RH** |
| **TRBC-T2A** gene block | 1 | TRBC1-Cys + furin-GSG-T2A (partial) | ← identical |

### Variable fragment (order per TCR)

One gene block per clonotype — **identical for both loci**:

```
[BbsI 5' / oh_1a] [TRBV bridge] [CDR3β] [TRBC bridge] [BsmBI 3' / oh_2a]
[BsmBI 5' / oh_2b] [T2A tail] [TRAV full] [CDR3α] [TRAC bridge] [BbsI 3' / oh_1b]
```

The bridges and tail are sequence segments derived from the adjacent constant gene (TRBV, TRBC, T2A, TRAC) that fall between the OH window and the true boundary, ensuring exact sequence reconstitution after ligation.

---

## Golden gate scheme

The scheme is identical for both loci.

### Step 1 — BbsI reaction

Assembles: **LH-TRBV** + **variable fragment** + **TRAC-RH**

| Junction | OH | Source |
|---|---|---|
| TRBV / CDR3β | **1a** (TRBV-specific) | Derived from within the TRBV NT sequence; unique per gene |
| CDR3α / TRAC | **1b** | Derived from within the TRAC-LVL-Cys NT sequence |

All 1a overhangs + 1b must satisfy **Hamming distance ≥ 2** from every other OH in this reaction (including RC comparisons). One TRBV gene (currently TRBV14) cannot satisfy this constraint and must be run in a **separate BbsI reaction pool**.

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
All sequences — LH arm, RH arm, leader (P2A or Kozak), TRBV, TRAV, TRBC, TRAC, T2A, CDR3 regions — must be free of BbsI (`GAAGAC`), BsmBI (`CGTCTC`), BsaI (`GGTCTC`), SapI, XbaI, EcoRI, XhoI, EcoRV, and NotI recognition sequences. The only permitted occurrences are the designed 1a, 1b, 2a, 2b golden gate junctions.

Codon-optimised sequences (TRBV, TRAV, TRBC, TRAC) are actively screened and resampled to remove any forbidden sites. Fixed sequences (LH, RH, P2A, Kozak, T2A) are verified clean at module load.

### 2. Overhangs must come from the flanking gene sequence
Every OH is a 4 nt window taken from within the sequence of the gene it flanks:

- **oh_1a** = `TRBV_nt[split_1a : split_1a+4]` → after ligation, full TRBV reconstituted
- **oh_1b** = `TRAC_nt[split_1b : split_1b+4]` → after ligation, full TRAC reconstituted
- **oh_2a** = `TRBC_nt[split_2a : split_2a+4]` → after ligation, full TRBC reconstituted
- **oh_2b** = `T2A_nt[split_2b : split_2b+4]` → after ligation, full T2A reconstituted

**This is the key rule that prevents frameshifts.** A synthetic 4-nt junction (4 mod 3 = 1) between T2A and TRAV would shift the TRAV reading frame by +1. Because oh_2b is derived from within T2A, the assembled product has no extra nucleotides at the T2A/TRAV boundary.

### 3. Overhang selection uses NEB ligation efficiency tables
OH candidates are ranked by NEB self-ligation efficiency (`df[OH, RC(OH)]`). The algorithm processes TRBV genes in order of fewest sequence-derived candidates first (hardest-first greedy assignment). Minimum Hamming distance = 2 between all OHs in the same reaction.

### 4. Codon optimisation rules
- CDR3 regions: weighted-random sampling from mouse codon table (freq > 0.15), resampled up to 30 times to remove forbidden RE sites and homopolymers (≥5 nt)
- Constant V/J/C gene bodies: deterministic best-codon selection, then cleaned for RE sites and homopolymers
- TRBC-Cys: codon-diverged from the endogenous TRAC-RH exonic sequence to reduce ectopic recombination risk (TRAC mode only; TRBC1 mode uses deterministic best-codon)
- P2A leader: **never codon-optimised** — used exactly as provided in `mm_trac_HDR_sequences.fa`
- Kozak context: **fixed sequence** (`GCCGCCACC`) — not codon-optimised (it is non-coding)

### 5. No homopolymers ≥ 5 nt in any synthesised sequence

---

## Reference sequences

All reference sequences are in `TCRAFT/references/mouse/`:

| File | Contents |
|---|---|
| `mm_trac_HDR_sequences.fa` | TRAC-LH, TRAC-RH, sgRNA, PAM, P2A leader, post-stop spacer, BGH polyA |
| `mm_trbc_HDR_sequences.fa` | TRBC1-LH, TRBC1-RH, sgRNA, Kozak context sequence |
| `imgt_genes_aa/IMGT_TRBV_aa.fa` + `IMGT_TRBV_extra_aa.fa` | Functional mouse TRBV AA sequences |
| `imgt_genes_aa/IMGT_TRAV_aa.fa` + `IMGT_TRAV_extra_aa.fa` | Mouse TRAV AA sequences |
| `imgt_genes_aa/IMGT_TRAC_LVL-Cys_aa.fa` | TRAC with LVL-Cys mutation |
| `imgt_genes_aa/IMGT_TRBC1_Cys_aa.fa` | TRBC1 with Cys mutation |
| `linkers_P2A_T2A.txt` | furin-GSG-T2A NT sequence |
| `mouse_codon_usage.json` | Mouse codon frequency table |
| `neb_tables/BbsI_ligation_efficiency.csv` | NEB BbsI overhang ligation counts |
| `neb_tables/BsmBI_ligation_efficiency.csv` | NEB BsmBI overhang ligation counts |

---

## TRBV gene coverage

Functional TRBV genes from IMGT + additional alleles from Cell Ranger GRCm38 VDJ v7.0.0. TRBV21 is excluded (only UTR sequence, no L-REGION+V-REGION in any reference). The extra Cell Ranger genes include TRBV10, TRBV26, TRBV30, and the TRBV12-2+TRBV13-2 alternate splice variant.

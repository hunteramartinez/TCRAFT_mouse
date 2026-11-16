# Mouse TRAC TCR HDR Template Assembly Protocol

**Purpose:** Assemble a full-length paired TCR HDR template via two sequential golden gate reactions, then PCR-amplify the product for Cas9-RNP electroporation into mouse T cells.

**Final product:** Linear dsDNA, ~2.5 kb:
```
5'-[LH 300 bp]-[P2A]-[TRBV]-[CDR3β]-[TRBC1-Cys]-[T2A]-[TRAV]-[CDR3α]-[TRAC-LVL-Cys]-[BGH polyA]-[RH 300 bp]-3'
```

**Assembly strategy:**
- **Step 1 (BbsI):** LH-TRBV + Variable gene block + TRAC-RH → ~2.1 kb intermediate
- **Step 2 (BsmBI):** Add BsmBI + TRBC-T2A directly to the heat-killed Step 1 reaction — **no cleanup between steps**
- **Step 3 (HiFi PCR):** Amplify assembled template to generate electroporation-grade DNA

---

## Materials

### Enzymes
| Reagent | Vendor | Cat # | Notes |
|---------|--------|-------|-------|
| BbsI | NEB | R0539 | 10,000 U/mL; active at 37°C |
| BsmBI-v2 | NEB | R0580 | 10,000 U/mL; active at 42°C |
| T4 DNA Ligase | NEB | M0202 | 400,000 U/mL |
| 10× T4 DNA Ligase Buffer | NEB | (with M0202) | Contains 10 mM ATP; keep frozen |
| Q5 Hot Start 2× Master Mix | NEB | M0494 | Or KAPA HiFi HotStart |
| Nuclease-free water | any | — | |

> **Critical:** Use 10× T4 Ligase Buffer (not CutSmart) for both golden gate steps — it supplies the ATP required for ligation. Do not substitute with restriction enzyme buffers.

### Cleanup
- Zymo DNA Clean & Concentrator-5 **or** AMPure XP beads (Beckman Coulter) — see notes below
- 1.5% agarose gel + SYBR Safe (or GelRed); TAE buffer
- Gel purification kit (e.g., Zymo Gel DNA Recovery) — only needed for troubleshooting

### Gene Blocks
Run TCRAFT (`TCRAFT-hdr your_airr.tsv`) to generate all sequences. Order as dsDNA gene blocks (IDT eBlocks or Twist Gene Fragments):

| Fragment | Length | Order |
|----------|--------|-------|
| LH-TRBV[N] | ~700–740 nt each | Once per TRBV gene needed |
| TRAC-RH | 954 nt | Once per lab |
| TRBC-T2A | 593 nt | Once per lab |
| Variable (CDR3β-TRAV-CDR3α) | ~490–530 nt | Per TCR clonotype |

> IDT eBlocks ship at ~250 fmol in 10 µL (25 nM). Twist fragments typically ship at ~20 ng/µL; convert to fmol: **fmol/µL = (ng/µL × 1000) / (length bp × 0.66)**.

---

## Part 0 — Stock Preparation

1. Briefly centrifuge all gene block tubes (2,000 × g, 30 s).
2. Resuspend or dilute each gene block to **25 nM** in nuclease-free water.
   - For IDT eBlocks (250 fmol / 10 µL = 25 nM): use as-is.
   - For Twist at 20 ng/µL (e.g., 700 nt fragment ≈ 43 fmol/µL = 43 nM): dilute ~1.7× with water.
3. Aliquot stocks into single-use tubes (5–10 µL each). Freeze at −20°C. Avoid repeated freeze-thaw of T4 Ligase Buffer (ATP degrades).
4. Thaw T4 Ligase Buffer on ice; vortex gently to dissolve any precipitate.

---

## Part 1 — Step 1: BbsI Golden Gate Assembly

**Joins:** LH-TRBV + Variable gene block + TRAC-RH  
**Junctions:** oh_1a (TRBV-specific) and oh_1b (GGAA)

> **TRBV14 note:** The LH-TRBV14 gene block cannot be pooled with other LH-TRBV blocks without risking junction cross-ligation. If assembling TRBV14-containing TCRs simultaneously with other TRBV genes in a pooled format, run TRBV14 in its own separate Step 1 reaction.

### 1.1 Reaction Setup (20 µL total, on ice)

| Component | Volume | Amount |
|-----------|--------|--------|
| Nuclease-free water | to 20 µL | — |
| 10× T4 DNA Ligase Buffer | 2.0 µL | 1× |
| LH-TRBV[N] gene block (25 nM) | 4.0 µL | **100 fmol** |
| Variable gene block (25 nM) | 4.0 µL | **100 fmol** |
| TRAC-RH gene block (25 nM) | 4.0 µL | **100 fmol** |
| BbsI (10,000 U/mL) | 1.0 µL | 10 U |
| T4 DNA Ligase (400,000 U/mL) | 2.5 µL | 1,000 U |

Mix gently by pipetting; spin briefly to collect. Keep on ice until thermocycler is ready.

### 1.2 Thermocycler Program

```
Step            Temp    Time    Cycles
───────────────────────────────────────
Digest          37°C    1 min   ┐
Ligate          16°C    1 min   ┘ × 25
Final ligate    16°C    10 min
Hold             4°C    ∞
```

> No heat inactivation after Step 1. BbsI has no recognition sites remaining in the assembled product (they were on the cut-away waste stubs) and is therefore inert in Step 2. T4 ligase and ATP remain active and carry directly into Step 2 — adding heat kill here would only destroy them unnecessarily. Proceed immediately to Step 2 setup.

### 1.3 Gel Verification (optional, recommended for first runs)

Remove **2 µL** and run on a 1.5% agarose gel alongside a 1 kb ladder.

- **Expected band:** ~2.1 kb (LH-TRBV ~720 nt + Variable ~510 nt + TRAC-RH 954 nt, minus ~80 nt in BbsI waste pieces).
- A faint smear of unligated fragments (500–950 nt) is normal; the 2.1 kb band should dominate.
- If no assembled band is visible, see Troubleshooting before proceeding.

Once you've validated Step 1 efficiency on your first assembly, you can skip the gel check in routine runs.

### 1.4 Proceed Directly to Step 2

No cleanup is needed. Why:

- BbsI is dead after 80°C/20 min and cannot re-cut anything.
- The assembled Step 1 product has **no remaining BbsI recognition sites** — the GAAGAC sequences were on the `ACAT` waste stubs that were cut off during cycling.
- The internal BsmBI sites on the Variable gene block are intact inside the assembled product and ready for Step 2.
- Unligated LH-TRBV and TRAC-RH fragments have **zero BsmBI sites** — BsmBI ignores them entirely.
- Any unligated Variable gene blocks do carry BsmBI sites, but truncated mis-ligation products they could form lack the LH/RH primer-binding sequences and are filtered out by the PCR in Part 3.
- ATP is preserved — no heat kill means no thermal degradation. T4 ligase and ATP carry over intact into Step 2; no replenishment needed.

Keep the tube on ice and proceed immediately.

---

## Part 2 — Step 2: BsmBI Golden Gate Assembly (one-pot)

**Joins:** Step 1 product (~2.1 kb) + TRBC-T2A gene block  
**Junctions:** oh_2a (GGAA) and oh_2b (AGGA)

### 2.1 Add to Existing Step 1 Reaction (on ice)

To the Step 1 reaction (cooled to 4°C), add directly:

| Component | Volume | Amount |
|-----------|--------|--------|
| TRBC-T2A gene block (25 nM) | 4.0 µL | **100 fmol** |
| BsmBI-v2 (10,000 U/mL) | 1.0 µL | 10 U |

**Total Step 2 volume: ~25 µL** (or ~23 µL if 2 µL was removed for the Step 1 gel check). No fresh T4 Ligase Buffer or T4 ligase is needed — both carry over intact from Step 1. BbsI from Step 1 is still present but has no recognition sites in the assembled product or TRBC-T2A and is inert.

Mix gently by pipetting. Spin briefly to collect.

> The Step 1 reaction starts with ~100 fmol of each input fragment. Assuming 50–80% ligation efficiency, the assembled 2.1 kb product is ~50–80 fmol — close to equimolar with the 100 fmol TRBC-T2A added. This slight TRBC-T2A excess is intentional and favors complete insertion.

### 2.2 Thermocycler Program

```
Step              Temp    Time    Cycles
──────────────────────────────────────────
Digest            42°C    1 min   ┐
Ligate            16°C    1 min   ┘ × 25
Final ligate      16°C    10 min
T4 ligase kill    65°C    10 min
BsmBI cleanup     55°C    20 min
Hold               4°C    ∞
```

> The T4 ligase kill (65°C) is the only heat inactivation required — it stops re-ligation so that BsmBI (still active, no sites remaining in the 2.5 kb full product) can cut any residual uninserted 2.1 kb intermediate at 55°C. BsmBI does not need a dedicated kill step; residual enzyme at 1:20 dilution is denatured by the 98°C PCR initial denaturation. 1 min/1 min cycling is sufficient — same rationale as Step 1.

### 2.3 Gel Verification

Remove **2 µL** and run on a 1% agarose gel. Load the saved Step 1 gel image (or a fresh 2 µL of a Step 1 reaction run in parallel) for size comparison.

- **Expected band:** ~2.5 kb — a clear upward shift from the 2.1 kb Step 1 product.
- **Do not proceed to PCR if any visible 2.1 kb band remains.** Both the 2.1 kb intermediate and the 2.5 kb full product carry the outer HDR_F/HDR_R primer sites and will be co-amplified exponentially. Shorter templates amplify more efficiently, so residual 2.1 kb can be enriched in the final PCR product — 2.1 kb template integrating at TRAC disrupts the locus but produces no functional TCR. If 2.1 kb remains visible, add an additional 50–100 fmol TRBC-T2A and run a further 10 cycles of the Step 2 program before re-checking.
- If no 2.5 kb band is visible at all, see Troubleshooting.


---

## Part 3 — Asymmetric PCR Amplification of Assembled Template

**Purpose:** Generate ≥ 2 µg of ssDNA-enriched HDR template for electroporation. A 20:1 molar excess of HDR_R over HDR_F is used: after the limiting HDR_F primer is consumed (typically within the first 15–20 cycles), subsequent cycles perform linear amplification of the antisense ssDNA strand using only HDR_R. The predominantly ssDNA product prevents annealing between amplicons from different clonotypes — which share substantial constant-region sequence — reducing the likelihood of heteroduplex formation during and after electroporation.

> **Alternatives:** **Linear PCR** — a single-primer reaction using only HDR_R directly synthesises ssDNA from a dsDNA template in a linear (non-exponential) fashion; simpler but requires a clean dsDNA starting template and more total reaction volume to reach µg quantities. **Emulsion PCR** — compartmentalises individual template molecules in water-in-oil droplets, completely preventing inter-amplicon annealing; preferred for complex pooled libraries where cross-annealing between clonotypes is a primary concern.

**Primers:** HDR_F (`5'-TGCCAGAAAGGACCCCGTGG-3'`) + HDR_R (`5'-TGCCTAATGATGAGGGTGAGATG-3'`)  
**Ratio:** 20:1 HDR_R : HDR_F  
**Expected product:** ~2.5 kb, predominantly antisense ssDNA (exact size varies by TRBV/TRAV gene and CDR3 length)

### 3.1 Primer Preparation

Pre-dilute HDR_F 1:10 in nuclease-free water to a **1 µM working stock**. Keep HDR_R at 10 µM. Prepare the working dilution fresh or aliquot and store at −20°C.

### 3.2 Reaction Setup (per 50 µL reaction)

| Component | Volume | Final concentration |
|-----------|--------|---------------------|
| Nuclease-free water | 14.5 µL | — |
| Q5 Hot Start 2× Master Mix | 25.0 µL | 1× |
| HDR_F (1 µM working stock) | 2.5 µL | **0.05 µM** (limiting) |
| HDR_R (10 µM) | 5.0 µL | **1.0 µM** (excess; 20:1 ratio) |
| Step 2 reaction (diluted 1:20 in water) | 3.0 µL | — |

Dilute the heat-killed Step 2 reaction 1:20 in nuclease-free water before use — no cleanup is needed. Set up **4–8 replicate reactions** in parallel (200–400 µL total).

### 3.3 Thermocycler Program

```
Step                 Temp    Time     Cycles
─────────────────────────────────────────────
Initial denaturation 98°C    30 s     1×
─────────────────────────────────────────────
Denature             98°C    10 s     ┐
Anneal               62°C    30 s     ├ × 35–40
Extend               72°C    90 s     ┘  (30 s/kb for ~2.5 kb)
─────────────────────────────────────────────
Final extension      72°C    2 min    1×
Hold                  4°C    ∞
```

> More cycles than standard PCR are used because the reaction transitions from exponential dsDNA synthesis (while both primers are present) to linear ssDNA synthesis (once HDR_F is exhausted) — the additional cycles accumulate the ssDNA product. 62°C annealing is a starting point; raise to 64°C if nonspecific bands appear or drop to 60°C if yield is low.

### 3.3 Gel Verification

Pool all replicate reactions. Remove **5 µL** for gel:

1. Run on 1% agarose alongside a 1 kb ladder. Expect a band at ~2.5 kb; the ssDNA product may appear as a faint, slightly diffuse band just above the expected dsDNA position due to secondary structure.
2. A doublet (tight dsDNA + diffuse ssDNA above) is normal and confirms asymmetric amplification.
3. **If a band at ~2.1 kb is present alongside the 2.5 kb product**, this indicates residual uninserted Step 1 intermediate was co-amplified (see Section 2.3). Gel-purify the 2.5 kb band before electroporation — do not use the mixed product.

> **Gel purification is recommended as a routine step** for the first assembly of any new clonotype to confirm complete Step 2 insertion. Once you have established that Step 2 runs to completion reliably with your reagent stocks, you may skip it in routine runs. Use a Zymo Gel DNA Recovery kit; elute in 15–20 µL nuclease-free water pre-warmed to 50°C. Requantify by Qubit ssDNA HS Assay after gel purification.

### 3.4 Cleanup — NEB Monarch 5 µg PCR Cleanup Kit

Pool all replicate reactions before cleanup. Follow the manufacturer's protocol (NEB T1030); key parameters:

1. Add **2 volumes** of DNA Binding Buffer to pooled PCR; mix by pipetting.
2. Load onto the Monarch column in batches; spin 1 min at 16,000 × g; discard flow-through.
3. Wash with **200 µL** DNA Wash Buffer; spin 1 min. Repeat wash once.
4. Transfer column to a clean 1.5 mL tube. Add **20 µL** nuclease-free water (pre-warmed to 50°C) directly to the centre of the membrane; incubate 1 min at room temperature; spin 1 min at 16,000 × g to elute.

> The 5 µg capacity accommodates 4–8 pooled 50 µL reactions (typical yield 0.5–2 µg each). If pooling ≥ 8 reactions, use two columns and combine eluates.

### 3.5 Quantification and Quality Check

**Quantify:** Use **Qubit ssDNA HS Assay** for accurate measurement of the predominantly ssDNA product — the dsDNA HS Assay will underestimate ssDNA concentration. Nanodrop 260/280 ≥ 1.8 confirms purity.

**Target final concentration:** ≥ 100 ng/µL in ≤ 20 µL (≥ 2 µg total).

Before freezing aliquots, run 100 ng on a 1% gel:

- Band at expected size (single or doublet): proceed.
- Visible degradation (smearing below 2 kb): repeat with fresh template.

**Aliquot and store:** 1–2 µg per tube at −20°C. Avoid repeated freeze-thaw; use a fresh aliquot for each electroporation.

---

## Part 3 (Alternative — Future) — IVT + Reverse Transcription for Pure ssDNA Template

> **Not yet implemented.** This section is a placeholder for a future alternative to asymmetric PCR. IVT + RT produces truly pure ssDNA rather than the predominantly ssDNA mixture from asymmetric PCR, and scales more efficiently to large quantities.

**Strategy:** A T7 promoter is appended to the 5′ end of HDR_F. Standard symmetric PCR generates a dsDNA template carrying the T7 promoter. T7 IVT transcribes the sense strand RNA, which is then reverse-transcribed with HDR_R to produce antisense ssDNA — the same strand orientation as asymmetric PCR.

**Key design note — T7 +1 G requirement:** T7 RNA polymerase requires a G at the +1 transcription start position. The current HDR_F begins with T (`TGCCAGAAAGGACCCCGTGG`), so the T7-modified primer must include a G spacer:

```
T7-HDR_F: 5'-TAATACGACTCACTATA-G-TGCCAGAAAGGACCCCGTGG-3'
```

The extra G is incorporated into the PCR product, transcribed into the RNA, and becomes a C at the 3′ tail of the final ssDNA — outside the LH arm and inert for HDR.

**Workflow outline (to be developed into a full protocol):**
1. Symmetric PCR with T7-HDR_F + HDR_R → dsDNA template
2. DNase I treatment → remove dsDNA template before IVT
3. T7 IVT → sense strand RNA
4. RNA cleanup (e.g., Zymo RNA Clean & Concentrator)
5. RT with HDR_R primer → antisense ssDNA (RNA:DNA hybrid)
6. RNase H or NaOH hydrolysis → degrade RNA template
7. ssDNA cleanup and Qubit ssDNA HS Assay quantification

---

## Electroporation Notes

- Typical input for mouse primary T cell RNP electroporation: **1–3 µg linear HDR template** in ≤ 3 µL.
- The template is co-electroporated with Cas9-RNP assembled with the TRAC sgRNA (`5'-GTTAAAAGATCCTCGGTCTC-3'`, PAM: AGG).
- The ssDNA-enriched template is less prone to inter-amplicon annealing than dsDNA but remains susceptible to cellular nucleases; use freshly thawed aliquots and electroporate promptly after thawing.

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|--------------|-----|
| No 2.1 kb band after Step 1 | Low fragment concentration; degraded ATP in buffer | Check stock concentrations by Qubit; use fresh T4 Ligase Buffer; increase cycling to 35 rounds |
| No upward shift to ~2.5 kb after Step 2 | Incomplete Step 1 product; or BsmBI not active | Verify Step 1 worked first; try 35 Step 2 cycles; confirm BsmBI-v2 stock is active (use a control digest) |
| Persistent 2.1 kb band with weak 2.5 kb | Partial TRBC-T2A insertion | Normal if minor; PCR selects for full-length product. If dominant, increase TRBC-T2A to 150–200 fmol (6–8 µL of 25 nM stock) relative to Step 1 product |
| PCR gives correct 2.5 kb band + extra lower bands | Residual unligated Variable gene blocks amplified off-target | Dilute Step 2 reaction further (1:50) before PCR; both HDR_F and HDR_R must bind for amplification so true off-targets are unlikely — more likely a primer-dimer artifact; check a no-template control |
| 2.1 kb band present in PCR product alongside 2.5 kb | Residual uninserted Step 1 intermediate co-amplified — both products share outer primer sites | Gel-purify the 2.5 kb band before electroporation; trace back to Step 2 and increase TRBC-T2A or cycling rounds |
| PCR yield < 50 ng/µL | Too little template; sub-optimal annealing; HDR_F working stock too dilute | Try 1:10 template dilution (more template); run gradient PCR 60–65°C; verify HDR_F 1 µM working stock is correctly prepared |
| Low 260/280 after Monarch cleanup | Residual binding buffer | Ensure two full wash steps; increase spin time to 2 min on final wash |

---

## Pooled One-Pot Library Assembly

Use this approach when assembling **multiple TCR HDR templates simultaneously** — e.g., a repertoire of 10–96 clonotypes for a polyclonal knock-in experiment.

### Concept

Because every oh_1a overhang is TRBV-specific and Hamming ≥ 2 from all others in the same reaction, BbsI ligation is directional and exclusive: only the LH-TRBV fragment with overhang X can join the Variable fragment that also carries overhang X on its 5′ end. All Variable fragments share the same oh_1b on their 3′ end, so they all join the single TRAC-RH fragment. The result is a pool of correctly paired full-length templates assembled in one tube.

Step 2 is even simpler — oh_2a and oh_2b are identical across all TCRs, so a single TRBC-T2A addition inserts the beta constant region into every assembled intermediate simultaneously.

### TRBV14 Separation Requirement

**TRBV14 (oh_1a = CTGA) must be assembled in its own BbsI reaction**, separate from the other 23 TRBV genes. Its overhang cannot satisfy Hamming ≥ 2 from all other oh_1a overhangs when pooled, risking mis-ligation at the 1a junction. After Step 1, the TRBV14 sub-pool and the main pool can be **combined before Step 2** — the BsmBI junctions (oh_2a, oh_2b) are identical across all TCRs and are not affected by which BbsI pool they came from.

### Fragment Loading Strategy

The goal is equimolar representation of every assembled TCR in the final pool. The key insight is that each LH-TRBV gene block must be loaded in proportion to the number of TCRs using that TRBV gene.

**Example:** 12 TCRs — TRBV5 (4 clones), TRBV12-1 (5 clones), TRBV19 (3 clones)

| Fragment | Relative moles | Rationale |
|----------|---------------|-----------|
| LH-TRBV5 | 4× | Must match the 4 Variable fragments it will pair with |
| LH-TRBV12-1 | 5× | Same logic |
| LH-TRBV19 | 3× | Same logic |
| Variable_TCR1…4 (TRBV5) | 1× each | Equal representation desired |
| Variable_TCR5…9 (TRBV12-1) | 1× each | |
| Variable_TCR10…12 (TRBV19) | 1× each | |
| TRAC-RH | 12–15× | At least 1× per Variable fragment; slight excess fine |

In practice, set a base amount of **50 fmol per Variable gene block** and scale the LH-TRBV and TRAC-RH amounts accordingly. Example for the 12-TCR pool above in a 50 µL Step 1 reaction:

| Component | Calculation | Volume (25 nM stock) |
|-----------|-------------|----------------------|
| LH-TRBV5 | 4 × 50 fmol = 200 fmol | 8 µL |
| LH-TRBV12-1 | 5 × 50 fmol = 250 fmol | 10 µL |
| LH-TRBV19 | 3 × 50 fmol = 150 fmol | 6 µL |
| Each Variable (×12) | 50 fmol each | 2 µL each = 24 µL total |
| TRAC-RH | 12 × 50 fmol = 600 fmol (or round to 500) | 20 µL |
| 10× T4 Ligase Buffer | — | 5 µL |
| BbsI | — | 2.5 µL |
| T4 Ligase | — | 5.0 µL |
| Water | to 50 µL | ~3.5 µL |

> Scale the reaction volume up (50–100 µL) proportionally to the number of TCRs and the total DNA input. DNA amounts above show the proportional logic; adjust to fit your available stocks and pipetting precision. For very large pools (> 48 TCRs), use a 100 µL reaction with the same proportional ratios.

### Step 1: Pooled BbsI Assembly

Use the same thermocycler program as Part 1.2 (37°C/16°C × 25 cycles + final ligation hold, no heat kill).

**For TRBV14 TCRs:** run a parallel, smaller Step 1 reaction containing only LH-TRBV14 + the TRBV14 Variable gene blocks + TRAC-RH (scaled as above). After both Step 1 reactions complete and cool to 4°C, **combine them in a single tube** before Step 2.

### Step 2: One-Pot BsmBI Insertion

Add directly to the combined Step 1 pool (cooled to 4°C, no heat kill):

| Component | Amount |
|-----------|--------|
| TRBC-T2A (25 nM) | 1× moles per Variable fragment total (e.g., 12 × 50 fmol = 600 fmol = 24 µL) |
| BsmBI-v2 | 1 µL per 20 µL reaction volume (scale proportionally) |

No fresh T4 Ligase Buffer or T4 ligase needed — both carry over intact from Step 1. Run the same Step 2 thermocycler program as Section 2.2 (42°C/16°C × 25 cycles + final hold, then 65°C T4 kill + 55°C BsmBI cleanup).

### PCR Amplification — Asymmetric PCR with Reduced Cycles to Preserve Pool Balance

Use the **same asymmetric primer ratio as Part 3** (20:1 HDR_R : HDR_F; HDR_F at 0.05 µM, HDR_R at 1.0 µM). Because all templates share the same primer-binding sites but differ in internal length (~50–150 nt across clonotypes), longer templates amplify slightly less efficiently per cycle during the exponential phase. Asymmetric PCR reduces this bias: once HDR_F is exhausted (early, given its low starting concentration), all templates switch to linear ssDNA amplification simultaneously regardless of length.

**Use 20–22 total cycles.** The limiting HDR_F is typically consumed within 10–15 cycles at 0.05 µM, so the exponential (length-biased) phase is brief. Set up more replicate reactions to compensate for lower per-reaction yield rather than adding cycles.

```
Step                 Temp    Time     Cycles
─────────────────────────────────────────────
Initial denaturation 98°C    30 s     1×
Denature             98°C    10 s     ┐
Anneal               62°C    30 s     ├ × 20–22
Extend               72°C    90 s     ┘
Final extension      72°C    2 min    1×
Hold                  4°C    ∞
```

Set up **8–16 × 50 µL reactions** in parallel using the Step 2 pool diluted 1:20 as template, then pool and clean up with the **NEB Monarch 5 µg PCR Cleanup Kit** as in Part 3.4 (use multiple columns if needed). Quantify with the Qubit ssDNA HS Assay.

### Pool QC Before Electroporation

After AMPure cleanup, run the following checks before use:

1. **Gel (1% agarose):** The pool should appear as a tight smear or a thick band centered at ~2.5 kb. A broad smear extending significantly below 2 kb suggests incomplete Step 2 assembly; repeat with more TRBC-T2A or more BsmBI cycling rounds.

2. **Bioanalyzer or TapeStation (optional, recommended):** A High Sensitivity D5000 or DNF-936 chip resolves individual template sizes and confirms the pool is not dominated by one or two clonotypes. Major size outliers suggest a pipetting error in loading.

3. **Qubit:** Target ≥ 100 ng/µL in ≤ 20 µL total (≥ 2 µg). For a 12-TCR pool each template contributes ~8% of the mass; Qubit measures total, not per-clonotype representation.

4. **Deep sequencing QC (optional):** For critical experiments, sequence the pooled PCR product (amplicon-seq with HDR_F/HDR_R) to verify per-clonotype representation. Aim for < 10× fold-difference between most and least abundant clonotypes.

> **If deep-sequencing reveals skewed representation:** normalize individual clonotypes by adding back the under-represented TCRs as separately assembled and PCR-amplified products (Part 1–3 single-TCR protocol), then remix to the desired ratio before electroporation.

---

## Overhang Reference (from TCRAFT output)

Overhangs used in the **current** design (from `custom_tcrs/assembly_report.txt`):

| Junction | Overhang | Reaction |
|----------|----------|----------|
| TRBV / CDR3β | oh_1a — TRBV-specific (see assembly_report.txt) | BbsI Step 1 |
| CDR3α / TRAC | oh_1b = **GGAA** | BbsI Step 1 |
| CDR3β / TRBC | oh_2a = **GGAA** | BsmBI Step 2 |
| T2A / TRAV | oh_2b = **AGGA** | BsmBI Step 2 |

> oh_1b and oh_2a both happen to be GGAA, but this is not a conflict — they are in separate reactions (BbsI vs BsmBI) and the enzymes do not cross-react.

**Outer PCR primers** (`custom_tcrs/ha_primers.tsv`):

| Primer | Sequence (5'→3') | Anchors at |
|--------|-----------------|------------|
| HDR_F | `TGCCAGAAAGGACCCCGTGG` | 5' end of LH arm |
| HDR_R | `TGCCTAATGATGAGGGTGAGATG` | 3' end of RH arm |

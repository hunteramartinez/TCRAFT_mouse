# TCRAFT — mouse implementation

TCRAFT (**T**CR **R**apid **A**ssembly for **F**unctional **T**esting) is a scalable, low-cost platform for pooled
assembly of synthetic TCR libraries by Golden Gate cloning. This repository is a **direct mouse implementation of
the human TCRAFT design**: the assembly logic, oligo layout, A/B pooling scheme and software interface are the same;
only the sequences (V/J genes, constant regions, codon usage, overhang sets) are mouse.

The software turns a table of TCRs into CDR3 oligonucleotide pools to order, and then simulates the Golden Gate
assembly of those oligos to check that every TCR assembles into exactly the intended protein.

> Original method: Gaglione *et al.*, "Scalable TCR synthesis and screening enable antigen reactivity mapping in
> vitiligo", *Immunity* 2026 (`docs/gaglione_et_al.pdf`). Original human code: Rachit Mukkamala,
> [birnbaumlab/TCRAFT](https://github.com/birnbaumlab/TCRAFT).
> The separate CRISPR-HDR knock-in design for mouse lives on the `mouse-trac-hdr` branch.

## The mouse cassette

```
TRBV ─ CDR3β ─ TRBJ ─ TRBC1(Cys) ─ GSG-P2A ─ TRAV ─ CDR3α ─ TRAJ(+J/C residue) ─ TRAC(Cys) ─ stop
```

- **Order:** β – P2A – α. The P2A includes the N-terminal GSG linker (`GSGATNFSLLKQAGDVEENPGP`).
- **Constants:** Cys-engineered mouse TRBC1 (S57C) and TRAC (T48C) for preferential α/β pairing. The LVL
  transmembrane mutation is **not** applied. All TCRs use TRBC1.
- **Oligo layout** (identical to human): `fwd primer | BbsI | TRBV-right | CDR3β-J | TRBC-left | BsmBI-or-BsaI middle
  site | TRAV-right | CDR3α-J | TRAC-left | BbsI | rev primer`.
- **Three Golden Gate steps:** (1) the oligo is ligated into TRBV–TRAC vectors with BbsI; (2) a TRBC1-P2A-TRAV insert is
  placed between the two CDR3s with BsmBI (TRAV pool A) or BsaI (pool B); (3) the finished cassette is cut out with SapI and moved into the pMSCV retroviral destination vector.
- **Overhangs:** each V-gene vector carries a unique pair of 4-nt overhangs written into the gene by synonymous
  recoding only (never changing the protein; Cys104 is always preserved).
- **J–C junction:** the residue between TRAJ and TRAC is encoded by `(last J nucleotide) + AC` and is N, D, Y or H
  depending on the J gene. It is carried by the TRAJ sequence and TRAC starts at `IQN…`, as in the human design.

### Gene coverage

The V genes are those observed in the clonotype table used to build the references
(HAM-07-007-012, 164,646 paired clones): **89 TRAV and 22 TRBV genes** (identical-protein genes share a vector:
85 unique TRAV blocks), and **55 J genes** (43 TRAJ, 12 TRBJ). TCRs that use any other gene are reported in
`Invalid_CDR3_table.csv`. V sequences are per-gene **consensus** proteins (leader → Cys104) taken from the Cell Ranger
`sequence_aa` of every clone using that gene, so they follow your animals' sequences, not IMGT. (IMGT is only the
fallback for genes with <30 clones; several IMGT "N"-gene records are also wrong or truncated.)

### Pools

| Reaction | Vectors | Pools | Model fidelity f1* |
|---|---|---|---|
| Step 1 TRBV–TRAC (BbsI) | 22 | 1 (A) | 1.000 |
| Step 2 TRBC1–TRAV | 85 | 2: A = BsmBI (43), B = BsaI (42) | 0.9987 / 0.9991 |

\*Worst-case predicted correct-assembly score from the NEB 4-nt ligation-fidelity tables, using the scoring function of
Gaglione *et al.* (human sets score 0.63–0.80 under the same function). These are model scores, not measured rates.

**The `Pool` column in the oligo files is the Step-1 (TRBV) pool, so every oligo is in pool A and the "Pool B" files
are empty.** The TRAV A/B split is built into each oligo (middle site BsmBI vs BsaI), so you order one pool and
Step 2 separates it by enzyme.

## Installation

Requires [uv](https://docs.astral.sh/uv/). Python 3.11 and the package versions used are pinned
(`.python-version`, `pyproject.toml`, `uv.lock`).

```bash
git clone https://github.com/hunteramartinez/TCRAFT_mouse.git
cd TCRAFT_mouse
git checkout mouse-tcraft-direct
uv sync
```

Run everything with `uv run` (e.g. `uv run TCRAFT-generate ...`). Dependencies: numpy, pandas, tqdm, biopython.

## Usage

### 1. Prepare the TCR table

CSV with these columns (extra columns are ignored):

| Column | Content |
|---|---|
| `V_alpha`, `V_beta` | mouse V gene, no allele, e.g. `TRAV6D-6`, `TRBV12-1` |
| `J_alpha`, `J_beta` | J gene, no allele, e.g. `TRAJ31`, `TRBJ1-1` |
| `CDR3_alpha`, `CDR3_beta` | junction amino acids, **C … F/W** (the conserved Cys through the J-anchor F/W) |

To convert a paired clonotype table (columns `tra_v_call`, `tra_j_call`, `tra_junction_aa`, and the `trb_*`
equivalents):

```bash
uv run python TCRAFT/scripts/prepare_input_from_clonotypes.py clonotypes.csv tcrs.csv [--sample 2000 --seed 1]
```

### 2. Generate the oligo pools

```bash
uv run TCRAFT-generate tcrs.csv --output_dir out/generate
```
```python
from TCRAFT.Generate import generate_cdr3_oligos
generate_cdr3_oligos("tcrs.csv", output_dir="out/generate")
```

Outputs (in `output_dir`, default `./TCRAFT-generate_<date>`): `All_<N>_CDR3_oligos.csv`,
`Under300_*`/`Over300_*` oligo files (≤300 bp for e.g. Twist, >300 bp for e.g. IDT), `Oligo_Pool_Metadata.csv`
and `Invalid_CDR3_table.csv` (TCRs with an unsupported gene or that failed generation). Codon choice is random, so
repeated runs give different (equivalent) DNA; seed `random`/`numpy` for reproducibility.

### 3. Validate by simulated assembly

```bash
uv run TCRAFT-validate out/generate/All_<N>_CDR3_oligos.csv --output_dir out/validate
```

Simulates all three Golden Gate steps for every oligo on the real cloning plasmids (Step 1 BbsI into the TRBV–TRAC
vector, Step 2 BsmBI/BsaI with the TRBC1-P2A-TRAV donor, Step 3 SapI into the pMSCV destination), translates the ORF in
the final plasmid and requires (a) the exact intended protein, (b) no stray BsmBI/BsaI/BbsI/SapI/XbaI/XhoI/EcoRI/EcoRV/NotI
sites in the cassette and (c) no SapI site left in the final plasmid. Outputs
`Reference_Assembled_TCR_Sequences.csv` (the ATG…TGA cassette of each passing TCR, plus the final plasmid length), `Failed_Assemblies.csv` and
`Assembly_Metadata.csv`.

> Note: the protein check uses the same reference tables as generation, so it tests the assembly mechanics, not
> the biology of the tables themselves. The tables were additionally checked against the source data (below).

## Reference tables (`TCRAFT/references/mouse/`)

The pipeline reads only this folder. The human tables in `TCRAFT/references/` are kept for provenance and are not used.
The files below mirror the human files one-to-one (same schemas):

`TRBV_TRAC_OHset.csv`, `TRBC_TRAV_OHset.csv` (overhang pairs + pool), `TRBV_TRAC_OH_frags.csv`,
`TRBC_TRAV_OH_frags.csv` (oligo-side/vector-side fragments), `TRBV_vector_inserts.csv`, `TRAV_vector_inserts.csv`
(gene blocks), `TRJ_seqs.csv`, `TRAC_protein.txt`, `TRBC1_protein.txt`, `TRAC_nt_conserved.txt`, `P2A.txt`,
`Oligo_Components.json`, `mouse_codon_usage.json`; plus the mouse-specific `V_reference_aa.csv`,
`mouse_V_consensus.csv`, `mouse_J_from_data.csv`, `overhang_sets.json` and `mouse_constants_cys_nt.fa`.
(`cellranger_genes*`, `imgt_genes_*`, `mm_*_HDR_sequences.fa`, `mouse_J_FR4.csv` and `tcr_antibodies/` are inputs or
leftovers from the HDR work and are not read by the pipeline.)

## Cloning vectors (`TCRAFT/references/mouse/vectors/`)

Complete plasmid sequences, intended to be synthesized directly (no platform part workflow needed):

| Step | Plasmids | Backbone / marker | Insert |
|---|---|---|---|
| 1 | `Step1_TRBV_plasmids.csv` (22) | Kan + pUC ori (1,763 bp from the platform KanR plasmid) | `GCTCTTC·a·ATG-TRBV…OH` – BbsI stuffer – `OH…TRAC…TGA·a·GAAGAGC` |
| 2 | `Step2_TRAV_plasmids.csv` (85 unique, 89 gene rows) | Cam + pUC ori (1,642 bp: cat + ori; EcoRV site removed) | pool-specific flank – `OH…TRBC1-P2A-TRAV…OH` – flank |
| 3 | `pMSCV_destination_SapI.fa/.gb` | pMSCV (AmpR) | `…GCCGCCACC·ATG·a·GAAGAGC`–[J23119-RBS-GFP dropout]–`GCTCTTC·a·TGA·GGAATT…` |

- Markers: Step 1 recipient = Kan, Step 2 donors = Cam, Step 3 destination = Amp, so each step's carry-over is counter-selected.
- Step 3 uses SapI with **ATG / TGA** as the 3-nt overhangs; the destination keeps its Kozak (`GCCGCCACC`) and the
  downstream IRES-CD8 marker. GFP-dropout colonies (green) are unrecombined vector.
- The destination differs from the supplied `pMSCV_TCR` vector by (i) removal of the SapI site in the pUC origin
  (`GCTCTTCC`→`GCTCTTTC`, as Gaglione *et al.* did by mutagenesis; **confirm copy number**) and (ii) replacement of the
  existing TCR ORF with the dropout. Grow LTR-containing plasmids in a recombination-deficient strain.
- GenBank maps are in `vectors/gb/`; `vector_QC.csv` lists restriction-site counts per plasmid; `vectors/source/` holds the
  supplied source plasmids.
- The E. coli promoter/RBS of the GFP dropout (BBa_J23119, BBa_B0034) were typed from standard parts, not from a supplied
  sequence: verify by sequencing.

### Rebuilding the tables

Run in this order (`clonotypes.csv` is the paired clonotype table):

```bash
uv run python TCRAFT/scripts/build_v_consensus.py clonotypes.csv     # per-gene V consensus (+ IMGT comparison)
uv run python TCRAFT/scripts/build_j_from_data.py clonotypes.csv     # J segments + TRAJ junction residue
uv run python TCRAFT/scripts/build_mouse_constants.py                # Cys constants, domesticated, checked vs IMGT
uv run python TCRAFT/scripts/design_overhangs.py                     # overhang-set search → overhang_sets.json
uv run python TCRAFT/scripts/build_mouse_tables.py                   # all CSV/TXT tables above
uv run python TCRAFT/scripts/build_mouse_vectors.py                  # complete plasmids, destination vector, GenBank maps
```

`TCRAFT/overhang_tools.py` holds the fidelity scoring (NEB tables in `references/neb_tables/`), the Metropolis
search and the synonymous-recoding/site-removal helpers. Overhang windows copy the human ranges: the first 10 codons of
each constant and the last 33 nt (TRBV) / 21 nt (TRAV) of each V block.

## Validation so far

Full run on all **164,646** paired clones of HAM-07-007-012 (97% singletons; every clone, expanded or not):

- **Generate:** 164,646 / 164,646 valid oligos, 0 invalid. 296 (0.18%) are >300 bp (mean 262 bp; mean of the >300 bp
  group 304 bp). About 650 oligos needed resampling to clear a junction restriction site.
- **Validate (simulated assembly, all three steps on the real cloning plasmids):** 164,646 / 164,646 assembled through
  Step 1 (BbsI, Kan vectors), Step 2 (90,379 via BsmBI pool A, 74,267 via BsaI pool B, Cam donors) and Step 3 (SapI into
  the pMSCV destination) to the intended protein, with no stray restriction sites and no SapI left in the final plasmid;
  0 failures. Final pMSCV-TCR plasmids average 7,326 bp. (Run time 223 s.)
- **Independent check against the source data:** the translated V+CDR3+J of each assembled TCR was compared with that
  clone's own `sequence_aa` in the input file (164,086 unique TCRs): identical for **99.69% of β chains and 98.31% of
  α chains**. The differences are mostly single residues where a clone deviates from its gene's consensus. This
  check is what exposed the J-gene-specific TRAJ/TRAC residue and two wrong J-table rows during development.

**No wet-lab validation yet.**

## Known placeholders and limitations

- The Step-2 **pool flanks** (XhoI/NotI ends, BsaI/BsmBI sites, `TAAT`/`AAAC`/`GTCA`/`GGGT`) and the Step-1 **BbsI stuffer
  cassette** are still copied unchanged from the human design (they are inside the new cloning plasmids). The plasmids
  themselves are untested; the GFP-dropout promoter/RBS and the origin SapI point mutation need experimental confirmation.
- Only genes seen in the source clonotype table are supported (10 further J genes in IMGT have no data and are omitted).
- Overhang fidelity numbers are model predictions.

## Version
1.1.0+mouse (human TCRAFT 1.1.0 → mouse implementation)

## Authors
Mouse implementation: Hunter Martinez. Original TCRAFT software: [Rachit Mukkamala](mailto:rsmukk@mit.edu).

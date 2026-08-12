# TCRAFT — Mouse

TCRAFT is a scalable and cost-effective platform for TCR library assembly. This mouse-specific fork extends the original human TCRAFT package with a second pipeline: **TCRAFT-hdr**, which designs all DNA fragments required for Cas9-RNP mediated knock-in of a paired TCRαβ into an endogenous mouse TCR locus.

---

## Installation

Requires Python 3.9+ and [uv](https://docs.astral.sh/uv/).

```bash
git clone <repo>
cd TCRAFT_mouse
uv sync
```

This creates a virtual environment and installs all dependencies from `uv.lock`. To install in editable mode for development:

```bash
uv sync --editable
```

---

## Two pipelines

### 1. TCRAFT-generate (human CDR3 oligo pools)

Converts a CSV of TCR sequences into CDR3 oligonucleotide pools for the human TCRAFT golden gate assembly.

**Input:** CSV with columns `V_alpha`, `V_beta`, `J_alpha`, `J_beta`, `CDR3_alpha`, `CDR3_beta`

```bash
TCRAFT-generate /path/to/tcr_list.csv --output_dir ./output/
```
```python
from TCRAFT.Generate import generate_cdr3_oligos
generate_cdr3_oligos("tcr_list.csv", output_dir="./output/")
```

**Outputs:** oligo CSVs split by pool (A/B) and length (≤300 bp / >300 bp), metadata table, invalid TCR log.

**Validate assembled sequences:**
```bash
TCRAFT-validate /path/to/oligos.csv
```

---

### 2. TCRAFT-hdr (mouse TCR knock-in design)

Designs all gene blocks and per-TCR variable fragments for a two-step BbsI → BsmBI golden gate assembly that produces an HDR template for Cas9-RNP knock-in at an endogenous mouse TCR locus.

Two target loci are supported:

| `--locus` | Target cells | Use when... |
|---|---|---|
| `TRAC` *(default)* | Mature T cells | Cells have already rearranged their TCR loci |
| `TRBC1` | DN3 double-negative thymocytes | Cells have **not** yet rearranged their TCRβ locus |

**Input:** AIRR-format rearrangement TSV (e.g., from Cell Ranger VDJ)

#### TRAC locus (mature T cells — default)

```bash
TCRAFT-hdr airr_rearrangement.tsv --output_dir ./output/hdr_trac/
```
```python
from TCRAFT.HDR_Design import generate_hdr_design
result = generate_hdr_design("airr_rearrangement.tsv", output_dir="./output/hdr_trac/")
```

HDR template structure:
```
[TRAC-LH] [P2A] [TRBV] [CDR3β] [TRBC-Cys] [T2A] [TRAV] [CDR3α] [TRAC-LVL-Cys] [BGH polyA] [TRAC-RH]
```

The P2A self-cleaving peptide rides in-frame with the upstream TRAC open reading frame. The insert is expressed as part of a polyprotein that is cleaved co-translationally.

#### TRBC1 locus (DN3 thymocytes)

```bash
TCRAFT-hdr airr_rearrangement.tsv --locus TRBC1 --output_dir ./output/hdr_trbc1/
```
```python
result = generate_hdr_design("airr_rearrangement.tsv", locus='TRBC1', output_dir="./output/hdr_trbc1/")
```

HDR template structure:
```
[TRBC1-LH] [Kozak] [TRBV] [CDR3β] [TRBC-Cys] [T2A] [TRAV] [CDR3α] [TRAC-LVL-Cys] [BGH polyA] [TRBC1-RH]
```

A Kozak context sequence (`GCCGCCACC`) is placed immediately 5' of each TRBV's native ATG start codon. The cassette is expressed from the endogenous Vβ promoter/Eβ enhancer at the TRBC1 locus. The golden gate assembly scheme — and the variable fragment — are **identical** to the TRAC design; only the LH-TRBV gene blocks and the right homology arm in the TRAC-RH gene block change.

#### Restrict to a single TRBV gene (testing / specific experiments)

```bash
TCRAFT-hdr airr_rearrangement.tsv --trbv TRBV5 --output_dir ./output/
TCRAFT-hdr airr_rearrangement.tsv --locus TRBC1 --trbv TRBV5 --output_dir ./output/
```

#### Outputs

```
output/hdr_*/
├── constant_gene_blocks.csv   — LH-TRBV + TRAC-RH + TRBC-T2A fragments (order once per lab)
├── variable_gene_blocks.csv   — 1 variable fragment per TCR clone (order per experiment)
├── overhangs.json             — all OH sequences and NEB efficiency scores
├── assembly_report.txt        — locus, OH table, pool assignments, fragment sizes
├── invalid_clones.csv         — skipped clones with error reasons
└── gb/
    ├── constants/             — annotated .gb file per constant fragment
    └── clonotypes/
        └── <clone_id>/
            ├── variable.gb    — gene block with GG site annotations
            └── assembled.gb   — final HDR template, all regions annotated
```

The `assembled.gb` files can be opened directly in SnapGene, Benchling, or Geneious to visually verify the design before ordering.

#### Reaction pools

All TRBV-specific LH-TRBV fragments, variable fragments, TRBC-T2A, and TRAC-RH assemble in a single BbsI reaction (step 1), followed by a single BsmBI reaction (step 2) to insert TRBC-T2A. One TRBV gene (currently TRBV14) requires a separate BbsI reaction pool — this is flagged in `assembly_report.txt`.

See `TCR_HDR_DESIGN.md` for full specifications.

---

## References

- Mouse TRBV/TRAV reference sequences: IMGT + Cell Ranger GRCm38 VDJ v7.0.0
- LVL-Cys constant region sequences: engineered for improved TCRαβ pairing
- Ligation efficiency tables: NEB BbsI and BsmBI (Pryor et al. 2020, PLOS ONE)
- TRAC sgRNA: `GTTAAAAGATCCTCGGTCTC` + PAM `AGG`
- TRBC1 sgRNA: `GAGAAATGTGACTCCACCCA`

## Version
1.3.0

## Authors
Original TCRAFT: Rachit Mukkamala  
Mouse HDR extension: Hunter Martinez

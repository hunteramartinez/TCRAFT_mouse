# TCRAFT — Mouse

TCRAFT is a scalable and cost-effective platform for TCR library assembly. This mouse-specific fork extends the original human TCRAFT package with a second pipeline: **TCRAFT-hdr**, which designs all DNA fragments required for Cas9-RNP mediated knock-in of a paired TCRαβ into the endogenous mouse TRAC locus.

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

### 2. TCRAFT-hdr (mouse TRAC knock-in design)

Designs all gene blocks and per-TCR variable fragments for a two-step BbsI → BsmBI golden gate assembly that produces a ~2.4 kb HDR template for Cas9-RNP knock-in at the endogenous mouse TRAC locus.

**Input:** AIRR-format rearrangement TSV (e.g., from Cell Ranger VDJ)

```bash
TCRAFT-hdr airr_rearrangement.tsv --output_dir ./output/hdr_design/
```
```python
from TCRAFT.HDR_Design import generate_hdr_design
result = generate_hdr_design("airr_rearrangement.tsv", output_dir="./output/hdr_design/")
```

The `--trbv` flag restricts output to clones using a specific TRBV gene (useful for testing):
```bash
TCRAFT-hdr airr_rearrangement.tsv --trbv TRBV5 --output_dir ./output/
```

#### Assembled HDR template structure

```
[Left Homology Arm] [P2A leader] [TRBV] [CDR3β] [TRBC-LVL-Cys] [T2A] [TRAV] [CDR3α] [TRAC-LVL-Cys] [Right Homology Arm]
```

All coding regions are in the correct reading frame relative to the upstream endogenous TRAC ORF. See `TCR_HDR_DESIGN.md` for full specifications.

#### Outputs

```
output/hdr_design/
├── constant_gene_blocks.csv   — 26 fragments (order once per lab)
├── variable_gene_blocks.csv   — 1 fragment per TCR clone (order per experiment)
├── overhangs.json             — all OH sequences and NEB efficiency scores
├── assembly_report.txt        — OH table, pool assignments, fragment sizes
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

All 24 TRBV-specific LH-TRBV fragments, variable fragments, TRBC-T2A, and TRAC-RH assemble in a single BbsI reaction (step 1), followed by a single BsmBI reaction (step 2) to insert TRBC-T2A. One TRBV gene (currently TRBV14) requires a separate BbsI reaction pool — this is flagged in `assembly_report.txt`.

---

## References

- Mouse TRBV/TRAV reference sequences: IMGT + Cell Ranger GRCm38 VDJ v7.0.0
- LVL-Cys constant region sequences: engineered for improved TCRαβ pairing
- Ligation efficiency tables: NEB BbsI and BsmBI (Pryor et al. 2020, PLOS ONE)
- sgRNA target: `GTTAAAAGATCCTCGGTCTC` + PAM `AGG` at mouse TRAC locus

## Version
1.2.0

## Authors
Original TCRAFT: Rachit Mukkamala  
Mouse HDR extension: Hunter Martinez

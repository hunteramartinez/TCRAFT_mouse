# -*- coding: utf-8 -*-
"""
Build Cell Ranger supplement entries for additional-genes.fasta.

Extracts L-REGION+V-REGION nucleotide sequences for mouse TCR V genes present in
the 10x Genomics Cell Ranger GRCm38 VDJ reference but absent from IMGT. These are
stored as ~VARIABLE entries in stitchr's additional-genes.fasta. When airrvector
stitches a clone annotated to one of these genes, it passes no_leader=True so the
combined L+V sequence serves as the full N-terminal germline sequence.

Excluded genes:
  TRBV21 — no L-REGION+V-REGION entry in Cell Ranger reference (UTR only)

Included alternate splice variants:
  TRBV12-2+TRBV13-2 — a very common alternate splice: exon 1 of TRBV12-2 combined with
    exon 2 of TRBV13-2. Named and added to the GRCm38 VDJ reference by 10x Genomics.
    The L+V sequence from that reference is stored under this name so clones annotated
    to this variant can be stitched directly.

Run from the stitchr repo root:
  python scripts/build_mouse_cellranger_supplement.py
      --regions path/to/regions.fa
      --output src/Data/additional-genes.fasta
"""

import argparse
import os
import re


EXTRA_GENES = [
    "TRAV3N-3",
    "TRAV4-4-DV10",
    "TRAV6-7-DV9",
    "TRAV7N-6",
    "TRAV13-4-DV7",
    "TRAV14D-3-DV8",
    "TRAV15-1-DV6-1",
    "TRAV15-2-DV6-2",
    "TRAV15D-1-DV6D-1",
    "TRAV16D-DV11",
    "TRAV21-DV12",
    "TRBV10",
    "TRBV26",
    "TRBV30",
    "TRBV12-2+TRBV13-2",   # alternate splice: TRBV12-2 exon 1 + TRBV13-2 exon 2 (named by 10x)
]

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def translate(nt_seq):
    aa = []
    for i in range(0, len(nt_seq) - (len(nt_seq) % 3), 3):
        codon = nt_seq[i:i+3].upper()
        aa.append(CODON_TABLE.get(codon, 'X'))
    return ''.join(aa)


def strip_stop(nt_seq):
    """Truncate NT sequence at the first in-frame stop codon."""
    aa = translate(nt_seq)
    if '*' in aa:
        stop_pos = aa.index('*')
        return nt_seq[:stop_pos * 3]
    return nt_seq


def parse_regions_fa(path):
    """
    Parse Cell Ranger regions.fa, return dict:
      {gene_name: {region_type: nt_sequence}}
    region_type is the 4th pipe-delimited field, e.g. 'L-REGION+V-REGION'.
    """
    genes = {}
    current_gene = None
    current_region = None
    current_seq = []

    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if current_gene and current_region:
                    genes.setdefault(current_gene, {})[current_region] = ''.join(current_seq)
                parts = line[1:].split('|')
                # format: index|gene_name ensid|gene_name|region_type|chain|locus|..
                current_gene = parts[2].strip() if len(parts) > 2 else None
                current_region = parts[3].strip() if len(parts) > 3 else None
                current_seq = []
            else:
                current_seq.append(line.upper())

    if current_gene and current_region:
        genes.setdefault(current_gene, {})[current_region] = ''.join(current_seq)

    return genes


def make_fasta_entry(gene_name, nt_seq, source='CellRanger-GRCm38-v7.0.0'):
    """Format a ~VARIABLE additional-genes.fasta entry."""
    allele = gene_name + '*01'
    header = f'>{source}|{allele}|Mus musculus|F|L-REGION+V-REGION||~VARIABLE'
    # wrap sequence at 60 chars
    seq_lines = [nt_seq[i:i+60] for i in range(0, len(nt_seq), 60)]
    return header + '\n' + '\n'.join(seq_lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--regions', required=True,
                        help='Path to Cell Ranger regions.fa')
    parser.add_argument('--output', required=True,
                        help='Path to additional-genes.fasta (entries will be appended)')
    args = parser.parse_args()

    if not os.path.isfile(args.regions):
        raise FileNotFoundError(f"regions.fa not found: {args.regions}")
    if not os.path.isfile(args.output):
        raise FileNotFoundError(f"output file not found: {args.output}")

    print(f"Parsing {args.regions}...")
    cr_genes = parse_regions_fa(args.regions)

    entries = []
    for gene in EXTRA_GENES:
        gene_data = cr_genes.get(gene, {})
        lv_seq = gene_data.get('L-REGION+V-REGION')
        if not lv_seq:
            print(f"  [SKIP] {gene}: no L-REGION+V-REGION entry found")
            continue

        cleaned = strip_stop(lv_seq)
        if len(cleaned) < len(lv_seq):
            removed = len(lv_seq) - len(cleaned)
            print(f"  [TRIM] {gene}: removed {removed} nt ({removed//3} codons) at in-frame stop")

        if len(cleaned) < 90:
            print(f"  [WARN] {gene}: sequence suspiciously short ({len(cleaned)} nt) after trimming — skipping")
            continue

        entries.append(make_fasta_entry(gene, cleaned))
        print(f"  [OK]   {gene}: {len(cleaned)} nt")

    print(f"\nAppending {len(entries)} entries to {args.output}...")
    with open(args.output, 'a') as fh:
        fh.write('\n')
        fh.write('\n'.join(entries))
        fh.write('\n')

    print("Done.")


if __name__ == '__main__':
    main()

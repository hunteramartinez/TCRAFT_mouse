# -*- coding: utf-8 -*-
"""
Build LVL-Cys constant region entries for additional-genes.fasta.

LVL-Cys TRAC and TRBC are full-length constant regions carrying an engineered
cysteine mutation near the transmembrane domain that forms an inter-chain disulfide
bond, forcing introduced TCR alpha and beta chains to pair with each other rather
than with endogenous chains in transduced T cells.

The source AA sequences begin with X — the J-C junction amino acid encoded at the
exon boundary whose codon is completed at runtime by stitchr's tidy_c_term. This X
is stripped before back-translation. The resulting synthetic codon-optimised NT
sequences are added as ~CONSTANT entries. The C-region start motif for frame-finding
is identical to the parent mTRAC*01 / mTRBC1*01 / mTRBC2*01 entries already in the
MOUSE C-region-motifs.tsv.

Run from the stitchr repo root:
  python scripts/build_lvlcys_constants.py
      --trac   path/to/IMGT_TRAC_LVL-Cys_aa.fa
      --trbc   path/to/IMGT_TRBC_LVL-Cys_aa.fa
      --codons src/Data/kazusa/MOUSE.txt
      --output src/Data/additional-genes.fasta
"""

import argparse
import collections
import os
import re


def parse_kazusa(path):
    """Parse Kazusa codon usage file, return {AA: best_codon} dict."""
    codon_table = {
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

    # freq[AA][codon] = frequency
    freq = collections.defaultdict(dict)
    with open(path) as fh:
        for line in fh:
            # strip parenthesized counts, split on whitespace
            cleaned = re.sub(r'\(.+?\)', '', line.rstrip()).upper().replace('U', 'T').split()
            if len(cleaned) % 2 != 0 or not cleaned:
                continue
            for i in range(0, len(cleaned), 2):
                codon, val = cleaned[i], float(cleaned[i+1])
                aa = codon_table.get(codon)
                if aa and aa != '*':
                    freq[aa][codon] = val

    return {aa: max(codons, key=codons.get) for aa, codons in freq.items()}


def rev_translate(aa_seq, codon_dict):
    result = []
    for aa in aa_seq:
        codon = codon_dict.get(aa)
        if codon is None:
            raise ValueError(f"No codon found for amino acid '{aa}'")
        result.append(codon)
    return ''.join(result)


def parse_fasta(path):
    """Return list of (header, sequence) tuples from a FASTA file."""
    records = []
    header, seq_parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(seq_parts)))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        records.append((header, ''.join(seq_parts)))
    return records


def make_constant_entry(gene_name, allele_suffix, nt_seq, source='LVL-Cys-synthetic'):
    """Format a ~CONSTANT additional-genes.fasta entry."""
    allele = f'{gene_name}*{allele_suffix}'
    header = f'>{source}|{allele}|Mus musculus|F|EX1+EX2+EX3+EX4||~CONSTANT'
    seq_lines = [nt_seq[i:i+60] for i in range(0, len(nt_seq), 60)]
    return header + '\n' + '\n'.join(seq_lines)


_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_DEFAULT_TRAC   = os.path.join(_REPO_ROOT, 'references', 'imgt_genes_aa', 'IMGT_TRAC_LVL-Cys_aa.fa')
_DEFAULT_TRBC   = os.path.join(_REPO_ROOT, 'references', 'imgt_genes_aa', 'IMGT_TRBC_LVL-Cys_aa.fa')
_DEFAULT_CODONS = os.path.join(_REPO_ROOT, 'src', 'Data', 'kazusa', 'MOUSE.txt')
_DEFAULT_OUTPUT = os.path.join(_REPO_ROOT, 'src', 'Data', 'additional-genes.fasta')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--trac', default=_DEFAULT_TRAC,
                        help='Path to IMGT_TRAC_LVL-Cys_aa.fa (default: references/imgt_genes_aa/)')
    parser.add_argument('--trbc', default=_DEFAULT_TRBC,
                        help='Path to IMGT_TRBC_LVL-Cys_aa.fa (default: references/imgt_genes_aa/)')
    parser.add_argument('--codons', default=_DEFAULT_CODONS,
                        help='Path to kazusa/MOUSE.txt codon usage file (default: src/Data/kazusa/)')
    parser.add_argument('--output', default=_DEFAULT_OUTPUT,
                        help='Path to additional-genes.fasta (entries will be appended; default: src/Data/)')
    args = parser.parse_args()

    for p in [args.trac, args.trbc, args.codons, args.output]:
        if not os.path.isfile(p):
            raise FileNotFoundError(f"File not found: {p}")

    print(f"Loading codon table from {args.codons}...")
    codon_dict = parse_kazusa(args.codons)

    entries = []

    # --- TRAC ---
    print(f"\nProcessing TRAC from {args.trac}...")
    trac_records = parse_fasta(args.trac)
    if not trac_records:
        raise ValueError("No records found in TRAC file")

    _, trac_aa = trac_records[0]
    if trac_aa.startswith('X'):
        trac_aa = trac_aa[1:]
        print("  Stripped leading X (J-C junction residue)")
    trac_aa = trac_aa.rstrip('-*')
    trac_nt = rev_translate(trac_aa, codon_dict)
    entries.append(make_constant_entry('mTRAC', 'LVLCys', trac_nt))
    print(f"  [OK] mTRAC*LVLCys: {len(trac_aa)} AA → {len(trac_nt)} nt")

    # --- TRBC ---
    # The TRBC LVL-Cys file contains one sequence covering both TRBC1 and TRBC2
    # (they differ only at a few positions in exon 1). We register it as both
    # mTRBC1*LVLCys and mTRBC2*LVLCys so either constant region call resolves.
    print(f"\nProcessing TRBC from {args.trbc}...")
    trbc_records = parse_fasta(args.trbc)
    if not trbc_records:
        raise ValueError("No records found in TRBC file")

    _, trbc_aa = trbc_records[0]
    if trbc_aa.startswith('X'):
        trbc_aa = trbc_aa[1:]
        print("  Stripped leading X (J-C junction residue)")
    trbc_aa = trbc_aa.rstrip('-*')
    trbc_nt = rev_translate(trbc_aa, codon_dict)
    entries.append(make_constant_entry('mTRBC1', 'LVLCys', trbc_nt))
    entries.append(make_constant_entry('mTRBC2', 'LVLCys', trbc_nt))
    print(f"  [OK] mTRBC1*LVLCys + mTRBC2*LVLCys: {len(trbc_aa)} AA → {len(trbc_nt)} nt")

    # Guard against accidental re-runs doubling the entries
    with open(args.output) as fh:
        existing = fh.read()
    already_present = [e for e in entries if e.split('\n')[0] in existing]
    if already_present:
        print(f"\nERROR: {len(already_present)} entr(ies) already present in {args.output}.")
        print("  Remove them first or use a fresh file. Aborting.")
        return

    print(f"\nAppending {len(entries)} entries to {args.output}...")
    with open(args.output, 'a') as fh:
        fh.write('\n')
        fh.write('\n'.join(entries))
        fh.write('\n')

    print("Done.")


if __name__ == '__main__':
    main()

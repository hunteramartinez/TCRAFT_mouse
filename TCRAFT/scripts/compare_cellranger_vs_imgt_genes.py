import os
from pathlib import Path
from Bio import SeqIO

def get_cellranger_genes(fasta_path, prefixes):
    """
    Parses regions.fa headers.
    Example: >644|TRAV13-1 GRCm38-release94|TRAV13-1|...
    Target: index 2 ('TRAV13-1')
    """
    genes = set()
    try:
        for record in SeqIO.parse(str(fasta_path), "fasta"):
            # Split by pipe
            parts = record.description.split('|')
            if len(parts) >= 3:
                # Based on your example, index 2 is the clean gene name
                gene_symbol = parts[2].strip()
                if gene_symbol.startswith(prefixes):
                    genes.add(gene_symbol)
    except Exception as e:
        print(f"Error reading Cell Ranger file {fasta_path}: {e}")
    return genes

def get_imgt_genes(fasta_path, prefixes):
    """
    Parses IMGT files.
    Example: >U04312|TRAV13-1*01|Mus musculus...
    Target: index 1 ('TRAV13-1*01'), then strip allele
    """
    genes = set()
    try:
        for record in SeqIO.parse(str(fasta_path), "fasta"):
            parts = record.description.split('|')
            if len(parts) >= 2:
                # Field 2 contains Gene*Allele; strip after '*'
                gene_symbol = parts[1].split('*')[0].strip()
                if gene_symbol.startswith(prefixes):
                    genes.add(gene_symbol)
    except Exception as e:
        print(f"Error reading IMGT file {fasta_path}: {e}")
    return genes

def main():
    # Defining the TCR subsets of interest
    tcr_prefixes = ('TRBV', 'TRAV', 'TRBJ', 'TRAJ')
    base_path = Path("references")
    
    if not base_path.exists():
        print(f"Error: '{base_path}' folder not found. Run script from the parent folder.")
        return

    # 1. Find all regions.fa files recursively in cellranger_genes
    cr_files = list(base_path.rglob("regions.fa"))
    
    # 2. Find all .fa files specifically in the imgt_genes_aa directory
    imgt_dir = base_path / "imgt_genes_aa"
    imgt_files = list(imgt_dir.glob("*.fa")) if imgt_dir.exists() else []

    if not cr_files:
        print("Warning: No 'regions.fa' files found.")
    if not imgt_files:
        print(f"Warning: No FASTA files found in {imgt_dir}")

    # Process Cell Ranger sources
    all_cr_genes = set()
    for f in cr_files:
        print(f"Processing Cell Ranger source: {f}")
        all_cr_genes.update(get_cellranger_genes(f, tcr_prefixes))
        
    # Process IMGT sources
    all_imgt_genes = set()
    for f in imgt_files:
        print(f"Processing IMGT source: {f.name}")
        all_imgt_genes.update(get_imgt_genes(f, tcr_prefixes))

    # Compare sets
    missing = sorted(list(all_cr_genes - all_imgt_genes))
    
    print("\n" + "="*60)
    print(f"Cell Ranger unique genes found: {len(all_cr_genes)}")
    print(f"IMGT unique genes found:        {len(all_imgt_genes)}")
    print("="*60)
    
    if missing:
        print(f"Found {len(missing)} genes in Cell Ranger NOT present in IMGT_aa files:")
        for gene in missing:
            print(f"  - {gene}")
    else:
        print("Everything matches! All specified TCR genes in Cell Ranger exist in IMGT_aa.")

if __name__ == "__main__":
    main()
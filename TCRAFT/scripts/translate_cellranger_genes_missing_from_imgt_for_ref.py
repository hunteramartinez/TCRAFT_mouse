import json
import os
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def get_translation_map(json_path):
    with open(json_path, 'r') as f:
        usage_data = json.load(f)
    translation_map = {}
    for aa, codon_list in usage_data.items():
        for codon_info in codon_list:
            codon = codon_info[0].upper()
            translation_map[codon] = aa
    return translation_map

def translate_custom(sequence, trans_map):
    protein = []
    for i in range(0, len(sequence) - (len(sequence) % 3), 3):
        codon = str(sequence[i:i+3]).upper()
        aa = trans_map.get(codon, "X")
        if aa == "*":
            break
        protein.append(aa)
    return "".join(protein)

def main():
    base_path = Path("references")
    codon_json = base_path / "codon_tables" / "mouse_codon_usage.json"
    output_dir = base_path / "imgt_genes_aa"
    
    # List of genes to extract and translate
    genes_to_process = [
        "TRAV3N-3","TRAV7N-6","TRAV13-4-DV7", "TRAV14D-3-DV8", "TRAV15-1-DV6-1", "TRAV15-2-DV6-2",
        "TRAV15D-1-DV6D-1", "TRAV16D-DV11", "TRAV21-DV12", "TRAV4-4-DV10",
        "TRAV6-7-DV9", "TRBV10", "TRBV12-2+TRBV13-2", "TRBV21", "TRBV26", "TRBV30"
    ]

    if not codon_json.exists():
        print(f"Error: Codon table not found.")
        return
    trans_map = get_translation_map(codon_json)

    cr_file = next(base_path.rglob("regions.fa"), None)
    if not cr_file:
        print("Error: regions.fa not found.")
        return

    # Dictionary to hold our records before writing
    # Key: filename (e.g., IMGT_TRAV_extra_aa.fa), Value: list of records
    files_to_write = {
        "IMGT_TRAV_extra_aa.fa": [],
        "IMGT_TRBV_extra_aa.fa": []
    }

    print(f"Scanning {cr_file} for target genes...")
    cr_records = {record.description.split('|')[2].strip(): record 
                  for record in SeqIO.parse(str(cr_file), "fasta") 
                  if len(record.description.split('|')) >= 3}

    for gene_name in genes_to_process:
        if gene_name in cr_records:
            nt_record = cr_records[gene_name]
            protein_str = translate_custom(nt_record.seq, trans_map)
            
            # Construct IMGT header
            full_header = (
                f"CUSTOM|{gene_name}*01|Mus musculus_Derived|F|L-PART1+V-REGION|"
                f"0..{len(nt_record.seq)}|{len(nt_record.seq)} nt|1||||"
                f"{len(protein_str)} AA||||"
            )
            
            new_record = SeqRecord(Seq(protein_str), id=full_header, description="")
            
            # Sort into TRAV or TRBV files
            if "TRAV" in gene_name:
                files_to_write["IMGT_TRAV_extra_aa.fa"].append(new_record)
            elif "TRBV" in gene_name:
                files_to_write["IMGT_TRBV_extra_aa.fa"].append(new_record)
            
            print(f"  [OK] Translated {gene_name}")
        else:
            print(f"  [!] Skipping {gene_name}: Not found in regions.fa")

    # Write the files
    os.makedirs(output_dir, exist_ok=True)
    for filename, records in files_to_write.items():
        if records:
            out_path = output_dir / filename
            with open(out_path, "w") as f:
                SeqIO.write(records, f, "fasta")
            print(f"Successfully created {out_path} with {len(records)} sequences.")

if __name__ == "__main__":
    main()
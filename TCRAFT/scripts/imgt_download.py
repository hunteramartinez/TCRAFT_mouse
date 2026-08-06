#!/usr/bin/env python3
"""
IMGT TRAV Downloader
Downloads T-cell Receptor Alpha Variable genes with CDR/FR annotations from IMGT

Usage:
    python imgt_download.py                           # Download human TRAV, TRAJ, TRAC
    python imgt_download.py --species mouse           # Download mouse genes
    python imgt_download.py --genes TRAV TRBV         # Download specific genes
    python imgt_download.py --output my_data          # Custom output directory
    python imgt_download.py --all-species             # Download all available species
"""

import requests
import os
import argparse
from pathlib import Path
from typing import Dict, List
import json
from datetime import datetime
import urllib3
import ssl
import certifi

# Disable SSL warnings (only use if you trust the source)
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

class IMGTDownloader:
    def __init__(self, output_dir: str = "imgt_data", verify_ssl: bool = False):
        """Initialize the IMGT downloader"""
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.verify_ssl = verify_ssl
        
        # IMGT reference directory URLs
        self.base_url = "https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory"
        
        # Alternative: try HTTP instead of HTTPS
        self.base_url_http = "http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory"
        
        # Species and gene options
        self.species_options = {
            "human": "Homo_sapiens",
            "mouse": "Mus_musculus",
            "rat": "Rattus_norvegicus"
        }
        
        # Setup session with retry logic
        self.session = requests.Session()
        retry_strategy = urllib3.util.retry.Retry(
            total=3,
            backoff_factor=1,
            status_forcelist=[429, 500, 502, 503, 504]
        )
        adapter = requests.adapters.HTTPAdapter(max_retries=retry_strategy)
        self.session.mount("http://", adapter)
        self.session.mount("https://", adapter)
        
    def download_file(self, url: str, use_http_fallback: bool = True) -> str:
        """
        Download file with SSL handling and HTTP fallback
        
        Args:
            url: URL to download
            use_http_fallback: If True, try HTTP if HTTPS fails
            
        Returns:
            File content as string
        """
        # Try HTTPS first (with SSL verification disabled for IMGT)
        try:
            print(f"  Attempting HTTPS download...")
            response = self.session.get(url, timeout=30, verify=self.verify_ssl)
            response.raise_for_status()
            return response.text
        except (requests.exceptions.SSLError, requests.exceptions.RequestException) as e:
            print(f"  HTTPS failed: {type(e).__name__}")
            
            # Try HTTP fallback
            if use_http_fallback:
                try:
                    http_url = url.replace("https://", "http://")
                    print(f"  Attempting HTTP fallback...")
                    response = self.session.get(http_url, timeout=30)
                    response.raise_for_status()
                    print(f"  ✓ HTTP fallback successful")
                    return response.text
                except requests.exceptions.RequestException as e2:
                    print(f"  HTTP fallback also failed: {type(e2).__name__}")
                    raise e2
            else:
                raise e
        
    def download_trav(self, species: str = "human", gene_types: List[str] = None):
        """
        Download TRAV genes for specified species
        
        Args:
            species: 'human', 'mouse', or 'rat'
            gene_types: List of gene types to download (default: ['TRAV', 'TRAJ', 'TRAC'])
        """
        if gene_types is None:
            gene_types = ['TRAV', 'TRAJ', 'TRAC']
        
        species_name = self.species_options.get(species.lower())
        if not species_name:
            raise ValueError(f"Species must be one of: {list(self.species_options.keys())}")
        
        results = {}
        
        for gene_type in gene_types:
            print(f"\n{'='*60}")
            print(f"Downloading {gene_type} for {species_name}...")
            print(f"{'='*60}")
            
            # Construct URL
            filename = f"{species_name}_{gene_type}.fasta"
            url = f"{self.base_url}/{filename}"
            
            try:
                # Download with fallback handling
                content = self.download_file(url)
                
                # Save FASTA file
                output_path = self.output_dir / filename
                with open(output_path, 'w') as f:
                    f.write(content)
                
                print(f"✓ Downloaded: {filename}")
                print(f"  Saved to: {output_path}")
                print(f"  Size: {len(content)} bytes")
                
                # Parse and extract CDR information
                sequences = self.parse_fasta(content)
                results[gene_type] = sequences
                
                print(f"  Found {len(sequences)} sequences")
                
                if len(sequences) > 0:
                    # Save parsed data as JSON
                    json_path = self.output_dir / f"{species_name}_{gene_type}_parsed.json"
                    with open(json_path, 'w') as f:
                        json.dump(sequences, f, indent=2)
                    print(f"  Parsed data saved to: {json_path}")
                else:
                    print(f"  ⚠ Warning: No sequences found in file")
                
            except requests.exceptions.RequestException as e:
                print(f"✗ Failed to download {filename}")
                print(f"  Error: {e}")
                print(f"  You may need to download manually from:")
                print(f"  {url}")
                results[gene_type] = None
        
        return results
    
    def parse_fasta(self, fasta_content: str) -> List[Dict]:
        """Parse FASTA content and extract sequence information"""
        sequences = []
        current_seq = None
        
        for line in fasta_content.split('\n'):
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    sequences.append(current_seq)
                
                # Parse header
                header = line[1:]
                current_seq = {
                    'header': header,
                    'sequence': '',
                    'metadata': self.parse_header(header)
                }
            elif line and current_seq is not None:
                current_seq['sequence'] += line
        
        if current_seq:
            sequences.append(current_seq)
        
        return sequences
    
    def parse_header(self, header: str) -> Dict:
        """Parse IMGT FASTA header to extract metadata"""
        parts = header.split('|')
        metadata = {
            'full_header': header,
            'gene_name': parts[1] if len(parts) > 1 else '',
            'allele': parts[1] if len(parts) > 1 else '',
        }
        
        # Extract additional information from header
        if len(parts) > 1:
            metadata['accession'] = parts[0]
            metadata['gene_allele'] = parts[1]
        
        return metadata
    
    def create_summary_report(self, results: Dict, species: str):
        """Create a summary report of downloaded data"""
        report_path = self.output_dir / f"download_report_{species}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
        
        with open(report_path, 'w') as f:
            f.write("IMGT Download Report\n")
            f.write("="*60 + "\n")
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Species: {species}\n")
            f.write(f"Output Directory: {self.output_dir}\n\n")
            
            total_success = 0
            total_failed = 0
            
            for gene_type, sequences in results.items():
                if sequences:
                    total_success += 1
                    f.write(f"\n{gene_type}: ✓ SUCCESS\n")
                    f.write(f"  Total sequences: {len(sequences)}\n")
                    
                    # List first few genes
                    if len(sequences) > 0:
                        f.write(f"  Sample genes:\n")
                        for seq in sequences[:5]:
                            gene_name = seq['metadata'].get('gene_allele', 'Unknown')
                            seq_len = len(seq['sequence'])
                            f.write(f"    - {gene_name} ({seq_len} bp)\n")
                        if len(sequences) > 5:
                            f.write(f"    ... and {len(sequences) - 5} more\n")
                else:
                    total_failed += 1
                    f.write(f"\n{gene_type}: ✗ FAILED\n")
            
            f.write(f"\n{'='*60}\n")
            f.write(f"Summary: {total_success} successful, {total_failed} failed\n")
        
        print(f"\n✓ Summary report saved to: {report_path}")
        return report_path


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Download TRAV genes with CDR annotations from IMGT',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s                                    # Download human TRAV, TRAJ, TRAC
  %(prog)s --species mouse                    # Download mouse genes
  %(prog)s --genes TRAV TRBV                  # Download specific genes
  %(prog)s --output my_data                   # Custom output directory
  %(prog)s --all-species                      # Download all species
  %(prog)s --species human --genes TRAV       # Human TRAV only
  %(prog)s --list-options                     # Show available options

Note: This script disables SSL verification for IMGT due to certificate issues.
Downloads are made over HTTP as fallback if HTTPS fails.

For detailed CDR annotations, visit:
  - IMGT/V-QUEST: http://www.imgt.org/IMGT_vquest/
  - IMGT/GENE-DB: http://www.imgt.org/genedb/
        """
    )
    
    parser.add_argument(
        '--species',
        type=str,
        default='human',
        choices=['human', 'mouse', 'rat'],
        help='Species to download (default: human)'
    )
    
    parser.add_argument(
        '--genes',
        nargs='+',
        default=['TRAV', 'TRAJ', 'TRAC'],
        help='Gene types to download (default: TRAV TRAJ TRAC)'
    )
    
    parser.add_argument(
        '--output', '-o',
        type=str,
        default='imgt_data',
        help='Output directory (default: imgt_data)'
    )
    
    parser.add_argument(
        '--all-species',
        action='store_true',
        help='Download for all available species'
    )
    
    parser.add_argument(
        '--list-options',
        action='store_true',
        help='List available species and gene types'
    )
    
    parser.add_argument(
        '--quiet', '-q',
        action='store_true',
        help='Minimal output'
    )
    
    parser.add_argument(
        '--verify-ssl',
        action='store_true',
        help='Enable SSL certificate verification (may fail with IMGT)'
    )
    
    return parser.parse_args()


def main():
    """Main execution function"""
    args = parse_arguments()
    
    # Handle list options
    if args.list_options:
        print("\n" + "="*60)
        print("Available Options")
        print("="*60)
        print("\nSpecies:")
        print("  - human (Homo sapiens)")
        print("  - mouse (Mus musculus)")
        print("  - rat (Rattus norvegicus)")
        print("\nCommon Gene Types:")
        print("  TRA: TRAV, TRAJ, TRAC")
        print("  TRB: TRBV, TRBD, TRBJ, TRBC")
        print("  TRG: TRGV, TRGJ, TRGC")
        print("  TRD: TRDV, TRDD, TRDJ, TRDC")
        print("\n")
        return
    
    if not args.quiet:
        print("="*60)
        print("IMGT TRAV Downloader")
        print("="*60)
        if not args.verify_ssl:
            print("⚠ SSL verification disabled (IMGT certificate issue)")
            print("  Files will be downloaded over HTTP if HTTPS fails")
        print()
    
    # Initialize downloader
    downloader = IMGTDownloader(output_dir=args.output, verify_ssl=args.verify_ssl)
    
    # Download for specified species or all species
    species_list = ['human', 'mouse', 'rat'] if args.all_species else [args.species]
    
    all_results = {}
    for species in species_list:
        if not args.quiet:
            print(f"\n{'='*60}")
            print(f"Processing {species.upper()}")
            print(f"{'='*60}")
        
        try:
            results = downloader.download_trav(
                species=species,
                gene_types=args.genes
            )
            all_results[species] = results
            
            # Create summary report
            downloader.create_summary_report(results, species)
            
        except Exception as e:
            print(f"✗ Error processing {species}: {e}")
            continue
    
    # Print final summary
    if not args.quiet:
        print("\n" + "="*60)
        print("Download Complete!")
        print("="*60)
        print(f"Output directory: {downloader.output_dir}")
        print("\nFiles downloaded:")
        total_sequences = 0
        for species, results in all_results.items():
            for gene_type, sequences in results.items():
                if sequences:
                    print(f"  ✓ {species} {gene_type}: {len(sequences)} sequences")
                    total_sequences += len(sequences)
                else:
                    print(f"  ✗ {species} {gene_type}: FAILED")
        
        print(f"\nTotal sequences downloaded: {total_sequences}")
        
        print("\n" + "="*60)
        print("Additional Resources:")
        print("="*60)
        print("For detailed CDR annotations:")
        print("  • IMGT/V-QUEST: http://www.imgt.org/IMGT_vquest/")
        print("  • IMGT/DomainGapAlign: http://www.imgt.org/3Dstructure-DB/cgi/DomainGapAlign.cgi")
        print("  • IMGT/GENE-DB: http://www.imgt.org/genedb/")
        print()


if __name__ == "__main__":
    main()
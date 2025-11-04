#!/usr/bin/env python3
"""
Extract sequences for each SD pair from assembly or reference FASTAs
"""

import pandas as pd
import sys
import os
import re
from pathlib import Path

def parse_region(region_str):
    """
    Parse region string like 'chr13:38070105-38072275'
    Returns (chrom, start, end)
    """
    if pd.isna(region_str) or region_str == "None" or region_str == "":
        return None, None, None

    match = re.match(r'(.+):(\d+)-(\d+)', str(region_str))
    if match:
        chrom = match.group(1)
        start = int(match.group(2))
        end = int(match.group(3))
        return chrom, start, end
    return None, None, None

def extract_sequence_samtools(fasta_file, chrom, start, end, output_file):
    """
    Extract sequence region using samtools faidx
    """
    region = f"{chrom}:{start}-{end}"

    # Index fasta if needed
    if not os.path.exists(f"{fasta_file}.fai"):
        print(f"  Indexing {fasta_file}...")
        os.system(f"samtools faidx {fasta_file}")

    # Extract sequence
    cmd = f"samtools faidx {fasta_file} {region} > {output_file}"
    ret = os.system(cmd)

    if ret != 0:
        print(f"  ERROR: Failed to extract {region} from {fasta_file}")
        return False

    # Check if output file was created and has content
    if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
        return True
    else:
        print(f"  ERROR: Output file {output_file} is empty or not created")
        return False

def extract_sequence_pysam(fasta_file, chrom, start, end, output_file):
    """
    Extract sequence region using pysam (fallback if samtools not available)
    """
    try:
        import pysam
    except ImportError:
        print("ERROR: Neither samtools nor pysam available for sequence extraction")
        print("Please install: conda install samtools or pip install pysam")
        sys.exit(1)

    try:
        fasta = pysam.FastaFile(fasta_file)
        sequence = fasta.fetch(chrom, start, end)

        with open(output_file, 'w') as f:
            f.write(f">{chrom}:{start}-{end}\n")
            # Write sequence in 60 character lines
            for i in range(0, len(sequence), 60):
                f.write(sequence[i:i+60] + "\n")

        return True
    except Exception as e:
        print(f"  ERROR: Failed to extract {chrom}:{start}-{end} from {fasta_file}")
        print(f"  {e}")
        return False

def find_assembly_file(sample, haplotype, assemblies_config):
    """
    Find assembly FASTA file for a given sample and haplotype
    """
    # Check if directly specified
    key = f"{sample}_{haplotype}"
    if key in assemblies_config:
        return assemblies_config[key]

    # Check if assembly_dir and pattern specified
    if 'assembly_dir' in assemblies_config and 'assembly_pattern' in assemblies_config:
        assembly_dir = assemblies_config['assembly_dir']
        pattern = assemblies_config['assembly_pattern']

        # Replace placeholders
        filename = pattern.replace('{sample}', sample).replace('{haplotype}', haplotype)
        filepath = os.path.join(assembly_dir, filename)

        if os.path.exists(filepath):
            return filepath

        # Try common extensions
        for ext in ['.fa', '.fasta', '.fa.gz', '.fasta.gz']:
            base = filename.rsplit('.', 1)[0] if '.' in filename else filename
            test_path = os.path.join(assembly_dir, base + ext)
            if os.path.exists(test_path):
                return test_path

    return None

def main(snakemake):
    """Main sequence extraction function"""

    # Load filtered SDs
    df = pd.read_csv(snakemake.input.filtered_sds, sep="\t")

    print(f"Extracting sequences for {len(df)} SD pairs...")

    # Get configuration
    assemblies_config = snakemake.params.assemblies
    reference_config = snakemake.params.reference
    coord_system = snakemake.params.coord_system

    # Create output directory
    os.makedirs(snakemake.output.sequence_dir, exist_ok=True)

    # Track which sequences were successfully extracted
    manifest_data = []

    # Check if we can use samtools or need pysam
    use_samtools = os.system("which samtools > /dev/null 2>&1") == 0

    if use_samtools:
        print("Using samtools for sequence extraction")
        extract_func = extract_sequence_samtools
    else:
        print("Samtools not found, using pysam")
        extract_func = extract_sequence_pysam

    # Process each SD pair
    for idx, row in df.iterrows():
        sd_id = row['SD_ID']
        sample = row['Sample']
        haplotype = row['Haplotype']
        coord_to_use = row.get('coord_to_use', 'contig')

        print(f"\nProcessing {sd_id} ({sample} {haplotype})...")
        print(f"  Using {coord_to_use} coordinates")

        # Determine which coordinates and FASTA to use
        if coord_to_use == 't2t':
            # Use T2T coordinates and reference FASTA
            sd1_region = row.get('SD1_T2T')
            sd2_region = row.get('SD2_T2T')
            fasta_file = reference_config.get('fasta')

            if not fasta_file or not os.path.exists(fasta_file):
                print(f"  ERROR: T2T reference not found: {fasta_file}")
                print(f"  Skipping {sd_id}")
                continue

        else:
            # Use contig coordinates and assembly FASTA
            sd1_region = row['SD1_original']
            sd2_region = row['SD2_original']
            fasta_file = find_assembly_file(sample, haplotype, assemblies_config)

            if not fasta_file:
                print(f"  ERROR: Assembly not found for {sample}_{haplotype}")
                print(f"  Skipping {sd_id}")
                continue

        print(f"  Using FASTA: {fasta_file}")
        print(f"  SD1: {sd1_region}")
        print(f"  SD2: {sd2_region}")

        # Parse coordinates
        sd1_chrom, sd1_start, sd1_end = parse_region(sd1_region)
        sd2_chrom, sd2_start, sd2_end = parse_region(sd2_region)

        if not sd1_chrom or not sd2_chrom:
            print(f"  ERROR: Could not parse coordinates")
            print(f"  Skipping {sd_id}")
            continue

        # Output files
        sd1_fasta = os.path.join(snakemake.output.sequence_dir, f"{sd_id}_SD1.fa")
        sd2_fasta = os.path.join(snakemake.output.sequence_dir, f"{sd_id}_SD2.fa")

        # Extract sequences
        print(f"  Extracting SD1...")
        success1 = extract_func(fasta_file, sd1_chrom, sd1_start, sd1_end, sd1_fasta)

        print(f"  Extracting SD2...")
        success2 = extract_func(fasta_file, sd2_chrom, sd2_start, sd2_end, sd2_fasta)

        if success1 and success2:
            print(f"  ✓ Successfully extracted both sequences")
            manifest_data.append({
                'SD_ID': sd_id,
                'Sample': sample,
                'Haplotype': haplotype,
                'SD1_region': sd1_region,
                'SD2_region': sd2_region,
                'SD1_fasta': sd1_fasta,
                'SD2_fasta': sd2_fasta,
                'coord_system': coord_to_use,
                'source_fasta': fasta_file
            })
        else:
            print(f"  ✗ Failed to extract sequences for {sd_id}")

    # Save manifest
    manifest_df = pd.DataFrame(manifest_data)
    manifest_df.to_csv(snakemake.output.manifest, sep="\t", index=False)

    print(f"\nSuccessfully extracted {len(manifest_df)} / {len(df)} SD pairs")
    print(f"Manifest saved to: {snakemake.output.manifest}")

    if len(manifest_df) == 0:
        print("\nERROR: No sequences were extracted!")
        print("Please check:")
        print("1. Assembly FASTA files exist and are in the correct location")
        print("2. Contig names in SD calls match those in FASTA files")
        print("3. samtools or pysam is installed")
        sys.exit(1)

if __name__ == "__main__":
    main(snakemake)

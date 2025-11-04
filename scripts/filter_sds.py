#!/usr/bin/env python3
"""
Filter SD calls to select which ones to visualize
"""

import pandas as pd
import sys
import re

def parse_region(region_str):
    """
    Parse region string like 'chr13:38070105-38072275' or 'h1tg000008l:61214599-61216770'
    Returns (chrom, start, end, size)
    """
    if pd.isna(region_str) or region_str == "None" or region_str == "":
        return None, None, None, 0

    match = re.match(r'(.+):(\d+)-(\d+)', str(region_str))
    if match:
        chrom = match.group(1)
        start = int(match.group(2))
        end = int(match.group(3))
        size = end - start
        return chrom, start, end, size
    return None, None, None, 0

def main(snakemake):
    """Main filtering function"""

    # Load configuration
    max_sds = snakemake.params.max_sds
    filter_opts = snakemake.params.filter_opts
    coord_system = snakemake.params.coord_system

    print(f"Loading SD calls from: {snakemake.input.sd_calls}")

    # Read SD calls
    df = pd.read_csv(snakemake.input.sd_calls, sep="\t")

    print(f"Total SD pairs in input: {len(df)}")

    # Create unique SD_ID for each pair
    df['SD_ID'] = [f"SD_{i:05d}" for i in range(len(df))]

    # Parse regions and calculate sizes
    print("Parsing coordinates...")

    # Parse SD1_original
    df[['SD1_chrom', 'SD1_start', 'SD1_end', 'SD1_size']] = df['SD1_original'].apply(
        lambda x: pd.Series(parse_region(x))
    )

    # Parse SD2_original
    df[['SD2_chrom', 'SD2_start', 'SD2_end', 'SD2_size']] = df['SD2_original'].apply(
        lambda x: pd.Series(parse_region(x))
    )

    # Parse T2T coordinates if available
    if 'SD1_T2T' in df.columns:
        df[['SD1_T2T_chrom', 'SD1_T2T_start', 'SD1_T2T_end', 'SD1_T2T_size']] = df['SD1_T2T'].apply(
            lambda x: pd.Series(parse_region(x))
        )
    if 'SD2_T2T' in df.columns:
        df[['SD2_T2T_chrom', 'SD2_T2T_start', 'SD2_T2T_end', 'SD2_T2T_size']] = df['SD2_T2T'].apply(
            lambda x: pd.Series(parse_region(x))
        )

    # Calculate average SD size
    df['avg_size'] = (df['SD1_size'] + df['SD2_size']) / 2

    # Add priority score
    df['priority_score'] = 0

    # Apply filters
    initial_count = len(df)

    # Min size filter
    min_size = filter_opts.get('min_sd_size', 0)
    if min_size > 0:
        df = df[df['avg_size'] >= min_size]
        print(f"After min_size filter ({min_size} bp): {len(df)} SDs")

    # Require genes filter
    if filter_opts.get('require_genes', False):
        if 'SD1_has_exon' in df.columns and 'SD2_has_exon' in df.columns:
            df = df[(df['SD1_has_exon'] == True) | (df['SD2_has_exon'] == True)]
            print(f"After require_genes filter: {len(df)} SDs")

    # Require T2T liftover filter
    if filter_opts.get('require_t2t_liftover', False):
        if 'SD1_T2T' in df.columns and 'SD2_T2T' in df.columns:
            df = df[(df['SD1_T2T'].notna()) & (df['SD2_T2T'].notna())]
            df = df[(df['SD1_T2T'] != "None") & (df['SD2_T2T'] != "None")]
            print(f"After T2T liftover filter: {len(df)} SDs")

    # Sample filter
    if filter_opts.get('samples'):
        samples = filter_opts['samples']
        df = df[df['Sample'].isin(samples)]
        print(f"After sample filter: {len(df)} SDs")

    # Haplotype filter
    if filter_opts.get('haplotypes'):
        haplotypes = filter_opts['haplotypes']
        df = df[df['Haplotype'].isin(haplotypes)]
        print(f"After haplotype filter: {len(df)} SDs")

    # Calculate priority scores
    priority = filter_opts.get('priority', 'genes')

    if priority == 'genes':
        # Prioritize SDs with genes
        if 'SD1_has_exon' in df.columns and 'SD2_has_exon' in df.columns:
            df['priority_score'] += df['SD1_has_exon'].astype(int) * 100
            df['priority_score'] += df['SD2_has_exon'].astype(int) * 100
        # Secondary: size
        df['priority_score'] += df['avg_size'] / 10000

    elif priority == 'size':
        # Prioritize larger SDs
        df['priority_score'] = df['avg_size']

    elif priority == 'random':
        # Random priority
        import numpy as np
        np.random.seed(42)
        df['priority_score'] = np.random.rand(len(df))

    # Sort by priority
    df = df.sort_values('priority_score', ascending=False)

    # Select top N
    if len(df) > max_sds:
        print(f"\nSelecting top {max_sds} SDs based on priority: {priority}")
        df = df.head(max_sds)

    print(f"\nFinal SD count: {len(df)}")

    # Determine which coordinate system to use for each SD
    df['coord_to_use'] = 'contig'
    if coord_system == 't2t':
        df['coord_to_use'] = 't2t'
    elif coord_system == 'auto':
        # Use T2T if available, otherwise contig
        if 'SD1_T2T' in df.columns and 'SD2_T2T' in df.columns:
            has_t2t = (df['SD1_T2T'].notna()) & (df['SD2_T2T'].notna())
            has_t2t &= (df['SD1_T2T'] != "None") & (df['SD2_T2T'] != "None")
            df.loc[has_t2t, 'coord_to_use'] = 't2t'

    # Save filtered SDs
    df.to_csv(snakemake.output.filtered, sep="\t", index=False)
    print(f"\nSaved filtered SDs to: {snakemake.output.filtered}")

    # Write stats
    with open(snakemake.output.stats, 'w') as f:
        f.write("SD Filtering Statistics\n")
        f.write("=" * 50 + "\n")
        f.write(f"Input SDs: {initial_count}\n")
        f.write(f"Filtered SDs: {len(df)}\n")
        f.write(f"Max SDs requested: {max_sds}\n")
        f.write(f"Priority: {priority}\n")
        f.write(f"Coordinate system: {coord_system}\n")
        f.write("\nFilters applied:\n")
        for key, value in filter_opts.items():
            f.write(f"  {key}: {value}\n")
        f.write("\nSelected SDs:\n")
        for idx, row in df.iterrows():
            genes = row.get('Gene_names', 'none')
            f.write(f"  {row['SD_ID']}: {row['Sample']} {row['Haplotype']}, "
                   f"size={row['avg_size']:.0f}bp, genes={genes}\n")

    print(f"Saved statistics to: {snakemake.output.stats}")

if __name__ == "__main__":
    main(snakemake)

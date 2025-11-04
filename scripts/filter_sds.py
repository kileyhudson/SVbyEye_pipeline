#!/usr/bin/env python3
"""Utility for filtering and prioritising segmental duplication (SD) pairs."""

from __future__ import annotations

import re
from typing import Optional, Tuple

import pandas as pd


def parse_region(region_str: Optional[str]) -> Tuple[Optional[str], Optional[int], Optional[int], int]:
    """Parse a genomic interval string into components.

    Parameters
    ----------
    region_str:
        String of the form ``chrom:start-end``. ``None`` and empty strings are
        treated as missing values.

    Returns
    -------
    tuple
        A tuple containing the chromosome name, start coordinate, end
        coordinate, and computed interval length. Missing values are returned as
        ``None`` with a length of ``0``.
    """

    if region_str is None or pd.isna(region_str) or region_str in {"", "None"}:
        return None, None, None, 0

    match = re.match(r"(.+):(\d+)-(\d+)", str(region_str))
    if not match:
        return None, None, None, 0

    chrom = match.group(1)
    start = int(match.group(2))
    end = int(match.group(3))
    return chrom, start, end, max(end - start, 0)


def main(snakemake) -> None:  # type: ignore[annotation-unchecked]
    """Filter SD calls and compute prioritisation metrics."""

    max_sds = snakemake.params.max_sds
    filter_opts = snakemake.params.filter_opts
    coord_system = snakemake.params.coord_system

    print(f"Loading SD callset: {snakemake.input.sd_calls}")
    df = pd.read_csv(snakemake.input.sd_calls, sep="\t")
    print(f"Total SD pairs provided: {len(df)}")

    df["SD_ID"] = [f"SD_{i:05d}" for i in range(len(df))]
    print("Deriving coordinate summaries...")

    df[["SD1_chrom", "SD1_start", "SD1_end", "SD1_size"]] = df["SD1_original"].apply(
        lambda value: pd.Series(parse_region(value))
    )
    df[["SD2_chrom", "SD2_start", "SD2_end", "SD2_size"]] = df["SD2_original"].apply(
        lambda value: pd.Series(parse_region(value))
    )

    if "SD1_T2T" in df.columns:
        df[["SD1_T2T_chrom", "SD1_T2T_start", "SD1_T2T_end", "SD1_T2T_size"]] = df[
            "SD1_T2T"
        ].apply(lambda value: pd.Series(parse_region(value)))
    if "SD2_T2T" in df.columns:
        df[["SD2_T2T_chrom", "SD2_T2T_start", "SD2_T2T_end", "SD2_T2T_size"]] = df[
            "SD2_T2T"
        ].apply(lambda value: pd.Series(parse_region(value)))

    df["avg_size"] = (df["SD1_size"].fillna(0) + df["SD2_size"].fillna(0)) / 2
    df["priority_score"] = 0.0

    initial_count = len(df)

    min_size = filter_opts.get("min_sd_size", 0)
    if min_size > 0:
        df = df[df["avg_size"] >= min_size]
        print(f"After applying minimum size threshold ({min_size} bp): {len(df)}")

    if filter_opts.get("require_genes", False) and {"SD1_has_exon", "SD2_has_exon"} <= set(df.columns):
        df = df[(df["SD1_has_exon"] == True) | (df["SD2_has_exon"] == True)]  # noqa: E712
        print(f"After requiring gene content: {len(df)}")

    if filter_opts.get("require_t2t_liftover", False) and {"SD1_T2T", "SD2_T2T"} <= set(df.columns):
        mask = df["SD1_T2T"].notna() & df["SD2_T2T"].notna()
        mask &= (df["SD1_T2T"] != "None") & (df["SD2_T2T"] != "None")
        df = df[mask]
        print(f"After requiring T2T liftover: {len(df)}")

    samples = filter_opts.get("samples") or []
    if samples:
        df = df[df["Sample"].isin(samples)]
        print(f"After restricting to selected samples: {len(df)}")

    haplotypes = filter_opts.get("haplotypes") or []
    if haplotypes:
        df = df[df["Haplotype"].isin(haplotypes)]
        print(f"After restricting to selected haplotypes: {len(df)}")

    priority = filter_opts.get("priority", "genes")
    if priority == "genes" and {"SD1_has_exon", "SD2_has_exon"} <= set(df.columns):
        df["priority_score"] += df["SD1_has_exon"].fillna(False).astype(int) * 100
        df["priority_score"] += df["SD2_has_exon"].fillna(False).astype(int) * 100
        df["priority_score"] += df["avg_size"] / 10000
    elif priority == "size":
        df["priority_score"] = df["avg_size"]
    elif priority == "random":
        import numpy as np

        np.random.seed(42)
        df["priority_score"] = np.random.rand(len(df))
    else:
        df["priority_score"] = df["avg_size"]

    df = df.sort_values("priority_score", ascending=False)

    if len(df) > max_sds:
        print(f"Selecting the top {max_sds} SD pairs using '{priority}' priority")
        df = df.head(max_sds)

    print(f"Final SD count: {len(df)}")

    df["coord_to_use"] = "contig"
    if coord_system == "t2t":
        df["coord_to_use"] = "t2t"
    elif coord_system == "auto" and {"SD1_T2T", "SD2_T2T"} <= set(df.columns):
        has_t2t = df["SD1_T2T"].notna() & df["SD2_T2T"].notna()
        has_t2t &= (df["SD1_T2T"] != "None") & (df["SD2_T2T"] != "None")
        df.loc[has_t2t, "coord_to_use"] = "t2t"

    df.to_csv(snakemake.output.filtered, sep="\t", index=False)
    print(f"Filtered SD list written to: {snakemake.output.filtered}")

    with open(snakemake.output.stats, "w", encoding="utf-8") as handle:
        handle.write("SD Filtering Statistics\n")
        handle.write("=" * 50 + "\n")
        handle.write(f"Input SDs: {initial_count}\n")
        handle.write(f"Filtered SDs: {len(df)}\n")
        handle.write(f"Max SDs requested: {max_sds}\n")
        handle.write(f"Priority: {priority}\n")
        handle.write(f"Coordinate system: {coord_system}\n")
        handle.write("\nFilters applied:\n")
        for key, value in filter_opts.items():
            handle.write(f"  {key}: {value}\n")
        handle.write("\nSelected SDs:\n")
        for _, row in df.iterrows():
            genes = row.get("Gene_names", "none")
            handle.write(
                f"  {row['SD_ID']}: {row['Sample']} {row['Haplotype']}, "
                f"size={row['avg_size']:.0f}bp, genes={genes}\n"
            )

    print(f"Saved statistics to: {snakemake.output.stats}")

if __name__ == "__main__":
    main(snakemake)

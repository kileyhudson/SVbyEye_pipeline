#!/usr/bin/env python3
"""Assign stable SD identifiers and derive helpful metadata without filtering."""

import csv
import re
from dataclasses import dataclass
from typing import Iterable, List, Optional, Tuple


@dataclass
class RegionInfo:
    chrom: Optional[str]
    start: Optional[int]
    end: Optional[int]
    length: int


REGION_PATTERN = re.compile(r"(.+):(\d+)-(\d+)")


def parse_region(value: Optional[str]) -> RegionInfo:
    """Parse ``chrom:start-end`` strings and return coordinate metadata."""

    if value in (None, "", "None"):
        return RegionInfo(None, None, None, 0)

    match = REGION_PATTERN.match(str(value))
    if not match:
        return RegionInfo(None, None, None, 0)

    chrom = match.group(1)
    start = int(match.group(2))
    end = int(match.group(3))
    length = max(end - start, 0)
    return RegionInfo(chrom, start, end, length)


def load_sd_calls(path: str) -> List[dict]:
    """Read the SD callset as a list of dictionaries."""

    with open(path, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [dict(row) for row in reader]


def write_sd_table(path: str, rows: Iterable[dict], fieldnames: List[str]) -> None:
    """Write the enriched SD table to ``path`` preserving tab separation."""

    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_summary(path: str, sd_rows: List[dict]) -> None:
    """Emit a lightweight summary describing the prepared SD pairs."""

    total = len(sd_rows)
    liftover_count = sum(
        1
        for row in sd_rows
        if row.get("SD1_T2T") not in (None, "", "None")
        and row.get("SD2_T2T") not in (None, "", "None")
    )

    with open(path, "w", encoding="utf-8") as handle:
        handle.write("SVbyEye SD Preparation Summary\n")
        handle.write("=" * 40 + "\n")
        handle.write(f"Total SD pairs: {total}\n")
        handle.write(f"Pairs with T2T coordinates: {liftover_count}\n")
        handle.write("\nSD identifiers:\n")
        for row in sd_rows:
            handle.write(f"  {row['SD_ID']}: {row['Sample']} {row['Haplotype']}\n")


def main(snakemake) -> None:  # type: ignore[annotation-unchecked]
    """Attach SD_IDs and derived columns while keeping every input row."""

    sd_rows = load_sd_calls(snakemake.input.sd_calls)
    if not sd_rows:
        raise SystemExit("No SD rows found; populate data/sd_calls.tsv before running.")

    enriched: List[dict] = []
    for idx, row in enumerate(sd_rows):
        sd_id = f"SD_{idx:05d}"
        sd1_info = parse_region(row.get("SD1_original"))
        sd2_info = parse_region(row.get("SD2_original"))

        row = dict(row)  # make a mutable copy
        row["SD_ID"] = sd_id
        row["SD1_size"] = sd1_info.length
        row["SD2_size"] = sd2_info.length
        row["avg_size"] = (sd1_info.length + sd2_info.length) / 2
        enriched.append(row)

    # Preserve original columns and append the derived ones at the end.
    base_fields = list(sd_rows[0].keys())
    extra_fields = [field for field in ("SD_ID", "SD1_size", "SD2_size", "avg_size") if field not in base_fields]
    fieldnames = base_fields + extra_fields

    write_sd_table(snakemake.output.table, enriched, fieldnames)
    write_summary(snakemake.output.summary, enriched)
    print(f"Prepared SD table with {len(enriched)} entries.")


if __name__ == "__main__":
    main(snakemake)

#!/usr/bin/env python3
"""Generate synthetic minimap2-like PAF alignments used by automated tests."""

from __future__ import annotations

import csv
import os
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


def parse_region(region: Optional[str]) -> Optional[Tuple[str, int, int]]:
    """Parse ``chrom:start-end`` strings into structured coordinates."""

    if region in (None, "", "None"):
        return None

    chrom, coords = region.split(":")
    start_str, end_str = coords.split("-")
    start = int(start_str)
    end = int(end_str)
    if end < start:
        start, end = end, start
    return chrom, start, end


def region_length(region: Optional[str]) -> int:
    """Return the length of a genomic interval string."""

    parsed = parse_region(region)
    if parsed is None:
        return 0
    _, start, end = parsed
    return max(end - start, 0)


def make_paf_record(
    sd_id: str,
    sample: str,
    haplotype: str,
    query_region: Optional[str],
    target_region: Optional[str],
    query_label: str,
    target_label: str,
) -> List[str]:
    """Construct a single PAF alignment line approximating minimap2 output."""

    qlen = region_length(query_region) or 1000
    tlen = region_length(target_region) or qlen

    qname = f"{sample}_{haplotype}_{sd_id}_{query_label}"
    tname = f"{sample}_{haplotype}_{sd_id}_{target_label}"

    aln_span = min(qlen, tlen)
    nmatch = max(aln_span - 10, int(aln_span * 0.95))

    cg = f"{nmatch}={aln_span - nmatch}X" if nmatch < aln_span else f"{aln_span}="

    return [
        qname,
        str(qlen),
        "0",
        str(aln_span),
        "+",
        tname,
        str(tlen),
        "0",
        str(aln_span),
        str(nmatch),
        str(aln_span),
        "60",
        f"cg:Z:{cg}",
    ]


def read_sd_calls(filepath: str) -> List[Dict[str, str]]:
    """Load the SD callset describing the synthetic test pairs."""

    with open(filepath, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def main() -> None:
    """Create one PAF per SD pair defined in ``data/sd_calls.tsv``."""

    sd_calls = read_sd_calls("data/sd_calls.tsv")
    align_dir = Path("data/alignments")
    align_dir.mkdir(parents=True, exist_ok=True)

    print("Generating synthetic PAF alignments for integration tests")
    for idx, row in enumerate(sd_calls):
        sd_id = f"SD_{idx:05d}"
        paf_path = align_dir / f"{sd_id}.paf"

        records: Iterable[List[str]] = [
            make_paf_record(sd_id, row["Sample"], row["Haplotype"], row["SD1_original"], row["SD2_original"], "SD1", "SD2"),
            make_paf_record(sd_id, row["Sample"], row["Haplotype"], row["SD2_original"], row["SD1_original"], "SD2", "SD1"),
        ]

        with open(paf_path, "w", encoding="utf-8") as handle:
            for record in records:
                handle.write("\t".join(record) + "\n")

        print(f"  Wrote {paf_path.relative_to(Path('.'))}")

    produced = sorted(str(p.relative_to(Path("."))) for p in align_dir.glob("*.paf"))
    print("\nSynthetic alignments:")
    for path in produced:
        size = os.path.getsize(path)
        print(f"  {path} ({size} bytes)")


if __name__ == "__main__":
    main()

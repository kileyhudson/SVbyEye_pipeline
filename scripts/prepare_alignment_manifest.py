#!/usr/bin/env python3
"""Validate precomputed PAF alignments for each SD pair and emit a manifest."""

import csv
from pathlib import Path
from typing import Dict, Iterable, List


def render_pattern(pattern: str, context: Dict[str, str]) -> str:
    """Render filename patterns, surfacing a helpful error on bad placeholders."""

    try:
        return pattern.format(**context)
    except KeyError as exc:
        missing = exc.args[0]
        raise ValueError(
            f"Pattern '{pattern}' references unknown placeholder '{missing}'. "
            "Allowed keys: {sd_id}, {sample}, {haplotype}, {index}."
        ) from exc


def load_sd_table(path: str) -> Iterable[Dict[str, str]]:
    """Read the SD table without depending on pandas."""

    with open(path, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            yield row


def main(snakemake) -> None:  # type: ignore[annotation-unchecked]
    """Ensure all SDs map to an existing PAF and persist the mapping."""

    align_cfg: Dict[str, str] = snakemake.params.alignments or {}
    align_dir = Path(align_cfg.get("dir", "data/alignments"))
    pattern = align_cfg.get("filename_pattern", "{sd_id}.paf")

    rows = list(load_sd_table(snakemake.input.sd_table))
    print(f"Validating alignments for {len(rows)} SD pairs")

    manifest: List[List[str]] = []
    missing: List[str] = []

    for idx, row in enumerate(rows):
        context = {
            "sd_id": row["SD_ID"],
            "sample": row["Sample"],
            "haplotype": row["Haplotype"],
            "index": f"{idx:05d}",
        }
        paf_path = align_dir / render_pattern(pattern, context)

        if not paf_path.exists():
            missing.append(f"{row['SD_ID']} -> {paf_path}")
            continue

        manifest.append([row["SD_ID"], row["Sample"], row["Haplotype"], str(paf_path)])

    if missing:
        print("ERROR: Missing PAF files for the following SD IDs:")
        for record in missing:
            print(f"  {record}")
        raise SystemExit(f"{len(missing)} PAF file(s) are missing. Populate data/alignments before running.")

    with open(snakemake.output.manifest, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["SD_ID", "Sample", "Haplotype", "paf_path"])
        writer.writerows(manifest)

    print(f"Manifest written to: {snakemake.output.manifest}")


if __name__ == "__main__":
    main(snakemake)

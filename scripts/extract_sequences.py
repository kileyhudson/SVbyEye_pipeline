#!/usr/bin/env python3
"""Extract FASTA segments for segmental duplication pairs."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Callable, Dict, Optional, Tuple

import pandas as pd


def parse_region(region_str: Optional[str]) -> Tuple[Optional[str], Optional[int], Optional[int]]:
    """Parse genomic interval strings of the form ``chrom:start-end``."""

    if region_str is None or pd.isna(region_str) or region_str in {"", "None"}:
        return None, None, None

    match = re.match(r"(.+):(\d+)-(\d+)", str(region_str))
    if match is None:
        return None, None, None

    chrom = match.group(1)
    start = int(match.group(2))
    end = int(match.group(3))
    return chrom, start, end


def _ensure_faidx(fasta_file: str) -> None:
    """Create a FASTA index using ``samtools faidx`` when required."""

    index_path = f"{fasta_file}.fai"
    if os.path.exists(index_path):
        return

    print(f"  Indexing {fasta_file} with samtools faidx")
    subprocess.check_call(["samtools", "faidx", fasta_file])


def extract_sequence_samtools(
    fasta_file: str,
    chrom: str,
    start: int,
    end: int,
    output_file: str,
) -> bool:
    """Extract a sequence interval using ``samtools faidx``."""

    _ensure_faidx(fasta_file)
    region = f"{chrom}:{start}-{end}"

    try:
        with open(output_file, "w", encoding="utf-8") as handle:
            subprocess.check_call(["samtools", "faidx", fasta_file, region], stdout=handle)
    except subprocess.CalledProcessError as error:
        print(f"  ERROR: samtools faidx failed for {region} ({error})")
        return False

    if Path(output_file).stat().st_size == 0:
        print(f"  ERROR: extracted file is empty: {output_file}")
        return False

    return True


def extract_sequence_pysam(
    fasta_file: str,
    chrom: str,
    start: int,
    end: int,
    output_file: str,
) -> bool:
    """Extract a sequence interval using ``pysam`` (fallback path)."""

    try:
        import pysam
    except ImportError as exc:  # pragma: no cover - environment dependent
        print("ERROR: samtools is unavailable and pysam is not installed.")
        print("Install samtools or run `pip install pysam` to enable sequence extraction.")
        raise SystemExit(1) from exc

    try:
        fasta = pysam.FastaFile(fasta_file)
        sequence = fasta.fetch(chrom, start, end)
    except Exception as error:  # pragma: no cover - defensive path
        print(f"  ERROR: failed to extract {chrom}:{start}-{end} from {fasta_file}")
        print(f"  {error}")
        return False

    with open(output_file, "w", encoding="utf-8") as handle:
        handle.write(f">{chrom}:{start}-{end}\n")
        for i in range(0, len(sequence), 60):
            handle.write(sequence[i : i + 60] + "\n")

    return True


def find_assembly_file(sample: str, haplotype: str, assemblies_config: Dict[str, str]) -> Optional[str]:
    """Resolve the appropriate assembly FASTA for a sample/haplotype pair."""

    explicit_key = f"{sample}_{haplotype}"
    if explicit_key in assemblies_config:
        return assemblies_config[explicit_key]

    assembly_dir = assemblies_config.get("assembly_dir")
    pattern = assemblies_config.get("assembly_pattern")
    if not assembly_dir or not pattern:
        return None

    filename = pattern.replace("{sample}", sample).replace("{haplotype}", haplotype)
    candidate = Path(assembly_dir, filename)
    if candidate.exists():
        return str(candidate)

    stem = candidate.stem
    for ext in (".fa", ".fasta", ".fa.gz", ".fasta.gz"):
        option = Path(assembly_dir, f"{stem}{ext}")
        if option.exists():
            return str(option)

    return None


def _select_extractor() -> Callable[[str, str, int, int, str], bool]:
    """Choose the sequence extraction implementation available on the host."""

    if shutil.which("samtools"):
        print("Using samtools for sequence extraction")
        return extract_sequence_samtools

    print("samtools not detected; falling back to pysam")
    return extract_sequence_pysam


def _extract_pair(
    sd_id: str,
    region_label: str,
    region: str,
    fasta_file: str,
    extractor: Callable[[str, str, int, int, str], bool],
    output_dir: str,
) -> Tuple[bool, Optional[str]]:
    """Extract a single SD region and report success state."""

    chrom, start, end = parse_region(region)
    if chrom is None or start is None or end is None:
        print(f"  ERROR: invalid coordinates for {region_label} ({region})")
        return False, None

    fasta_path = os.path.join(output_dir, f"{sd_id}_{region_label}.fa")
    success = extractor(fasta_file, chrom, start, end, fasta_path)
    return success, fasta_path if success else None


def main(snakemake) -> None:  # type: ignore[annotation-unchecked]
    """Extract SD sequences required by downstream alignment and visualisation rules."""

    df = pd.read_csv(snakemake.input.filtered_sds, sep="\t")
    print(f"Extracting sequences for {len(df)} SD pairs")

    assemblies_config = snakemake.params.assemblies
    reference_config = snakemake.params.reference

    output_dir = snakemake.output.sequence_dir
    os.makedirs(output_dir, exist_ok=True)

    extractor = _select_extractor()

    manifest_records = []
    for _, row in df.iterrows():
        sd_id = row["SD_ID"]
        sample = row["Sample"]
        haplotype = row["Haplotype"]
        coord_to_use = row.get("coord_to_use", "contig")

        print(f"\nProcessing {sd_id} ({sample} {haplotype})")
        print(f"  Coordinate system: {coord_to_use}")

        if coord_to_use == "t2t":
            sd1_region = row.get("SD1_T2T")
            sd2_region = row.get("SD2_T2T")
            fasta_file = reference_config.get("fasta")
            if not fasta_file or not os.path.exists(fasta_file):
                print(f"  ERROR: reference FASTA not found for T2T extraction ({fasta_file})")
                continue
        else:
            sd1_region = row["SD1_original"]
            sd2_region = row["SD2_original"]
            fasta_file = find_assembly_file(sample, haplotype, assemblies_config)
            if not fasta_file or not os.path.exists(fasta_file):
                print(f"  ERROR: assembly not found for {sample}_{haplotype}")
                continue

        print(f"  FASTA source: {fasta_file}")
        print(f"  SD1 region: {sd1_region}")
        print(f"  SD2 region: {sd2_region}")

        success1, sd1_path = _extract_pair(sd_id, "SD1", sd1_region, fasta_file, extractor, output_dir)
        success2, sd2_path = _extract_pair(sd_id, "SD2", sd2_region, fasta_file, extractor, output_dir)

        if not (success1 and success2 and sd1_path and sd2_path):
            print(f"  ERROR: extraction failed for {sd_id}")
            continue

        manifest_records.append(
            {
                "SD_ID": sd_id,
                "Sample": sample,
                "Haplotype": haplotype,
                "SD1_region": sd1_region,
                "SD2_region": sd2_region,
                "SD1_fasta": sd1_path,
                "SD2_fasta": sd2_path,
                "coord_system": coord_to_use,
                "source_fasta": fasta_file,
            }
        )

    manifest_df = pd.DataFrame(manifest_records)
    manifest_df.to_csv(snakemake.output.manifest, sep="\t", index=False)

    print(f"\nSuccessfully extracted {len(manifest_df)} of {len(df)} SD pairs")
    print(f"Manifest written to: {snakemake.output.manifest}")

    if manifest_df.empty:
        print("No sequences were extracted. Verify assembly paths, coordinate definitions, and tool availability.")
        raise SystemExit(1)


if __name__ == "__main__":
    main(snakemake)

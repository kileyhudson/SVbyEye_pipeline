"""Snakemake workflow for SVbyEye-based segmental duplication visualisation."""

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

configfile: "config.yaml"

import os
import pandas as pd

# ---------------------------------------------------------------------------
# Global variables
# ---------------------------------------------------------------------------

OUTPUT_DIR = config.get("output_dir", "results")

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def get_all_sd_ids(wildcards):
    """Return all SD identifiers recorded in the filtered list."""
    filtered_file = f"{OUTPUT_DIR}/filtered_sds.tsv"
    if os.path.exists(filtered_file):
        df = pd.read_csv(filtered_file, sep="\t")
        return df['SD_ID'].tolist()
    return []

# ---------------------------------------------------------------------------
# Target rules
# ---------------------------------------------------------------------------

rule all:
    """Default target that generates filtered tables and the summary report."""
    input:
        # Filtered SD list
        filtered_sds = f"{OUTPUT_DIR}/filtered_sds.tsv",
        # Summary report
        report = f"{OUTPUT_DIR}/summary_report.html" if config.get("generate_report", True) else []

rule all_plots:
    """Generate every individual SVbyEye plot."""
    input:
        f"{OUTPUT_DIR}/plots_complete.txt"

# ---------------------------------------------------------------------------
# Filtering rules
# ---------------------------------------------------------------------------

rule filter_sds:
    """Filter the SD callset and compute prioritisation metrics."""
    input:
        sd_calls = config["sd_calls"]
    output:
        filtered = f"{OUTPUT_DIR}/filtered_sds.tsv",
        stats = f"{OUTPUT_DIR}/filtering_stats.txt"
    params:
        max_sds = config.get("max_sds", 20),
        filter_opts = config.get("filter", {}),
        coord_system = config.get("coordinate_system", "auto")
    log:
        f"{OUTPUT_DIR}/logs/filter_sds.log"
    script:
        "scripts/filter_sds.py"

# ---------------------------------------------------------------------------
# Sequence extraction
# ---------------------------------------------------------------------------

checkpoint extract_all_sequences:
    """Extract sequences for all filtered SD pairs."""
    input:
        filtered_sds = rules.filter_sds.output.filtered,
        assemblies = []  # Assembly files will be found dynamically
    output:
        sequence_dir = directory(f"{OUTPUT_DIR}/sequences"),
        manifest = f"{OUTPUT_DIR}/sequence_manifest.txt"
    params:
        assemblies = config.get("assemblies", {}),
        reference = config.get("reference", {}),
        coord_system = config.get("coordinate_system", "auto")
    log:
        f"{OUTPUT_DIR}/logs/extract_sequences.log"
    script:
        "scripts/extract_sequences.py"

rule get_sd_ids:
    """Write SD identifiers to a manifest supporting downstream rules."""
    input:
        rules.extract_all_sequences.output.manifest
    output:
        f"{OUTPUT_DIR}/sd_ids.txt"
    shell:
        """
        cut -f1 {input} | tail -n +2 > {output}
        """

# ---------------------------------------------------------------------------
# Alignment
# ---------------------------------------------------------------------------

def get_sd_sequences(wildcards):
    """Return the FASTA paths for a specific SD pair."""
    manifest = checkpoints.extract_all_sequences.get().output.manifest
    df = pd.read_csv(manifest, sep="\t")
    row = df[df['SD_ID'] == wildcards.sd_id].iloc[0]
    return {
        "seq1": row['SD1_fasta'],
        "seq2": row['SD2_fasta']
    }

rule align_sd_pair:
    """Align the two SD sequences using minimap2."""
    input:
        unpack(get_sd_sequences)
    output:
        paf = f"{OUTPUT_DIR}/alignments/{{sd_id}}.paf"
    params:
        preset = config.get("minimap2", {}).get("preset", "asm20"),
        extra = config.get("minimap2", {}).get("extra_params", "--eqx --secondary=no")
    threads:
        config.get("threads", {}).get("minimap2", 4)
    log:
        f"{OUTPUT_DIR}/logs/align_{{sd_id}}.log"
    shell:
        """
        minimap2 -x {params.preset} -c {params.extra} -t {threads} \
            {input.seq1} {input.seq2} > {output.paf} 2> {log}
        """

# ---------------------------------------------------------------------------
# Visualisation
# ---------------------------------------------------------------------------

rule visualize_sd:
    """Generate a plot for a single SD pair."""
    input:
        paf = rules.align_sd_pair.output.paf,
        filtered_sds = rules.filter_sds.output.filtered
    output:
        plot = f"{OUTPUT_DIR}/plots/{{sd_id}}.png"
    params:
        sd_id = "{sd_id}",
        svbyeye_opts = config.get("svbyeye", {}),
        add_genes = config.get("add_gene_annotations", True)
    log:
        f"{OUTPUT_DIR}/logs/visualize_{{sd_id}}.log"
    script:
        "scripts/visualize_sd.R"

def get_all_plots(wildcards):
    """Return all expected plot paths based on the sequence manifest."""
    manifest = checkpoints.extract_all_sequences.get().output.manifest
    if os.path.exists(manifest):
        df = pd.read_csv(manifest, sep="\t")
        sd_ids = df['SD_ID'].tolist()
        return expand(f"{OUTPUT_DIR}/plots/{{sd_id}}.png", sd_id=sd_ids)
    return []

rule collect_plots:
    """Produce a marker file once all plots have been generated."""
    input:
        get_all_plots
    output:
        f"{OUTPUT_DIR}/plots_complete.txt"
    shell:
        """
        echo "Generated {input} plots" > {output}
        ls {OUTPUT_DIR}/plots/*.png >> {output}
        """

# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

rule generate_report:
    """Assemble the final HTML report and ensure dependencies are complete."""
    input:
        filtered_sds = rules.filter_sds.output.filtered,
        plots_done = rules.collect_plots.output
    output:
        report = f"{OUTPUT_DIR}/summary_report.html"
    params:
        output_dir = OUTPUT_DIR
    log:
        f"{OUTPUT_DIR}/logs/generate_report.log"
    script:
        "scripts/generate_report.py"

# ---------------------------------------------------------------------------
# Cleanup
# ---------------------------------------------------------------------------

rule clean:
    """Remove every generated artefact."""
    shell:
        f"""
        rm -rf {OUTPUT_DIR}
        echo "Cleaned all output files"
        """

rule clean_intermediate:
    """Remove intermediate artefacts while preserving final deliverables."""
    shell:
        f"""
        rm -rf {OUTPUT_DIR}/sequences
        rm -rf {OUTPUT_DIR}/alignments
        rm -f {OUTPUT_DIR}/sequence_manifest.txt
        rm -f {OUTPUT_DIR}/sd_ids.txt
        echo "Cleaned intermediate files"
        """

# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------

rule test_config:
    """Perform basic validation of the configuration file."""
    run:
        import sys
        print("Testing configuration...")

        # Check SD calls file exists
        if not os.path.exists(config["sd_calls"]):
            print(f"ERROR: SD calls file not found: {config['sd_calls']}")
            sys.exit(1)

        print(f"✓ SD calls file found: {config['sd_calls']}")

        # Check assemblies configuration
        assemblies = config.get("assemblies", {})
        if "assembly_dir" in assemblies:
            if os.path.exists(assemblies["assembly_dir"]):
                print(f"✓ Assembly directory found: {assemblies['assembly_dir']}")
            else:
                print(f"WARNING: Assembly directory not found: {assemblies['assembly_dir']}")

        # Check reference if using T2T coordinates
        coord_system = config.get("coordinate_system", "auto")
        reference = config.get("reference", {})
        if coord_system in ["t2t", "auto"] and reference.get("fasta"):
            if os.path.exists(reference["fasta"]):
                print(f"✓ Reference genome found: {reference['fasta']}")
            else:
                print(f"WARNING: Reference genome not found: {reference['fasta']}")

        print("\nConfiguration validation completed")
        print(f"Maximum SD pairs to process: {config.get('max_sds', 20)}")
        print(f"Output directory: {OUTPUT_DIR}")

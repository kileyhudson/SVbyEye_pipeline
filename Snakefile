"""Snakemake workflow for SVbyEye-based segmental duplication visualisation."""

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

configfile: "config.yaml"

import csv
import os

# ---------------------------------------------------------------------------
# Global variables
# ---------------------------------------------------------------------------

OUTPUT_DIR = config.get("output_dir", "results")
SD_TABLE = f"{OUTPUT_DIR}/sd_table.tsv"

# ---------------------------------------------------------------------------
# Target rules
# ---------------------------------------------------------------------------

rule all:
    """Default target that generates the SD table and summary report."""
    input:
        sd_table = SD_TABLE,
        report = f"{OUTPUT_DIR}/summary_report.html" if config.get("generate_report", True) else []

rule all_plots:
    """Generate every individual SVbyEye plot."""
    input:
        f"{OUTPUT_DIR}/plots_complete.txt"

# ---------------------------------------------------------------------------
# Filtering rules
# ---------------------------------------------------------------------------

rule prepare_sd_table:
    """Attach SD identifiers and derived columns without filtering."""
    input:
        sd_calls = config["sd_calls"]
    output:
        table = SD_TABLE,
        summary = f"{OUTPUT_DIR}/sd_summary.txt"
    log:
        f"{OUTPUT_DIR}/logs/prepare_sd_table.log"
    script:
        "scripts/prepare_sd_table.py"

# ---------------------------------------------------------------------------
# Alignment manifest
# ---------------------------------------------------------------------------

checkpoint prepare_alignment_manifest:
    """Validate PAF availability and record file paths per SD."""
    input:
        sd_table = SD_TABLE
    output:
        manifest = f"{OUTPUT_DIR}/alignment_manifest.tsv"
    params:
        alignments = config.get("alignments", {})
    log:
        f"{OUTPUT_DIR}/logs/prepare_alignment_manifest.log"
    script:
        "scripts/prepare_alignment_manifest.py"

def get_paf_for_sd(wildcards):
    """Look up the precomputed PAF path for a given SD identifier."""
    manifest = checkpoints.prepare_alignment_manifest.get().output.manifest
    with open(manifest, encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["SD_ID"] == wildcards.sd_id:
                return row["paf_path"]
    raise ValueError(f"PAF path not found for {wildcards.sd_id}")

# ---------------------------------------------------------------------------
# Visualisation
# ---------------------------------------------------------------------------

rule visualize_sd:
    """Generate a plot for a single SD pair."""
    input:
        paf = get_paf_for_sd,
        sd_table = SD_TABLE
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
    """Return all expected plot paths based on the alignment manifest."""
    manifest = checkpoints.prepare_alignment_manifest.get().output.manifest
    if os.path.exists(manifest):
        with open(manifest, encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            sd_ids = [row["SD_ID"] for row in reader]
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
        sd_table = SD_TABLE,
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
        rm -f {OUTPUT_DIR}/alignment_manifest.tsv
        rm -f {OUTPUT_DIR}/plots_complete.txt
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

        # Check alignments configuration
        alignments = config.get("alignments", {})
        align_dir = alignments.get("dir", "data/alignments")
        if os.path.exists(align_dir):
            print(f"✓ Alignment directory found: {align_dir}")
        else:
            print(f"WARNING: Alignment directory not found: {align_dir}")
        pattern = alignments.get("filename_pattern", "{sd_id}.paf")
        print(f"Alignment filename pattern: {pattern}")

        print("\nConfiguration validation completed")
        print(f"Output directory: {OUTPUT_DIR}")

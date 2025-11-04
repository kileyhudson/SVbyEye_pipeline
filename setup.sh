#!/bin/bash
#
# Provisioning script for the SVbyEye pipeline.
# Installs the dependencies required to execute the workflow locally.
#

set -e

echo "========================================="
echo "SVbyEye Pipeline Setup"
echo "========================================="
echo ""

# Ensure conda (or mamba) is available before continuing.
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda not found"
    echo "Install conda or mamba before running this script:"
    echo "  https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

echo "Creating conda environment for SVbyEye pipeline..."
echo ""

# Create a dedicated environment with the required dependencies.
conda create -n svbyeye_pipeline -y \
    python=3.9 \
    snakemake \
    minimap2 \
    samtools \
    pandas \
    r-base \
    bioconductor-genomicranges \
    bioconductor-iranges \
    r-ggplot2

echo ""
# Install the SVbyEye R package and its dependencies.
echo ""

# Activate the environment and install SVbyEye.
conda activate svbyeye_pipeline || source activate svbyeye_pipeline

# Install SVbyEye from Bioconductor.
Rscript -e 'if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("SVbyEye")'

echo ""
echo "========================================="
echo "Setup complete"
echo "========================================="
echo ""
echo "To use the pipeline:"
echo "  1. Activate the environment:"
echo "     conda activate svbyeye_pipeline"
echo ""
echo "  2. Execute the workflow as needed, for example:"
echo "     snakemake --cores 4"
echo ""
echo "  3. Deactivate the environment when finished:"
echo "     conda deactivate"
echo ""

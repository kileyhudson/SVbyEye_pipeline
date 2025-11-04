#!/bin/bash
#
# Setup script for SVbyEye pipeline
# This installs all required dependencies
#

set -e

echo "========================================="
echo "SVbyEye Pipeline Setup"
echo "========================================="
echo ""

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda not found!"
    echo "Please install conda/mamba first:"
    echo "  https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

echo "Creating conda environment for SVbyEye pipeline..."
echo ""

# Create conda environment
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
echo "Installing SVbyEye R package..."
echo ""

# Activate environment and install SVbyEye
conda activate svbyeye_pipeline || source activate svbyeye_pipeline

# Install SVbyEye from Bioconductor
Rscript -e 'if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("SVbyEye")'

echo ""
echo "========================================="
echo "Setup complete!"
echo "========================================="
echo ""
echo "To use the pipeline:"
echo "  1. Activate the environment:"
echo "     conda activate svbyeye_pipeline"
echo ""
echo "  2. Run the pipeline:"
echo "     snakemake --cores 4"
echo ""
echo "  3. To deactivate when done:"
echo "     conda deactivate"
echo ""

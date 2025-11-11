#!/bin/bash
#
# Provisioning script for the SVbyEye pipeline.
# Installs the dependencies required to execute the workflow locally.
#

set -euo pipefail

echo "========================================="
echo "SVbyEye Pipeline Setup"
echo "========================================="
echo ""

# Ensure conda (or mamba) is available before continuing.
if ! command -v conda >/dev/null 2>&1 && ! command -v mamba >/dev/null 2>&1; then
    echo "ERROR: neither conda nor mamba was found on PATH."
    echo "Install Miniconda/Mambaforge before running this script."
    exit 1
fi

SOLVER="conda"
if command -v mamba >/dev/null 2>&1; then
    SOLVER="mamba"
fi

echo "Creating conda environment for SVbyEye pipeline using ${SOLVER}..."
echo ""

${SOLVER} env remove -n svbyeye_pipeline -y >/dev/null 2>&1 || true

if [ ! -f "environment.yml" ]; then
    echo "ERROR: environment.yml not found in the current directory."
    exit 1
fi

${SOLVER} env create -f environment.yml

echo ""
# Install the SVbyEye R package and its dependencies.
echo ""

CONDA_BASE="$(conda info --base 2>/dev/null || true)"
if [ -n "${CONDA_BASE}" ] && [ -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]; then
    # shellcheck disable=SC1090
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate svbyeye_pipeline
    Rscript install_svbyeye_github.R
    conda deactivate
else
    echo "WARNING: Unable to locate conda.sh to activate the environment automatically."
    echo "Please run the following manually to finish installing SVbyEye:"
    echo "  conda activate svbyeye_pipeline"
    echo "  Rscript install_svbyeye_github.R"
fi

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
echo "     # Ensure data/alignments contains PAF files generated upstream with minimap2 -c --eqx"
echo ""
echo "  3. Deactivate the environment when finished:"
echo "     conda deactivate"
echo ""

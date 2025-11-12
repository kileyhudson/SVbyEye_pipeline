#!/bin/bash
#
# Provisioning script for the SVbyEye pipeline.
# Installs dependencies and sets up the workflow environment.
# Works on local machines and HPC clusters.
#

set -euo pipefail

echo "========================================="
echo "SVbyEye Pipeline Setup"
echo "========================================="
echo ""

# Ensure conda (or mamba) is available before continuing.
if ! command -v conda >/dev/null 2>&1 && ! command -v mamba >/dev/null 2>&1; then
    echo "ERROR: neither conda nor mamba was found on PATH."
    echo ""
    echo "Local installation:"
    echo "  Install Miniconda/Mambaforge before running this script."
    echo ""
    echo "HPC cluster:"
    echo "  Try: module load conda"
    echo "  Or:  module load miniforge"
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
echo "Installing SVbyEye R package and dependencies..."
echo ""

# Try to activate the environment and install SVbyEye
CONDA_BASE="$(conda info --base 2>/dev/null || ${SOLVER} info --base 2>/dev/null || true)"
if [ -n "${CONDA_BASE}" ] && [ -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]; then
    # shellcheck disable=SC1090
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate svbyeye_pipeline
    
    # Install SVbyEye R package
    Rscript install_svbyeye_github.R
    
    # Check if snakemake installation is functional by testing the actual import
    echo ""
    echo "Verifying snakemake installation..."
    
    # Test the Python import that actually fails on broken installations
    if ! python -c "from snakemake.cli import main" >/dev/null 2>&1; then
        echo "⚠️  WARNING: Conda snakemake is broken (Python import failed)"
        echo "This can happen on some HPC systems."
        echo ""
        
        # Find the environment path and remove broken snakemake
        ENV_PATH=$(conda env list | grep "^svbyeye_pipeline " | awk '{print $NF}')
        if [ -f "${ENV_PATH}/bin/snakemake" ]; then
            echo "Removing broken snakemake binary..."
            rm -f "${ENV_PATH}/bin/snakemake"
            echo "✓ Removed"
        fi
        
        echo ""
        echo "Use system snakemake instead:"
        echo "  - HPC cluster: module load snakemake"
        echo "  - Local: mamba install snakemake -c conda-forge -c bioconda"
        SNAKEMAKE_BROKEN=true
    else
        echo "✓ Snakemake is working"
        SNAKEMAKE_BROKEN=false
    fi
    
    conda deactivate
else
    echo "WARNING: Unable to locate conda.sh to activate the environment automatically."
    echo "Please run the following manually to finish setup:"
    echo "  conda activate svbyeye_pipeline"
    echo "  Rscript install_svbyeye_github.R"
    echo ""
    SNAKEMAKE_BROKEN=unknown
fi

# Create activation helper script
echo ""
echo "Creating activation helper script..."

cat > activate_pipeline.sh << 'ACTIVATE_SCRIPT'
#!/bin/bash
# Activation helper for SVbyEye pipeline
# Usage: source activate_pipeline.sh

# Find the conda environment path
ENV_PATH=$(conda env list 2>/dev/null | grep "^svbyeye_pipeline " | awk '{print $NF}')
if [ -z "$ENV_PATH" ]; then
    ENV_PATH=$(mamba env list 2>/dev/null | grep "^svbyeye_pipeline " | awk '{print $NF}')
fi

if [ -z "$ENV_PATH" ]; then
    echo "ERROR: svbyeye_pipeline environment not found"
    echo "Run: ./setup.sh"
    return 1 2>/dev/null || exit 1
fi

# Check if svbyeye_pipeline has a working snakemake
SNAKEMAKE_WORKS=false
if [ -f "${ENV_PATH}/bin/snakemake" ]; then
    # Test if it actually works (not just if the binary exists)
    if "${ENV_PATH}/bin/python" -c "from snakemake.cli import main" >/dev/null 2>&1; then
        SNAKEMAKE_WORKS=true
    fi
fi

# If snakemake is broken, try to find snakemake_plus (common on HPC clusters)
if [ "$SNAKEMAKE_WORKS" = "false" ]; then
    SNAKEMAKE_PLUS=$(conda env list 2>/dev/null | grep "^snakemake_plus " | awk '{print $NF}')
    if [ -z "$SNAKEMAKE_PLUS" ]; then
        SNAKEMAKE_PLUS=$(mamba env list 2>/dev/null | grep "^snakemake_plus " | awk '{print $NF}')
    fi
    
    # Add snakemake_plus to PATH first if it exists
    if [ -n "$SNAKEMAKE_PLUS" ] && [ -d "${SNAKEMAKE_PLUS}/bin" ]; then
        export PATH="${SNAKEMAKE_PLUS}/bin:$PATH"
    fi
fi

# Add svbyeye_pipeline environment to PATH (for R and Python)
export PATH="${ENV_PATH}/bin:$PATH"

# Clear bash command cache
hash -r 2>/dev/null || true

# Verify setup
echo "✓ SVbyEye pipeline environment activated"
echo ""
echo "Active tools:"
echo "  Python: $(which python)"
echo "  R: $(which R)"

if command -v snakemake >/dev/null 2>&1; then
    SNAKEMAKE_VERSION=$(snakemake --version 2>/dev/null || echo 'unknown')
    echo "  Snakemake: $(which snakemake) (${SNAKEMAKE_VERSION})"
else
    echo "  Snakemake: NOT FOUND"
    echo ""
    echo "Load snakemake:"
    echo "  - HPC: module load snakemake"
    echo "  - Local: mamba install snakemake -c conda-forge -c bioconda"
fi

echo ""
echo "Ready to run: snakemake --cores N"
ACTIVATE_SCRIPT

chmod +x activate_pipeline.sh

echo ""
echo "========================================="
echo "Setup Complete!"
echo "========================================="
echo ""

if [ "${SNAKEMAKE_BROKEN:-unknown}" = "true" ]; then
    echo "⚠️  Note: Conda snakemake was removed (broken installation)"
    echo ""
    echo "To use the pipeline:"
    echo "  1. Load system snakemake:"
    echo "     module load snakemake  # on HPC"
    echo ""
    echo "  2. Activate this environment:"
    echo "     source activate_pipeline.sh"
    echo ""
    echo "  3. Run the pipeline:"
    echo "     snakemake --cores 4"
    echo ""
else
    echo "To use the pipeline:"
    echo "  1. Activate the environment:"
    echo "     source activate_pipeline.sh"
    echo "     # OR: conda activate svbyeye_pipeline"
    echo ""
    echo "  2. Run the pipeline:"
    echo "     snakemake --cores 4"
    echo ""
    echo "  3. Deactivate when finished:"
    echo "     conda deactivate"
    echo ""
fi

echo "Note: Ensure data/alignments/ contains PAF files"
echo "      generated with: minimap2 -c --eqx"
echo ""

#!/usr/bin/env Rscript
# Intelligent SVbyEye installation script
# Automatically handles library setup, dependency detection, and installation

cat("SVbyEye Installation Beginning...\n")
cat("========================================\n\n")

# ---------------------------------------------------------------------------
# Step 1: Setup library path (personal library if system is not writable)
# ---------------------------------------------------------------------------

cat("Step 1: Checking R library configuration...\n")

default_lib <- .libPaths()[1]
can_write <- file.access(default_lib, 2) == 0

if (!can_write) {
    cat("  System library not writable: ", default_lib, "\n", sep = "")
    
    # Determine personal library path based on R version
    r_version <- paste(R.version$major, 
                      strsplit(R.version$minor, "\\.")[[1]][1], 
                      sep = ".")
    personal_lib <- path.expand(paste0("~/R/x86_64-pc-linux-gnu-library/", r_version))
    
    # Create personal library if it doesn't exist
    if (!dir.exists(personal_lib)) {
        cat("  Creating personal library: ", personal_lib, "\n", sep = "")
        dir.create(personal_lib, recursive = TRUE, showWarnings = FALSE)
    } else {
        cat("  Personal library exists: ", personal_lib, "\n", sep = "")
    }
    
    # Set personal library as primary
    .libPaths(c(personal_lib, .libPaths()))
    cat("  ✓ Using personal library\n")
} else {
    cat("  ✓ System library is writable\n")
}

cat("  Active library: ", .libPaths()[1], "\n\n", sep = "")

# ---------------------------------------------------------------------------
# Step 2: Check if SVbyEye is already installed
# ---------------------------------------------------------------------------

cat("Step 2: Checking for existing SVbyEye installation...\n")

if (requireNamespace("SVbyEye", quietly = TRUE)) {
    cat("  ✓ SVbyEye is already installed\n")
    library(SVbyEye)
    cat("  Version: ", as.character(packageVersion("SVbyEye")), "\n", sep = "")
    cat("\n========================================\n")
    cat("Installation check complete!\n")
    cat("========================================\n")
    quit(save = "no", status = 0)
}

cat("  SVbyEye not found, will install\n\n")

# ---------------------------------------------------------------------------
# Step 3: Install remotes if needed
# ---------------------------------------------------------------------------

cat("Step 3: Checking remotes package...\n")

if (!requireNamespace("remotes", quietly = TRUE)) {
    cat("  Installing remotes...\n")
    install.packages("remotes", repos = "https://cloud.r-project.org", quiet = TRUE)
    cat("  ✓ remotes installed\n")
} else {
    cat("  ✓ remotes already available\n")
}
cat("\n")

# ---------------------------------------------------------------------------
# Step 4: Check and install CRAN dependencies
# ---------------------------------------------------------------------------

cat("Step 4: Checking CRAN dependencies...\n")

cran_deps <- c("ggplot2", "ggnewscale", "ggforce", "gggenes", "scales",
               "dplyr", "stringr", "magrittr", "data.table", "tibble",
               "wesanderson", "randomcoloR")

missing_cran <- character(0)
for (pkg in cran_deps) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        missing_cran <- c(missing_cran, pkg)
    }
}

if (length(missing_cran) > 0) {
    cat("  Installing", length(missing_cran), "CRAN packages:", 
        paste(missing_cran, collapse = ", "), "\n")
    cat("  This may take 2-3 minutes...\n")
    install.packages(missing_cran, repos = "https://cloud.r-project.org", quiet = TRUE)
    cat("  ✓ CRAN dependencies installed\n")
} else {
    cat("  ✓ All CRAN dependencies already installed\n")
}
cat("\n")

# ---------------------------------------------------------------------------
# Step 5: Check and install Bioconductor dependencies
# ---------------------------------------------------------------------------

cat("Step 5: Checking Bioconductor dependencies...\n")

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    cat("  Installing BiocManager...\n")
    install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
    cat("  ✓ BiocManager installed\n")
} else {
    cat("  ✓ BiocManager already available\n")
}

bioc_deps <- c("S4Vectors", "BSgenome", "GenomicRanges", "GenomicAlignments",
               "Biostrings", "IRanges", "GenomeInfoDb", "Rsamtools")

missing_bioc <- character(0)
for (pkg in bioc_deps) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        missing_bioc <- c(missing_bioc, pkg)
    }
}

if (length(missing_bioc) > 0) {
    cat("  Installing", length(missing_bioc), "Bioconductor packages:", 
        paste(missing_bioc, collapse = ", "), "\n")
    cat("  This may take 5-10 minutes as packages compile...\n")
    BiocManager::install(missing_bioc, update = FALSE, ask = FALSE, quiet = TRUE)
    cat("  ✓ Bioconductor dependencies installed\n")
} else {
    cat("  ✓ All Bioconductor dependencies already installed\n")
}
cat("\n")

# ---------------------------------------------------------------------------
# Step 6: Install SVbyEye from GitHub
# ---------------------------------------------------------------------------

cat("Step 6: Installing SVbyEye from GitHub...\n")
cat("  Source: https://github.com/daewoooo/SVbyEye\n")
cat("  This may take 1-2 minutes...\n")

remotes::install_github(
    "daewoooo/SVbyEye",
    dependencies = FALSE,  # We handled dependencies already
    upgrade = "never",
    build_vignettes = FALSE,
    quiet = TRUE
)

cat("  ✓ SVbyEye installed\n\n")

# ---------------------------------------------------------------------------
# Step 7: Verify installation
# ---------------------------------------------------------------------------

cat("Step 7: Verifying installation...\n")

if (!requireNamespace("SVbyEye", quietly = TRUE)) {
    cat("  ✗ ERROR: SVbyEye installation failed\n")
    cat("  Check error messages above for details\n")
    quit(save = "no", status = 1)
}

library(SVbyEye)
version <- as.character(packageVersion("SVbyEye"))
install_path <- find.package("SVbyEye")

cat("  ✓ SVbyEye loaded successfully\n")
cat("  Version: ", version, "\n", sep = "")
cat("  Location: ", install_path, "\n", sep = "")

cat("Installation complete!\n")

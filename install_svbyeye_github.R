#!/usr/bin/env Rscript
# Install SVbyEye package from GitHub

cat("Ensuring remotes is available...\n")
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("SVbyEye", quietly = TRUE)) {
    cat("Installing SVbyEye from GitHub (dependencies handled via Conda)...\n")
    remotes::install_github(
        "daewoooo/SVbyEye",
        dependencies = FALSE,
        upgrade = "never",
        build_vignettes = FALSE
    )
} else {
    cat("SVbyEye already installed, skipping download.\n")
}

cat("Checking installation...\n")
library(SVbyEye)
cat("SVbyEye version:", as.character(packageVersion("SVbyEye")), "\n")
cat("Installation successful!\n")

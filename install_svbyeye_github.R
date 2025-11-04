#!/usr/bin/env Rscript
# Install SVbyEye package from GitHub

cat("Installing remotes package...\n")
if (!require("remotes", quietly = TRUE)) {
    install.packages("remotes", repos="https://cloud.r-project.org")
}

cat("Installing dependencies...\n")
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos="https://cloud.r-project.org")
}

# Install Bioconductor dependencies
BiocManager::install(c("GenomicRanges", "IRanges", "Biostrings"), update=FALSE, ask=FALSE)

cat("Installing ggplot2...\n")
install.packages("ggplot2", repos="https://cloud.r-project.org", quiet=TRUE)

cat("Installing SVbyEye from GitHub...\n")
remotes::install_github("daewoooo/SVbyEye")

cat("Checking installation...\n")
library(SVbyEye)
cat("SVbyEye version:", as.character(packageVersion("SVbyEye")), "\n")
cat("Installation successful!\n")

#!/usr/bin/env Rscript
# Install SVbyEye package

cat("Installing BiocManager...\n")
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos="https://cloud.r-project.org")
}

cat("Installing SVbyEye...\n")
BiocManager::install("SVbyEye", update=FALSE, ask=FALSE)

cat("Checking installation...\n")
library(SVbyEye)
cat("SVbyEye version:", as.character(packageVersion("SVbyEye")), "\n")

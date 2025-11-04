#!/usr/bin/env Rscript
#
# Visualize SD pair using SVbyEye
#

# Suppress warnings
options(warn = -1)

# Load required libraries
suppressPackageStartupMessages({
    library(SVbyEye)
    library(ggplot2)
})

# Get snakemake parameters
paf_file <- snakemake@input[["paf"]]
filtered_sds_file <- snakemake@input[["filtered_sds"]]
output_plot <- snakemake@output[["plot"]]
sd_id <- snakemake@params[["sd_id"]]
svbyeye_opts <- snakemake@params[["svbyeye_opts"]]
add_genes <- snakemake@params[["add_genes"]]

# Log
cat("Visualizing", sd_id, "\n")
cat("PAF file:", paf_file, "\n")
cat("Output plot:", output_plot, "\n")

# Read PAF file
cat("Reading PAF alignments...\n")
paf_table <- tryCatch({
    readPaf(paf.file = paf_file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
}, error = function(e) {
    cat("ERROR reading PAF file:", conditionMessage(e), "\n")
    # Create empty PAF for plotting (indicates no alignment)
    data.frame()
})

# Check if PAF is empty
if (nrow(paf_table) == 0) {
    cat("WARNING: No alignments found in PAF file!\n")
    cat("Creating placeholder plot...\n")

    # Create a simple plot indicating no alignment
    p <- ggplot() +
        annotate("text",
            x = 0.5, y = 0.5,
            label = paste0(sd_id, "\nNo alignments found"),
            size = 8, hjust = 0.5, vjust = 0.5
        ) +
        theme_void() +
        labs(title = paste("SD Pair:", sd_id))

    # Save plot
    ggsave(output_plot, plot = p, width = 12, height = 8, dpi = 300)
    cat("Saved placeholder plot to:", output_plot, "\n")
    quit(save = "no", status = 0)
}

cat("Found", nrow(paf_table), "alignment(s)\n")

# Get visualization options
color_by <- svbyeye_opts$color_by %||% "direction"
binsize <- svbyeye_opts$binsize %||% 0
highlight_indels <- svbyeye_opts$highlight_indels %||% FALSE
min_del_size <- svbyeye_opts$min_deletion_size %||% 50
min_ins_size <- svbyeye_opts$min_insertion_size %||% 50
plot_width <- svbyeye_opts$width %||% 12
plot_height <- svbyeye_opts$height %||% 8
plot_dpi <- svbyeye_opts$dpi %||% 300

# Create base plot
cat("Creating SVbyEye plot...\n")

plot_params <- list(
    paf.table = paf_table,
    color.by = color_by
)

# Add binning if requested
if (binsize > 0) {
    plot_params$binsize <- binsize
}

# Add indel highlighting if requested
if (highlight_indels) {
    plot_params$min.deletion.size <- min_del_size
    plot_params$min.insertion.size <- min_ins_size
    plot_params$highlight.sv <- "outline"
}

# Create plot
plt <- tryCatch({
    do.call(plotMiro, plot_params)
}, error = function(e) {
    cat("ERROR creating plot:", conditionMessage(e), "\n")
    cat("Trying simplified plot...\n")
    plotMiro(paf.table = paf_table)
})

# Add title with SD info
filtered_sds <- read.table(filtered_sds_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sd_info <- filtered_sds[filtered_sds$SD_ID == sd_id, ]

if (nrow(sd_info) > 0) {
    title_text <- paste0(
        sd_id, ": ", sd_info$Sample, " ", sd_info$Haplotype,
        "\nSD1: ", sd_info$SD1_original,
        "\nSD2: ", sd_info$SD2_original
    )

    # Add gene info if available
    if ("Gene_names" %in% colnames(sd_info) && !is.na(sd_info$Gene_names) && sd_info$Gene_names != "None") {
        title_text <- paste0(title_text, "\nGenes: ", sd_info$Gene_names)
    }

    plt <- plt + labs(title = title_text)
}

# Add gene annotations if requested
if (add_genes && nrow(sd_info) > 0) {
    # Check if either SD has genes
    has_genes <- FALSE
    if ("SD1_has_exon" %in% colnames(sd_info) && !is.na(sd_info$SD1_has_exon)) {
        has_genes <- has_genes || sd_info$SD1_has_exon
    }
    if ("SD2_has_exon" %in% colnames(sd_info) && !is.na(sd_info$SD2_has_exon)) {
        has_genes <- has_genes || sd_info$SD2_has_exon
    }

    if (has_genes) {
        cat("Note: Gene annotations requested but not yet implemented\n")
        cat("      To add gene annotations, provide gene coordinates as GenomicRanges\n")
    }
}

# Apply theme adjustments
plt <- plt + theme(
    plot.title = element_text(size = 10, hjust = 0),
    legend.position = "bottom"
)

# Save plot
cat("Saving plot to:", output_plot, "\n")
ggsave(
    output_plot,
    plot = plt,
    width = plot_width,
    height = plot_height,
    dpi = plot_dpi
)

cat("Done!\n")

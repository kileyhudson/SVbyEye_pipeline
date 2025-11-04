#!/usr/bin/env Rscript
# Visualise a segmental duplication pair using SVbyEye.

# Suppress non-critical warnings to keep logs concise.
options(warn = -1)

# Load required libraries.
suppressPackageStartupMessages({
    library(SVbyEye)
    library(ggplot2)
})

# Retrieve snakemake parameters supplied by the workflow.
paf_file <- snakemake@input[["paf"]]
filtered_sds_file <- snakemake@input[["filtered_sds"]]
output_plot <- snakemake@output[["plot"]]
sd_id <- snakemake@params[["sd_id"]]
svbyeye_opts <- snakemake@params[["svbyeye_opts"]]
add_genes <- snakemake@params[["add_genes"]]

# Log the primary inputs for traceability.
cat("Rendering visualisation for", sd_id, "\n")
cat("PAF file:", paf_file, "\n")
cat("Output plot:", output_plot, "\n")

# Read the PAF alignments produced by minimap2.
cat("Reading PAF alignments...\n")
paf_table <- tryCatch({
    readPaf(paf.file = paf_file, include.paf.tags = TRUE, restrict.paf.tags = "cg")
}, error = function(e) {
    cat("ERROR reading PAF file:", conditionMessage(e), "\n")
    # Create empty structure to trigger placeholder plot generation.
    data.frame()
})

# Handle the case where no alignments were detected.
if (nrow(paf_table) == 0) {
    cat("WARNING: No alignments found in PAF file!\n")
    cat("Creating placeholder plot...\n")

    # Create a placeholder plot highlighting the absence of alignments.
    p <- ggplot() +
        annotate("text",
            x = 0.5, y = 0.5,
            label = paste0(sd_id, "\nNo alignments found"),
            size = 8, hjust = 0.5, vjust = 0.5
        ) +
        theme_void() +
        labs(title = paste("SD Pair:", sd_id))

    # Persist the placeholder plot before exiting.
    ggsave(output_plot, plot = p, width = 12, height = 8, dpi = 300)
    cat("Saved placeholder plot to:", output_plot, "\n")
    quit(save = "no", status = 0)
}

cat("Found", nrow(paf_table), "alignment(s)\n")

# Extract visualisation parameters from the configuration block.
color_by <- svbyeye_opts$color_by %||% "direction"
binsize <- svbyeye_opts$binsize %||% 0
highlight_indels <- svbyeye_opts$highlight_indels %||% FALSE
min_del_size <- svbyeye_opts$min_deletion_size %||% 50
min_ins_size <- svbyeye_opts$min_insertion_size %||% 50
plot_width <- svbyeye_opts$width %||% 12
plot_height <- svbyeye_opts$height %||% 8
plot_dpi <- svbyeye_opts$dpi %||% 300

# Initialise the parameter list passed to SVbyEye.
cat("Creating SVbyEye plot...\n")

plot_params <- list(
    paf.table = paf_table,
    color.by = color_by
)

# Add binning when explicitly configured.
if (binsize > 0) {
    plot_params$binsize <- binsize
}

# Enable structural variant highlighting when requested.
if (highlight_indels) {
    plot_params$min.deletion.size <- min_del_size
    plot_params$min.insertion.size <- min_ins_size
    plot_params$highlight.sv <- "outline"
}

# Generate the SVbyEye plot and fall back to defaults if necessary.
plt <- tryCatch({
    do.call(plotMiro, plot_params)
}, error = function(e) {
    cat("ERROR creating plot:", conditionMessage(e), "\n")
    cat("Trying simplified plot...\n")
    plotMiro(paf.table = paf_table)
})

# Annotate the plot with metadata describing the SD pair.
filtered_sds <- read.table(filtered_sds_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
sd_info <- filtered_sds[filtered_sds$SD_ID == sd_id, ]

if (nrow(sd_info) > 0) {
    title_text <- paste0(
        sd_id, ": ", sd_info$Sample, " ", sd_info$Haplotype,
        "\nSD1: ", sd_info$SD1_original,
        "\nSD2: ", sd_info$SD2_original
    )

    # Append gene information when present in the metadata.
    if ("Gene_names" %in% colnames(sd_info) && !is.na(sd_info$Gene_names) && sd_info$Gene_names != "None") {
        title_text <- paste0(title_text, "\nGenes: ", sd_info$Gene_names)
    }

    plt <- plt + labs(title = title_text)
}

# Report planned gene-annotation support when requested.
if (add_genes && nrow(sd_info) > 0) {
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

# Apply minimal theming for consistency across outputs.
plt <- plt + theme(
    plot.title = element_text(size = 10, hjust = 0),
    legend.position = "bottom"
)

# Write the plot to disk.
cat("Saving plot to:", output_plot, "\n")
ggsave(
    output_plot,
    plot = plt,
    width = plot_width,
    height = plot_height,
    dpi = plot_dpi
)

cat("Visualisation complete.\n")

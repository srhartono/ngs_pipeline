#!/usr/bin/env Rscript

# DESeq2 Differential Expression Analysis for NGS Pipeline
# This script performs comprehensive differential expression analysis

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(pheatmap)
  library(EnhancedVolcano)
  library(VennDiagram)
  library(RColorBrewer)
  library(reshape2)
  library(scales)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  cat("Usage: Rscript deseq2.R <count_matrix> <sample_metadata> <output_dir> [design_formula]\n")
  cat("  count_matrix: Tab-separated count matrix with genes as rows, samples as columns\n")
  cat("  sample_metadata: Tab-separated metadata with sample_id, condition, batch columns\n")
  cat("  output_dir: Output directory for results\n")
  cat("  design_formula: DESeq2 design formula (default: '~ batch + condition')\n")
  quit(status = 1)
}

count_file <- args[1]
metadata_file <- args[2]
output_dir <- args[3]
design_formula <- ifelse(length(args) >= 4, args[4], "~ batch + condition")

cat("NGS Pipeline - DESeq2 Analysis\n")
cat("=========================================\n")
cat("Count matrix:", count_file, "\n")
cat("Metadata:", metadata_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("Design formula:", design_formula, "\n\n")

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Function to safely load data
safe_read <- function(file, ...) {
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }
  read_tsv(file, ...)
}

# Load count matrix
cat("Loading count matrix...\n")
count_data <- safe_read(count_file)

# Prepare count matrix
gene_ids <- count_data[[1]]
count_matrix <- as.matrix(count_data[, -1])
rownames(count_matrix) <- gene_ids

cat("Count matrix dimensions:", nrow(count_matrix), "genes x", ncol(count_matrix), "samples\n")

# Load sample metadata
cat("Loading sample metadata...\n")
metadata <- safe_read(metadata_file)

# Ensure sample order matches
sample_order <- colnames(count_matrix)
metadata <- metadata[match(sample_order, metadata$sample_id), ]

if (any(is.na(metadata$sample_id))) {
  missing_samples <- sample_order[is.na(metadata$sample_id)]
  stop("Missing metadata for samples: ", paste(missing_samples, collapse = ", "))
}

cat("Samples matched successfully\n")

# Prepare metadata for DESeq2
metadata$condition <- factor(metadata$condition, 
                           levels = c("Untreated", "siXRN2", "siXRN2_XRN2WT_OE", "siXRN2_XRN2MT_OE"))

if ("batch" %in% colnames(metadata)) {
  metadata$batch <- factor(metadata$batch)
} else {
  cat("No batch column found, creating dummy batch\n")
  metadata$batch <- factor(rep(1, nrow(metadata)))
}

cat("Conditions found:", paste(levels(metadata$condition), collapse = ", "), "\n")
cat("Batches found:", paste(levels(metadata$batch), collapse = ", "), "\n")

# Filter low-count genes
cat("Filtering low-count genes...\n")
keep <- rowSums(count_matrix >= 10) >= 3  # At least 10 counts in at least 3 samples
count_matrix <- count_matrix[keep, ]
cat("Retained", nrow(count_matrix), "genes after filtering\n")

# Create DESeq2 dataset
cat("Creating DESeq2 dataset...\n")
if (design_formula == "~ condition") {
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = metadata,
    design = ~ condition
  )
} else {
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = metadata,
    design = ~ batch + condition
  )
}

# Pre-filtering
dds <- dds[rowSums(counts(dds)) > 1, ]

# Run DESeq2 analysis
cat("Running DESeq2 analysis...\n")
dds <- DESeq(dds)

# Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
vst_counts <- vst(dds, blind = FALSE)

# Save normalized counts
write_tsv(
  as.data.frame(normalized_counts) %>% 
    rownames_to_column("gene_id"),
  file.path(output_dir, "normalized_counts.tsv")
)

write_tsv(
  as.data.frame(assay(vst_counts)) %>% 
    rownames_to_column("gene_id"),
  file.path(output_dir, "vst_counts.tsv")
)

cat("Normalized counts saved\n")

# Principal Component Analysis
cat("Generating PCA plot...\n")
pca_data <- plotPCA(vst_counts, intgroup = c("condition", "batch"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = condition, shape = batch)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("Principal Component Analysis") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "pca_plot.png"), pca_plot, width = 8, height = 6, dpi = 300)

# Sample correlation heatmap
cat("Generating sample correlation heatmap...\n")
sample_cor <- cor(assay(vst_counts))
pheatmap(
  sample_cor,
  filename = file.path(output_dir, "sample_correlation_heatmap.png"),
  width = 8, height = 6,
  annotation_col = metadata[, c("condition", "batch"), drop = FALSE],
  main = "Sample Correlation Heatmap"
)

# Define contrasts of interest
contrasts <- list(
  "siXRN2_vs_Untreated" = c("condition", "siXRN2", "Untreated"),
  "siXRN2_XRN2WT_OE_vs_siXRN2" = c("condition", "siXRN2_XRN2WT_OE", "siXRN2"),
  "siXRN2_XRN2MT_OE_vs_siXRN2" = c("condition", "siXRN2_XRN2MT_OE", "siXRN2")
)

# Perform differential expression analysis for each contrast
de_results <- list()
significant_genes <- list()

for (contrast_name in names(contrasts)) {
  cat("Analyzing contrast:", contrast_name, "\n")
  
  contrast <- contrasts[[contrast_name]]
  res <- results(dds, contrast = contrast, alpha = 0.05)
  
  # Order by adjusted p-value
  res <- res[order(res$padj), ]
  
  # Add gene symbols and descriptions (if available)
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    filter(!is.na(padj))
  
  # Save results
  output_file <- file.path(output_dir, paste0(contrast_name, "_results.csv"))
  write_csv(res_df, output_file)
  
  # Store for downstream analysis
  de_results[[contrast_name]] <- res_df
  
  # Get significant genes
  sig_genes <- res_df %>%
    filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
    pull(gene_id)
  
  significant_genes[[contrast_name]] <- sig_genes
  
  cat("  Significant genes (FDR<0.05, |log2FC|>1):", length(sig_genes), "\n")
  
  # MA plot
  ma_plot <- ggplot(res_df, aes(x = baseMean, y = log2FoldChange)) +
    geom_point(aes(color = padj < 0.05), alpha = 0.6) +
    scale_x_log10() +
    scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red")) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "blue") +
    ggtitle(paste("MA Plot -", contrast_name)) +
    theme_bw() +
    labs(color = "FDR < 0.05")
  
  ggsave(
    file.path(output_dir, paste0(contrast_name, "_ma_plot.png")),
    ma_plot, width = 8, height = 6, dpi = 300
  )
  
  # Volcano plot
  volcano_plot <- EnhancedVolcano(
    res_df,
    lab = res_df$gene_id,
    x = 'log2FoldChange',
    y = 'padj',
    title = paste("Volcano Plot -", contrast_name),
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 2.0,
    labSize = 3.0,
    maxoverlapsConnectors = 10,
    drawConnectors = TRUE
  )
  
  ggsave(
    file.path(output_dir, paste0(contrast_name, "_volcano_plot.png")),
    volcano_plot, width = 10, height = 8, dpi = 300
  )
}

# Summary statistics
cat("\nSummary of Differential Expression Analysis:\n")
cat("============================================\n")

summary_stats <- data.frame(
  Contrast = names(significant_genes),
  Significant_Genes = sapply(significant_genes, length),
  Upregulated = sapply(names(significant_genes), function(x) {
    sum(de_results[[x]]$padj < 0.05 & de_results[[x]]$log2FoldChange > 1)
  }),
  Downregulated = sapply(names(significant_genes), function(x) {
    sum(de_results[[x]]$padj < 0.05 & de_results[[x]]$log2FoldChange < -1)
  })
)

print(summary_stats)
write_csv(summary_stats, file.path(output_dir, "de_summary.csv"))

# Venn diagram of significant genes (if multiple contrasts)
if (length(significant_genes) > 1) {
  cat("Generating Venn diagram of significant genes...\n")
  
  # Limit to first 3 contrasts for visualization
  sig_genes_for_venn <- significant_genes[1:min(3, length(significant_genes))]
  
  venn_plot <- venn.diagram(
    sig_genes_for_venn,
    filename = file.path(output_dir, "significant_genes_venn.png"),
    category.names = names(sig_genes_for_venn),
    output = TRUE,
    imagetype = "png",
    height = 480,
    width = 480,
    resolution = 300,
    compression = "lzw",
    fill = brewer.pal(length(sig_genes_for_venn), "Set2"),
    cex = 0.6,
    fontface = "bold",
    fontfamily = "sans"
  )
}

# Top variable genes heatmap
cat("Generating heatmap of top variable genes...\n")
top_var_genes <- order(rowVars(assay(vst_counts)), decreasing = TRUE)[1:50]
top_var_matrix <- assay(vst_counts)[top_var_genes, ]

# Scale for visualization
top_var_matrix_scaled <- t(scale(t(top_var_matrix)))

pheatmap(
  top_var_matrix_scaled,
  filename = file.path(output_dir, "top_variable_genes_heatmap.png"),
  width = 10, height = 12,
  annotation_col = metadata[, c("condition", "batch"), drop = FALSE],
  show_rownames = FALSE,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  main = "Top 50 Most Variable Genes"
)

# XRN2 rescue analysis
cat("Performing XRN2 rescue analysis...\n")

# Find genes that are:
# 1. Significantly different in siXRN2 vs Untreated
# 2. Rescued by WT but not by MT overexpression

rescue_analysis <- function() {
  if (!"siXRN2_vs_Untreated" %in% names(de_results) ||
      !"siXRN2_XRN2WT_OE_vs_siXRN2" %in% names(de_results) ||
      !"siXRN2_XRN2MT_OE_vs_siXRN2" %in% names(de_results)) {
    cat("Not all required contrasts available for rescue analysis\n")
    return(NULL)
  }
  
  # Get significantly affected genes in siXRN2 vs Untreated
  sirna_affected <- de_results$siXRN2_vs_Untreated %>%
    filter(padj < 0.05, abs(log2FoldChange) > 1)
  
  # Check rescue by WT
  wt_rescue <- de_results$siXRN2_XRN2WT_OE_vs_siXRN2 %>%
    filter(gene_id %in% sirna_affected$gene_id)
  
  # Check rescue by MT
  mt_rescue <- de_results$siXRN2_XRN2MT_OE_vs_siXRN2 %>%
    filter(gene_id %in% sirna_affected$gene_id)
  
  # Combine results
  rescue_df <- sirna_affected %>%
    select(gene_id, sirna_log2FC = log2FoldChange, sirna_padj = padj) %>%
    left_join(
      wt_rescue %>% select(gene_id, wt_rescue_log2FC = log2FoldChange, wt_rescue_padj = padj),
      by = "gene_id"
    ) %>%
    left_join(
      mt_rescue %>% select(gene_id, mt_rescue_log2FC = log2FoldChange, mt_rescue_padj = padj),
      by = "gene_id"
    ) %>%
    mutate(
      wt_rescues = (sign(sirna_log2FC) != sign(wt_rescue_log2FC)) & (wt_rescue_padj < 0.05),
      mt_rescues = (sign(sirna_log2FC) != sign(mt_rescue_log2FC)) & (mt_rescue_padj < 0.05),
      wt_specific_rescue = wt_rescues & !mt_rescues
    ) %>%
    replace_na(list(wt_rescues = FALSE, mt_rescues = FALSE, wt_specific_rescue = FALSE))
  
  # Save rescue analysis
  write_csv(rescue_df, file.path(output_dir, "xrn2_rescue_analysis.csv"))
  
  # Summary of rescue
  rescue_summary <- rescue_df %>%
    summarise(
      total_sirna_affected = n(),
      wt_rescued = sum(wt_rescues, na.rm = TRUE),
      mt_rescued = sum(mt_rescues, na.rm = TRUE),
      wt_specific_rescued = sum(wt_specific_rescue, na.rm = TRUE)
    )
  
  cat("XRN2 Rescue Analysis Summary:\n")
  cat("  Total siRNA-affected genes:", rescue_summary$total_sirna_affected, "\n")
  cat("  WT-rescued genes:", rescue_summary$wt_rescued, "\n")
  cat("  MT-rescued genes:", rescue_summary$mt_rescued, "\n")
  cat("  WT-specific rescued genes:", rescue_summary$wt_specific_rescued, "\n")
  
  return(rescue_df)
}

rescue_results <- rescue_analysis()

# Generate biological insights summary
cat("Generating biological insights...\n")

insights <- c(
  "# Top 5 XRN2 Biology Insights",
  "## Based on Differential Expression Analysis",
  "",
  paste("1. Total genes significantly affected by siXRN2:", 
        length(significant_genes$siXRN2_vs_Untreated)),
  paste("2. Genes rescued by XRN2-WT overexpression:", 
        ifelse(is.null(rescue_results), "N/A", sum(rescue_results$wt_rescues, na.rm = TRUE))),
  paste("3. Genes rescued by XRN2-MT overexpression:", 
        ifelse(is.null(rescue_results), "N/A", sum(rescue_results$mt_rescues, na.rm = TRUE))),
  paste("4. WT-specific rescue genes (potential XRN2 targets):", 
        ifelse(is.null(rescue_results), "N/A", sum(rescue_results$wt_specific_rescue, na.rm = TRUE))),
  "5. See pathway analysis results for enriched biological processes"
)

writeLines(insights, file.path(output_dir, "biological_insights.txt"))

# Session info
cat("Saving session info...\n")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))

cat("\nDESeq2 analysis completed successfully!\n")
cat("Results saved to:", output_dir, "\n")

# Summary of output files
cat("\nOutput files generated:\n")
output_files <- c(
  "normalized_counts.tsv",
  "vst_counts.tsv", 
  "de_summary.csv",
  "pca_plot.png",
  "sample_correlation_heatmap.png",
  "top_variable_genes_heatmap.png",
  "biological_insights.txt",
  "session_info.txt"
)

for (contrast_name in names(contrasts)) {
  output_files <- c(
    output_files,
    paste0(contrast_name, "_results.csv"),
    paste0(contrast_name, "_ma_plot.png"),
    paste0(contrast_name, "_volcano_plot.png")
  )
}

if (!is.null(rescue_results)) {
  output_files <- c(output_files, "xrn2_rescue_analysis.csv")
}

if (length(significant_genes) > 1) {
  output_files <- c(output_files, "significant_genes_venn.png")
}

for (file in output_files) {
  file_path <- file.path(output_dir, file)
  if (file.exists(file_path)) {
    cat("  ✓", file, "\n")
  } else {
    cat("  ✗", file, "(not generated)\n")
  }
}

cat("\nAnalysis complete!\n")
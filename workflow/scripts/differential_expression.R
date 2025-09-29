#!/usr/bin/env Rscript
# Differential Expression Analysis using DESeq2

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(dplyr)

# Get parameters from snakemake
count_file <- snakemake@input[["counts"]]
metadata_file <- snakemake@input[["metadata"]]
output_results <- snakemake@output[["results"]]
output_plots <- snakemake@output[["plots"]]
comparison_name <- snakemake@params[["comparison"]]
output_dir <- snakemake@params[["outdir"]]

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
count_matrix <- read.table(count_file, header = TRUE, row.names = 1, sep = "\t")
metadata <- read.table(metadata_file, header = TRUE, sep = "\t")

# Parse comparison from config
# This is a simplified version - in practice, you'd parse the YAML config
# For now, assume comparison_name follows pattern: "Treat1_vs_Control1"
comparison_parts <- strsplit(comparison_name, "_vs_")[[1]]
case_condition <- comparison_parts[1]
control_condition <- comparison_parts[2]

# Filter samples for this comparison
case_samples <- metadata[metadata$condition == case_condition, "sample_id"]
control_samples <- metadata[metadata$condition == control_condition, "sample_id"]
comparison_samples <- c(case_samples, control_samples)

# Filter count matrix and metadata
count_subset <- count_matrix[, comparison_samples]
metadata_subset <- metadata[metadata$sample_id %in% comparison_samples, ]

# Ensure sample order matches
metadata_subset <- metadata_subset[match(colnames(count_subset), metadata_subset$sample_id), ]

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = count_subset,
  colData = metadata_subset,
  design = ~ batch + condition
)

# Filter low count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Set reference level
dds$condition <- relevel(dds$condition, ref = control_condition)

# Run DESeq2
dds <- DESeq(dds)

# Get results
res <- results(dds, contrast = c("condition", case_condition, control_condition))
res <- lfcShrink(dds, contrast = c("condition", case_condition, control_condition), res = res, type = "normal")

# Convert to dataframe and add gene names
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df <- res_df[order(res_df$pvalue), ]

# Save results
write.table(res_df, output_results, sep = "\t", quote = FALSE, row.names = FALSE)

# Generate plots
pdf(output_plots, width = 12, height = 8)

# MA plot
plotMA(res, main = paste("MA Plot:", comparison_name))

# Volcano plot
if (nrow(res_df) > 0) {
  volcano_plot <- EnhancedVolcano(
    res_df,
    lab = res_df$gene_id,
    x = 'log2FoldChange',
    y = 'pvalue',
    title = paste('Volcano Plot:', comparison_name),
    subtitle = paste('Case:', case_condition, 'vs Control:', control_condition),
    pCutoff = 0.05,
    FCcutoff = 1.0,
    labSize = 3.0
  )
  print(volcano_plot)
}

# PCA plot
vsd <- vst(dds, blind = FALSE)
pca_plot <- plotPCA(vsd, intgroup = c("condition", "batch")) +
  ggtitle(paste("PCA Plot:", comparison_name))
print(pca_plot)

# Heatmap of top genes
if (sum(!is.na(res_df$padj) & res_df$padj < 0.05) > 1) {
  top_genes <- head(res_df[!is.na(res_df$padj) & res_df$padj < 0.05, ], 50)
  
  if (nrow(top_genes) > 1) {
    vsd_matrix <- assay(vsd)[rownames(top_genes), ]
    pheatmap(
      vsd_matrix,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = FALSE,
      annotation_col = data.frame(
        condition = metadata_subset$condition,
        batch = metadata_subset$batch,
        row.names = metadata_subset$sample_id
      ),
      main = paste("Top 50 DE Genes:", comparison_name)
    )
  }
}

dev.off()

# Print summary
cat("Differential expression analysis completed for:", comparison_name, "\n")
cat("Total genes tested:", nrow(res_df), "\n")
cat("Significant genes (padj < 0.05):", sum(!is.na(res_df$padj) & res_df$padj < 0.05), "\n")
cat("Upregulated genes (LFC > 1, padj < 0.05):", sum(!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange > 1), "\n")
cat("Downregulated genes (LFC < -1, padj < 0.05):", sum(!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange < -1), "\n")
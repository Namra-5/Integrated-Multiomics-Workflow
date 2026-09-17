if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2", "EnhancedVolcano", "pheatmap"))
install.packages(c("ggplot2", "dplyr", "RColorBrewer"))
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
# Load count matrix (FeatureCounts output)
counts <- read.table("FeatureCounts_Mod_clean.txt", 
                     header = TRUE, 
                     row.names = 1,    # Gene IDs as row names
                     sep = "\t",
                     check.names = FALSE)

# Load sample metadata
coldata <- read.csv("coldata.csv", 
                    header = TRUE, 
                    row.names = 1)     # Sample names as row names

# Quick check
head(counts)
head(coldata)
dim(counts)
# Keep only numeric count columns (drop annotation columns if present)
counts <- counts[, sapply(counts, is.numeric)]

# Make sure column names of counts match row names of coldata
# (They MUST match exactly — same order)
all(colnames(counts) == rownames(coldata))  # Should print TRUE

# If not matching, reorder:
counts <- counts[, rownames(coldata)]

# Make sure condition column is a factor
coldata$condition <- factor(coldata$condition, levels = c("Normal", "COVID"))

# Construct DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = coldata,
                              design    = ~ condition)

# Pre-filter: remove genes with very low counts
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

# Get results: COVID vs Normal
res <- results(dds, contrast = c("condition", "COVID", "Normal"))

# Summary
summary(res)

# Convert to data frame
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# Filter significant DEGs: padj < 0.01
sig_genes <- res_df[!is.na(res_df$padj) & res_df$padj < 0.01, ]

cat("Number of significant DEGs:", nrow(sig_genes), "\n")

# Save results
write.csv(res_df, "DEG_results_all.csv", row.names = TRUE)
write.csv(sig_genes, "DEG_results_significant.csv", row.names = TRUE)


png("volcano_plot.png", width = 2400, height = 2000, res = 300)

EnhancedVolcano(res_df,
                lab        = res_df$gene,
                x          = 'log2FoldChange',
                y          = 'padj',
                title      = 'COVID vs Normal',
                subtitle   = 'DESeq2 | FDR < 0.01',
                pCutoff    = 0.01,
                FCcutoff   = 1.5,
                pointSize  = 2.0,
                labSize    = 3.0,
                col        = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                legendPosition = 'right')

dev.off()
cat("Volcano plot saved!\n")



# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)[rownames(sig_genes), ]

# Scale each gene (row)
mat_scaled <- t(scale(t(mat)))

# Color palette
colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)

# Annotation bar showing Normal/COVID
anno <- as.data.frame(colData(dds)[, "condition", drop = FALSE])

png("heatmap.png", width = 2400, height = 3000, res = 300)

pheatmap(mat_scaled,
         annotation_col  = anno,
         color           = colors,
         show_rownames   = FALSE,
         show_colnames   = TRUE,
         cluster_rows    = TRUE,
         cluster_cols    = TRUE,
         fontsize        = 8,
         main            = paste0("Heatmap of Significant DEGs (n=", nrow(sig_genes), ")"))

dev.off()
cat("Heatmap saved!\n")



# Sort by log2FoldChange
sig_sorted <- sig_genes[order(sig_genes$log2FoldChange, decreasing = TRUE), ]

# Top 20 Upregulated
top20_up <- head(sig_sorted, 20)

# Top 20 Downregulated
top20_down <- tail(sig_sorted, 20)

# Save both
write.csv(top20_up,   "Top20_Upregulated.csv",   row.names = TRUE)
write.csv(top20_down, "Top20_Downregulated.csv", row.names = TRUE)

# View in console
cat("\n--- TOP 20 UPREGULATED ---\n")
print(top20_up[, c("gene", "log2FoldChange", "pvalue", "padj")])

cat("\n--- TOP 20 DOWNREGULATED ---\n")
print(top20_down[, c("gene", "log2FoldChange", "pvalue", "padj")])

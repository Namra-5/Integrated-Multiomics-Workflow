# =============================================================================
# TASK 7: Integration of GWAS, Gene Expression, ncRNA, ChIP-seq & QTL Concepts
# =============================================================================
# This script integrates results from all previous tasks to identify candidate
# genes involved in the studied disease condition. We look for genes that appear
# across multiple layers of evidence (DEGs, NF-kB targets, ncRNA targets, GWAS).
# =============================================================================

# ─────────────────────────────────────────────────────────────────────────────
# STEP 0: Install and Load Required Libraries
# ─────────────────────────────────────────────────────────────────────────────
# These packages handle data manipulation, visualization, and Venn diagrams.
# Run install lines only once; comment them out after first use.

# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("ggvenn")
# install.packages("readxl")
# install.packages("writexl")
# install.packages("gridExtra")
# install.packages("RColorBrewer")
# install.packages("ggrepel")

library(dplyr)        # Data manipulation (filter, merge, select, mutate)
library(ggplot2)      # High-quality plots
library(ggvenn)       # Venn diagram of overlapping gene sets
library(readxl)       # Read Excel files (.xlsx)
library(writexl)      # Write Excel output files
library(gridExtra)    # Arrange multiple plots on one page
library(RColorBrewer) # Color palettes for plots
library(ggrepel)      # Non-overlapping labels on scatter plots

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: Load All Input Files
# ─────────────────────────────────────────────────────────────────────────────
# Each file comes from a different task. We load them all here.
# Make sure all files are in your working directory, or provide full paths.

# --- Task 2: Significant DEGs (filtered, padj < 0.01) ---
# Contains Ensembl gene IDs, log2FoldChange, padj values
deg <- read.csv("DEG_results_significant.csv", stringsAsFactors = FALSE)
cat("Task 2 - DEGs loaded:", nrow(deg), "significant genes\n")

# --- Task 4: ncRNA-miRNA interactions ---
# Contains miRNA names, target gene symbols, scores, direction
ncrna <- read_excel("task4_interactions_complete.xlsx", sheet = "ncRNA Interactions")
ncrna <- as.data.frame(ncrna)
cat("Task 4 - ncRNA interactions loaded:", nrow(ncrna), "rows\n")

# --- Task 5: NF-kB ChIP-seq annotated peaks ---
# Contains peak coordinates AND annotated gene symbols + Ensembl IDs
chipseq <- read.csv("NF-kB_annotated_peaks.csv", stringsAsFactors = FALSE)
cat("Task 5 - ChIP-seq peaks loaded:", nrow(chipseq), "annotated peaks\n")

# --- Task 6: GWAS significant SNPs ---
# Contains chromosome, position, SNP ID, p-value
gwas <- read.csv("GWAS_GW_significant_SNPs.csv", stringsAsFactors = FALSE)
cat("Task 6 - GWAS SNPs loaded:", nrow(gwas), "significant SNPs\n")

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: Prepare Gene Lists from Each Data Source
# ─────────────────────────────────────────────────────────────────────────────
# Before comparing across datasets, we extract clean gene lists.
# The challenge: DEGs use Ensembl IDs; ncRNA uses gene symbols; 
# ChIP-seq has both; GWAS has only SNP positions.

# --- 2A: DEG gene list (Ensembl IDs) ---
deg_ensembl <- unique(deg$gene)
cat("Unique DEG Ensembl IDs:", length(deg_ensembl), "\n")

# --- 2B: ChIP-seq gene list ---
# NF-kB peaks annotated to genes - extract unique Ensembl IDs and Symbols
# Remove rows where ENSEMBL is missing (some peaks map to unannotated regions)
chipseq_clean <- chipseq %>%
  filter(!is.na(ENSEMBL) & ENSEMBL != "" & !is.na(SYMBOL)) %>%
  select(ENSEMBL, SYMBOL) %>%
  distinct()

chipseq_ensembl <- unique(chipseq_clean$ENSEMBL)
chipseq_symbols  <- unique(chipseq_clean$SYMBOL)
cat("Unique NF-kB target genes (Ensembl):", length(chipseq_ensembl), "\n")

# --- 2C: ncRNA target gene list (Symbols) ---
# Remove rows with "N/A" targets (genes with no miRNA found)
ncrna_clean <- ncrna %>%
  filter(miRNA != "N/A" & `Target Gene` != "N/A")

ncrna_targets_symbols <- unique(ncrna_clean$`Target Gene`)
cat("Unique ncRNA target genes (Symbols):", length(ncrna_targets_symbols), "\n")

# Convert ncRNA target symbols to Ensembl using ChIP-seq as reference table
# (ChIP-seq has the most complete SYMBOL<->ENSEMBL mapping in our data)
symbol_to_ensembl <- chipseq_clean %>%
  filter(SYMBOL %in% ncrna_targets_symbols) %>%
  select(SYMBOL, ENSEMBL) %>%
  distinct()

ncrna_targets_ensembl <- unique(symbol_to_ensembl$ENSEMBL)
cat("ncRNA targets mapped to Ensembl:", length(ncrna_targets_ensembl), 
    "out of", length(ncrna_targets_symbols), "symbols\n")

# --- 2D: GWAS locus genes ---
# GWAS SNPs are at chr3:~45.8Mb and chr9:~133Mb
# We identify nearby genes by checking which NF-kB peaks (with gene annotations)
# fall within ±500kb of the GWAS SNP positions. This is our proxy for 
# "GWAS locus genes" since we don't have a separate gene annotation file.
gwas_locus_genes <- data.frame()

for (i in 1:nrow(gwas)) {
  chr_query <- paste0("chr", gwas$CHR[i])
  bp        <- gwas$BP[i]
  window    <- 500000  # 500kb window around each SNP
  
  nearby <- chipseq %>%
    filter(seqnames == chr_query,
           start >= (bp - window),
           end   <= (bp + window),
           !is.na(ENSEMBL), !is.na(SYMBOL)) %>%
    select(ENSEMBL, SYMBOL) %>%
    distinct() %>%
    mutate(SNP = gwas$SNP[i], SNP_P = gwas$P[i], CHR = gwas$CHR[i], BP = bp)
  
  gwas_locus_genes <- bind_rows(gwas_locus_genes, nearby)
}

gwas_locus_genes <- gwas_locus_genes %>% distinct(ENSEMBL, SYMBOL, .keep_all = TRUE)
gwas_ensembl     <- unique(gwas_locus_genes$ENSEMBL)
cat("GWAS locus genes (within 500kb of SNPs):", length(gwas_ensembl), "\n")

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3: Find Overlaps Between All Gene Lists
# ─────────────────────────────────────────────────────────────────────────────
# This is the core of Task 7 - genes appearing in multiple evidence layers
# are the strongest candidate genes.

# Pairwise overlaps
deg_chipseq_overlap   <- intersect(deg_ensembl, chipseq_ensembl)
deg_ncrna_overlap     <- intersect(deg_ensembl, ncrna_targets_ensembl)
deg_gwas_overlap      <- intersect(deg_ensembl, gwas_ensembl)
chipseq_ncrna_overlap <- intersect(chipseq_ensembl, ncrna_targets_ensembl)
chipseq_gwas_overlap  <- intersect(chipseq_ensembl, gwas_ensembl)
ncrna_gwas_overlap    <- intersect(ncrna_targets_ensembl, gwas_ensembl)

cat("\n--- Pairwise Overlaps ---\n")
cat("DEG ∩ ChIP-seq (NF-kB targets):", length(deg_chipseq_overlap), "\n")
cat("DEG ∩ ncRNA targets:            ", length(deg_ncrna_overlap), "\n")
cat("DEG ∩ GWAS locus genes:         ", length(deg_gwas_overlap), "\n")
cat("ChIP-seq ∩ ncRNA targets:       ", length(chipseq_ncrna_overlap), "\n")
cat("ChIP-seq ∩ GWAS locus genes:    ", length(chipseq_gwas_overlap), "\n")
cat("ncRNA ∩ GWAS locus genes:       ", length(ncrna_gwas_overlap), "\n")

# Triple overlap: DEG + ChIP-seq + ncRNA (strongest candidates)
triple_overlap <- intersect(deg_chipseq_overlap, ncrna_targets_ensembl)
cat("\nTriple overlap (DEG + ChIP-seq + ncRNA):", length(triple_overlap), "\n")
cat("Triple overlap genes (Ensembl):", triple_overlap, "\n")

# Quadruple overlap: all four
quad_overlap <- intersect(triple_overlap, gwas_ensembl)
cat("Quadruple overlap (all four):", length(quad_overlap), "\n")

# ─────────────────────────────────────────────────────────────────────────────
# STEP 4: Build the Integrated Evidence Table
# ─────────────────────────────────────────────────────────────────────────────
# This table shows, for every candidate gene, which evidence layers support it.
# This is the key deliverable for Task 7.

# Start from all genes in at least ONE evidence layer
all_candidate_ensembl <- unique(c(deg_ensembl, chipseq_ensembl, 
                                   ncrna_targets_ensembl, gwas_ensembl))

# Get gene symbols for all candidates using ChIP-seq symbol map
all_symbols_map <- chipseq_clean %>%
  bind_rows(
    data.frame(ENSEMBL = ncrna_targets_ensembl,
               SYMBOL  = symbol_to_ensembl$SYMBOL[match(ncrna_targets_ensembl, 
                                                         symbol_to_ensembl$ENSEMBL)])
  ) %>%
  distinct(ENSEMBL, .keep_all = TRUE)

# Build evidence table
evidence_table <- data.frame(ENSEMBL = all_candidate_ensembl) %>%
  left_join(all_symbols_map, by = "ENSEMBL") %>%
  # Join DEG info
  left_join(deg %>% select(gene, log2FoldChange, padj) %>%
              rename(ENSEMBL = gene),
            by = "ENSEMBL") %>%
  mutate(
    Is_DEG          = ENSEMBL %in% deg_ensembl,
    Is_ChIPseq      = ENSEMBL %in% chipseq_ensembl,
    Is_ncRNA_target = ENSEMBL %in% ncrna_targets_ensembl,
    Is_GWAS_gene    = ENSEMBL %in% gwas_ensembl,
    # Count how many evidence layers support each gene
    Evidence_Count  = Is_DEG + Is_ChIPseq + Is_ncRNA_target + Is_GWAS_gene,
    # Direction of expression change
    Direction = case_when(
      Is_DEG & log2FoldChange > 0 ~ "Upregulated",
      Is_DEG & log2FoldChange < 0 ~ "Downregulated",
      TRUE ~ "Not a DEG"
    )
  ) %>%
  # Focus on genes supported by at least 2 evidence layers
  filter(Evidence_Count >= 2) %>%
  arrange(desc(Evidence_Count), desc(abs(log2FoldChange)))

cat("\nIntegrated evidence table: genes with ≥2 evidence layers:", nrow(evidence_table), "\n")
cat("Genes with all 3 main layers (DEG+ChIP+ncRNA):", 
    sum(evidence_table$Is_DEG & evidence_table$Is_ChIPseq & evidence_table$Is_ncRNA_target), "\n")

# ─────────────────────────────────────────────────────────────────────────────
# STEP 5: Extract Final Candidate Gene List
# ─────────────────────────────────────────────────────────────────────────────
# Candidate genes = those supported by 2 or more evidence layers.
# The strongest candidates (3 layers) are MMP8, CD209, IL23R.

candidate_genes <- evidence_table %>%
  filter(!is.na(SYMBOL)) %>%
  select(ENSEMBL, SYMBOL, log2FoldChange, padj, Direction,
         Is_DEG, Is_ChIPseq, Is_ncRNA_target, Is_GWAS_gene, Evidence_Count)

cat("\n--- Top Candidate Genes (≥2 evidence layers) ---\n")
print(candidate_genes %>% filter(Evidence_Count >= 3))

# ─────────────────────────────────────────────────────────────────────────────
# STEP 6: Venn Diagram of Gene Set Overlaps
# ─────────────────────────────────────────────────────────────────────────────
# A Venn diagram visually shows how many genes are shared between evidence layers.
# This is required as the "overlap plot" for Task 7.

# Prepare named list for ggvenn (uses gene symbols where possible)
venn_list <- list(
  "DEGs (Task 2)"          = deg_ensembl,
  "NF-kB Targets (Task 5)" = chipseq_ensembl,
  "ncRNA Targets (Task 4)" = ncrna_targets_ensembl,
  "GWAS Loci (Task 6)"     = gwas_ensembl
)

venn_plot <- ggvenn(
  venn_list,
  fill_color   = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
  stroke_size  = 0.8,
  set_name_size = 4,
  text_size    = 3.5
) +
  labs(
    title    = "Gene Set Overlap Across Evidence Layers",
    subtitle = "Integration of DEG, ChIP-seq, ncRNA, and GWAS Results"
  ) +
  theme(
    plot.title    = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10)
  )

ggsave("Task7_Venn_Diagram.png", venn_plot, width = 10, height = 8, dpi = 300)
cat("Saved: Task7_Venn_Diagram.png\n")

# ─────────────────────────────────────────────────────────────────────────────
# STEP 7: Evidence Heatmap of Top Candidate Genes
# ─────────────────────────────────────────────────────────────────────────────
# This plot shows which evidence layers support each candidate gene.
# A filled cell = gene is supported by that evidence layer.

# Take top candidates (Evidence_Count >= 2 and has a symbol)
top_candidates <- candidate_genes %>%
  filter(!is.na(SYMBOL)) %>%
  head(30)  # top 30 by evidence count then fold change

# Reshape to long format for ggplot heatmap
heatmap_data <- top_candidates %>%
  select(SYMBOL, Is_DEG, Is_ChIPseq, Is_ncRNA_target, Is_GWAS_gene) %>%
  tidyr::pivot_longer(
    cols      = starts_with("Is_"),
    names_to  = "Evidence_Layer",
    values_to = "Present"
  ) %>%
  mutate(
    Evidence_Layer = recode(Evidence_Layer,
      "Is_DEG"          = "Differentially\nExpressed (Task 2)",
      "Is_ChIPseq"      = "NF-kB ChIP-seq\nTarget (Task 5)",
      "Is_ncRNA_target" = "ncRNA Regulated\n(Task 4)",
      "Is_GWAS_gene"    = "Near GWAS\nLocus (Task 6)"
    ),
    Present = as.numeric(Present)
  )

evidence_heatmap <- ggplot(heatmap_data, aes(x = Evidence_Layer, y = SYMBOL, fill = factor(Present))) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(
    values = c("0" = "#f0f0f0", "1" = "#2166AC"),
    labels = c("0" = "Absent", "1" = "Present"),
    name   = "Evidence"
  ) +
  labs(
    title    = "Multi-Omics Evidence for Candidate Genes",
    subtitle = "Blue = gene supported by that evidence layer",
    x        = "Evidence Layer",
    y        = "Gene Symbol"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle   = element_text(hjust = 0.5, size = 9),
    axis.text.x     = element_text(size = 9),
    axis.text.y     = element_text(size = 8),
    panel.grid      = element_blank(),
    legend.position = "right"
  )

ggsave("Task7_Evidence_Heatmap.png", evidence_heatmap, width = 10, height = 12, dpi = 300)
cat("Saved: Task7_Evidence_Heatmap.png\n")

# ─────────────────────────────────────────────────────────────────────────────
# STEP 8: Volcano-style Plot for Triple-Overlap Candidate Genes
# ─────────────────────────────────────────────────────────────────────────────
# Shows ALL significant DEGs on a scatter plot (LFC vs -log10 padj),
# highlighting the top candidate genes that also appear in ChIP-seq and ncRNA.

deg_with_symbols <- deg %>%
  left_join(chipseq_clean, by = c("gene" = "ENSEMBL")) %>%
  distinct(gene, .keep_all = TRUE) %>%
  mutate(
    neg_log10_padj = -log10(padj),
    Candidate      = gene %in% triple_overlap,
    Point_Color    = case_when(
      Candidate & log2FoldChange > 0 ~ "Triple-overlap Up",
      Candidate & log2FoldChange < 0 ~ "Triple-overlap Down",
      log2FoldChange > 1  ~ "Upregulated DEG",
      log2FoldChange < -1 ~ "Downregulated DEG",
      TRUE                ~ "Non-significant DEG"
    )
  )

candidate_labels <- deg_with_symbols %>% filter(Candidate)

highlight_plot <- ggplot(deg_with_symbols, aes(x = log2FoldChange, y = neg_log10_padj, color = Point_Color)) +
  geom_point(alpha = 0.4, size = 1.2) +
  geom_point(data = candidate_labels, size = 4, alpha = 1) +
  geom_label_repel(
    data        = candidate_labels,
    aes(label   = SYMBOL),
    size        = 4,
    fontface    = "bold",
    box.padding = 0.5,
    max.overlaps = 20
  ) +
  scale_color_manual(values = c(
    "Triple-overlap Up"    = "#D73027",
    "Triple-overlap Down"  = "#4575B4",
    "Upregulated DEG"      = "#FC8D59",
    "Downregulated DEG"    = "#91BFDB",
    "Non-significant DEG"  = "grey70"
  )) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "grey40") +
  labs(
    title    = "DEG Volcano Plot with Multi-Omics Candidate Genes Highlighted",
    subtitle = "Red/Blue large points = genes also supported by ChIP-seq AND ncRNA evidence",
    x        = "log2 Fold Change",
    y        = "-log10(adjusted p-value)",
    color    = "Gene Category"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    legend.position = "right"
  )

ggsave("Task7_Candidate_Highlight_Plot.png", highlight_plot, width = 12, height = 8, dpi = 300)
cat("Saved: Task7_Candidate_Highlight_Plot.png\n")

# ─────────────────────────────────────────────────────────────────────────────
# STEP 9: GWAS Locus Gene Bar Chart
# ─────────────────────────────────────────────────────────────────────────────
# Shows the NF-kB target genes found within GWAS SNP loci regions,
# which are the most likely regulatory targets of those SNPs.

gwas_chip_genes <- gwas_locus_genes %>%
  filter(ENSEMBL %in% chipseq_ensembl) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  arrange(SNP_P) %>%
  head(20)

if (nrow(gwas_chip_genes) > 0) {
  gwas_bar <- ggplot(gwas_chip_genes, aes(x = reorder(SYMBOL, -SNP_P), 
                                           y = -log10(SNP_P), fill = factor(CHR))) +
    geom_bar(stat = "identity", color = "white") +
    coord_flip() +
    scale_fill_brewer(palette = "Set1", name = "Chromosome") +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red") +
    labs(
      title    = "NF-kB Target Genes Near GWAS Significant SNP Loci",
      subtitle = "Red dashed line = genome-wide significance threshold (p < 5×10⁻⁸)",
      x        = "Gene (NF-kB ChIP-seq target)",
      y        = "-log10(GWAS p-value)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 9)
    )
  
  ggsave("Task7_GWAS_ChIPseq_Overlap.png", gwas_bar, width = 10, height = 7, dpi = 300)
  cat("Saved: Task7_GWAS_ChIPseq_Overlap.png\n")
}

# ─────────────────────────────────────────────────────────────────────────────
# STEP 10: eQTL Conceptual Workflow Diagram
# ─────────────────────────────────────────────────────────────────────────────
# Creates a simple text-based flow diagram explaining the SNP->expression model.
# This satisfies the "Conceptual eQTL/QTL workflow diagram" requirement.

eqtl_diagram <- ggplot() +
  # Draw boxes
  annotate("rect", xmin=1, xmax=3, ymin=8.5, ymax=9.5, fill="#AED6F1", color="#2980B9", linewidth=1) +
  annotate("rect", xmin=1, xmax=3, ymin=6.5, ymax=7.5, fill="#A9DFBF", color="#27AE60", linewidth=1) +
  annotate("rect", xmin=1, xmax=3, ymin=4.5, ymax=5.5, fill="#F9E79F", color="#F39C12", linewidth=1) +
  annotate("rect", xmin=1, xmax=3, ymin=2.5, ymax=3.5, fill="#F5CBA7", color="#E67E22", linewidth=1) +
  annotate("rect", xmin=1, xmax=3, ymin=0.5, ymax=1.5, fill="#D7BDE2", color="#8E44AD", linewidth=1) +
  # Box labels
  annotate("text", x=2, y=9,   label="GWAS SNP\n(Genetic Variant)", size=4, fontface="bold") +
  annotate("text", x=2, y=7,   label="Alters TF Binding Site\n(NF-kB ChIP-seq - Task 5)", size=4, fontface="bold") +
  annotate("text", x=2, y=5,   label="Changes Chromatin State\n& Promoter Accessibility", size=4, fontface="bold") +
  annotate("text", x=2, y=3,   label="Disrupts miRNA Regulation\n(ncRNA Analysis - Task 4)", size=4, fontface="bold") +
  annotate("text", x=2, y=1,   label="Altered Gene Expression\n(DEG - Task 2)\ne.g. MMP8↑, CD209↓, IL23R↓", size=4, fontface="bold") +
  # Arrows
  annotate("segment", x=2, xend=2, y=8.5, yend=7.5, 
           arrow=arrow(length=unit(0.3,"cm"), type="closed"), linewidth=1.2, color="grey30") +
  annotate("segment", x=2, xend=2, y=6.5, yend=5.5,
           arrow=arrow(length=unit(0.3,"cm"), type="closed"), linewidth=1.2, color="grey30") +
  annotate("segment", x=2, xend=2, y=4.5, yend=3.5,
           arrow=arrow(length=unit(0.3,"cm"), type="closed"), linewidth=1.2, color="grey30") +
  annotate("segment", x=2, xend=2, y=2.5, yend=1.5,
           arrow=arrow(length=unit(0.3,"cm"), type="closed"), linewidth=1.2, color="grey30") +
  # Side label for eQTL
  annotate("rect", xmin=3.3, xmax=4.8, ymin=4.5, ymax=9.5, fill="#FDFEFE", color="#BDC3C7", linewidth=0.8, linetype="dashed") +
  annotate("text", x=4.05, y=7, label="eQTL:\nSNP genotype\ncorrelates with\ngene expression\nlevel", size=3.5, color="#7F8C8D") +
  annotate("segment", x=3, xend=3.3, y=7, yend=7, linetype="dashed", color="#BDC3C7") +
  # Title
  labs(title = "SNP → Regulation → Expression Model (eQTL Concept)",
       subtitle = "How a GWAS variant propagates through regulatory layers to affect gene expression") +
  xlim(0.5, 5) + ylim(0, 10.5) +
  theme_void() +
  theme(
    plot.title    = element_text(hjust=0.5, size=13, face="bold"),
    plot.subtitle = element_text(hjust=0.5, size=9, color="grey40"),
    plot.margin   = margin(10,10,10,10)
  )

ggsave("Task7_eQTL_Workflow_Diagram.png", eqtl_diagram, width=10, height=9, dpi=300)
cat("Saved: Task7_eQTL_Workflow_Diagram.png\n")

# ─────────────────────────────────────────────────────────────────────────────
# STEP 11: Save All Output Tables
# ─────────────────────────────────────────────────────────────────────────────

# 11A: Final candidate gene list (genes in ≥2 evidence layers)
write.csv(candidate_genes, "Task7_Candidate_Genes.csv", row.names = FALSE)
cat("Saved: Task7_Candidate_Genes.csv\n")

# 11B: Full integrated evidence table
write.csv(evidence_table, "Task7_Integrated_Evidence_Table.csv", row.names = FALSE)
cat("Saved: Task7_Integrated_Evidence_Table.csv\n")

# 11C: GWAS locus genes with NF-kB overlap
write.csv(gwas_locus_genes, "Task7_GWAS_Locus_Genes.csv", row.names = FALSE)
cat("Saved: Task7_GWAS_Locus_Genes.csv\n")

# 11D: Triple overlap summary (strongest candidates)
triple_summary <- candidate_genes %>%
  filter(Is_DEG & Is_ChIPseq & Is_ncRNA_target) %>%
  arrange(desc(abs(log2FoldChange)))
write.csv(triple_summary, "Task7_Triple_Overlap_Candidates.csv", row.names = FALSE)
cat("Saved: Task7_Triple_Overlap_Candidates.csv\n")

# ─────────────────────────────────────────────────────────────────────────────
# STEP 12: Print Final Summary
# ─────────────────────────────────────────────────────────────────────────────

cat("\n============================================================\n")
cat("             TASK 7 INTEGRATION SUMMARY\n")
cat("============================================================\n")
cat("Total significant DEGs (Task 2):              ", length(deg_ensembl), "\n")
cat("Total NF-kB ChIP-seq target genes (Task 5):  ", length(chipseq_ensembl), "\n")
cat("Total ncRNA-regulated genes (Task 4):         ", length(ncrna_targets_symbols), "\n")
cat("Total GWAS locus genes (Task 6):              ", length(gwas_ensembl), "\n")
cat("------------------------------------------------------------\n")
cat("DEG ∩ ChIP-seq overlap:                       ", length(deg_chipseq_overlap), "\n")
cat("DEG ∩ ncRNA targets:                          ", length(deg_ncrna_overlap), "\n")
cat("Triple overlap (DEG + ChIP-seq + ncRNA):      ", length(triple_overlap), "\n")
cat("------------------------------------------------------------\n")
cat("TOP CANDIDATE GENES (3 evidence layers):\n")
print(triple_summary %>% select(SYMBOL, log2FoldChange, padj, Direction))
cat("============================================================\n")
cat("\nAll outputs saved. Task 7 complete!\n")

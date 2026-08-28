# =============================================================
# TASK 3: GO and KEGG Enrichment Analysis
# COVID-19 vs Normal Samples
# Tools: biomaRt + enrichR + ggplot2
# =============================================================

# -------------------------------------------------------------
# STEP 1: Install Required Packages
# -------------------------------------------------------------
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
install.packages("enrichR")
install.packages("ggplot2")

# -------------------------------------------------------------
# STEP 2: Load Libraries
# -------------------------------------------------------------
library(biomaRt)
library(enrichR)
library(ggplot2)

# -------------------------------------------------------------
# STEP 3: Load DEG Results from Task 2
# -------------------------------------------------------------
setwd("C:/Users/ACG/Downloads/Genomics Task 2 DEG Analysis")

sig_genes   <- read.csv("DEG_results_significant.csv", row.names = 1)
ensembl_ids <- sig_genes$gene

cat("Total Ensembl IDs to convert:", length(ensembl_ids), "\n")

# -------------------------------------------------------------
# STEP 4: Convert Ensembl IDs to Gene Symbols via biomaRt
# -------------------------------------------------------------
mart <- useMart(biomart = "ensembl",
                dataset = "hsapiens_gene_ensembl",
                host    = "https://sep2025.archive.ensembl.org")

cat("Querying Ensembl... this may take 2-3 minutes...\n")

converted <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                   filters    = "ensembl_gene_id",
                   values     = ensembl_ids,
                   mart       = mart)

converted     <- converted[converted$hgnc_symbol != "", ]
valid_symbols <- unique(converted$hgnc_symbol)

cat("Successfully converted:", length(valid_symbols), "gene symbols\n")
cat("Sample symbols:", paste(head(valid_symbols, 5), collapse = ", "), "\n")

# -------------------------------------------------------------
# STEP 5: Run GO and KEGG Enrichment Analysis
# -------------------------------------------------------------
setEnrichrSite("Enrichr")

dbs <- c("GO_Biological_Process_2023",
         "GO_Molecular_Function_2023",
         "GO_Cellular_Component_2023",
         "KEGG_2021_Human")

cat("Running enrichment analysis...\n")
enriched <- enrichr(valid_symbols, dbs)

cat("GO BP terms found:",   nrow(enriched[["GO_Biological_Process_2023"]]),  "\n")
cat("GO MF terms found:",   nrow(enriched[["GO_Molecular_Function_2023"]]),  "\n")
cat("GO CC terms found:",   nrow(enriched[["GO_Cellular_Component_2023"]]),  "\n")
cat("KEGG pathways found:", nrow(enriched[["KEGG_2021_Human"]]),             "\n")

# -------------------------------------------------------------
# STEP 6: Filter Significant Results (adjusted p < 0.05)
# -------------------------------------------------------------
go_bp <- enriched[["GO_Biological_Process_2023"]]
go_bp <- go_bp[go_bp$Adjusted.P.value < 0.05, ]
go_bp <- go_bp[order(go_bp$Adjusted.P.value), ]

go_mf <- enriched[["GO_Molecular_Function_2023"]]
go_mf <- go_mf[go_mf$Adjusted.P.value < 0.05, ]
go_mf <- go_mf[order(go_mf$Adjusted.P.value), ]

go_cc <- enriched[["GO_Cellular_Component_2023"]]
go_cc <- go_cc[go_cc$Adjusted.P.value < 0.05, ]
go_cc <- go_cc[order(go_cc$Adjusted.P.value), ]

kegg <- enriched[["KEGG_2021_Human"]]
kegg <- kegg[kegg$Adjusted.P.value < 0.05, ]
kegg <- kegg[order(kegg$Adjusted.P.value), ]

cat("Significant GO BP terms:",  nrow(go_bp), "\n")
cat("Significant GO MF terms:",  nrow(go_mf), "\n")
cat("Significant GO CC terms:",  nrow(go_cc), "\n")
cat("Significant KEGG pathways:", nrow(kegg), "\n")

# -------------------------------------------------------------
# STEP 7: Save Result Tables
# -------------------------------------------------------------
write.csv(go_bp,          "GO_BP_results.csv",  row.names = FALSE)
write.csv(go_mf,          "GO_MF_results.csv",  row.names = FALSE)
write.csv(go_cc,          "GO_CC_results.csv",  row.names = FALSE)
write.csv(kegg,           "KEGG_results.csv",   row.names = FALSE)
write.csv(head(go_bp,20), "GO_BP_top20.csv",    row.names = FALSE)
write.csv(head(kegg, 20), "KEGG_top20.csv",     row.names = FALSE)

cat("All result tables saved!\n")

# -------------------------------------------------------------
# STEP 8: GO Biological Process Bar Plot
# -------------------------------------------------------------
top_bp       <- head(go_bp, 15)
top_bp$Term  <- gsub("\\s*\\(GO:\\d+\\)", "", top_bp$Term)

png("GO_BP_barplot.png", width = 3200, height = 2400, res = 300)
ggplot(top_bp, aes(x = reorder(Term, -Adjusted.P.value),
                   y = -log10(Adjusted.P.value),
                   fill = Adjusted.P.value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "steelblue", high = "lightblue") +
  labs(title    = "Top 15 GO Biological Processes",
       subtitle = "COVID-19 vs Normal",
       x = NULL, y = "-log10(Adjusted P-value)", fill = "Adj. P-value") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))
dev.off()
cat("GO BP bar plot saved!\n")

# -------------------------------------------------------------
# STEP 9: KEGG Bar Plot
# -------------------------------------------------------------
top_kegg <- head(kegg, 15)

png("KEGG_barplot.png", width = 3200, height = 2400, res = 300)
ggplot(top_kegg, aes(x = reorder(Term, -Adjusted.P.value),
                     y = -log10(Adjusted.P.value),
                     fill = Adjusted.P.value)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "darkgreen", high = "lightgreen") +
  labs(title    = "Top 15 KEGG Pathways",
       subtitle = "COVID-19 vs Normal",
       x = NULL, y = "-log10(Adjusted P-value)", fill = "Adj. P-value") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))
dev.off()
cat("KEGG bar plot saved!\n")

# -------------------------------------------------------------
# STEP 10: GO Dot Plot
# -------------------------------------------------------------
top_bp$GeneCount <- sapply(strsplit(top_bp$Genes, ";"), length)

png("GO_dotplot.png", width = 3200, height = 2400, res = 300)
ggplot(top_bp, aes(x = GeneCount,
                   y = reorder(Term, GeneCount),
                   color = Adjusted.P.value,
                   size  = GeneCount)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  labs(title    = "GO Biological Process Dot Plot",
       subtitle = "COVID-19 vs Normal",
       x = "Gene Count", y = NULL,
       color = "Adj. P-value", size = "Gene Count") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))
dev.off()
cat("GO dot plot saved!\n")

# -------------------------------------------------------------
# STEP 11: KEGG Dot Plot
# -------------------------------------------------------------
top_kegg$GeneCount <- sapply(strsplit(top_kegg$Genes, ";"), length)

png("KEGG_dotplot.png", width = 3200, height = 2400, res = 300)
ggplot(top_kegg, aes(x = GeneCount,
                     y = reorder(Term, GeneCount),
                     color = Adjusted.P.value,
                     size  = GeneCount)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  labs(title    = "KEGG Pathway Dot Plot",
       subtitle = "COVID-19 vs Normal",
       x = "Gene Count", y = NULL,
       color = "Adj. P-value", size = "Gene Count") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))
dev.off()
cat("KEGG dot plot saved!\n")

# -------------------------------------------------------------
# STEP 12: Print Top Results to Console
# -------------------------------------------------------------
cat("\n--- TOP 10 GO BIOLOGICAL PROCESSES ---\n")
print(head(go_bp[, c("Term", "Adjusted.P.value", "Overlap")], 10))

cat("\n--- TOP 10 KEGG PATHWAYS ---\n")
print(head(kegg[, c("Term", "Adjusted.P.value", "Overlap")], 10))

# -------------------------------------------------------------
# DONE
# -------------------------------------------------------------
cat("\n========================================\n")
cat("TASK 3 ANALYSIS COMPLETE\n")
cat("========================================\n")
cat("Working directory:", getwd(), "\n\n")
cat("Files saved:\n")
cat("  GO_BP_results.csv\n  GO_MF_results.csv\n  GO_CC_results.csv\n")
cat("  KEGG_results.csv\n  GO_BP_top20.csv\n  KEGG_top20.csv\n")
cat("  GO_BP_barplot.png\n  KEGG_barplot.png\n")
cat("  GO_dotplot.png\n  KEGG_dotplot.png\n")

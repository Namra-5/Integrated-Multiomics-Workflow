# ============================================================
# 8_cleanup.R — Count Matrix Cleaning
# GSE152418: COVID-19 Severe vs Healthy PBMC RNA-seq
# SINES-NUST Submission | Task 1
#
# Purpose:
#   1. Load FeatureCounts_Raw.txt (raw featureCounts output)
#   2. Rename columns from BAM file paths to clean sample names
#   3. Filter out genes with zero counts across all 10 samples
#   4. Save cleaned matrix as FeatureCounts_Mod_clean.txt
#
# Note: coldata creation is NOT included here (Task 2).
#
# Run from terminal:
#   cd ~/rnaseq_project
#   Rscript 8_cleanup.R
#
# Or from within R / RStudio:
#   setwd("~/rnaseq_project")
#   source("8_cleanup.R")
# ============================================================

cat("============================================================\n")
cat(" SINES-NUST | GSE152418 RNA-seq Pipeline\n")
cat(" Step 8: R — Count Matrix Cleanup\n")
cat("============================================================\n\n")

# ------------------------------------------------------------
# SECTION A: Set Working Directory
# Resolves relative paths used throughout this script.
# Adjust if running from a different working directory.
# ------------------------------------------------------------

setwd("~/rnaseq_project")

# ------------------------------------------------------------
# SECTION B: Load Raw featureCounts Output
#
# FeatureCounts_Raw.txt format:
#   Row 1  : Comment line (# Program:featureCounts ...)
#   Row 2  : Column headers (Geneid, Chr, Start, End, Strand, Length, BAM1 ... BAM10)
#   Rows 3+: Gene data
#
# read.table parameters:
#   header    = TRUE  : row 2 is the header
#   row.names = 1     : use column 1 (Geneid) as row names
#   sep       = "\t"  : tab-delimited
#   skip      = 1     : skip the comment line (row 1)
#   check.names= FALSE: preserve exact column names (BAM paths with slashes)
# ------------------------------------------------------------

cat("[A] Loading FeatureCounts_Raw.txt...\n")

RAW_FILE <- "Results/Counts/FeatureCounts_Raw.txt"

if (!file.exists(RAW_FILE)) {
  stop(paste("ERROR: File not found:", RAW_FILE,
             "\n       Run 7_featurecounts.sh first."))
}

counts_raw <- read.table(
  RAW_FILE,
  header     = TRUE,
  row.names  = 1,
  sep        = "\t",
  skip       = 1,          # Skip the '# Program:featureCounts ...' comment line
  check.names = FALSE
)

cat("  Loaded successfully.\n")
cat("  Raw dimensions:", nrow(counts_raw), "genes x", ncol(counts_raw), "columns\n")
cat("  (Expected: ~60,000 genes × 16 columns including annotation)\n\n")

# ------------------------------------------------------------
# SECTION C: Extract Count Columns Only
#
# featureCounts columns layout:
#   Col 1   : Geneid (used as row names above, already removed)
#   Cols 2–5: Chr, Start, End, Strand (annotation — discard)
#   Col 6   : Length (annotation — discard)
#   Cols 7–16: Raw counts for each BAM (10 samples)
#
# We keep only the 10 count columns (columns 6–15 in 0-indexed,
# or columns 6:15 in R after Geneid is set as rownames).
# ------------------------------------------------------------

cat("[B] Extracting 10 count columns (dropping annotation columns 1–5)...\n")

# Columns 1–5 in the loaded dataframe are: Chr, Start, End, Strand, Length
# Columns 6–15 are the 10 count columns (BAM paths as headers)
counts <- counts_raw[ , 6:15]

cat("  Count matrix dimensions:", nrow(counts), "genes x", ncol(counts), "samples\n\n")

# ------------------------------------------------------------
# SECTION D: Rename Columns to Clean Sample Names
#
# Raw column names are full BAM paths, e.g.:
#   Results/Alignments/Normal_1_sorted.bam
#
# Target names (matching featureCounts BAM input order):
#   Normal_1, Normal_2, Normal_3, Normal_4, Normal_5,
#   COVID_1,  COVID_2,  COVID_3,  COVID_4,  COVID_5
#
# The order MUST match the BAM order used in 7_featurecounts.sh
# ------------------------------------------------------------

cat("[C] Renaming columns from BAM paths to clean sample names...\n")

cat("  Current column names:\n")
for (i in seq_along(colnames(counts))) {
  cat(sprintf("    [%2d] %s\n", i, colnames(counts)[i]))
}

# Assign clean names — order mirrors the BAM order in 7_featurecounts.sh
colnames(counts) <- c(
  "Normal_1", "Normal_2", "Normal_3", "Normal_4", "Normal_5",
  "COVID_1",  "COVID_2",  "COVID_3",  "COVID_4",  "COVID_5"
)

cat("\n  New column names:\n")
cat(" ", paste(colnames(counts), collapse = ", "), "\n\n")

# ------------------------------------------------------------
# SECTION E: Filter Zero-Count Genes
#
# Remove genes where the sum across ALL 10 samples equals zero.
# These genes are not expressed in any sample and add noise.
# rowSums(counts) > 0 keeps any gene with at least 1 count
# in at least 1 sample.
# ------------------------------------------------------------

cat("[D] Filtering genes with zero counts across all samples...\n")

n_before <- nrow(counts)
counts   <- counts[rowSums(counts) > 0, ]
n_after  <- nrow(counts)
n_removed <- n_before - n_after

cat(sprintf("  Genes before filter : %d\n", n_before))
cat(sprintf("  Genes removed       : %d  (all-zero rows)\n", n_removed))
cat(sprintf("  Genes after filter  : %d\n\n", n_after))

# ------------------------------------------------------------
# SECTION F: Save Cleaned Matrix
#
# Output: Results/Counts/FeatureCounts_Mod_clean.txt
# Format : Tab-delimited, with Geneid as row names in first column
# write.table parameters:
#   sep       = "\t"  : tab-delimited
#   quote     = FALSE : no quotes around strings
#   col.names = NA    : leave top-left cell blank (standard R behaviour)
# ------------------------------------------------------------

cat("[E] Saving cleaned count matrix...\n")

OUT_FILE <- "Results/Counts/FeatureCounts_Mod_clean.txt"

write.table(
  counts,
  file      = OUT_FILE,
  sep       = "\t",
  quote     = FALSE,
  col.names = NA      # Blank top-left corner, matching standard count matrix format
)

cat(sprintf("  Saved: %s\n\n", OUT_FILE))

# ------------------------------------------------------------
# SECTION G: Final Summary & Matrix Preview
# ------------------------------------------------------------

cat("============================================================\n")
cat(" Cleanup Complete — Summary\n")
cat("============================================================\n\n")

cat(sprintf("  Output file  : %s\n", OUT_FILE))
cat(sprintf("  Dimensions   : %d genes x %d samples\n", nrow(counts), ncol(counts)))
cat(sprintf("  Samples      : %s\n", paste(colnames(counts), collapse = ", ")))
cat("\n")

cat("  Column totals (library sizes — all should be in millions):\n")
col_sums <- colSums(counts)
for (nm in names(col_sums)) {
  cat(sprintf("    %-10s : %s reads\n", nm, format(col_sums[nm], big.mark = ",")))
}

cat("\n  First 5 genes preview:\n")
print(head(counts, 5))

cat("\nNext step: Task 2 — DESeq2 differential expression analysis\n")

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

# Load the peak file
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peak <- readPeakFile("~/Documents/genomics/Project B/ENCFF190KNC.bed.gz")

# Annotate peaks to nearest genes
peakAnno <- annotatePeak(peak, 
                          tssRegion = c(-3000, 3000),
                          TxDb = txdb,
                          annoDb = "org.Hs.eg.db")

# Save annotated results as a table
anno_df <- as.data.frame(peakAnno)
write.csv(anno_df, 
          "~/Documents/genomics/Project B/CTCF_annotated_peaks.csv",
          row.names = FALSE)

# Make the annotation plot 
pdf("~/Documents/genomics/Project B/CTCF_anno_pie.pdf")
plotAnnoPie(peakAnno)
dev.off()

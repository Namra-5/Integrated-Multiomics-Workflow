#!/bin/bash
# ============================================================
# 2_reference_and_index.sh — Reference Genome Download & HISAT2 Index
#
# Reference : GRCh38 (hg38) — Ensembl Release 109
# Genome    : Primary Assembly FASTA (~3.1 GB decompressed)
# Annotation: Ensembl GTF Release 109 (~52 MB decompressed)
# Aligner   : HISAT2 2.2.1 (chosen for 12 GB RAM constraint)
# Time      : Index build takes ~45–90 minutes — do NOT interrupt
# ============================================================

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq_env

cd ~/rnaseq_project

echo "============================================================"
echo " Step 2: Reference Genome Download & HISAT2 Index Build"
echo "============================================================"
echo ""

# ------------------------------------------------------------
# SECTION A: Download Reference FASTA
# Source: Ensembl FTP — Release 109
# File  : Primary assembly (no haplotype contigs)
#         Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz (~900 MB)
# ------------------------------------------------------------

GENOME_DIR="./reference/genome"
ANNOT_DIR="./reference/annotation"
INDEX_DIR="./reference/hisat2_index"

mkdir -p "$GENOME_DIR" "$ANNOT_DIR" "$INDEX_DIR"

GENOME_FA="${GENOME_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GENOME_GZ="${GENOME_FA}.gz"
GTF="${ANNOT_DIR}/Homo_sapiens.GRCh38.109.gtf"
GTF_GZ="${GTF}.gz"

echo "[A] Downloading GRCh38 Primary Assembly FASTA (Ensembl Release 109)..."
echo "    Size: ~900 MB compressed | ~3.1 GB decompressed"
echo ""

if [[ -f "$GENOME_FA" ]]; then
    echo "  [SKIP] FASTA already exists: ${GENOME_FA}"
else
    wget -P "$GENOME_DIR" \
        "https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

    echo "  Decompressing FASTA..."
    gunzip "$GENOME_GZ"
    echo "  Done: ${GENOME_FA}  ($(du -sh "$GENOME_FA" | cut -f1))"
fi

echo ""

# ------------------------------------------------------------
# SECTION B: Download GTF Annotation
# Source: Ensembl FTP — Release 109
# File  : Homo_sapiens.GRCh38.109.gtf.gz (~52 MB compressed)
# NOTE  : Genome and GTF versions must always match (both 109)
# ------------------------------------------------------------

echo "[B] Downloading GRCh38 GTF Annotation (Ensembl Release 109)..."
echo "    Size: ~52 MB compressed"
echo ""

if [[ -f "$GTF" ]]; then
    echo "  [SKIP] GTF already exists: ${GTF}"
else
    wget -P "$ANNOT_DIR" \
        "https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz"

    echo "  Decompressing GTF..."
    gunzip "$GTF_GZ"
    echo "  Done: ${GTF}  ($(du -sh "$GTF" | cut -f1))"
fi

echo ""

# ------------------------------------------------------------
# SECTION C: Build HISAT2 Genome Index
# HISAT2 was selected over STAR due to the 12 GB RAM constraint.
# STAR requires ≥ 32 GB RAM for human genome indexing.
# HISAT2 indexes the same GRCh38 genome in ~8–12 GB RAM.
#
# Index prefix : reference/hisat2_index/GRCh38_109
# Time         : 45–90 minutes depending on CPU/disk speed
# Output files : GRCh38_109.1.ht2 through GRCh38_109.8.ht2
# ------------------------------------------------------------

INDEX_PREFIX="${INDEX_DIR}/GRCh38_109"

echo "[C] Building HISAT2 genome index..."
echo "    This will take 45–90 minutes — do NOT close the terminal."
echo "    RAM usage peak: ~8–12 GB (within 12 GB limit)"
echo ""

# Check if index already built (all 8 .ht2 files present)
if ls "${INDEX_PREFIX}".*.ht2 1>/dev/null 2>&1; then
    echo "  [SKIP] HISAT2 index already exists at: ${INDEX_PREFIX}"
else
    echo "  Running hisat2-build..."
    echo "  Log: ./logs/hisat2_build.log"
    echo ""

    hisat2-build \
        -p 4 \
        "$GENOME_FA" \
        "$INDEX_PREFIX" \
        2>&1 | tee ./logs/hisat2_build.log

    echo ""
    echo "  Index build complete."
fi

# ------------------------------------------------------------
# SECTION D: Verification
# ------------------------------------------------------------

echo ""
echo "============================================================"
echo " Reference Prep Complete — Verification"
echo "============================================================"
echo ""
echo "Genome FASTA:"
ls -lh "${GENOME_FA}"
echo ""
echo "GTF Annotation:"
ls -lh "${GTF}"
echo ""
echo "HISAT2 index files (expected: 8 .ht2 files):"
ls -lh "${INDEX_PREFIX}".*.ht2 2>/dev/null || echo "  ERROR: Index files not found."
echo ""
echo "Index file count: $(ls "${INDEX_PREFIX}".*.ht2 2>/dev/null | wc -l)  (expected: 8)"
echo ""

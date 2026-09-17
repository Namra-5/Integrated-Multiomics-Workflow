#!/bin/bash
# ============================================================
# 5_mapping.sh — HISAT2 Alignment
#
# Aligner   : HISAT2 2.2.1
# Index     : GRCh38 Ensembl Release 109 (built in 2_reference_prep.sh)
# Input     : data/trimmed/SAMPLE_trimmed.fq.gz (single-end)
# Output    : Results/Alignments/SAMPLE.bam (coordinate-sorted BAM)
#
# Key flag:
#   --rna-strandness R : Library = TruSeq Stranded mRNA (reverse)
#                        'R' = single-end reverse-stranded
#                        'RF' would be paired-end equivalent
#
# Pipeline per sample:
#   hisat2  →  samtools view (SAM→BAM)  →  samtools sort  →  SAMPLE.bam
# ============================================================

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq_env

cd ~/rnaseq_project

echo "============================================================"
echo " Step 5: HISAT2 Alignment"
echo "============================================================"
echo ""

# ------------------------------------------------------------
# SECTION A: Path Configuration
# ------------------------------------------------------------

TRIM_DIR="./data/trimmed"
ALIGN_DIR="./Results/Alignments"
LOG_DIR="./logs"
INDEX_PREFIX="./reference/hisat2_index/GRCh38_109"

mkdir -p "$ALIGN_DIR" "$LOG_DIR"

SAMPLES=(Normal_1 Normal_2 Normal_3 Normal_4 Normal_5
         COVID_1  COVID_2  COVID_3  COVID_4  COVID_5)

# ------------------------------------------------------------
# SECTION B: Validate HISAT2 index
# ------------------------------------------------------------

echo "[A] Validating HISAT2 index..."
if ! ls "${INDEX_PREFIX}".*.ht2 1>/dev/null 2>&1; then
    echo "  ERROR: HISAT2 index not found at ${INDEX_PREFIX}"
    echo "         Run 2_reference_prep.sh first."
    exit 1
fi
echo "  Index OK: $(ls "${INDEX_PREFIX}".*.ht2 | wc -l) .ht2 files found"
echo ""

# ------------------------------------------------------------
# SECTION C: Alignment loop
#
# hisat2 flags used:
#   -p 4                  : 4 alignment threads
#   --dta                 : downstream transcriptome assembly mode
#                           (recommended for featureCounts use)
#   --rna-strandness R    : reverse-stranded single-end library
#   -x INDEX_PREFIX       : path to HISAT2 index (without .ht2 ext)
#   -U                    : single-end input FASTQ
#   --no-unal             : suppress unaligned reads from SAM output
#
# samtools flags:
#   view -bS              : convert SAM (from stdin) to BAM
#   sort -@ 4             : coordinate sort with 4 threads
#   -o SAMPLE.bam         : write coordinate-sorted BAM
# ------------------------------------------------------------

echo "[B] Running alignment for all 10 samples..."
echo "    --rna-strandness R (TruSeq Stranded mRNA = reverse)"
echo ""

for SAMPLE in "${SAMPLES[@]}"; do

    INPUT="${TRIM_DIR}/${SAMPLE}_trimmed.fq.gz"
    OUTPUT="${ALIGN_DIR}/${SAMPLE}.bam"
    HISAT2_LOG="${LOG_DIR}/${SAMPLE}.hisat2.log"

    echo "----------------------------------------------------"
    echo " Aligning: ${SAMPLE}"
    echo "----------------------------------------------------"

    if [[ -f "$OUTPUT" ]]; then
        echo "  [SKIP] BAM already exists: ${OUTPUT}"
        echo ""
        continue
    fi

    if [[ ! -f "$INPUT" ]]; then
        echo "  ERROR: Trimmed file not found: ${INPUT}"
        echo "         Run 4_trimming.sh first."
        exit 1
    fi

    # Core alignment pipeline:
    # HISAT2 → SAM → BAM conversion → coordinate sort
    hisat2 \
        -p 4 \
        --dta \
        --rna-strandness R \
        -x "$INDEX_PREFIX" \
        -U "$INPUT" \
        --no-unal \
        2> "$HISAT2_LOG" \
    | samtools view -bS - \
    | samtools sort \
        -@ 4 \
        -o "$OUTPUT"

    echo ""
    echo "  Done: ${OUTPUT}  ($(du -sh "$OUTPUT" | cut -f1))"
    echo ""
    echo "  Alignment summary for ${SAMPLE}:"
    cat "$HISAT2_LOG"
    echo ""

done

# ------------------------------------------------------------
# SECTION D: Alignment Rate Summary
# All samples should achieve ≥ 70% overall alignment rate
# If any sample is below 70%: check trimmed input and index path
# ------------------------------------------------------------

echo "========================================================"
echo " Alignment Complete — Summary (all must be ≥ 70%)"
echo "========================================================"
echo ""

for SAMPLE in "${SAMPLES[@]}"; do
    LOG="${LOG_DIR}/${SAMPLE}.hisat2.log"
    if [[ -f "$LOG" ]]; then
        RATE=$(grep "overall alignment rate" "$LOG" | head -1)
        echo "  ${SAMPLE}: ${RATE}"
    else
        echo "  ${SAMPLE}: log not found"
    fi
done

echo ""
echo "BAM files in ${ALIGN_DIR}/ (expected: 10):"
ls -lh "${ALIGN_DIR}"/*.bam 2>/dev/null || echo "  No BAM files found."
echo ""
echo "Next step: bash 6_sort_and_index.sh"

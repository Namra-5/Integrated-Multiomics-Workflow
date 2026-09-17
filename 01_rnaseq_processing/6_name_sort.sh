#!/bin/bash
# ============================================================
# 6_name_sort.sh — BAM Sorting by Read Name
#
# Input  : Results/Alignments/SAMPLE.bam   (coordinate-sorted)
# Output : Results/Alignments/SAMPLE_sorted.bam  (name-sorted)
#
# Why name-sort?
#   featureCounts for single-end data accepts coordinate-sorted BAMs,
#   but this project spec requires name-sorted output. Name-sorted
#   BAMs are also required for paired-end featureCounts (-p flag),
#   establishing a consistent convention for Task 2 extension.
#
# samtools sort -n : sort reads alphabetically by QNAME (read name)
#                    rather than by chromosomal coordinates (default)
# Output naming    : SAMPLE_sorted.bam  (used in 7_featurecounts.sh)
# ============================================================

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq_env

cd ~/rnaseq_project

echo "============================================================"
echo " Step 6: Samtools — Name-Sort BAM Files"
echo "============================================================"
echo ""

# ------------------------------------------------------------
# SECTION A: Setup
# ------------------------------------------------------------

ALIGN_DIR="./Results/Alignments"
LOG_DIR="./logs"

mkdir -p "$LOG_DIR"

SAMPLES=(Normal_1 Normal_2 Normal_3 Normal_4 Normal_5
         COVID_1  COVID_2  COVID_3  COVID_4  COVID_5)

echo "[A] Sorting BAMs by read name (-n flag)..."
echo "    Input  : SAMPLE.bam"
echo "    Output : SAMPLE_sorted.bam"
echo ""

# ------------------------------------------------------------
# SECTION B: Name-sort each BAM
#
# samtools sort flags:
#   -n        : sort by read name (QNAME) — NOT by coordinate
#   -@ 4      : 4 compression/sort threads
#   -O BAM    : output format BAM
#   -o FILE   : write output to FILE
# ------------------------------------------------------------

for SAMPLE in "${SAMPLES[@]}"; do

    INPUT="${ALIGN_DIR}/${SAMPLE}.bam"
    OUTPUT="${ALIGN_DIR}/${SAMPLE}_sorted.bam"

    echo "----------------------------------------------------"
    echo " Name-sorting: ${SAMPLE}"
    echo "----------------------------------------------------"

    if [[ -f "$OUTPUT" ]]; then
        echo "  [SKIP] Name-sorted BAM already exists: ${OUTPUT}"
        echo ""
        continue
    fi

    if [[ ! -f "$INPUT" ]]; then
        echo "  ERROR: Coordinate-sorted BAM not found: ${INPUT}"
        echo "         Run 5_mapping.sh first."
        exit 1
    fi

    samtools sort \
        -n \
        -@ 4 \
        -O BAM \
        -o "$OUTPUT" \
        "$INPUT"

    echo "  Done: ${OUTPUT}  ($(du -sh "$OUTPUT" | cut -f1))"
    echo ""

done

# ------------------------------------------------------------
# SECTION C: Verification — flagstat on name-sorted BAMs
# Note: samtools index is NOT created here because index files
#       require coordinate-sorted BAMs. Name-sorted BAMs cannot
#       be indexed with samtools index.
# ------------------------------------------------------------

echo "========================================================"
echo " Sort Complete — Flagstat Summary"
echo "========================================================"
echo ""

for SAMPLE in "${SAMPLES[@]}"; do
    SORTED_BAM="${ALIGN_DIR}/${SAMPLE}_sorted.bam"
    if [[ -f "$SORTED_BAM" ]]; then
        MAPPED=$(samtools flagstat "$SORTED_BAM" | grep "mapped (" | head -1)
        echo "  ${SAMPLE}: ${MAPPED}"
    else
        echo "  ${SAMPLE}: MISSING — ${SORTED_BAM}"
    fi
done

echo ""
echo "Name-sorted BAM files (expected: 10):"
ls -lh "${ALIGN_DIR}"/*_sorted.bam 2>/dev/null || echo "  None found."
echo ""
echo "Total: $(ls "${ALIGN_DIR}"/*_sorted.bam 2>/dev/null | wc -l)  (expected: 10)"
echo ""


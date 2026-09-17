#!/bin/bash
# ============================================================
# 7_featurecounts.sh — Gene-Level Quantification
#
# Tool   : featureCounts (Subread 2.0.6)
# Input  : Results/Alignments/SAMPLE_sorted.bam (name-sorted, ×10)
# Output : Results/Counts/FeatureCounts_Raw.txt
#
# CRITICAL flags for this dataset:
#   -s 2   : Reverse-stranded (TruSeq Stranded mRNA = reverse)
#            '-s 1' = forward, '-s 2' = reverse, '-s 0' = unstranded
#   -T 4   : 4 threads
#   NO -p  : This is SINGLE-END data. The -p flag is for paired-end
#            only and will cause incorrect counts if used here.
#
# BAM sample order (preserved in output columns 7–16):
#   Normal_1, Normal_2, Normal_3, Normal_4, Normal_5,
#   COVID_1,  COVID_2,  COVID_3,  COVID_4,  COVID_5
# ============================================================

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq_env

cd ~/rnaseq_project

echo "============================================================"
echo " Step 7: featureCounts — Gene-Level Quantification"
echo "============================================================"
echo ""

# ------------------------------------------------------------
# SECTION A: Path Configuration
# ------------------------------------------------------------

ALIGN_DIR="./Results/Alignments"
COUNT_DIR="./Results/Counts"
GTF="./reference/annotation/Homo_sapiens.GRCh38.109.gtf"
OUTPUT="${COUNT_DIR}/FeatureCounts_Raw.txt"

mkdir -p "$COUNT_DIR"

# ------------------------------------------------------------
# SECTION B: Input Validation
# ------------------------------------------------------------

echo "[A] Validating input BAMs..."
echo ""

SAMPLES=(Normal_1 Normal_2 Normal_3 Normal_4 Normal_5
         COVID_1  COVID_2  COVID_3  COVID_4  COVID_5)

MISSING=0
for SAMPLE in "${SAMPLES[@]}"; do
    BAM="${ALIGN_DIR}/${SAMPLE}_sorted.bam"
    if [[ -f "$BAM" ]]; then
        echo "  OK: ${BAM}  ($(du -sh "$BAM" | cut -f1))"
    else
        echo "  MISSING: ${BAM}"
        MISSING=$((MISSING + 1))
    fi
done

echo ""
if [[ $MISSING -gt 0 ]]; then
    echo "ERROR: ${MISSING} BAM file(s) missing."
    echo "       Run 6_sort_and_index.sh first."
    exit 1
fi

if [[ ! -f "$GTF" ]]; then
    echo "ERROR: GTF not found: ${GTF}"
    echo "       Run 2_reference_prep.sh first."
    exit 1
fi
echo "  GTF OK: ${GTF}"
echo ""

# ------------------------------------------------------------
# SECTION C: Run featureCounts
#
# All 10 BAMs are passed in a single featureCounts call.
# This produces one output file with a count column per sample.
# BAM order determines column order in the output matrix.
#
# Output file: FeatureCounts_Raw.txt
#   Columns 1–6  : GeneID, Chr, Start, End, Strand, Length
#   Columns 7–16 : Raw counts for each sample (full BAM path as header)
# ------------------------------------------------------------

echo "[B] Running featureCounts..."
echo "    -s 2  (reverse-stranded: TruSeq Stranded mRNA)"
echo "    -T 4  (4 threads)"
echo "    SINGLE-END mode (no -p flag)"
echo ""

featureCounts \
    -s 2 \
    -T 4 \
    -a "$GTF" \
    -o "$OUTPUT" \
    "${ALIGN_DIR}/Normal_1_sorted.bam" \
    "${ALIGN_DIR}/Normal_2_sorted.bam" \
    "${ALIGN_DIR}/Normal_3_sorted.bam" \
    "${ALIGN_DIR}/Normal_4_sorted.bam" \
    "${ALIGN_DIR}/Normal_5_sorted.bam" \
    "${ALIGN_DIR}/COVID_1_sorted.bam" \
    "${ALIGN_DIR}/COVID_2_sorted.bam" \
    "${ALIGN_DIR}/COVID_3_sorted.bam" \
    "${ALIGN_DIR}/COVID_4_sorted.bam" \
    "${ALIGN_DIR}/COVID_5_sorted.bam"

# ------------------------------------------------------------
# SECTION D: Verification
# ------------------------------------------------------------

echo ""
echo "========================================================"
echo " featureCounts Complete — Summary"
echo "========================================================"
echo ""

echo "Assignment rate summary (Assigned = reads counted):"
echo "  A low Assigned % with a high Unassigned_Strand % means"
echo "  the -s flag is wrong. For this library, -s 2 is correct."
echo ""
cat "${OUTPUT}.summary"

echo ""
echo "Output file:"
ls -lh "$OUTPUT"
echo ""
echo "Matrix dimensions:"
echo "  Rows: $(wc -l < "$OUTPUT")  (includes 2 header rows; gene rows ≈ 60,000)"
echo "  Cols: $(awk 'NR==2{print NF}' "$OUTPUT")  (expected: 16 = 6 annotation + 10 count)"
echo ""
echo "First 2 rows preview:"
head -2 "$OUTPUT" | cut -f1-3,7-8
echo ""

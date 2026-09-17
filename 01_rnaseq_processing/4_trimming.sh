#!/bin/bash
# ============================================================
# 4_trimming.sh — Adapter & Quality Trimming
#
# Tool    : Trim Galore 0.6.10 (wraps Cutadapt)
# Library : Single-end, TruSeq Stranded mRNA (Illumina Universal Adapter)
# Input   : raw_data/SAMPLE.fastq.gz
# Output  : data/trimmed/SAMPLE_trimmed.fq.gz
#
# Parameters:
#   --quality 20  : Remove bases with Phred score < 20 from 3' end
#   --length  36  : Discard reads shorter than 36 bp after trimming
#   --cores   4   : Parallel Cutadapt threads
#
# Output naming (Trim Galore default for single-end):
#   Normal_1.fastq.gz  →  Normal_1_trimmed.fq.gz
#   COVID_1.fastq.gz   →  COVID_1_trimmed.fq.gz
# ============================================================

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq_env

cd ~/rnaseq_project

echo "============================================================"
echo " Step 4: Trim Galore — Adapter & Quality Trimming"
echo "============================================================"
echo ""

# ------------------------------------------------------------
# SECTION A: Setup
# ------------------------------------------------------------

RAW_DIR="./raw_data"
TRIM_DIR="./data/trimmed"
LOG_DIR="./logs"

mkdir -p "$TRIM_DIR" "$LOG_DIR"

SAMPLES=(Normal_1 Normal_2 Normal_3 Normal_4 Normal_5
         COVID_1  COVID_2  COVID_3  COVID_4  COVID_5)

echo "[A] Trimming parameters:"
echo "    --quality 20  (remove 3' bases with Phred < 20)"
echo "    --length  36  (discard reads < 36 bp after trimming)"
echo "    --cores   4"
echo "    Library  : Single-end (no --paired flag)"
echo "    Adapter  : Auto-detected (Illumina Universal — TruSeq)"
echo ""

# ------------------------------------------------------------
# SECTION B: Trim each sample
# Trim Galore single-end output naming convention:
#   {SAMPLE}.fastq.gz  →  {SAMPLE}_trimmed.fq.gz
# The '_trimmed.fq.gz' suffix is used downstream in 5_mapping.sh
# ------------------------------------------------------------

for SAMPLE in "${SAMPLES[@]}"; do

    INPUT="${RAW_DIR}/${SAMPLE}.fastq.gz"
    EXPECTED_OUT="${TRIM_DIR}/${SAMPLE}_trimmed.fq.gz"

    echo "--------------------------------------------------------"
    echo " Trimming: ${SAMPLE}"
    echo "--------------------------------------------------------"

    # Skip if already trimmed
    if [[ -f "$EXPECTED_OUT" ]]; then
        echo "  [SKIP] Output already exists: ${EXPECTED_OUT}"
        echo ""
        continue
    fi

    # Confirm input exists
    if [[ ! -f "$INPUT" ]]; then
        echo "  ERROR: Input not found: ${INPUT}"
        echo "         Run 1_download_and_merge.sh first."
        exit 1
    fi

    trim_galore \
        --quality 20 \
        --length 36 \
        --cores 4 \
        --output_dir "$TRIM_DIR" \
        "$INPUT" \
        2>&1 | tee "${LOG_DIR}/${SAMPLE}_trimming.log"

    echo "  Done: ${EXPECTED_OUT}  ($(du -sh "$EXPECTED_OUT" | cut -f1))"
    echo ""

done

# ------------------------------------------------------------
# SECTION C: Verification
# ------------------------------------------------------------

echo "========================================================"
echo " Trimming Complete — Verification"
echo "========================================================"
echo ""
echo "Trimmed files in ${TRIM_DIR}/ (expected: 10):"
ls -lh "${TRIM_DIR}"/*_trimmed.fq.gz 2>/dev/null || echo "  No trimmed files found."
echo ""
echo "Total: $(ls "${TRIM_DIR}"/*_trimmed.fq.gz 2>/dev/null | wc -l)  (expected: 10)"
echo ""


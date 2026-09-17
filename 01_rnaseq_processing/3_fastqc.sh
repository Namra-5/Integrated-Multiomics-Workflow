#!/bin/bash
# ============================================================
# 3_fastqc.sh — Raw Read Quality Control
# GSE152418: COVID-19 Severe vs Healthy PBMC RNA-seq
#
# Input  : 10 merged FASTQ files in raw_data/
#          (Normal_1.fastq.gz ... COVID_5.fastq.gz)
# Output : HTML + ZIP QC reports in Results/fastqc/
# Run    : bash 3_fastqc.sh
# After  : Open .html files in a browser to inspect quality
# ============================================================

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq_env

cd ~/rnaseq_project

echo "============================================================"
echo " Step 3: FastQC — Raw Read Quality Control"
echo "============================================================"
echo ""

# ------------------------------------------------------------
# SECTION A: Input Validation
# Confirm all 10 merged FASTQ files are present before starting
# ------------------------------------------------------------

RAW_DIR="./raw_data"
QC_DIR="./Results/fastqc"
LOG_DIR="./logs"

mkdir -p "$QC_DIR" "$LOG_DIR"

SAMPLES=(Normal_1 Normal_2 Normal_3 Normal_4 Normal_5
         COVID_1  COVID_2  COVID_3  COVID_4  COVID_5)

echo "[A] Checking input files..."
echo ""

MISSING=0
for SAMPLE in "${SAMPLES[@]}"; do
    FILE="${RAW_DIR}/${SAMPLE}.fastq.gz"
    if [[ -f "$FILE" ]]; then
        echo "  OK: ${FILE}  ($(du -sh "$FILE" | cut -f1))"
    else
        echo "  MISSING: ${FILE}"
        MISSING=$((MISSING + 1))
    fi
done

echo ""
if [[ $MISSING -gt 0 ]]; then
    echo "ERROR: ${MISSING} file(s) missing. Run 1_download_and_merge.sh first."
    exit 1
fi

echo "  All 10 input files confirmed."
echo ""

# ------------------------------------------------------------
# SECTION B: Run FastQC
# --outdir : write HTML and ZIP reports to Results/fastqc/
# --threads : 8 parallel threads (one thread per file)
# Input    : all 10 .fastq.gz files passed at once
# Logs     : timestamped stdout/stderr in logs/
# ------------------------------------------------------------

TIMESTAMP=$(date "+%Y%m%d_%H%M%S")

echo "[B] Running FastQC on all 10 merged raw FASTQ files..."
echo "    Output → ${QC_DIR}/"
echo "    Log    → ${LOG_DIR}/fastqc_${TIMESTAMP}.log"
echo ""

fastqc \
    --outdir "$QC_DIR" \
    --threads 8 \
    "${RAW_DIR}"/Normal_1.fastq.gz \
    "${RAW_DIR}"/Normal_2.fastq.gz \
    "${RAW_DIR}"/Normal_3.fastq.gz \
    "${RAW_DIR}"/Normal_4.fastq.gz \
    "${RAW_DIR}"/Normal_5.fastq.gz \
    "${RAW_DIR}"/COVID_1.fastq.gz \
    "${RAW_DIR}"/COVID_2.fastq.gz \
    "${RAW_DIR}"/COVID_3.fastq.gz \
    "${RAW_DIR}"/COVID_4.fastq.gz \
    "${RAW_DIR}"/COVID_5.fastq.gz \
    2>&1 | tee "${LOG_DIR}/fastqc_${TIMESTAMP}.log"

# ------------------------------------------------------------
# SECTION C: Output Summary
# ------------------------------------------------------------

echo ""
echo "============================================================"
echo " FastQC Complete — Output Files"
echo "============================================================"
echo ""
echo "HTML reports (open in browser):"
ls "${QC_DIR}"/*.html 2>/dev/null || echo "  No HTML files found."
echo ""
echo "Report count: $(ls "${QC_DIR}"/*.html 2>/dev/null | wc -l)  (expected: 10)"
echo ""


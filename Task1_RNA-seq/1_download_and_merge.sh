#!/bin/bash
# ============================================================
# 1_download_and_merge.sh — Data Download & Lane Merging
#
# Dataset   : GSE152418 | BioProject: PRJNA639275
# Samples   : 10 biological samples × 2 SRR lane splits = 20 runs
# IMPORTANT : These are NOT paired-end reads. Each sample was
#             sequenced across 2 NovaSeq flow-cell lanes.
#             The two SRR files are merged with 'cat' BEFORE FastQC.
#             Do NOT use --split-files with fasterq-dump.
# Output    : 10 merged .fastq.gz files in raw_data/
# ============================================================

set -euo pipefail

source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq_env

cd ~/rnaseq_project
mkdir -p raw_data

echo "============================================================"
echo " Step 1: Download SRR Runs & Merge Lane Splits"
echo "============================================================"
echo ""

# ------------------------------------------------------------
# SECTION A: Sample-to-SRR Mapping
# Source: SRA metadata, verified from PRJNA639275
#
# | Label    | Lane 1 SRR   | Lane 2 SRR   | Condition      |
# |----------|--------------|--------------|----------------|
# | Normal_1 | SRR12007855  | SRR12007856  | Healthy        |
# | Normal_2 | SRR12007857  | SRR12007858  | Healthy        |
# | Normal_3 | SRR12007859  | SRR12007860  | Healthy        |
# | Normal_4 | SRR12007861  | SRR12007862  | Healthy        |
# | Normal_5 | SRR12007863  | SRR12007864  | Healthy        |
# | COVID_1  | SRR12007825  | SRR12007826  | COVID-19 Sev.  |
# | COVID_2  | SRR12007827  | SRR12007828  | COVID-19 Sev.  |
# | COVID_3  | SRR12007831  | SRR12007832  | COVID-19 Sev.  |
# | COVID_4  | SRR12007835  | SRR12007836  | COVID-19 Sev.  |
# | COVID_5  | SRR12007839  | SRR12007840  | COVID-19 Sev.  |
# ------------------------------------------------------------

# Associative arrays: LANE1[sample]=SRR_run1, LANE2[sample]=SRR_run2
declare -A LANE1
declare -A LANE2

# Healthy controls
LANE1[Normal_1]=SRR12007855 ; LANE2[Normal_1]=SRR12007856
LANE1[Normal_2]=SRR12007857 ; LANE2[Normal_2]=SRR12007858
LANE1[Normal_3]=SRR12007859 ; LANE2[Normal_3]=SRR12007860
LANE1[Normal_4]=SRR12007861 ; LANE2[Normal_4]=SRR12007862
LANE1[Normal_5]=SRR12007863 ; LANE2[Normal_5]=SRR12007864

# COVID-19 Severe cases
LANE1[COVID_1]=SRR12007825  ; LANE2[COVID_1]=SRR12007826
LANE1[COVID_2]=SRR12007827  ; LANE2[COVID_2]=SRR12007828
LANE1[COVID_3]=SRR12007831  ; LANE2[COVID_3]=SRR12007832
LANE1[COVID_4]=SRR12007835  ; LANE2[COVID_4]=SRR12007836
LANE1[COVID_5]=SRR12007839  ; LANE2[COVID_5]=SRR12007840

# Ordered sample list (Normal first, then COVID)
SAMPLES=(Normal_1 Normal_2 Normal_3 Normal_4 Normal_5
         COVID_1  COVID_2  COVID_3  COVID_4  COVID_5)

echo "[INFO] Disk space check (need ≥ 25 GB free):"
df -h ~/rnaseq_project
echo ""

# ------------------------------------------------------------
# SECTION B: Download, Convert & Merge per Sample
# For each biological sample:
#   1. prefetch  → downloads .sra cache file
#   2. fasterq-dump → converts .sra to .fastq (NO --split-files)
#   3. gzip → compress .fastq to .fastq.gz
#   4. cat  → merge the two lane files into one
#   5. Cleanup per-lane intermediates to save disk space
# ------------------------------------------------------------

for NAME in "${SAMPLES[@]}"; do

    SRR1=${LANE1[$NAME]}
    SRR2=${LANE2[$NAME]}
    MERGED="./raw_data/${NAME}.fastq.gz"

    echo "========================================================"
    echo " Processing: ${NAME}  (${SRR1} + ${SRR2})"
    echo "========================================================"

    # Skip if merged file already exists (safe to re-run after crash)
    if [[ -f "$MERGED" ]]; then
        echo "  [SKIP] ${MERGED} already exists. Skipping..."
        echo ""
        continue
    fi

    # --- Lane 1 ---
    echo "  [1/4] Downloading lane 1: ${SRR1}..."
    prefetch "${SRR1}" --output-directory ./raw_data/

    echo "  [2/4] Converting ${SRR1}.sra to FASTQ..."
    fasterq-dump "./raw_data/${SRR1}/${SRR1}.sra" \
        --outdir ./raw_data/ \
        --threads 8 \
        --progress
    gzip "./raw_data/${SRR1}.fastq"

    # --- Lane 2 ---
    echo "  [1/4] Downloading lane 2: ${SRR2}..."
    prefetch "${SRR2}" --output-directory ./raw_data/

    echo "  [2/4] Converting ${SRR2}.sra to FASTQ..."
    fasterq-dump "./raw_data/${SRR2}/${SRR2}.sra" \
        --outdir ./raw_data/ \
        --threads 8 \
        --progress
    gzip "./raw_data/${SRR2}.fastq"

    # --- Merge lanes ---
    echo "  [3/4] Merging ${SRR1}.fastq.gz + ${SRR2}.fastq.gz → ${NAME}.fastq.gz..."
    cat "./raw_data/${SRR1}.fastq.gz" \
        "./raw_data/${SRR2}.fastq.gz" \
        > "$MERGED"

    # --- Cleanup intermediates ---
    echo "  [4/4] Cleaning up per-lane files to save disk space..."
    rm -f "./raw_data/${SRR1}.fastq.gz" "./raw_data/${SRR2}.fastq.gz"
    rm -rf "./raw_data/${SRR1}/" "./raw_data/${SRR2}/"

    echo "  Done: ${NAME}.fastq.gz  ($(du -sh "$MERGED" | cut -f1))"
    echo ""

done

# ------------------------------------------------------------
# SECTION C: Final Verification
# ------------------------------------------------------------

echo "========================================================"
echo " Download & Merge Complete — Verification"
echo "========================================================"
echo ""
echo "Files in raw_data/ (expected: 10):"
ls -lh ./raw_data/*.fastq.gz
echo ""
echo "Total files: $(ls ./raw_data/*.fastq.gz | wc -l)  (expected: 10)"
echo ""

#!/bin/bash
# ============================================================
# 0_setup.sh — Project Setup
#
# Purpose : Create directory structure and install all
#           required bioinformatics tools via conda.
# Run once: bash 0_setup.sh
# ============================================================

set -euo pipefail

echo "============================================================"
echo " Step 0: Environment & Directory Setup"
echo "============================================================"
echo ""

# ------------------------------------------------------------
# SECTION A: Directory Structure
# Creates the full project tree under ~/rnaseq_project/
# ------------------------------------------------------------

echo "[A] Creating project directory structure..."

mkdir -p ~/rnaseq_project/{raw_data,data/trimmed,reference/{genome,annotation},Results/{Alignments,Counts,fastqc},logs}

echo ""
echo "    Directory tree established at ~/rnaseq_project/:"
ls ~/rnaseq_project/
echo ""
echo "    Subdirectories:"
find ~/rnaseq_project -type d | sort
echo ""

# ------------------------------------------------------------
# SECTION B: Conda Environment Creation
# Environment name : rnaseq_env
# Channels         : bioconda, conda-forge, defaults
# Tools installed  : samtools, subread (featureCounts),
#                    fastqc, trim-galore, hisat2, sra-tools
# ------------------------------------------------------------

echo "[B] Creating conda environment: rnaseq_env"
echo "    This may take 5-15 minutes on first run..."
echo ""

# Create environment with Python 3.10 as base
conda create -n rnaseq_env -y python=3.10

echo ""
echo "[C] Installing bioinformatics tools into rnaseq_env..."
echo "    Tools: samtools, subread, fastqc, trim-galore, hisat2, sra-tools"
echo ""

conda install -n rnaseq_env -c bioconda -c conda-forge -y \
    samtools=1.18 \
    subread=2.0.6 \
    fastqc=0.12.1 \
    trim-galore=0.6.10 \
    hisat2=2.2.1 \
    sra-tools=3.1.1

# ------------------------------------------------------------
# SECTION D: Tool Verification
# ------------------------------------------------------------

echo ""
echo "[D] Verifying installations..."
echo ""

source ~/miniconda3/etc/profile.d/conda.sh
conda activate rnaseq_env

echo -n "    samtools    : "; samtools --version 2>/dev/null | head -1
echo -n "    featureCounts: "; featureCounts -v 2>&1 | head -1
echo -n "    fastqc      : "; fastqc --version 2>/dev/null
echo -n "    trim_galore : "; trim_galore --version 2>/dev/null | grep version | head -1
echo -n "    hisat2      : "; hisat2 --version 2>/dev/null | head -1
echo -n "    fasterq-dump: "; fasterq-dump --version 2>/dev/null | head -1

echo ""
echo "============================================================"
echo " Setup complete."
echo " Activate your environment before running any script:"
echo "   conda activate rnaseq_env"
echo "============================================================"

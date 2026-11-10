#!/bin/bash
# =============================================================================
# 00_setup.sh
# Install all tools for the ARV haplotype phasing / read-pair linkage project
#
# Tools installed:
#   - SRA Toolkit        (download data)
#   - Entrez Direct      (NCBI API queries)
#   - fastp              (read trimming)
#   - BWA-MEM2           (read mapping)
#   - SPAdes             (de novo assembly of divergent strain)
#   - BLAST+             (classify assembled contigs by strain)
#   - samtools           (BAM processing)
#   - bcftools           (variant calling)
#   - lofreq             (low-frequency variant calling)
#   - HapCUT2            (haplotype phasing from reads)
#   - Python packages    (pysam, pandas, matplotlib, networkx, biopython)
# =============================================================================

set -euo pipefail

echo "========================================"
echo " ARV Haplotype Project — Setup"
echo "========================================"

if ! command -v conda &> /dev/null; then
    echo "[ERROR] conda not found."
    echo "Install from: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

ENV_NAME="arv_haplo"

# Create environment
if conda env list | grep -q "^${ENV_NAME}"; then
    echo "[INFO] Environment '${ENV_NAME}' already exists."
else
    echo "[INFO] Creating conda environment: ${ENV_NAME}"
    conda create -y -n "${ENV_NAME}" python=3.10
fi

echo "[INFO] Installing bioinformatics tools..."
conda install -y -n "${ENV_NAME}" \
    -c bioconda -c conda-forge \
    sra-tools \
    entrez-direct \
    fastp \
    bwa-mem2 \
    spades \
    blast \
    samtools \
    bcftools \
    lofreq \
    hapcut2 \
    bedtools \
    seqkit

echo "[INFO] Installing Python packages..."
conda run -n "${ENV_NAME}" pip install \
    pysam \
    pandas \
    numpy \
    matplotlib \
    seaborn \
    scipy \
    networkx \
    biopython

# Verify
echo ""
echo "[INFO] Verifying tool installations..."
TOOLS=(prefetch fasterq-dump fastp bwa-mem2 spades.py makeblastdb blastn samtools bcftools lofreq hapcut2 esearch efetch)
for tool in "${TOOLS[@]}"; do
    if conda run -n "${ENV_NAME}" which "${tool}" &>/dev/null; then
        echo "  [OK] ${tool}"
    else
        echo "  [MISSING] ${tool}"
    fi
done

echo ""
echo "[SUCCESS] Setup complete."
echo "Activate with: conda activate ${ENV_NAME}"

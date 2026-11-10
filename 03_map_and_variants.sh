#!/bin/bash
# =============================================================================
# 03_map_and_variants.sh
# Map ALL 6 samples to the dual reference (S1133 + divergent combined)
# Then call variants with lofreq (designed for low-frequency viral variants)
#
# WHY LOFREQ:
#   Standard variant callers assume diploid organisms and filter out
#   variants <20% frequency. ARV is a mixed viral population — we need
#   to detect variants at 1–99% frequency. LoFreq is purpose-built for this.
#
# OUTPUTS per sample:
#   data/bam/{SAMPLE}_dual.bam        — reads mapped to dual reference
#   data/bam/{SAMPLE}_dual.bam.bai    — index
#   data/stats/{SAMPLE}_dual_flagstat.txt
#   data/stats/{SAMPLE}_dual_idxstats.txt  — per-segment read counts (KEY!)
#   data/variants/{SAMPLE}.vcf        — variants called by lofreq
# =============================================================================

set -euo pipefail
cd "${PBS_O_WORKDIR:-$PWD}"
source /apps/x86-64/apps/anaconda_3-2025.12/etc/profile.d/conda.sh
conda activate arv_haplo

echo "========================================"
echo " Step 3: Map to Dual Reference + Call Variants"
echo "========================================"

DUAL_REF="data/reference/dual_reference.fasta"
# Raw reads already downloaded by week12/01_download_sra.sh — reuse them
RAW_DIR="/home/aubnxn001/week12/data/raw"
TRIM_DIR="data/trimmed"
BAM_DIR="data/bam"
STATS_DIR="data/stats"
VCF_DIR="data/variants"
THREADS=4

mkdir -p "${BAM_DIR}" "${STATS_DIR}" "${VCF_DIR}" "${TRIM_DIR}"

SAMPLES=("AL" "S1133" "P1" "PP3" "PPP" "EP")

for SAMPLE in "${SAMPLES[@]}"; do
    R1="${RAW_DIR}/${SAMPLE}/${SAMPLE}_R1.fastq.gz"
    R2="${RAW_DIR}/${SAMPLE}/${SAMPLE}_R2.fastq.gz"

    if [ ! -f "${R1}" ]; then
        echo "[SKIP] ${SAMPLE}: raw reads not found (${R1})"
        continue
    fi

    echo ""
    echo "-------- Processing ${SAMPLE} --------"

    # ---- Trim (skip if already done) ----
    TRIM_R1="${TRIM_DIR}/${SAMPLE}/${SAMPLE}_R1.fastq.gz"
    if [ ! -f "${TRIM_R1}" ]; then
        mkdir -p "${TRIM_DIR}/${SAMPLE}"
        fastp \
            -i "${R1}" -I "${R2}" \
            -o "${TRIM_DIR}/${SAMPLE}/${SAMPLE}_R1.fastq.gz" \
            -O "${TRIM_DIR}/${SAMPLE}/${SAMPLE}_R2.fastq.gz" \
            --json "${STATS_DIR}/${SAMPLE}_fastp.json" \
            --thread "${THREADS}" \
            --detect_adapter_for_pe \
            --qualified_quality_phred 20 \
            --length_required 50 \
            2>&1 | grep -E "reads|passed"
    else
        echo "[INFO] ${SAMPLE}: trimmed reads already exist, skipping fastp"
    fi

    # ---- Map to dual reference ----
    echo "[INFO] ${SAMPLE}: Mapping to dual reference..."
    bwa-mem2 mem \
        -t "${THREADS}" \
        -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
        "${DUAL_REF}" \
        "${TRIM_DIR}/${SAMPLE}/${SAMPLE}_R1.fastq.gz" \
        "${TRIM_DIR}/${SAMPLE}/${SAMPLE}_R2.fastq.gz" \
    | samtools sort -@ "${THREADS}" -o "${BAM_DIR}/${SAMPLE}_dual.bam"

    samtools index "${BAM_DIR}/${SAMPLE}_dual.bam"

    # ---- Collect statistics ----
    echo "[INFO] ${SAMPLE}: Collecting stats..."

    # Overall mapping rate
    samtools flagstat "${BAM_DIR}/${SAMPLE}_dual.bam" \
        > "${STATS_DIR}/${SAMPLE}_dual_flagstat.txt"

    # Per-segment read counts — THE KEY DATA for haplotype analysis
    # Shows how many reads mapped to S1133_* segments vs DIV_* contigs
    samtools idxstats "${BAM_DIR}/${SAMPLE}_dual.bam" \
        > "${STATS_DIR}/${SAMPLE}_dual_idxstats.txt"

    echo "[INFO] ${SAMPLE}: Read distribution across S1133 vs DIV segments:"
    echo "  S1133 segments:"
    grep "^S1133" "${STATS_DIR}/${SAMPLE}_dual_idxstats.txt" | \
        awk '{printf "    %-30s %d reads\n", $1, $3}'
    echo "  Divergent contigs:"
    grep "^DIV" "${STATS_DIR}/${SAMPLE}_dual_idxstats.txt" | \
        awk '{printf "    %-30s %d reads\n", $1, $3}' | head -15

    # ---- Variant calling with LoFreq ----
    echo "[INFO] ${SAMPLE}: Calling variants with lofreq..."

    # lofreq requires base quality scores to be recalibrated
    # First, fill in missing base qualities for indel calling
    lofreq indelqual \
        --dindel \
        -f "${DUAL_REF}" \
        -o "${BAM_DIR}/${SAMPLE}_iq.bam" \
        "${BAM_DIR}/${SAMPLE}_dual.bam"

    samtools index "${BAM_DIR}/${SAMPLE}_iq.bam"

    # Call variants (SNPs + indels)
    # --call-indels: also call insertions/deletions
    # --sig 0.01: significance threshold
    lofreq call-parallel \
        --pp-threads "${THREADS}" \
        --call-indels \
        -f "${DUAL_REF}" \
        -o "${VCF_DIR}/${SAMPLE}.vcf" \
        "${BAM_DIR}/${SAMPLE}_iq.bam" \
        2>&1 | grep -v "^#" | head -5 || true

    # Count variants
    SNP_COUNT=$(grep -v "^#" "${VCF_DIR}/${SAMPLE}.vcf" 2>/dev/null | wc -l || echo 0)
    echo "[INFO] ${SAMPLE}: ${SNP_COUNT} variants called"
    echo "[DONE] ${SAMPLE}"
done

echo ""
echo "========================================"
echo " Mapping and variant calling complete."
echo ""
echo " KEY FILE — per-segment read counts:"
echo "   data/stats/{SAMPLE}_dual_idxstats.txt"
echo ""
echo " KEY FILE — variants:"
echo "   data/variants/{SAMPLE}.vcf"
echo "========================================"

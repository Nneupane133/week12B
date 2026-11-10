#!/bin/bash
# =============================================================================
# run_all.sh
# Master script — runs the complete ARV haplotype / reassortment pipeline
#
# QUESTION: Is the AL field isolate ONE reassortant virus or TWO co-infecting viruses?
#
# Steps:
#   00 — Set up conda environment
#   01 — Download SRA data (6 samples from PRJNA1217582)
#   02 — Build dual reference (S1133 + divergent field strain)
#   03 — Map all samples to dual reference; call variants with LoFreq
#   04 — Read-pair linkage analysis (cross-strain pairs = key signal)
#   05 — Haplotype phasing with HapCUT2
#   06 — Visualize all results
#
# Runtime: ~60–90 min (download-dependent), ~8 GB RAM, ~4 CPU cores
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${SCRIPT_DIR}/logs"
mkdir -p "${LOG_DIR}"

START_TIME=$(date +%s)

echo ""
echo "╔══════════════════════════════════════════════════════════╗"
echo "║   ARV Haplotype / Reassortment Analysis Pipeline        ║"
echo "║   BioProject PRJNA1217582                               ║"
echo "╚══════════════════════════════════════════════════════════╝"
echo ""
echo "Working directory: ${SCRIPT_DIR}"
echo "Started:           $(date)"
echo ""

# ──────────────────────────────────────────────────────────────
# Helper: run step with logging and timing
# ──────────────────────────────────────────────────────────────
run_step() {
    local STEP_NAME="$1"
    local STEP_CMD="$2"
    local LOG="${LOG_DIR}/${STEP_NAME}.log"

    echo "▶ ${STEP_NAME}"
    echo "  Log: ${LOG}"

    STEP_START=$(date +%s)

    if eval "${STEP_CMD}" > "${LOG}" 2>&1; then
        STEP_END=$(date +%s)
        ELAPSED=$(( STEP_END - STEP_START ))
        echo "  ✓ Done in ${ELAPSED}s"
    else
        echo "  ✗ FAILED — check ${LOG}"
        echo ""
        echo "Last 20 lines of log:"
        tail -20 "${LOG}"
        exit 1
    fi
    echo ""
}

# ──────────────────────────────────────────────────────────────
# Parse flags
# ──────────────────────────────────────────────────────────────
SKIP_SETUP=0
SKIP_DOWNLOAD=0
START_FROM=0

while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-setup)     SKIP_SETUP=1;    shift ;;
        --skip-download)  SKIP_DOWNLOAD=1; shift ;;
        --from)           START_FROM="$2"; shift 2 ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --skip-setup      Skip conda env setup (env already exists)"
            echo "  --skip-download   Skip SRA download (data already downloaded)"
            echo "  --from N          Start from step N (0–6)"
            echo ""
            echo "Examples:"
            echo "  $0                       # Full pipeline from scratch"
            echo "  $0 --skip-setup --from 2 # Rebuild reference and run analysis"
            echo "  $0 --skip-setup --from 4 # Re-run analysis only"
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ──────────────────────────────────────────────────────────────
# Step 00 — Environment setup
# ──────────────────────────────────────────────────────────────
if [[ "${SKIP_SETUP}" -eq 0 ]] && [[ "${START_FROM}" -le 0 ]]; then
    run_step "00_setup" "bash ${SCRIPT_DIR}/00_setup.sh"
else
    echo "— Skipping step 00 (setup)"
fi

# ──────────────────────────────────────────────────────────────
# Step 01 — Download SRA data
# ──────────────────────────────────────────────────────────────
if [[ "${SKIP_DOWNLOAD}" -eq 0 ]] && [[ "${START_FROM}" -le 1 ]]; then
    run_step "01_download" "bash ${SCRIPT_DIR}/01_download_sra.sh"
else
    echo "— Skipping step 01 (download)"
fi

# ──────────────────────────────────────────────────────────────
# Step 02 — Build dual reference
# ──────────────────────────────────────────────────────────────
if [[ "${START_FROM}" -le 2 ]]; then
    run_step "02_dual_ref" "bash ${SCRIPT_DIR}/02_build_dual_reference.sh"
else
    echo "— Skipping step 02 (dual reference)"
fi

# ──────────────────────────────────────────────────────────────
# Step 03 — Map and call variants
# ──────────────────────────────────────────────────────────────
if [[ "${START_FROM}" -le 3 ]]; then
    run_step "03_map_variants" "bash ${SCRIPT_DIR}/03_map_and_variants.sh"
else
    echo "— Skipping step 03 (mapping + variants)"
fi

# ──────────────────────────────────────────────────────────────
# Step 04 — Read-pair linkage analysis (Python)
# ──────────────────────────────────────────────────────────────
if [[ "${START_FROM}" -le 4 ]]; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate arv_haplo
    run_step "04_linkage" "python3 ${SCRIPT_DIR}/04_read_pair_linkage.py"
else
    echo "— Skipping step 04 (linkage analysis)"
fi

# ──────────────────────────────────────────────────────────────
# Step 05 — Haplotype phasing
# ──────────────────────────────────────────────────────────────
if [[ "${START_FROM}" -le 5 ]]; then
    run_step "05_phasing" "bash ${SCRIPT_DIR}/05_phase_haplotypes.sh"
else
    echo "— Skipping step 05 (haplotype phasing)"
fi

# ──────────────────────────────────────────────────────────────
# Step 06 — Visualize
# ──────────────────────────────────────────────────────────────
if [[ "${START_FROM}" -le 6 ]]; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate arv_haplo
    run_step "06_visualize" "python3 ${SCRIPT_DIR}/06_visualize.py"
else
    echo "— Skipping step 06 (visualization)"
fi

# ──────────────────────────────────────────────────────────────
# Summary
# ──────────────────────────────────────────────────────────────
END_TIME=$(date +%s)
TOTAL=$(( END_TIME - START_TIME ))
MINS=$(( TOTAL / 60 ))
SECS=$(( TOTAL % 60 ))

echo "╔══════════════════════════════════════════════════════════╗"
echo "║   Pipeline Complete                                     ║"
echo "╚══════════════════════════════════════════════════════════╝"
echo ""
echo "Total runtime: ${MINS}m ${SECS}s"
echo "Finished:      $(date)"
echo ""
echo "KEY RESULTS:"
echo "  Cross-strain pair report: results/linkage_report.txt"
echo "  Cross-strain pair list:   results/cross_strain_pairs.tsv"
echo "  Haplotype blocks:         data/phasing/AL_haplotypes"
echo ""
echo "FIGURES:"
echo "  results/figures/fig1_cross_strain_rates.png"
echo "  results/figures/fig2_segment_genome_map.png"
echo "  results/figures/fig3_haplotype_blocks.png"
echo "  results/figures/fig4_summary_poster.png"
echo ""
echo "INTERPRETATION GUIDE: See README.md"
echo ""
echo "QUICK VERDICT:"
if [ -f "results/linkage_report.txt" ]; then
    echo "  AL isolate cross-strain pair rate:"
    grep -A2 "  AL (" results/linkage_report.txt | grep "CROSS-STRAIN" || \
        echo "  (run step 04 to generate)"
fi

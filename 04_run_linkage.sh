#!/bin/bash
# =============================================================================
# 04_run_linkage.sh
# PBS wrapper: run the read-pair linkage analysis (04_read_pair_linkage.py)
# Reads:   data/bam/*_dual.bam  (from 03_map_and_variants.sh)
# Writes:  results/linkage_report.txt, results/cross_strain_pairs.tsv,
#          results/figures/
# =============================================================================

set -euo pipefail
cd "${PBS_O_WORKDIR:-$PWD}"
source /apps/x86-64/apps/anaconda_3-2025.12/etc/profile.d/conda.sh
conda activate arv_haplo

echo "========================================"
echo " Step 4: Read-Pair Linkage Analysis"
echo "========================================"

python 04_read_pair_linkage.py

echo ""
echo "[INFO] Outputs:"
ls -l results/ results/figures/

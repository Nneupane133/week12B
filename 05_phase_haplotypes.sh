#!/bin/bash
# =============================================================================
# 05_phase_haplotypes.sh
# Phase variants within each ARV segment using HapCUT2
#
# WHY HAPLOTYPE PHASING?
#   The read-pair linkage analysis (script 04) looks for cross-strain signals
#   at the whole-segment level. Phasing goes one level deeper: within a single
#   segment, it reconstructs which set of SNPs travel together on the same
#   molecule — i.e., which variants define each of the two haplotypes.
#
#   If two viruses co-exist:
#     → Each segment should phase into two clean haplotypes:
#       Haplotype 1 = S1133-like SNPs
#       Haplotype 2 = Divergent SNPs
#
#   If only one virus (highly diverse or reassortant):
#     → Haplotypes may be blurry or recombinant within a segment
#
# TOOL: HapCUT2
#   - Designed for phasing heterozygous variants from short reads
#   - Works perfectly for a mixed viral population (treat the two strains
#     as "two haplotypes" on the same reference)
#   - Input: BAM + VCF → Output: phased haplotype blocks
#
# We run on the AL isolate mapped to the DUAL reference, focusing on
# the segments where both S1133 and DIV reads were observed.
# =============================================================================

set -euo pipefail
cd "${PBS_O_WORKDIR:-$PWD}"
source /apps/x86-64/apps/anaconda_3-2025.12/etc/profile.d/conda.sh
conda activate arv_haplo

echo "========================================"
echo " Step 5: Haplotype Phasing with HapCUT2"
echo "========================================"

DUAL_REF="data/reference/dual_reference.fasta"
BAM_DIR="data/bam"
VCF_DIR="data/variants"
PHASE_DIR="data/phasing"
mkdir -p "${PHASE_DIR}"

# Focus on the AL isolate — this is the sample we're investigating
SAMPLE="AL"
BAM="${BAM_DIR}/${SAMPLE}_iq.bam"
VCF="${VCF_DIR}/${SAMPLE}.vcf"

if [ ! -f "${BAM}" ]; then
    echo "[ERROR] BAM not found: ${BAM}"
    echo "Run 03_map_and_variants.sh first."
    exit 1
fi

echo "[INFO] Processing sample: ${SAMPLE}"
echo "[INFO] BAM:  ${BAM}"
echo "[INFO] VCF:  ${VCF}"

# ============================================================
# STEP 1 — Filter VCF: keep only high-quality biallelic SNPs
#           in regions with both S1133 and DIV reads (the
#           informative segments for phasing)
# ============================================================
echo ""
echo "[STEP 1] Filtering VCF for phasing-informative variants..."

# LoFreq VCFs have no FORMAT/GT column (only INFO/AF), but HapCUT2 phases
# heterozygous (0/1) genotypes. So we build a single-sample VCF ourselves:
#   - PASS, biallelic SNVs only (no indels; HapCUT2 works best with SNPs)
#   - mixed sites only: MIN_AF <= AF <= MAX_AF  → GT = 0/1
#     (sites near 0% or 100% have nothing to phase)
#   - contig lines added from the reference .fai (LoFreq omits them)
MIN_AF=0.05
MAX_AF=0.95

python3 - "${VCF}" "${DUAL_REF}.fai" "${PHASE_DIR}/${SAMPLE}_snps.vcf" \
    "${SAMPLE}" "${MIN_AF}" "${MAX_AF}" <<'PYEOF'
import sys
vcf_in, fai, vcf_out, sample, min_af, max_af = sys.argv[1:]
min_af, max_af = float(min_af), float(max_af)

meta, records, seen = [], [], set()
with open(vcf_in) as f:
    for line in f:
        if line.startswith("##"):
            if not line.startswith(("##contig", "##fileformat")):
                meta.append(line)
            continue
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        chrom, pos, _id, ref, alt, qual, filt, info = c[:8]
        if filt != "PASS" or "INDEL" in info.split(";"):
            continue
        if len(ref) != 1 or len(alt) != 1 or "," in alt:
            continue
        af = next((float(x[3:]) for x in info.split(";") if x.startswith("AF=")), None)
        if af is None or not (min_af <= af <= max_af):
            continue
        if (chrom, pos) in seen:          # multiple ALTs at one site → not biallelic
            records = [r for r in records if (r[0], r[1]) != (chrom, pos)]
            continue
        seen.add((chrom, pos))
        records.append([chrom, pos, _id, ref, alt, qual, filt, info, "GT", "0/1"])

with open(vcf_out, "w") as out:
    out.write("##fileformat=VCFv4.2\n")
    with open(fai) as f:
        for line in f:
            name, length = line.split("\t")[:2]
            out.write(f"##contig=<ID={name},length={length}>\n")
    out.writelines(meta)
    out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n")
    for r in records:
        out.write("\t".join(r) + "\n")
PYEOF

SNP_COUNT=$(grep -vc "^#" "${PHASE_DIR}/${SAMPLE}_snps.vcf" || true)
echo "[INFO] Mixed-frequency SNPs retained for phasing (AF ${MIN_AF}-${MAX_AF}): ${SNP_COUNT}"
if [ "${SNP_COUNT}" -eq 0 ]; then
    echo "[ERROR] No phasing-informative SNPs found in ${VCF}" >&2
    exit 1
fi

# ============================================================
# STEP 2 — HapCUT2: Extract fragment information from reads
#           This creates a "fragment file" that records which
#           SNPs each read (or read pair) covers.
# ============================================================
echo ""
echo "[STEP 2] Extracting haplotype fragments from BAM..."

extractHAIRS \
    --bam "${BAM}" \
    --VCF "${PHASE_DIR}/${SAMPLE}_snps.vcf" \
    --out "${PHASE_DIR}/${SAMPLE}_fragments.txt" \
    --ref "${DUAL_REF}" \
    --10X 0

FRAG_COUNT=$(wc -l < "${PHASE_DIR}/${SAMPLE}_fragments.txt")
echo "[INFO] Fragment file: ${FRAG_COUNT} read fragments"

# ============================================================
# STEP 3 — HapCUT2: Phase the variants
# ============================================================
echo ""
echo "[STEP 3] Running HapCUT2 phasing..."

HAPCUT2 \
    --fragments "${PHASE_DIR}/${SAMPLE}_fragments.txt" \
    --VCF "${PHASE_DIR}/${SAMPLE}_snps.vcf" \
    --output "${PHASE_DIR}/${SAMPLE}_haplotypes" \
    --verbose 1

echo "[INFO] Haplotype blocks written to: ${PHASE_DIR}/${SAMPLE}_haplotypes"

# ============================================================
# STEP 4 — Parse and summarize haplotype blocks
# ============================================================
echo ""
echo "[STEP 4] Summarizing haplotype blocks..."

python3 - <<'PYEOF'
import os

haplo_file = "data/phasing/AL_haplotypes"
if not os.path.exists(haplo_file):
    print("[WARN] Haplotype file not found.")
    exit()

blocks = []
current_block = {"variants": [], "block_id": None}

with open(haplo_file) as f:
    for line in f:
        line = line.strip()
        if line.startswith("BLOCK:"):
            if current_block["variants"]:
                blocks.append(current_block)
            current_block = {"variants": [], "block_id": line}
        elif line and not line.startswith("#"):
            parts = line.split("\t")
            if len(parts) >= 8:
                current_block["variants"].append({
                    "phase_set":  parts[0],
                    "chrom":      parts[3],
                    "pos":        int(parts[4]),
                    "allele0":    parts[5],
                    "allele1":    parts[6],
                    "genotype":   parts[7],
                })

if current_block["variants"]:
    blocks.append(current_block)

print(f"\nHaplotype Summary for AL Isolate:")
print(f"  Total phased blocks: {len(blocks)}")
print(f"  Total phased SNPs:   {sum(len(b['variants']) for b in blocks)}")
print()

for i, block in enumerate(blocks, 1):
    if not block["variants"]:
        continue
    chrom     = block["variants"][0]["chrom"]
    start_pos = block["variants"][0]["pos"]
    end_pos   = block["variants"][-1]["pos"]
    n_snps    = len(block["variants"])
    span      = end_pos - start_pos

    # Count S1133 vs DIV variants in this block
    s1133_vars = sum(1 for v in block["variants"] if "S1133" in v["chrom"])
    div_vars   = sum(1 for v in block["variants"] if "DIV" in v["chrom"])

    print(f"  Block {i:>3}: {chrom}")
    print(f"    Position: {start_pos:,} – {end_pos:,}  (span: {span:,} bp)")
    print(f"    SNPs phased: {n_snps}")
    print(f"    S1133 reference variants: {s1133_vars}")
    print(f"    DIV contig variants: {div_vars}")

    if n_snps >= 5:
        # Sample the first few variants to show the haplotype structure
        print(f"    First 5 phased SNPs:")
        for v in block["variants"][:5]:
            print(f"      pos {v['pos']:>8}  hap0={v['allele0']}  hap1={v['allele1']}")
    print()

PYEOF

echo ""
echo "========================================"
echo " Haplotype phasing complete."
echo " Key file: ${PHASE_DIR}/${SAMPLE}_haplotypes"
echo " Run 06_visualize.py to see haplotype block plots"
echo "========================================"

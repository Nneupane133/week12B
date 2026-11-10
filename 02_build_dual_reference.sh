#!/bin/bash
# =============================================================================
# 02_build_dual_reference.sh
# Build a "dual reference" containing BOTH the S1133 reference genome AND
# the divergent field strain assembled from AL isolate reads.
#
# STRATEGY (read subtraction approach):
#   Step A. Download ARV S1133 reference from NCBI (10 known segments)
#   Step B. Trim the AL isolate reads (fastp)
#   Step C. Map AL reads to S1133 → extract UNMAPPED reads
#           These unmapped reads represent the divergent strain
#   Step D. De novo assemble the unmapped reads with SPAdes
#   Step E. BLAST assembled contigs against S1133 to remove any S1133 survivors
#   Step F. Combine: S1133 (NCBI) + divergent contigs (assembled) → dual reference
#   Step G. Index the dual reference
#
# The dual reference has two "worlds" of sequences:
#   - Segments named "S1133_L1", "S1133_L2", ..., "S1133_S4"  (reference strain)
#   - Segments named "DIV_contig_1", "DIV_contig_2", ...       (divergent strain)
#
# This lets us map ALL reads to both and see which world each read prefers.
# =============================================================================

set -euo pipefail
cd "${PBS_O_WORKDIR:-$PWD}"
source /apps/x86-64/apps/anaconda_3-2025.12/etc/profile.d/conda.sh
conda activate arv_haplo

echo "========================================"
echo " Step 2: Build Dual Reference Genome"
echo "========================================"

REF_DIR="data/reference"
# Raw reads already downloaded by week12/01_download_sra.sh — reuse them
RAW_DIR="/home/aubnxn001/week12/data/raw"
AL_RAW_R1="${RAW_DIR}/AL/AL_R1.fastq.gz"
AL_RAW_R2="${RAW_DIR}/AL/AL_R2.fastq.gz"
THREADS=4

for F in "${AL_RAW_R1}" "${AL_RAW_R2}"; do
    if [ ! -f "${F}" ]; then
        echo "[ERROR] Raw read file not found: ${F}" >&2
        exit 1
    fi
done
echo "[INFO] Using raw reads from: ${RAW_DIR}"

mkdir -p "${REF_DIR}"
mkdir -p "data/trimmed/AL"
mkdir -p "data/assembly"

# ============================================================
# STEP A — Download ARV S1133 reference from NCBI
# ============================================================
echo ""
echo "[STEP A] Downloading ARV S1133 reference sequences from NCBI..."

S1133_FA="${REF_DIR}/S1133_genome.fasta"

# S1133 segments: GenBank KF741756-KF741765 (one complete sequence each)
ACCESSIONS=(
    KF741756  # L1
    KF741757  # L2
    KF741758  # L3
    KF741759  # M1
    KF741760  # M2
    KF741761  # M3
    KF741762  # S1
    KF741763  # S2
    KF741764  # S3
    KF741765  # S4
)
EXPECTED=${#ACCESSIONS[@]}
ID_LIST=$(IFS=,; echo "${ACCESSIONS[*]}")

# NCBI occasionally drops connections, so retry a few times
for ATTEMPT in 1 2 3; do
    if efetch -db nuccore -id "${ID_LIST}" -format fasta > "${S1133_FA}" \
        && [ "$(grep -c '^>' "${S1133_FA}" || true)" -eq "${EXPECTED}" ]; then
        break
    fi
    echo "[WARN] Attempt ${ATTEMPT} failed or incomplete; retrying in 10s..."
    sleep 10
done

# Stop if any segment is missing — a partial S1133 reference makes reads
# from the missing segments look "divergent" and corrupts the dual reference
SEQ_COUNT=$(grep -c "^>" "${S1133_FA}" || true)
if [ "${SEQ_COUNT}" -ne "${EXPECTED}" ]; then
    echo "[ERROR] Expected ${EXPECTED} S1133 segments but got ${SEQ_COUNT}." >&2
    exit 1
fi
for ACC in "${ACCESSIONS[@]}"; do
    if ! grep -q "^>${ACC}" "${S1133_FA}"; then
        echo "[ERROR] ${ACC} missing from download." >&2
        exit 1
    fi
done

# Rename headers to include "S1133_" prefix
python3 - <<'PYEOF'
from Bio import SeqIO

seen = set()
records = []
for rec in SeqIO.parse("data/reference/S1133_genome.fasta", "fasta"):
    acc = rec.id.split(".")[0]
    if acc not in seen:
        seen.add(acc)
        # Rename so they're clearly labeled as S1133
        rec.id = f"S1133_{rec.id}"
        rec.description = f"S1133 reference | {rec.description}"
        records.append(rec)

SeqIO.write(records, "data/reference/S1133_genome.fasta", "fasta")
print(f"[INFO] S1133 reference: {len(records)} segments saved.")
PYEOF

echo "[INFO] S1133 segments:"
grep "^>" "${S1133_FA}"

# ============================================================
# STEP B — Trim AL isolate reads
# ============================================================
echo ""
echo "[STEP B] Trimming AL isolate reads..."

fastp \
    -i "${AL_RAW_R1}" -I "${AL_RAW_R2}" \
    -o "data/trimmed/AL/AL_R1.fastq.gz" \
    -O "data/trimmed/AL/AL_R2.fastq.gz" \
    --json "data/trimmed/AL/fastp.json" \
    --html "data/trimmed/AL/fastp.html" \
    --thread "${THREADS}" \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 50 \
    2>&1 | grep -E "reads|passed|rate"

# ============================================================
# STEP C — Map AL reads to S1133; extract UNMAPPED reads
# ============================================================
echo ""
echo "[STEP C] Mapping AL reads to S1133 to isolate divergent reads..."

# Index S1133 for BWA-MEM2
bwa-mem2 index "${S1133_FA}"
samtools faidx "${S1133_FA}"

# Map
bwa-mem2 mem -t "${THREADS}" \
    -R "@RG\tID:AL\tSM:AL\tPL:ILLUMINA" \
    "${S1133_FA}" \
    "data/trimmed/AL/AL_R1.fastq.gz" \
    "data/trimmed/AL/AL_R2.fastq.gz" \
| samtools sort -@ "${THREADS}" -o "data/assembly/AL_vs_S1133.bam"

samtools index "data/assembly/AL_vs_S1133.bam"

# Report mapping rate
echo "[INFO] AL isolate mapping to S1133:"
samtools flagstat "data/assembly/AL_vs_S1133.bam" | grep -E "mapped|total"

# Extract UNMAPPED read pairs (both mates unmapped: flag -f 12)
# -f 12 = read unmapped AND mate unmapped
samtools view -f 12 -b "data/assembly/AL_vs_S1133.bam" \
| samtools sort -n -@ "${THREADS}" -o "data/assembly/AL_unmapped.bam"

# Convert unmapped BAM to paired FASTQ
samtools fastq \
    -1 "data/assembly/AL_divergent_R1.fastq.gz" \
    -2 "data/assembly/AL_divergent_R2.fastq.gz" \
    -0 /dev/null -s /dev/null \
    "data/assembly/AL_unmapped.bam"

UNMAPPED_COUNT=$(zcat "data/assembly/AL_divergent_R1.fastq.gz" | wc -l | awk '{print $1/4}')
echo "[INFO] Unmapped (divergent) reads extracted: ${UNMAPPED_COUNT}"

# ============================================================
# STEP D — De novo assemble unmapped (divergent) reads with SPAdes
# ============================================================
echo ""
echo "[STEP D] Assembling divergent reads with SPAdes..."

# --rnaviral: optimized for viral RNA genomes
# If rnaviral mode not available, use --rna or regular mode
spades.py \
    --rnaviral \
    -1 "data/assembly/AL_divergent_R1.fastq.gz" \
    -2 "data/assembly/AL_divergent_R2.fastq.gz" \
    -o "data/assembly/divergent_spades" \
    --threads "${THREADS}" \
    --memory 8 \
    || \
spades.py \
    --rna \
    -1 "data/assembly/AL_divergent_R1.fastq.gz" \
    -2 "data/assembly/AL_divergent_R2.fastq.gz" \
    -o "data/assembly/divergent_spades" \
    --threads "${THREADS}" \
    --memory 8

echo "[INFO] Assembly complete."
ls -lh "data/assembly/divergent_spades/"*.fasta 2>/dev/null || \
    ls -lh "data/assembly/divergent_spades/"

# ============================================================
# STEP E — BLAST contigs against S1133; remove S1133-like ones
# ============================================================
echo ""
echo "[STEP E] Filtering assembled contigs to keep only divergent ones..."

ASSEMBLY_FA="data/assembly/divergent_spades/transcripts.fasta"
[ ! -f "${ASSEMBLY_FA}" ] && ASSEMBLY_FA="data/assembly/divergent_spades/contigs.fasta"

# Count total contigs
TOTAL_CONTIGS=$(grep -c "^>" "${ASSEMBLY_FA}" || echo 0)
echo "[INFO] Total assembled contigs: ${TOTAL_CONTIGS}"

# Build BLAST database from S1133
makeblastdb -in "${S1133_FA}" -dbtype nucl -out "data/reference/S1133_blastdb"

# BLAST assembled contigs against S1133
blastn \
    -query "${ASSEMBLY_FA}" \
    -db "data/reference/S1133_blastdb" \
    -out "data/assembly/contigs_vs_S1133.blast" \
    -outfmt "6 qseqid sseqid pident length qcovs evalue" \
    -perc_identity 90 \
    -qcov_hsp_perc 50 \
    -num_threads "${THREADS}"

# Remove contigs that match S1133 at ≥90% identity over ≥50% of their length
# Keep only the truly divergent ones
python3 - <<'PYEOF'
from Bio import SeqIO

# Read BLAST hits - which contigs are S1133-like?
s1133_like = set()
blast_file = "data/assembly/contigs_vs_S1133.blast"
try:
    with open(blast_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 5:
                qid, sid, pident, length, qcov = parts[:5]
                if float(pident) >= 90 and float(qcov) >= 50:
                    s1133_like.add(qid)
except FileNotFoundError:
    pass

print(f"[INFO] Contigs classified as S1133-like (removed): {len(s1133_like)}")

# Keep only divergent contigs, minimum length 300 bp
assembly = "data/assembly/divergent_spades/transcripts.fasta"
import os
if not os.path.exists(assembly):
    assembly = "data/assembly/divergent_spades/contigs.fasta"

divergent = []
for rec in SeqIO.parse(assembly, "fasta"):
    if rec.id not in s1133_like and len(rec.seq) >= 300:
        # Rename with DIV_ prefix
        rec.id = f"DIV_{rec.id}"
        rec.description = f"Divergent ARV field strain | {rec.description}"
        divergent.append(rec)

SeqIO.write(divergent, "data/reference/divergent_genome.fasta", "fasta")
print(f"[INFO] Divergent contigs kept: {len(divergent)}")
for rec in divergent:
    print(f"  {rec.id}: {len(rec.seq)} bp")
PYEOF

# ============================================================
# STEP F — Combine S1133 + divergent → dual reference
# ============================================================
echo ""
echo "[STEP F] Building dual reference (S1133 + divergent)..."

DUAL_REF="${REF_DIR}/dual_reference.fasta"
cat "${S1133_FA}" "data/reference/divergent_genome.fasta" > "${DUAL_REF}"

TOTAL_SEQS=$(grep -c "^>" "${DUAL_REF}")
echo "[INFO] Dual reference contains ${TOTAL_SEQS} sequences:"
grep "^>" "${DUAL_REF}"

# ============================================================
# STEP G — Index the dual reference
# ============================================================
echo ""
echo "[STEP G] Indexing dual reference..."
bwa-mem2 index "${DUAL_REF}"
samtools faidx "${DUAL_REF}"

echo ""
echo "========================================"
echo " Dual reference built: ${DUAL_REF}"
echo " S1133 segments: $(grep -c 'S1133_' ${DUAL_REF} | head -1 || grep -c '^>S1133' ${DUAL_REF})"
echo " Divergent contigs: $(grep -c '^>DIV_' ${DUAL_REF})"
echo "========================================"

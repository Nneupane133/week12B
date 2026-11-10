# ARV Haplotype / Reassortment Analysis
### BioProject PRJNA1217582 — Data Re-use Project (Week 12)

---

## Biological Question

> **Is the Alabama (AL) field isolate of Avian Orthoreovirus (ARV) a single reassortant virus,
> or two separate co-infecting viruses that failed to separate during plaque purification?**

The original paper (Avian Diseases 70:23, 2026) reported that the AL field isolate contained
sequences divergent from the vaccine strain S1133, but plaque purification did not cleanly
separate the two populations. Two competing explanations exist:

| Hypothesis | Prediction |
|---|---|
| **One reassortant virus** | Its 10-segment genome is mosaic — some segments S1133-like, others from a divergent lineage. All 10 segments travel together in the same particle. |
| **Two separate viruses** | Two distinct virions co-infect the same cells. Plaque purification failed because they depend on each other (complementation) or have identical plaque morphology. |

---

## Strategy: Dual Reference + Read-Pair Linkage

### Why a dual reference?

Standard mapping to S1133 alone forces every AL read onto the nearest S1133 segment, masking
divergence. Instead, we build a **combined reference** with two "worlds":

```
S1133_L1, S1133_L2, ..., S1133_S4   ← vaccine strain (10 segments from NCBI)
DIV_contig_1, DIV_contig_2, ...      ← divergent field strain (assembled from unmapped AL reads)
```

Each read competes for its best match across both worlds. The winning alignment reveals which
strain that read "belongs to" on a per-segment basis.

### The key test: cross-strain read pairs

An Illumina paired-end read pair comes from a **single DNA molecule**. If R1 maps to `S1133_L2`
and R2 maps to `DIV_contig_7`, that molecule must span both strains — which is only possible if
both sequences exist on the **same physical virion**.

| Cross-strain pair rate | Interpretation |
|---|---|
| > 10% | ONE virus — segments physically co-packaged |
| 1–10% | Ambiguous / mixed population |
| < 1% | TWO viruses — strains cleanly segregated |

### Haplotype phasing (confirmatory)

HapCUT2 phases SNPs within each segment using read-pair overlap. If two strains co-exist as
separate virions, each segment should phase into two clean haplotypes — one S1133-like, one
divergent. If only one reassortant exists, phasing within segments may be uninformative (each
molecule carries only one haplotype already).

---

## Pipeline Overview

```
00_setup.sh              Install conda environment (arv_haplo)
01_download_sra.sh       Download 6 SRA runs (AL, S1133, P1, PP3, PPP, EP)
02_build_dual_reference.sh  Build S1133 + divergent combined reference
03_map_and_variants.sh   Map all samples to dual reference; call variants with LoFreq
04_read_pair_linkage.py  Cross-strain read pair analysis (the key test)
05_phase_haplotypes.sh   HapCUT2 haplotype phasing within AL isolate
06_visualize.py          Generate all publication-ready figures
```

---

## Sample Table

| Sample ID | SRA Accession  | Description |
|-----------|---------------|-------------|
| AL        | SRX27528194   | Alabama field isolate — KEY SAMPLE |
| S1133     | SRX30741467   | S1133 vaccine strain control |
| P1        | SRX30741463   | 1st CELi passage of AL |
| PP3       | SRX30741464   | 3× plaque-purified AL |
| PPP       | SRX30741465   | Post-purification CELi passage |
| EP        | SRX30741466   | Embryo passage |

---

## Script Details

### `00_setup.sh`
Creates a conda environment `arv_haplo` and installs:
- **SRA-toolkit** — SRA download
- **fastp** — read trimming
- **BWA-MEM2** — fast short-read mapper
- **SPAdes** — de novo assembly (`--rnaviral` mode for viral RNA)
- **BLAST+** — classify assembled contigs
- **samtools / bcftools** — BAM/VCF manipulation
- **LoFreq** — low-frequency viral variant caller (detects variants 1–99%)
- **HapCUT2** (`extractHAIRS`, `HAPCUT2`) — read-based haplotype phasing
- **bedtools, seqkit** — sequence utilities
- **Python**: pysam, pandas, numpy, matplotlib, seaborn, scipy, networkx, biopython

### `01_download_sra.sh`
Downloads all 6 SRX accessions using `prefetch` + `fasterq-dump --split-files`.
Renames files to `{SAMPLE}_R1/R2.fastq.gz` for consistent naming downstream.

### `02_build_dual_reference.sh`
Seven-step pipeline to construct the combined reference:

1. Download S1133 genome from NCBI (esearch + efetch + known accessions)
2. Add `S1133_` prefix to all S1133 sequence headers
3. Trim AL reads with fastp (Q20, length ≥ 50 bp)
4. Map AL reads to S1133; extract **unmapped** read pairs (`samtools view -f 12`)
5. De novo assemble unmapped reads with **SPAdes `--rnaviral`**
6. BLAST assembled contigs against S1133; remove any with ≥90% identity over ≥50% coverage
7. Rename surviving contigs with `DIV_` prefix; concatenate with S1133 → `dual_reference.fasta`

### `03_map_and_variants.sh`
For each of the 6 samples:
- Trims with fastp (or reuses trimmed reads from step 02)
- Maps to `dual_reference.fasta` with BWA-MEM2
- Collects `samtools flagstat` (overall mapping rate) and `samtools idxstats` (per-segment read counts)
- Runs `lofreq indelqual` + `lofreq call-parallel` to call variants at any frequency

### `04_read_pair_linkage.py`
Core analysis. For each sample:
1. Reads the BAM file; groups reads by `query_name` (read pair)
2. Classifies each read as mapping to `S1133`, `DIV`, or `unmapped`
3. Classifies each **pair** as:
   - `S1133_S1133` — both mates map to S1133 segments
   - `DIV_DIV` — both mates map to divergent contigs
   - `CROSS_STRAIN` — mates map to different strains ← **the key signal**
   - `unmapped` — one or both mates unmapped
4. Computes cross-strain rate = CROSS_STRAIN / (all classified pairs)
5. Also reads `idxstats` to compute per-segment % S1133 vs DIV reads

**Output files:**
- `results/cross_strain_pairs.tsv` — one row per cross-strain pair with read names
- `results/linkage_report.txt` — per-sample summary

### `05_phase_haplotypes.sh`
HapCUT2 haplotype phasing on the AL isolate:
1. Filter VCF to biallelic SNPs
2. `extractHAIRS` — build fragment file (which SNPs each read covers)
3. `HAPCUT2` — phase SNPs into haplotype blocks
4. Inline Python summary: block count, SNP count, span, and first 5 phased variants per block

### `06_visualize.py`
Generates four publication-ready figures:

| Figure | Content |
|--------|---------|
| `fig1_cross_strain_rates.png` | Bar chart of cross-strain pair rate for each sample |
| `fig2_segment_genome_map.png` | Per-segment S1133% vs DIV% reads, all samples stacked |
| `fig3_haplotype_blocks.png` | HapCUT2 block swimlanes showing phased SNP positions |
| `fig4_summary_poster.png` | 4-panel poster: verdict panel + interpretation guide |

---

## Running the Pipeline

```bash
# Full pipeline (takes ~1 hour with download):
bash run_all.sh

# Or step by step:
bash 00_setup.sh
bash 01_download_sra.sh
bash 02_build_dual_reference.sh
bash 03_map_and_variants.sh
python3 04_read_pair_linkage.py
bash 05_phase_haplotypes.sh
python3 06_visualize.py
```

---

## Interpreting the Results

### Decision tree

```
Cross-strain rate in AL isolate
         |
     > 5%  ──────────────────────────────► ONE REASSORTANT VIRUS
         |                                   • Segments A and B travel together
         |                                   • Look at fig2 to see which segments
         |                                     are S1133-like vs. divergent
         |                                   • Plaque purification failed because
         |                                     there was only ONE virus to separate
         |
     1–5% ──────────────────────────────► AMBIGUOUS
         |                                   • Mixed population possible
         |                                   • Check per-segment map (fig2):
         |                                     if one strain dominates ALL segments,
         |                                     lean toward two viruses
         |
     < 1% ──────────────────────────────► TWO SEPARATE CO-INFECTING VIRUSES
                                            • Strains cleanly segregated
                                            • Plaque purification failed because
                                              of complementation or similar plaques
                                            • Per-segment map should show ALL
                                              segments assigned to same strain
```

### Haplotype blocks (fig3)

- **Few large blocks (>500 bp, >10 SNPs)** → clean haplotype structure → two strains
- **Many small, fragmentary blocks** → no clean haplotypes → one complex virus

### Per-segment map (fig2)

For the AL isolate, look for the "mosaic" pattern:
- If some segments are 100% S1133 and others are 100% DIV → reassortment (one mosaic virus)
- If ALL segments show the same strain → consistent with two separate virions

---

## Tools Reference

| Tool | Purpose | Citation |
|------|---------|---------|
| BWA-MEM2 | Fast short-read alignment | Vasimuddin et al. 2019 |
| SPAdes `--rnaviral` | De novo viral RNA assembly | Bankevich et al. 2012 |
| LoFreq | Low-frequency variant calling | Wilm et al. 2012 |
| HapCUT2 | Read-based haplotype phasing | Edge et al. 2017 |
| samtools | BAM manipulation | Li et al. 2009 |
| bcftools | VCF manipulation | Danecek et al. 2021 |
| fastp | Read trimming + QC | Chen et al. 2018 |
| BLAST+ | Sequence similarity search | Camacho et al. 2009 |

---

## Project Structure

```
arv_haplotype_project/
├── 00_setup.sh
├── 01_download_sra.sh
├── 02_build_dual_reference.sh
├── 03_map_and_variants.sh
├── 04_read_pair_linkage.py
├── 05_phase_haplotypes.sh
├── 06_visualize.py
├── run_all.sh
├── README.md
└── data/                     (created at runtime)
    ├── raw/                  SRA FASTQ downloads
    ├── trimmed/              fastp output
    ├── reference/            S1133 + dual reference FASTA
    ├── assembly/             SPAdes output
    ├── bam/                  mapped BAM files
    ├── stats/                flagstat, idxstats
    ├── variants/             LoFreq VCF files
    └── phasing/              HapCUT2 output
└── results/                  (created at runtime)
    ├── linkage_report.txt
    ├── cross_strain_pairs.tsv
    └── figures/
        ├── fig1_cross_strain_rates.png
        ├── fig2_segment_genome_map.png
        ├── fig3_haplotype_blocks.png
        └── fig4_summary_poster.png
```

---

*Data re-use project — Week 12 | Poultry Pathology, Auburn University*
*BioProject: PRJNA1217582 | Original paper: Avian Diseases 70:23 (2026)*

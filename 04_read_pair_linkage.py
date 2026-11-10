#!/usr/bin/env python3
"""
04_read_pair_linkage.py
=======================
THE CORE ANALYSIS — determines whether the AL field isolate contains
one virus (a reassortant) or two completely separate co-infecting viruses.

HOW IT WORKS:
  For every sequenced read PAIR (R1 + R2 from the same DNA fragment):
    - If both mates map to S1133 segments → pair is "S1133-S1133"
    - If both mates map to divergent (DIV) contigs → pair is "DIV-DIV"
    - If one mate maps to S1133, the other to DIV → pair is "CROSS-STRAIN"

  CROSS-STRAIN pairs are the key signal:
  ┌─────────────────────────────────────────────────────────────┐
  │ HIGH cross-strain rate (>5%)                                │
  │   → Segments from both strains co-occur on the same DNA    │
  │     fragment → they were PACKAGED TOGETHER in the same     │
  │     particle, or recombination occurred WITHIN a segment   │
  │   → Suggests ONE virus (reassortant) not TWO separate ones │
  │                                                             │
  │ LOW cross-strain rate (<1%)                                 │
  │   → Each read pair comes cleanly from one strain           │
  │   → Suggests TWO separate, independently replicating       │
  │     viruses that simply co-exist in the sample             │
  └─────────────────────────────────────────────────────────────┘

ADDITIONAL ANALYSIS — Per-segment strain assignment:
  For ARV's 10 segments (L1-L3, M1-M3, S1-S4), we ask:
  "What fraction of reads for this segment came from S1133 vs DIV?"

  If truly reassorted → some segments will be 100% S1133-like,
  others 100% DIV-like (clean segment-level separation)

  If single co-infection → all segments show mixed mapping

Outputs:
  results/linkage_report.txt
  results/figures/read_pair_linkage.png
  results/figures/segment_strain_map.png
  results/cross_strain_pairs.tsv  (the actual cross-strain read pairs)
"""

import os
import sys
import pysam
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from collections import defaultdict, Counter

# ---------- Paths ----------
BAM_DIR   = "data/bam"
STATS_DIR = "data/stats"
OUT_DIR   = "results"
FIG_DIR   = os.path.join(OUT_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

SAMPLES = {
    "AL":    "AL Isolate (field) ← KEY",
    "S1133": "Pure S1133 Control",
    "P1":    "1st CELi Passage",
    "PP3":   "3× Plaque Purification",
    "PPP":   "Post-Purif. CELi",
    "EP":    "Embryo Passage",
}

# ARV genome segments in order
SEGMENT_ORDER = ["L1", "L2", "L3", "M1", "M2", "M3", "S1", "S2", "S3", "S4"]

# ============================================================
# Helper: classify reference name as S1133 or divergent
# ============================================================
def classify_ref(ref_name):
    """Return 'S1133', 'DIV', or 'unknown'."""
    if ref_name is None:
        return "unmapped"
    if ref_name.startswith("S1133_"):
        return "S1133"
    if ref_name.startswith("DIV_"):
        return "DIV"
    return "unknown"

def get_segment_name(ref_name):
    """Extract segment name from reference sequence name."""
    for seg in SEGMENT_ORDER:
        if seg in ref_name.upper():
            return seg
    return "unknown"

# ============================================================
# ANALYSIS 1 — Cross-strain read pair detection
# ============================================================
def analyze_read_pairs(sample, bam_path, report, all_cross_pairs):
    """Scan every read pair in the BAM and classify by strain."""
    if not os.path.exists(bam_path):
        report.append(f"  {sample}: [BAM not found — run 03_map_and_variants.sh]")
        return None

    print(f"[INFO] Analyzing read pairs: {sample}...")

    bam = pysam.AlignmentFile(bam_path, "rb")
    refs = bam.references  # list of all reference sequence names

    # Collect read pairs by query name
    # We iterate through the BAM and group reads by name
    # For large BAMs use name-sorted BAM; for these small files, load all into memory
    reads_by_name = defaultdict(list)
    for read in bam.fetch(until_eof=True):
        if read.is_secondary or read.is_supplementary:
            continue  # only primary alignments
        reads_by_name[read.query_name].append(read)
    bam.close()

    counts = Counter()
    cross_pairs = []

    for name, reads in reads_by_name.items():
        if len(reads) == 1:
            counts["singleton"] += 1
            continue
        if len(reads) > 2:
            # Keep only read1 and read2
            r1_list = [r for r in reads if r.is_read1]
            r2_list = [r for r in reads if r.is_read2]
            if not r1_list or not r2_list:
                counts["singleton"] += 1
                continue
            r1, r2 = r1_list[0], r2_list[0]
        else:
            r1, r2 = reads[0], reads[1]

        # Get reference name for each mate
        ref1 = refs[r1.reference_id] if not r1.is_unmapped else None
        ref2 = refs[r2.reference_id] if not r2.is_unmapped else None

        strain1 = classify_ref(ref1)
        strain2 = classify_ref(ref2)
        seg1    = get_segment_name(ref1 or "")
        seg2    = get_segment_name(ref2 or "")

        # Classify the pair
        if strain1 == "unmapped" or strain2 == "unmapped":
            counts["unmapped"] += 1
        elif strain1 == "S1133" and strain2 == "S1133":
            counts["S1133_S1133"] += 1
        elif strain1 == "DIV" and strain2 == "DIV":
            counts["DIV_DIV"] += 1
        elif (strain1 == "S1133" and strain2 == "DIV") or \
             (strain1 == "DIV" and strain2 == "S1133"):
            counts["CROSS_STRAIN"] += 1
            # Save details of cross-strain pairs
            cross_pairs.append({
                "sample":   sample,
                "read_name": name,
                "ref_R1":   ref1,
                "ref_R2":   ref2,
                "strain_R1": strain1,
                "strain_R2": strain2,
                "segment_R1": seg1,
                "segment_R2": seg2,
            })
        else:
            counts["other"] += 1

    # Calculate rates
    total_paired = sum(counts[k] for k in ["S1133_S1133", "DIV_DIV", "CROSS_STRAIN"])
    cross_rate   = 100 * counts["CROSS_STRAIN"] / total_paired if total_paired > 0 else 0

    report.append(f"\n  {sample} ({SAMPLES.get(sample, '')})")
    report.append(f"    Total pairs classified:   {total_paired:>10,}")
    report.append(f"    Both S1133:               {counts['S1133_S1133']:>10,}")
    report.append(f"    Both DIV (divergent):     {counts['DIV_DIV']:>10,}")
    report.append(f"    CROSS-STRAIN pairs:       {counts['CROSS_STRAIN']:>10,}  ({cross_rate:.2f}%)")
    report.append(f"    Unmapped / singletons:    {counts['unmapped'] + counts['singleton']:>10,}")

    # Interpretation
    if cross_rate > 10:
        report.append(f"    → HIGH cross-strain rate ({cross_rate:.1f}%) — ONE VIRUS (reassortant) likely")
    elif cross_rate > 2:
        report.append(f"    → MODERATE cross-strain rate ({cross_rate:.1f}%) — ambiguous, check per-segment")
    else:
        report.append(f"    → LOW cross-strain rate ({cross_rate:.1f}%) — TWO SEPARATE VIRUSES likely")

    all_cross_pairs.extend(cross_pairs)
    counts["cross_rate"] = cross_rate
    return dict(counts)

# ============================================================
# ANALYSIS 2 — Per-segment strain assignment
# ============================================================
def analyze_per_segment(report):
    """
    For each ARV segment, calculate what fraction of reads mapped to
    the S1133 version vs the DIV version across all samples.
    """
    report.append("\n" + "="*60)
    report.append("ANALYSIS 2: Per-Segment Strain Assignment")
    report.append("="*60)
    report.append("For each of ARV's 10 genome segments, what fraction of reads")
    report.append("mapped to the S1133 reference vs the divergent (DIV) assembly?")
    report.append("")

    seg_data = {}  # {sample: {segment: {'S1133': n, 'DIV': n}}}

    for sample in SAMPLES:
        idx_path = os.path.join(STATS_DIR, f"{sample}_dual_idxstats.txt")
        if not os.path.exists(idx_path):
            continue

        seg_data[sample] = defaultdict(lambda: {"S1133": 0, "DIV": 0})
        with open(idx_path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 3:
                    continue
                ref_name  = parts[0]
                read_count = int(parts[2])
                strain     = classify_ref(ref_name)
                segment    = get_segment_name(ref_name)
                if strain in ("S1133", "DIV") and segment != "unknown":
                    seg_data[sample][segment][strain] += read_count

    # Print table for AL isolate (the most important sample)
    if "AL" in seg_data:
        report.append(f"{'Segment':<10} {'S1133 reads':>14} {'DIV reads':>14} {'% S1133':>10}  Verdict")
        report.append("-"*65)
        for seg in SEGMENT_ORDER:
            if seg in seg_data["AL"]:
                s = seg_data["AL"][seg]["S1133"]
                d = seg_data["AL"][seg]["DIV"]
                total = s + d
                pct_s = 100 * s / total if total > 0 else 0
                if pct_s > 80:
                    verdict = "S1133-like"
                elif pct_s < 20:
                    verdict = "DIVERGENT"
                else:
                    verdict = "MIXED ← interesting"
                report.append(f"  {seg:<8} {s:>14,} {d:>14,} {pct_s:>9.1f}%  {verdict}")

        # Reassortment conclusion
        s1133_segs = [seg for seg in SEGMENT_ORDER
                      if seg in seg_data["AL"] and
                      seg_data["AL"][seg]["S1133"] > seg_data["AL"][seg]["DIV"]]
        div_segs   = [seg for seg in SEGMENT_ORDER
                      if seg in seg_data["AL"] and
                      seg_data["AL"][seg]["DIV"] > seg_data["AL"][seg]["S1133"]]
        report.append("")
        report.append(f"Segments predominantly S1133-like: {', '.join(s1133_segs) or 'none'}")
        report.append(f"Segments predominantly divergent:  {', '.join(div_segs) or 'none'}")
        if s1133_segs and div_segs:
            report.append("")
            report.append("→ REASSORTMENT DETECTED: Different segments assigned to")
            report.append("  different strains. This is consistent with a SINGLE")
            report.append("  reassortant virus, not two separate co-infecting strains.")
        elif len(div_segs) == len(SEGMENT_ORDER):
            report.append("→ ALL segments divergent: original paper's S1133-like contigs")
            report.append("  may be an artifact (see contamination analysis).")
        elif len(s1133_segs) == len(SEGMENT_ORDER):
            report.append("→ ALL segments S1133-like: the AL isolate is S1133, not a field strain.")
        else:
            report.append("→ Unclear pattern. Check cross-strain pair rate.")

    return seg_data

# ============================================================
# FIGURES
# ============================================================
def make_linkage_figure(pair_counts, seg_data):
    """Create two figures: (1) cross-strain pair bar chart, (2) segment genome map."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })

    # --- Figure 1: Cross-strain pair rates ---
    if pair_counts:
        fig, ax = plt.subplots(figsize=(9, 5))
        samples  = list(pair_counts.keys())
        rates    = [pair_counts[s].get("cross_rate", 0) for s in samples]
        colors   = ["#e63946" if s == "AL" else "#adb5bd" for s in samples]
        bars     = ax.bar(samples, rates, color=colors, edgecolor="white")
        ax.axhline(1.0, ls="--", color="gray", lw=1.2,
                   label="1% — typical noise threshold")
        ax.axhline(5.0, ls="--", color="#2a9d8f", lw=1.2,
                   label="5% — reassortment signal threshold")
        for bar, rate in zip(bars, rates):
            ax.text(bar.get_x() + bar.get_width()/2,
                    bar.get_height() + 0.1,
                    f"{rate:.2f}%", ha="center", va="bottom", fontsize=9)
        ax.set_ylabel("Cross-Strain Read Pair Rate (%)")
        ax.set_title(
            "Cross-Strain Read Pair Rate per Sample\n"
            "High rate → segments co-packaged → ONE virus (reassortant)\n"
            "Low rate  → segments segregated  → TWO separate viruses",
            fontsize=10
        )
        ax.legend(fontsize=9)
        ax.set_ylim(0, max(rates + [6]) * 1.15)
        red_patch   = mpatches.Patch(color="#e63946", label="AL Isolate (key sample)")
        gray_patch  = mpatches.Patch(color="#adb5bd", label="Other samples")
        ax.legend(handles=[red_patch, gray_patch] +
                  [ax.get_lines()[0], ax.get_lines()[1]], fontsize=8)
        fig.tight_layout()
        path = os.path.join(FIG_DIR, "cross_strain_read_pairs.png")
        fig.savefig(path, dpi=180, bbox_inches="tight")
        plt.close()
        print(f"[SAVED] {path}")

    # --- Figure 2: Per-segment genome map (AL isolate) ---
    if "AL" in seg_data:
        fig, ax = plt.subplots(figsize=(12, 4))
        ax.set_title(
            "ARV Genome Segment Strain Assignment — AL Isolate\n"
            "Each bar shows the fraction of reads mapping to S1133 (blue) vs Divergent (red)",
            fontsize=11
        )
        x = np.arange(len(SEGMENT_ORDER))
        width = 0.6

        s1133_pct = []
        div_pct   = []
        for seg in SEGMENT_ORDER:
            if seg in seg_data["AL"]:
                s = seg_data["AL"][seg]["S1133"]
                d = seg_data["AL"][seg]["DIV"]
                total = s + d if (s + d) > 0 else 1
                s1133_pct.append(100 * s / total)
                div_pct.append(100 * d / total)
            else:
                s1133_pct.append(0)
                div_pct.append(0)

        ax.bar(x, s1133_pct, width, label="S1133-like", color="#457b9d")
        ax.bar(x, div_pct,   width, label="Divergent",  color="#e63946",
               bottom=s1133_pct)

        ax.set_xticks(x)
        ax.set_xticklabels(SEGMENT_ORDER, fontsize=11)
        ax.set_ylabel("% of Reads")
        ax.set_ylim(0, 115)
        ax.set_xlabel("ARV Genome Segment")
        ax.axhline(50, color="k", lw=0.8, ls="--", alpha=0.5)
        ax.legend(loc="upper right", fontsize=9)

        # Annotate each bar with the assignment
        for i, (s, d) in enumerate(zip(s1133_pct, div_pct)):
            label = "S1133" if s > 80 else ("DIV" if d > 80 else "MIXED")
            ax.text(i, 102, label, ha="center", va="bottom",
                    fontsize=8, fontweight="bold",
                    color="#457b9d" if label == "S1133" else
                          "#e63946" if label == "DIV" else "darkorange")

        fig.tight_layout()
        path = os.path.join(FIG_DIR, "segment_strain_map.png")
        fig.savefig(path, dpi=180, bbox_inches="tight")
        plt.close()
        print(f"[SAVED] {path}")

# ============================================================
# MAIN
# ============================================================
def main():
    report = []
    report.append("=" * 60)
    report.append("ARV HAPLOTYPE / READ-PAIR LINKAGE ANALYSIS REPORT")
    report.append("BioProject: PRJNA1217582")
    report.append("Question: Is the AL isolate ONE reassortant or TWO separate viruses?")
    report.append("=" * 60)

    # --- Analysis 1: cross-strain read pairs ---
    report.append("\n" + "="*60)
    report.append("ANALYSIS 1: Cross-Strain Read Pair Detection")
    report.append("="*60)
    report.append("Read pairs where one mate maps to S1133 and the other to DIV")
    report.append("indicate the two strain sequences co-occur on the same molecule.")
    report.append("")

    pair_counts   = {}
    all_cross     = []
    for sample in SAMPLES:
        bam_path = os.path.join(BAM_DIR, f"{sample}_dual.bam")
        result   = analyze_read_pairs(sample, bam_path, report, all_cross)
        if result:
            pair_counts[sample] = result

    # Save cross-strain pair details
    if all_cross:
        df_cross = pd.DataFrame(all_cross)
        cross_path = os.path.join(OUT_DIR, "cross_strain_pairs.tsv")
        df_cross.to_csv(cross_path, sep="\t", index=False)
        report.append(f"\n[INFO] Cross-strain pair details saved: {cross_path}")
        report.append(f"[INFO] Total cross-strain pairs across all samples: {len(all_cross)}")

    # --- Analysis 2: per-segment strain assignment ---
    seg_data = analyze_per_segment(report)

    # --- Overall verdict ---
    report.append("\n" + "="*60)
    report.append("OVERALL CONCLUSION")
    report.append("="*60)
    report.append("")
    report.append("ONE REASSORTANT VIRUS is likely if:")
    report.append("  ✓ Cross-strain pair rate in AL isolate > 5%")
    report.append("  ✓ Different segments assigned to different strains (S1133 vs DIV)")
    report.append("  ✓ Cross-strain pairs share the same segment in R1 and R2")
    report.append("    (indicating within-segment recombination, not just co-packaging)")
    report.append("")
    report.append("TWO SEPARATE VIRUSES are likely if:")
    report.append("  ✓ Cross-strain pair rate < 1% (noise level)")
    report.append("  ✓ All 10 segments map cleanly to ONE strain only")
    report.append("  ✓ No cross-strain pairs share the same segment")
    report.append("")
    report.append("IMPLICATION:")
    report.append("  If ONE reassortant → original paper may have assembled a")
    report.append("    single virus's natural genetic diversity into 'two strains'.")
    report.append("    Plaque purification failed because there was only one virus.")
    report.append("")
    report.append("  If TWO separate → original paper's conclusion is correct,")
    report.append("    and the follow-up question becomes: why did purification fail?")

    report_text = "\n".join(report)
    print(report_text)
    report_path = os.path.join(OUT_DIR, "linkage_report.txt")
    with open(report_path, "w") as f:
        f.write(report_text)
    print(f"\n[SAVED] Report: {report_path}")

    # Make figures
    make_linkage_figure(pair_counts, seg_data)
    print("\n[DONE] Linkage analysis complete.")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
06_visualize.py
===============
Final visualization — creates all publication-ready figures for the
haplotype / reassortment analysis.

Figures produced:
  Figure 1 — Cross-strain read pair rates per sample (bar chart)
  Figure 2 — ARV genome segment strain map for AL isolate (stacked bar)
  Figure 3 — Haplotype block structure within each segment (swimlane plot)
  Figure 4 — Summary poster: "One virus or two?" with verdict

Reads from:
  results/linkage_report.txt
  results/cross_strain_pairs.tsv
  data/phasing/AL_haplotypes
  data/stats/*_dual_idxstats.txt
"""

import os
import re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from collections import defaultdict

OUT_DIR   = "results"
FIG_DIR   = os.path.join(OUT_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

SEGMENT_ORDER = ["L1", "L2", "L3", "M1", "M2", "M3", "S1", "S2", "S3", "S4"]

# ARV segment approximate lengths (bp) — used for proportional drawing
SEG_LENGTHS = {
    "L1": 3959, "L2": 3828, "L3": 3900,
    "M1": 2283, "M2": 2159, "M3": 1996,
    "S1": 1643, "S2": 1324, "S3": 1198, "S4": 1096,
}

COLORS = {
    "S1133":  "#457b9d",
    "DIV":    "#e63946",
    "MIXED":  "#f4a261",
    "AL":     "#e63946",
    "other":  "#adb5bd",
}

SAMPLES_ORDER = ["AL", "S1133", "P1", "PP3", "PPP", "EP"]
SAMPLES_LABEL = {
    "AL":    "AL Isolate\n(field) ★",
    "S1133": "S1133\n(control)",
    "P1":    "1st CELi\nPassage",
    "PP3":   "3× Plaque\nPurif.",
    "PPP":   "Post-Purif.\nCELi",
    "EP":    "Embryo\nPassage",
}

# ============================================================
# Load data
# ============================================================

# S1133 reference accessions (GenBank KF741756-KF741765) → segment
S1133_ACC_TO_SEG = {
    "KF741756": "L1", "KF741757": "L2", "KF741758": "L3",
    "KF741759": "M1", "KF741760": "M2", "KF741761": "M3",
    "KF741762": "S1", "KF741763": "S2", "KF741764": "S3", "KF741765": "S4",
}

DUAL_REF     = "data/reference/dual_reference.fasta"
DIV_SEG_FILE = "data/reference/div_contig_segments.tsv"


def s1133_segment(ref_name):
    """'S1133_KF741756.1' → 'L1' (None if not a known accession)."""
    acc = ref_name.replace("S1133_", "", 1).split(".")[0]
    return S1133_ACC_TO_SEG.get(acc)


def assign_div_contigs():
    """
    Assign each DIV contig to an ARV segment by its best tblastx hit
    against the S1133 segments in the dual reference (translated search,
    so it works even for divergent strains). Cached in DIV_SEG_FILE.
    """
    if os.path.exists(DIV_SEG_FILE) and \
       os.path.getmtime(DIV_SEG_FILE) >= os.path.getmtime(DUAL_REF):
        with open(DIV_SEG_FILE) as f:
            return dict(line.rstrip("\n").split("\t")[:2] for line in f if line.strip())

    import subprocess, tempfile
    from Bio import SeqIO
    recs = list(SeqIO.parse(DUAL_REF, "fasta"))
    s1133 = [r for r in recs if r.id.startswith("S1133_")]
    div   = [r for r in recs if r.id.startswith("DIV_")]
    if not s1133 or not div:
        return {}

    with tempfile.TemporaryDirectory() as tmp:
        sub, qry = os.path.join(tmp, "s1133.fa"), os.path.join(tmp, "div.fa")
        SeqIO.write(s1133, sub, "fasta")
        SeqIO.write(div, qry, "fasta")
        out = subprocess.run(
            ["tblastx", "-query", qry, "-subject", sub, "-evalue", "1e-5",
             "-outfmt", "6 qseqid sseqid bitscore"],
            capture_output=True, text=True, check=True).stdout

    best = {}
    for line in out.splitlines():
        q, s, bits = line.split("\t")
        if q not in best or float(bits) > best[q][1]:
            best[q] = (s, float(bits))
    mapping = {q: (s1133_segment(s) or "unknown") for q, (s, _) in best.items()}
    for r in div:
        mapping.setdefault(r.id, "unknown")

    with open(DIV_SEG_FILE, "w") as f:
        for q, seg in sorted(mapping.items()):
            f.write(f"{q}\t{seg}\n")
    print(f"[INFO] Assigned {len(mapping)} DIV contigs to segments → {DIV_SEG_FILE}")
    return mapping


DIV_TO_SEG = None


def load_idxstats(sample):
    """Load per-segment read counts from idxstats."""
    global DIV_TO_SEG
    path = f"data/stats/{sample}_dual_idxstats.txt"
    if not os.path.exists(path):
        return {}
    if DIV_TO_SEG is None:
        DIV_TO_SEG = assign_div_contigs()
    seg = defaultdict(lambda: {"S1133": 0, "DIV": 0})
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            ref  = parts[0]
            cnt  = int(parts[2])
            if ref.startswith("S1133_"):
                seg[s1133_segment(ref) or "unknown"]["S1133"] += cnt
            elif ref.startswith("DIV_"):
                seg[DIV_TO_SEG.get(ref, "unknown")]["DIV"] += cnt
    return seg


def load_cross_pairs():
    """Load cross-strain pair details if available."""
    path = os.path.join(OUT_DIR, "cross_strain_pairs.tsv")
    if os.path.exists(path):
        return pd.read_csv(path, sep="\t")
    return pd.DataFrame()


def parse_linkage_report():
    """Extract cross-strain rates from the linkage report."""
    path = os.path.join(OUT_DIR, "linkage_report.txt")
    rates = {}
    if not os.path.exists(path):
        return rates
    with open(path) as f:
        current = None
        for line in f:
            for s in SAMPLES_ORDER:
                if f"  {s} (" in line:
                    current = s
            if current and "CROSS-STRAIN pairs:" in line:
                m = re.search(r"\((\d+\.\d+)%\)", line)
                if m:
                    rates[current] = float(m.group(1))
    return rates


def load_haplotype_blocks():
    """Parse HapCUT2 output."""
    path = "data/phasing/AL_haplotypes"
    if not os.path.exists(path):
        return []
    blocks = []
    current = {"variants": [], "block_id": ""}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("BLOCK:"):
                if current["variants"]:
                    blocks.append(current)
                current = {"variants": [], "block_id": line}
            elif line and not line.startswith("#"):
                parts = line.split("\t")
                if len(parts) >= 7:
                    current["variants"].append({
                        "chrom": parts[3],
                        "pos":   int(parts[4]),
                        "h0":    parts[5],
                        "h1":    parts[6],
                    })
    if current["variants"]:
        blocks.append(current)
    return blocks


# ============================================================
# Figure 1 — Cross-strain rates
# ============================================================

def fig_cross_strain_rates(rates):
    fig, ax = plt.subplots(figsize=(10, 5))
    samps  = [s for s in SAMPLES_ORDER if s in rates]
    vals   = [rates[s] for s in samps]
    labels = [SAMPLES_LABEL.get(s, s) for s in samps]
    colors = [COLORS["AL"] if s == "AL" else COLORS["other"] for s in samps]

    bars = ax.bar(range(len(samps)), vals, color=colors, edgecolor="white", width=0.6)
    ax.axhline(1.0, ls="--", color="#2d6a4f", lw=1.5, label="1% — noise ceiling")
    ax.axhline(5.0, ls="--", color="#b5838d", lw=1.5, label="5% — reassortment signal")

    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2,
                bar.get_height() + 0.1,
                f"{val:.2f}%", ha="center", va="bottom", fontsize=10, fontweight="bold")

    ax.set_xticks(range(len(samps)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Cross-Strain Read Pair Rate (%)", fontsize=11)
    ax.set_title(
        "Cross-Strain Read Pair Rate per Sample\n"
        "One mate → S1133; other mate → Divergent strain",
        fontsize=12
    )
    ax.set_ylim(0, max(vals + [6.5]) * 1.2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    al_patch   = mpatches.Patch(color=COLORS["AL"],    label="AL Isolate (key sample)")
    oth_patch  = mpatches.Patch(color=COLORS["other"], label="Other samples")
    ax.legend(handles=[al_patch, oth_patch,
                        plt.Line2D([0],[0], color="#2d6a4f", ls="--", label="1% noise ceiling"),
                        plt.Line2D([0],[0], color="#b5838d", ls="--", label="5% reassortment signal")],
              fontsize=9, loc="upper right")

    path = os.path.join(FIG_DIR, "fig1_cross_strain_rates.png")
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"[SAVED] {path}")


# ============================================================
# Figure 2 — Segment genome map
# ============================================================

def fig_segment_map(all_seg_data):
    n_samples = len([s for s in SAMPLES_ORDER if s in all_seg_data])
    fig, axes = plt.subplots(n_samples, 1, figsize=(12, 2 * n_samples + 1), sharex=True)
    if n_samples == 1:
        axes = [axes]

    fig.suptitle(
        "ARV Genome: Per-Segment Strain Assignment Across Samples\n"
        "Blue = S1133-like reads | Red = Divergent reads",
        fontsize=12, fontweight="bold"
    )

    ax_idx = 0
    for s in SAMPLES_ORDER:
        if s not in all_seg_data:
            continue
        ax = axes[ax_idx]
        seg_data = all_seg_data[s]

        x     = np.arange(len(SEGMENT_ORDER))
        s1_v  = []
        div_v = []
        for seg in SEGMENT_ORDER:
            sv = seg_data.get(seg, {"S1133": 0, "DIV": 0})
            total = sv["S1133"] + sv["DIV"] if (sv["S1133"] + sv["DIV"]) > 0 else 1
            s1_v.append(100 * sv["S1133"] / total)
            div_v.append(100 * sv["DIV"] / total)

        ax.bar(x, s1_v,  0.7, color=COLORS["S1133"], label="S1133")
        ax.bar(x, div_v, 0.7, color=COLORS["DIV"],   label="DIV", bottom=s1_v)
        ax.set_yticks([0, 50, 100])
        ax.set_ylabel(SAMPLES_LABEL.get(s, s), fontsize=9, rotation=0, labelpad=55, va="center")
        ax.set_ylim(0, 115)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.axhline(50, color="k", lw=0.5, ls="--", alpha=0.4)

        for i, (sv, dv) in enumerate(zip(s1_v, div_v)):
            label = "S1133" if sv > 80 else ("DIV" if dv > 80 else "MIX")
            color = COLORS["S1133"] if label == "S1133" else \
                    COLORS["DIV"]   if label == "DIV"   else COLORS["MIXED"]
            ax.text(i, 103, label, ha="center", va="bottom",
                    fontsize=7, color=color, fontweight="bold")
        ax_idx += 1

    axes[-1].set_xticks(range(len(SEGMENT_ORDER)))
    axes[-1].set_xticklabels(SEGMENT_ORDER, fontsize=10)
    axes[-1].set_xlabel("ARV Genome Segment", fontsize=11)

    path = os.path.join(FIG_DIR, "fig2_segment_genome_map.png")
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"[SAVED] {path}")


# ============================================================
# Figure 3 — Haplotype block swimlanes
# ============================================================

def fig_haplotype_blocks(blocks):
    if not blocks:
        print("[SKIP] No haplotype blocks to plot.")
        return

    n_blocks = sum(1 for b in blocks if b["variants"])
    fig, ax = plt.subplots(figsize=(14, max(5, 0.45 * n_blocks + 1.5)))
    ax.set_title(
        "Haplotype Blocks in AL Isolate — HapCUT2 Phasing\n"
        "Each block = a region where two distinct haplotypes were resolved",
        fontsize=11
    )

    y_pos   = 0
    colors  = ["#457b9d", "#e63946", "#2a9d8f", "#e9c46a", "#f4a261"]
    legend_handles = []
    max_end = 0

    for i, block in enumerate(blocks):
        if not block["variants"]:
            continue
        chrom  = block["variants"][0]["chrom"]
        start  = block["variants"][0]["pos"]
        end    = block["variants"][-1]["pos"]
        n_snps = len(block["variants"])
        col    = colors[i % len(colors)]
        max_end = max(max_end, end)

        # Plain rectangle (FancyBboxPatch pad is in data units, which
        # inflated each box by 5 rows); min width keeps tiny blocks visible
        rect = mpatches.Rectangle(
            (start, y_pos - 0.4), max(end - start, 20), 0.8,
            fc=col, ec="white", lw=0.5, alpha=0.85
        )
        ax.add_patch(rect)
        mid = (start + end) / 2
        if end - start >= 250:   # label only blocks wide enough to hold it
            ax.text(mid, y_pos, f"{n_snps} SNPs", ha="center", va="center",
                    fontsize=8, color="white", fontweight="bold")

        label = f"Block {i+1}: {chrom} ({n_snps} SNPs)"
        legend_handles.append(mpatches.Patch(fc=col, label=label))
        y_pos += 1.2

    # add_patch() does not autoscale the axes, so set limits explicitly;
    # without this the x-axis stays 0-1 and savefig(bbox_inches="tight")
    # tries to build an enormous image and runs out of memory
    ax.set_xlim(0, max_end * 1.05 + 50)
    ax.set_ylim(-0.8, y_pos + 0.5)
    ax.invert_yaxis()
    ax.set_xlabel("Position within segment/contig (bp)", fontsize=10)
    ax.set_yticks([])
    ax.legend(handles=legend_handles, fontsize=8,
              loc="upper left", bbox_to_anchor=(1.01, 1), framealpha=0.7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    path = os.path.join(FIG_DIR, "fig3_haplotype_blocks.png")
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"[SAVED] {path}")


# ============================================================
# Figure 4 — Summary poster
# ============================================================

def fig_summary_poster(rates, all_seg_data):
    fig = plt.figure(figsize=(14, 8))
    fig.patch.set_facecolor("#f8f9fa")
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.5, wspace=0.4)

    fig.suptitle(
        "ARV Haplotype Analysis — BioProject PRJNA1217582\n"
        "\"Is the AL isolate ONE reassortant or TWO separate viruses?\"",
        fontsize=14, fontweight="bold", y=1.0
    )

    # Panel A: Cross-strain rate for AL only
    ax_a = fig.add_subplot(gs[0, 0])
    if rates:
        al_rate  = rates.get("AL", 0)
        s1_rate  = rates.get("S1133", 0)
        ax_a.bar(["AL Isolate", "S1133\nControl"], [al_rate, s1_rate],
                 color=[COLORS["AL"], COLORS["S1133"]], edgecolor="white")
        ax_a.axhline(1.0, color="gray", ls="--", lw=1, label="1% noise")
        ax_a.axhline(5.0, color="green", ls="--", lw=1, label="5% signal")
        ax_a.set_ylabel("Cross-Strain Pair Rate (%)")
        ax_a.set_title("A  Cross-Strain Read Pairs", fontsize=10)
        ax_a.legend(fontsize=7)
        verdict = "One virus?" if al_rate > 5 else "Two viruses?"
        ax_a.text(0, al_rate + 0.2, verdict, ha="center", fontsize=9,
                  color="darkgreen" if al_rate > 5 else "darkred", fontweight="bold")

    # Panel B: Segment map for AL
    ax_b = fig.add_subplot(gs[0, 1:])
    if "AL" in all_seg_data:
        seg_data = all_seg_data["AL"]
        x = np.arange(len(SEGMENT_ORDER))
        s1_v  = []
        div_v = []
        for seg in SEGMENT_ORDER:
            sv = seg_data.get(seg, {"S1133": 0, "DIV": 0})
            total = sv["S1133"] + sv["DIV"] if (sv["S1133"] + sv["DIV"]) > 0 else 1
            s1_v.append(100 * sv["S1133"] / total)
            div_v.append(100 * sv["DIV"] / total)
        ax_b.bar(x, s1_v,  0.7, color=COLORS["S1133"], label="S1133-like")
        ax_b.bar(x, div_v, 0.7, color=COLORS["DIV"],   label="Divergent", bottom=s1_v)
        ax_b.set_xticks(x)
        ax_b.set_xticklabels(SEGMENT_ORDER)
        ax_b.set_ylabel("% of reads")
        ax_b.set_title("B  AL Isolate: Per-Segment Strain Assignment", fontsize=10)
        ax_b.legend(fontsize=8)
        ax_b.set_ylim(0, 120)
        for i, (sv, dv) in enumerate(zip(s1_v, div_v)):
            lbl = "S" if sv > 80 else ("D" if dv > 80 else "M")
            ax_b.text(i, 103, lbl, ha="center", fontsize=8, fontweight="bold",
                      color=COLORS["S1133"] if lbl == "S" else
                            COLORS["DIV"] if lbl == "D" else COLORS["MIXED"])

    # Panel C: Interpretation box
    ax_c = fig.add_subplot(gs[1, :])
    ax_c.axis("off")

    interpretation = """
ONE REASSORTANT VIRUS is supported when:
  • Cross-strain pair rate > 5%  (segments co-packaged in same particle)
  • Different genome segments map to different strains (true reassortment)
  • HapCUT2 finds only 1-2 clean haplotype blocks per segment

TWO SEPARATE CO-INFECTING VIRUSES is supported when:
  • Cross-strain pair rate < 1%  (clean strain segregation)
  • All 10 segments map to the same strain
  • HapCUT2 produces many small, fragmentary haplotype blocks

IMPLICATION FOR THE ORIGINAL PAPER:
  If ONE reassortant → Plaque purification failed because there was only one virus to purify.
                       The "two strains" were two divergent regions of a single mosaic genome.
  If TWO viruses     → Original conclusion stands. Future question: why did purification fail
                       (similar plaque size? mutual growth advantage? packaging interference?)
    """
    ax_c.text(0.02, 0.95, interpretation, transform=ax_c.transAxes,
              fontsize=9, va="top", fontfamily="monospace",
              bbox=dict(boxstyle="round", facecolor="white", edgecolor="#dee2e6", alpha=0.9))
    ax_c.set_title("C  Interpretation Guide", fontsize=10, loc="left")

    path = os.path.join(FIG_DIR, "fig4_summary_poster.png")
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"[SAVED] {path}")


# ============================================================
# MAIN
# ============================================================
def main():
    print("[INFO] Loading analysis results...")

    rates       = parse_linkage_report()
    cross_pairs = load_cross_pairs()
    blocks      = load_haplotype_blocks()
    all_seg     = {s: load_idxstats(s) for s in SAMPLES_ORDER}
    all_seg     = {s: v for s, v in all_seg.items() if v}

    print(f"  Samples with cross-strain rate data: {list(rates.keys())}")
    print(f"  Samples with segment data:            {list(all_seg.keys())}")
    print(f"  Haplotype blocks found:               {len(blocks)}")

    print("\n[INFO] Generating figures...")
    if rates:
        fig_cross_strain_rates(rates)
    if all_seg:
        fig_segment_map(all_seg)
    if blocks:
        fig_haplotype_blocks(blocks)
    fig_summary_poster(rates, all_seg)

    print("\n[DONE] All figures saved to:", FIG_DIR)
    print("  fig1_cross_strain_rates.png   — bar chart of cross-strain pair rates")
    print("  fig2_segment_genome_map.png   — per-segment S1133 vs DIV reads, all samples")
    print("  fig3_haplotype_blocks.png     — HapCUT2 phasing blocks in AL isolate")
    print("  fig4_summary_poster.png       — summary panel for presentation")


if __name__ == "__main__":
    main()

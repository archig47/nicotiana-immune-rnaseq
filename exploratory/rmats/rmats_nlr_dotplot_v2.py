#!/usr/bin/env python3
"""NLR gene-level ΔPSI dot plot — genes sig in ≥2 comparisons."""
import os, subprocess
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import font_manager as fm

BASE    = "/users/fyp/fyp5/project/PRJNA945175/deseq2_clean"
RMATS   = "/users/fyp/fyp5/project/PRJNA945175/rmats_v2"
OUT_DIR = os.path.join(BASE, "figures")
FONT    = "/users/fyp/fyp5/project/genome/Arial.ttf"

if os.path.exists(FONT):
    fm.fontManager.addfont(FONT)
    plt.rcParams["font.family"] = "Arial"
else:
    plt.rcParams["font.family"] = "sans-serif"

NLR_CLS = os.path.join(BASE, "nlr_prr/NLR_phylogeneti_classification.xlsx")
nlr_cls_df = pd.read_excel(NLR_CLS)
nlr_subclass = dict(zip(
    nlr_cls_df["ID"].astype(str).str.replace(r'-mRNA.*$', '', regex=True),
    nlr_cls_df["NLR_class"].astype(str)
))

NLR_COLS = {
    "CC-clade1":               "#4E79A7",
    "CC-clade2":               "#F28E2B",
    "CC-clade2-NRC-helper":    "#E15759",
    "CC-clade2-NRC-sensor-01": "#76B7B2",
    "CC-clade2-NRC-sensor-02": "#59A14F",
    "CC-clade2-NRC-sensor-03": "#EDC948",
    "TIR":                     "#B07AA1",
    "TNP":                     "#FF9DA7",
    "Unknown":                 "#BAB0AC",
}

COMPARISONS = [
    "HCRV_1dpi","HCRV_7dpi","HCRV_14dpi",
    "TSWV_1dpi","TSWV_7dpi","TSWV_14dpi",
    "TH_1dpi",  "TH_7dpi",  "TH_14dpi"
]
COMP_LABELS = [
    "1 dpi","7 dpi","14 dpi",
    "1 dpi","7 dpi","14 dpi",
    "1 dpi","7 dpi","14 dpi",
]

sig_df = pd.read_csv(os.path.join(OUT_DIR, "rmats_nlr_prr_sig_events.csv"))
nlr_sig = sig_df[sig_df["gene_class"] == "NLR"].copy()

# Only genes significant in ≥2 comparisons
gene_counts = nlr_sig.groupby("gene_id")["comparison"].nunique()
nlr_genes = gene_counts[gene_counts >= 2].index.tolist()
print(f"NLR genes with sig AS in ≥2 comparisons: {len(nlr_genes)}")
print(gene_counts[gene_counts >= 2].sort_values(ascending=False))

# Sort: most comparisons first, then by subclass
def sort_key(g):
    sc = nlr_subclass.get(g, "Unknown")
    sc_order = list(NLR_COLS.keys())
    return (-gene_counts.get(g, 0), sc_order.index(sc) if sc in sc_order else 99)
nlr_genes_sorted = sorted(nlr_genes, key=sort_key)

# Build ΔPSI and FDR matrices — for each gene×comparison, pick best SE event
# (highest |ΔPSI|); if no event at all, use NaN
dpsi_mat = {}
fdr_mat  = {}

for gene in nlr_genes_sorted:
    dpsi_mat[gene] = {}
    fdr_mat[gene]  = {}
    for comp in COMPARISONS:
        # First check significant events
        rows = nlr_sig[(nlr_sig["gene_id"] == gene) & (nlr_sig["comparison"] == comp)]
        if not rows.empty:
            best = rows.loc[rows["DPSI"].abs().idxmax()]
            dpsi_mat[gene][comp] = float(best["DPSI"])
            fdr_mat[gene][comp]  = float(best["FDR"])
        else:
            # Load raw SE file for this comparison to get the event even if not sig
            se_f = os.path.join(RMATS, comp, "SE.MATS.JC.txt")
            raw  = pd.read_csv(se_f, sep="\t")
            raw["gene_id"] = (raw["GeneID"].astype(str).str.strip('"')
                              .str.replace(r'^Nb(\d)', r'NbT2T\1', regex=True))
            raw_g = raw[raw["gene_id"] == gene]
            if not raw_g.empty:
                best_raw = raw_g.loc[raw_g["IncLevelDifference"].abs().idxmax()]
                dpsi_mat[gene][comp] = float(best_raw["IncLevelDifference"]) if pd.notna(best_raw["IncLevelDifference"]) else np.nan
                fdr_mat[gene][comp]  = float(best_raw["FDR"]) if pd.notna(best_raw["FDR"]) else np.nan
            else:
                dpsi_mat[gene][comp] = np.nan
                fdr_mat[gene][comp]  = np.nan

n_genes = len(nlr_genes_sorted)
fig_h   = 1.0 + n_genes * 0.65
fig_w   = 8.5

fig, ax = plt.subplots(figsize=(fig_w, fig_h))
fig.patch.set_facecolor("white")

x_pos = {c: i for i, c in enumerate(COMPARISONS)}

for gi, gene in enumerate(nlr_genes_sorted):
    y = n_genes - 1 - gi   # top gene at highest y
    sc  = nlr_subclass.get(gene, "Unknown")
    col = NLR_COLS.get(sc, NLR_COLS["Unknown"])
    ax.axhline(y, color="#F0F0F0", lw=0.5, zorder=0)

    for comp in COMPARISONS:
        x    = x_pos[comp]
        dpsi = dpsi_mat[gene][comp]
        fdr  = fdr_mat[gene][comp]
        if pd.isna(dpsi):
            continue
        sig = (not pd.isna(fdr)) and fdr < 0.05 and abs(dpsi) > 0.1
        size   = max(20, abs(dpsi) * 400)
        marker = "^" if dpsi > 0 else "v"
        alpha  = 0.9 if sig else 0.2
        ec     = "black" if sig else "none"
        ew     = 0.5 if sig else 0
        ax.scatter(x, y, s=size, c=col, marker=marker,
                   alpha=alpha, edgecolors=ec, linewidths=ew, zorder=3)

# Y-axis labels: gene_id [subclass] (n=X comparisons sig)
ax.set_yticks(range(n_genes))
ylabels = []
for gene in nlr_genes_sorted[::-1]:
    sc = nlr_subclass.get(gene, "Unknown")
    nc = gene_counts[gene]
    ylabels.append(f"{gene}  [{sc}]  (sig in {nc}/9)")
ax.set_yticklabels(ylabels, fontsize=8.5)

# X-axis
ax.set_xticks(range(9))
ax.set_xticklabels(COMP_LABELS, fontsize=8.5)
ax.tick_params(axis="x", length=0)

# Virus group dividers and labels
for xsep in [2.5, 5.5]:
    ax.axvline(xsep, color="#CCCCCC", lw=0.8, ls="--", zorder=1)

virus_label_y = n_genes - 0.25
for label, cx in [("HCRV", 1), ("TSWV", 4), ("TH", 7)]:
    ax.text(cx, virus_label_y, label, ha="center", va="bottom",
            fontsize=10, fontweight="bold", color="#333333")

# Legend — subclass colours
legend_patches = []
present_sc = {nlr_subclass.get(g, "Unknown") for g in nlr_genes_sorted}
for sc, col in NLR_COLS.items():
    if sc in present_sc:
        legend_patches.append(mpatches.Patch(facecolor=col, label=sc, linewidth=0))

# Size legend
size_legend = [
    plt.scatter([], [], s=0.1*400, c="grey", marker="^", label="|ΔPSI| = 0.1"),
    plt.scatter([], [], s=0.2*400, c="grey", marker="^", label="|ΔPSI| = 0.2"),
    plt.scatter([], [], s=0.3*400, c="grey", marker="^", label="|ΔPSI| = 0.3"),
]
legend_patches += size_legend

ax.legend(handles=legend_patches, fontsize=7.5, loc="lower right",
          frameon=True, framealpha=0.9, ncol=2,
          title="Subclass  /  |ΔPSI|", title_fontsize=8)

ax.set_xlim(-0.7, 8.7)
ax.set_ylim(-0.6, n_genes + 0.1)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.spines["left"].set_linewidth(0.5)
ax.spines["bottom"].set_linewidth(0.5)
ax.grid(axis="x", color="#F5F5F5", lw=0.3, zorder=0)

ax.set_title(
    "Alternative splicing events in NLR genes across viral infections\n"
    r"$\it{N. benthamiana}$ | filled = FDR < 0.05, |ΔPSI| > 0.1  |  ▲ exon inclusion  ▼ exon skipping",
    fontsize=10, pad=28, loc="center"
)

plt.subplots_adjust(left=0.38, right=0.97, top=0.86, bottom=0.1)

out_pdf = os.path.join(OUT_DIR, "rmats_nlr_dotplot_v2.pdf")
out_png = os.path.join(OUT_DIR, "rmats_nlr_dotplot_v2.png")
fig.savefig(out_pdf, dpi=300)
subprocess.run(["gs","-dNOPAUSE","-dBATCH","-sDEVICE=png16m","-r300",
                "-dGraphicsAlphaBits=4","-dTextAlphaBits=4",
                f"-sOutputFile={out_png}", out_pdf], check=True)
print(f"Saved {out_pdf}")
print(f"Saved {out_png}")

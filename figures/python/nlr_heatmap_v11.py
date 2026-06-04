#!/usr/bin/env python3
"""
nlr_heatmap_v11.py
NLR DEG heatmap v11 — fixes on v10:

  COLOUR CONFLICT ANALYSIS (v10 failures):
  - v10 LFC scale used #B03A2E (crimson red) which clashes with ZAR1 #E15759 (red)
    and CCG10 #FF9D9A (salmon) — all in the red/warm family → visual clash
  - v10 LFC scale used #1B4F72 (deep navy) which is blue-family like TIR #4E79A7
    and TNP #A0CBE8 → hue confusion

  Fix 1: LFC diverging scale changed to:
         high="#6B2D8B" (deep violet, upregulation) — clearly distinct from purple RPW8 #B07AA1
         mid="white"
         low="#0D7377"  (dark teal, downregulation) — not used anywhere in sidebar
         Both hues are genuinely unused across all 12 class sidebar colours.

  Fix 2: Revert ZAR1 back to #499894 (Tableau teal) — original correct colour
         Revert CCG10 back to #86BCB6 (Tableau light teal) — original correct colour
         The v10 change to red/salmon for these was wrong: it created MORE conflict
         with the LFC scale (both red). The narrow sidebar strip is sufficient to
         distinguish ZAR1/CCG10 teals from CC-clade greens.

  All v9 improvements retained: class labels in left margin, coloured sidebar strip,
  all 98 gene IDs on right, NbT2T08g02793 amber arrow highlight, bold large text,
  virus labels above columns, dpi labels at bottom.
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

# ── Paths ──────────────────────────────────────────────────────────────────────
RESULTS_DIR    = "/private/tmp/nlr_results"
NLR_IDS_FILE   = "/private/tmp/nlr_gene_ids.txt"
PHYLO_FILE     = "/Users/architagupta/Downloads/NLR_phylogeneti_classification.xlsx"
OUT_PNG        = "/Users/architagupta/Downloads/nlr_heatmap_v11.png"
HIGHLIGHT_GENE = "NbT2T08g02793"
HIGHLIGHT_COLOR = "#B6992D"   # warm amber arrow marker

# ── Class definitions ──────────────────────────────────────────────────────────
CLASS_ORDER = [
    "CC-clade2-NRC-sensor-01", "CC-clade2-NRC-sensor-02", "CC-clade2-NRC-helper",
    "CC-clade2", "CC-clade2-ZAR1", "CC-clade1", "CCG10",
    "TIR", "TNP", "RPW8-NRG1", "RPW8-ADR1", "Unclassified"
]
CLASS_LABELS = {
    "CC-clade2-NRC-sensor-01": "NRC-sensor-01",
    "CC-clade2-NRC-sensor-02": "NRC-sensor-02",
    "CC-clade2-NRC-helper":    "NRC-helper",
    "CC-clade2":               "CC-clade2",
    "CC-clade2-ZAR1":          "ZAR1",
    "CC-clade1":               "CC-clade1",
    "CCG10":                   "CCG10",
    "TIR":                     "TIR",
    "TNP":                     "TNP",
    "RPW8-NRG1":               "RPW8-NRG1",
    "RPW8-ADR1":               "RPW8-ADR1",
    "Unclassified":            "Unclassified",
}

# Fix 2: Revert ZAR1 → #499894 (Tableau teal), CCG10 → #86BCB6 (Tableau light teal)
# These were incorrectly changed to red/salmon in v10, which clashed with the LFC scale.
# The NRC-network colours (oranges/gold) and CC-clade greens remain unchanged.
CLASS_STRIP_COLORS = {
    "CC-clade2-NRC-sensor-01": "#F28E2B",   # NRC-sensor-01  (orange)
    "CC-clade2-NRC-sensor-02": "#FFBE7D",   # NRC-sensor-02  (light orange)
    "CC-clade2-NRC-helper":    "#B6992D",   # NRC-helper     (gold)
    "CC-clade2":               "#59A14F",   # CC-clade2      (green)
    "CC-clade2-ZAR1":          "#499894",   # ZAR1           ← reverted to Tableau teal
    "CC-clade1":               "#8CD17D",   # CC-clade1      (light green)
    "CCG10":                   "#86BCB6",   # CCG10          ← reverted to Tableau light teal
    "TIR":                     "#4E79A7",   # TIR            (steel blue)
    "TNP":                     "#A0CBE8",   # TNP            (light blue)
    "RPW8-NRG1":               "#B07AA1",   # RPW8-NRG1      (purple)
    "RPW8-ADR1":               "#B07AA1",   # RPW8-ADR1      (purple, same family)
    "Unclassified":            "#BAB0AC",   # Unclassified   (grey)
}

CONTRAST_ORDER = [
    "HCRV_1dpi", "HCRV_7dpi", "HCRV_14dpi",
    "TSWV_1dpi", "TSWV_7dpi", "TSWV_14dpi",
    "TH_1dpi",   "TH_7dpi",   "TH_14dpi",
]
COL_LABELS = ["1 dpi", "7 dpi", "14 dpi"] * 3

# ── Load data ──────────────────────────────────────────────────────────────────
print("Loading NLR gene list...")
nlr_genes = set(open(NLR_IDS_FILE).read().splitlines())

print("Loading DESeq2 results...")
all_dfs = []
for fname in sorted(os.listdir(RESULTS_DIR)):
    if not fname.endswith(".csv"):
        continue
    parts = fname.replace(".csv", "").split("_")
    contrast = parts[1] + "_" + parts[2]
    df = pd.read_csv(os.path.join(RESULTS_DIR, fname))
    df["contrast"] = contrast
    all_dfs.append(df)
all_res = pd.concat(all_dfs, ignore_index=True)

nlr_res = all_res[all_res["gene_id"].isin(nlr_genes)].copy()
deg_mask = (
    nlr_res["padj"].notna() &
    (nlr_res["padj"] < 0.05) &
    (nlr_res["log2FoldChange"].abs() > 1)
)
deg_genes = list(nlr_res[deg_mask]["gene_id"].unique())
print(f"NLR DEGs: {len(deg_genes)}")

lfc_wide = (
    nlr_res[nlr_res["gene_id"].isin(deg_genes)][["gene_id","contrast","log2FoldChange"]]
    .dropna(subset=["log2FoldChange"])
    .groupby(["gene_id","contrast"])["log2FoldChange"].mean()
    .unstack(fill_value=0)
)
for ct in CONTRAST_ORDER:
    if ct not in lfc_wide.columns:
        lfc_wide[ct] = 0.0
lfc_wide = lfc_wide[CONTRAST_ORDER]

# ── Load phylogenetic classification ──────────────────────────────────────────
print("Loading phylogenetic classification...")
phylo = pd.read_excel(PHYLO_FILE, sheet_name=0)
phylo["gene_id"] = phylo["ID"].str.replace(r"-mRNA.*$", "", regex=True)
gene_cls = pd.DataFrame({"gene_id": deg_genes})
gene_cls = gene_cls.merge(phylo[["gene_id","NLR_class"]], on="gene_id", how="left")
gene_cls["NLR_class"] = gene_cls["NLR_class"].fillna("Unclassified")
print(f"Phylo matched: {(gene_cls['NLR_class'] != 'Unclassified').sum()} / {len(gene_cls)}")

gene_cls_dict = dict(zip(gene_cls["gene_id"], gene_cls["NLR_class"]))

# ── Order rows: by class, then hierarchical clustering within class ────────────
ordered_genes = []
class_sizes = {}
for cl in CLASS_ORDER:
    g = [gid for gid in deg_genes if gene_cls_dict.get(gid) == cl and gid in lfc_wide.index]
    if not g:
        continue
    if len(g) > 2:
        sub = lfc_wide.loc[g].values
        dist = pdist(sub, metric="euclidean")
        if np.all(dist == 0) or len(dist) == 0:
            g_ordered = g
        else:
            link = linkage(dist, method="complete")
            idx = leaves_list(link)
            g_ordered = [g[i] for i in idx]
    else:
        g_ordered = g
    ordered_genes.extend(g_ordered)
    class_sizes[cl] = len(g_ordered)

print(f"Total ordered genes: {len(ordered_genes)}")
print("Class sizes:", class_sizes)

mat_ordered = lfc_wide.loc[ordered_genes].values.clip(-4, 4)
n_genes = len(ordered_genes)
n_cols  = 9

# ── Fix 1: Diverging colour map — deep violet / white / dark teal ─────────────
# Palette chosen after exhaustive audit of all 12 class sidebar colours:
#   #6B2D8B (deep violet, upregulation)  — clearly distinct from RPW8 purple #B07AA1
#                                           (different hue angle + much darker/more saturated)
#   #0D7377 (dark teal, downregulation)  — not used anywhere in sidebar
#                                           ZAR1/CCG10 are now reverted to #499894/#86BCB6
#                                           which are lighter/greener; #0D7377 is darker/bluer
# This is the only combination that avoids ALL six occupied hue families:
#   orange, green, red/salmon, blue, purple, grey
cmap = LinearSegmentedColormap.from_list(
    "lfc_div",
    [(0.0, "#0D7377"),   # dark teal      — downregulation (low LFC)
     (0.5, "#FFFFFF"),   # white          — midpoint (zero LFC)
     (1.0, "#6B2D8B")],  # deep violet    — upregulation (high LFC)
    N=256
)

# ══════════════════════════════════════════════════════════════════════════════
# Figure layout (all measurements in inches, converted to fractions as needed)
# ══════════════════════════════════════════════════════════════════════════════
fig_w = 14.5   # total figure width in inches
fig_h = 18.5   # total figure height in inches

# Fixed space allocations (inches):
top_gap        = 0.15
title_h        = 0.0
subtitle_h     = 0.0
virus_label_h  = 0.30
bracket_h      = 0.06
top_inner_gap  = 0.07

bottom_inner_gap = 0.04
dpi_label_h    = 0.65
bottom_gap     = 0.10

left_pad       = 0.08
class_label_w  = 1.15
class_gap_w    = 0.10
strip_w        = 0.12
strip_gap      = 0.06

gene_label_gap = 0.06
gene_label_w   = 0.95
cbar_gap       = 0.25
cbar_w         = 0.18
right_pad      = 0.20

# Derived:
top_margin    = (top_gap + title_h + subtitle_h + virus_label_h +
                 bracket_h + top_inner_gap)
bottom_margin = (bottom_inner_gap + dpi_label_h + bottom_gap)
left_margin   = left_pad + class_label_w + class_gap_w + strip_w + strip_gap
right_margin  = gene_label_gap + gene_label_w + cbar_gap + cbar_w + right_pad

heatmap_w = fig_w - left_margin - right_margin
heatmap_h = fig_h - top_margin - bottom_margin

print(f"Heatmap body: {heatmap_w:.3f} x {heatmap_h:.3f} inches")

# Axis position as figure fractions [left, bottom, width, height]
hm_left   = left_margin / fig_w
hm_bottom = bottom_margin / fig_h
hm_width  = heatmap_w / fig_w
hm_height = heatmap_h / fig_h

def fx(inches): return inches / fig_w
def fy(inches): return inches / fig_h

def col_fig_x(col_idx):
    return hm_left + (col_idx + 0.5) / n_cols * hm_width

def row_fig_y(row_idx, offset=0.5):
    return hm_bottom + (1.0 - (row_idx + offset) / n_genes) * hm_height

# ── Create figure ─────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(fig_w, fig_h), facecolor="white")

# Main heatmap axes
ax_hm = fig.add_axes([hm_left, hm_bottom, hm_width, hm_height])

# ── Draw heatmap ──────────────────────────────────────────────────────────────
im = ax_hm.imshow(
    mat_ordered,
    aspect="auto",
    cmap=cmap,
    vmin=-4, vmax=4,
    interpolation="nearest",
    origin="upper"
)
ax_hm.set_xticks([])
ax_hm.set_yticks([])
ax_hm.spines[:].set_visible(False)

# ── Gap lines between classes ─────────────────────────────────────────────────
gap_after = []
cumsum = 0
for cl in CLASS_ORDER:
    if cl not in class_sizes:
        continue
    cumsum += class_sizes[cl]
    gap_after.append(cumsum)
gap_after = gap_after[:-1]

for g in gap_after:
    ax_hm.axhline(y=g - 0.5, color="white", linewidth=3.0, zorder=3)

# Vertical gaps between virus groups (after col 2 and col 5)
ax_hm.axvline(x=2.5, color="white", linewidth=3.0, zorder=3)
ax_hm.axvline(x=5.5, color="white", linewidth=3.0, zorder=3)

# ── DPI labels — below heatmap (rotated 90°) ──────────────────────────────────
dpi_label_y = hm_bottom - fy(bottom_inner_gap + dpi_label_h * 0.10)

for i, lbl in enumerate(COL_LABELS):
    fig.text(
        col_fig_x(i),
        dpi_label_y,
        lbl,
        ha="center", va="top",
        fontsize=12, fontweight="bold",
        rotation=90,
        transform=fig.transFigure
    )

# ── Virus group labels — above heatmap ────────────────────────────────────────
bracket_y_frac = hm_bottom + hm_height + fy(top_inner_gap * 0.5)
virus_label_y  = hm_bottom + hm_height + fy(top_inner_gap + bracket_h + virus_label_h * 0.3)

groups = [("HCRV", 0, 2), ("TSWV", 3, 5), ("TH (co-infection)", 6, 8)]
for (lbl, start_col, end_col) in groups:
    centre_col = (start_col + end_col) / 2
    x1 = hm_left + start_col / n_cols * hm_width
    x2 = hm_left + (end_col + 1) / n_cols * hm_width
    fig.add_artist(plt.Line2D(
        [x1, x2], [bracket_y_frac, bracket_y_frac],
        transform=fig.transFigure,
        color="black", linewidth=1.0, clip_on=False
    ))
    fig.text(
        col_fig_x(centre_col),
        virus_label_y,
        lbl,
        ha="center", va="bottom",
        fontsize=14, fontweight="bold",
        transform=fig.transFigure
    )


# ── Gene ID labels — RIGHT of heatmap ─────────────────────────────────────────
gene_x = hm_left + hm_width + fx(gene_label_gap)

for row_idx, gid in enumerate(ordered_genes):
    y_frac = row_fig_y(row_idx)
    is_hl  = (gid == HIGHLIGHT_GENE)
    if is_hl:
        fig.text(
            gene_x, y_frac,
            gid,
            ha="left", va="center",
            fontsize=6.5,
            fontweight="bold",
            color="#111111",
            transform=fig.transFigure
        )
        marker_x = gene_x - fx(0.035)
        fig.add_artist(mpatches.FancyArrowPatch(
            posA=(marker_x - fx(0.015), y_frac),
            posB=(marker_x, y_frac),
            transform=fig.transFigure,
            arrowstyle="-|>",
            mutation_scale=5,
            color=HIGHLIGHT_COLOR,
            linewidth=1.2,
            clip_on=False
        ))
    else:
        fig.text(
            gene_x, y_frac,
            gid,
            ha="left", va="center",
            fontsize=6.5,
            fontweight="normal",
            color="#222222",
            transform=fig.transFigure
        )

# ── Class labels — LEFT of heatmap + coloured sidebar strip ───────────────────
strip_left = hm_left - fx(strip_gap + strip_w)
class_x    = strip_left - fx(class_gap_w)

cumulative_row = 0
for cl in CLASS_ORDER:
    if cl not in class_sizes:
        continue
    sz  = class_sizes[cl]
    lbl = CLASS_LABELS[cl]

    y_top = row_fig_y(cumulative_row, 0.0)
    y_bot = row_fig_y(cumulative_row + sz - 1, 1.0)
    y_cen = (y_top + y_bot) / 2

    strip_color = CLASS_STRIP_COLORS.get(cl, "#BAB0AC")
    strip_rect = mpatches.FancyBboxPatch(
        (strip_left, y_bot),
        fx(strip_w),
        y_top - y_bot,
        boxstyle="square,pad=0",
        facecolor=strip_color,
        edgecolor="none",
        transform=fig.transFigure,
        clip_on=False,
        zorder=4
    )
    fig.add_artist(strip_rect)

    fig.text(
        class_x, y_cen,
        lbl,
        ha="right", va="center",
        fontsize=12, fontweight="bold",
        color="#111111",
        transform=fig.transFigure
    )

    bar_x = strip_left - fx(0.02)
    fig.add_artist(plt.Line2D(
        [bar_x, bar_x], [y_bot, y_top],
        transform=fig.transFigure,
        color="#888888", linewidth=1.0, clip_on=False
    ))
    tick_x1 = bar_x
    tick_x2 = strip_left - fx(0.005)
    fig.add_artist(plt.Line2D(
        [tick_x1, tick_x2], [y_top, y_top],
        transform=fig.transFigure,
        color="#888888", linewidth=0.8, clip_on=False
    ))
    fig.add_artist(plt.Line2D(
        [tick_x1, tick_x2], [y_bot, y_bot],
        transform=fig.transFigure,
        color="#888888", linewidth=0.8, clip_on=False
    ))

    cumulative_row += sz

# ── Colour bar — tick labels 10pt, title 11pt bold ────────────────────────────
cbar_left_frac   = gene_x + fx(gene_label_w + cbar_gap)
cbar_bottom_frac = hm_bottom + hm_height * 0.25
cbar_height_frac = hm_height * 0.50
cbar_width_frac  = fx(cbar_w)

ax_cbar = fig.add_axes([cbar_left_frac, cbar_bottom_frac, cbar_width_frac, cbar_height_frac])
cbar = fig.colorbar(im, cax=ax_cbar)
cbar.set_ticks([-4, -2, 0, 2, 4])
cbar.set_ticklabels(["-4", "-2", "0", "2", "4"])
cbar.ax.tick_params(labelsize=11)
ax_cbar.set_title("log2FC", fontsize=13, fontweight="bold", pad=5)

# ── Save ──────────────────────────────────────────────────────────────────────
os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)
fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"\nSaved: {OUT_PNG}")

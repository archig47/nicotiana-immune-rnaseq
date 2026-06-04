import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.lines import Line2D
from matplotlib.text import Text
import matplotlib.transforms as transforms
import os
from adjustText import adjust_text

base  = "/users/fyp/fyp5/project/PRJNA945175/deseq2_clean"
res_d = os.path.join(base, "results")
nlr_d = os.path.join(base, "nlr_prr")

# ── Fonts ──────────────────────────────────────────────────────────────────
fm.fontManager.addfont("/users/fyp/fyp5/project/genome/Arial.ttf")
fm.fontManager.addfont("/users/fyp/fyp5/project/genome/Arial Bold.ttf")
fm.fontManager.addfont("/users/fyp/fyp5/project/genome/Arial Italic.ttf")
plt.rcParams["font.family"] = "Arial"
plt.rcParams["mathtext.fontset"] = "dejavusans"

# ── Colour palette (exact from nlr_upset_custom.R) ─────────────────────────
CLASS_ORDER = [
    "NRC-sensor-01", "NRC-sensor-02", "NRC-helper",
    "CC-clade2", "CC-clade1",
    "ZAR1", "CCG10",
    "TIR", "TNP",
    "RPW8-NRG1", "RPW8-ADR1",
    "Unclassified",
]
CLASS_COLS = {
    "NRC-sensor-01": "#F28E2B",
    "NRC-sensor-02": "#FFBE7D",
    "NRC-helper":    "#B6992D",
    "CC-clade2":     "#59A14F",
    "CC-clade1":     "#8CD17D",
    "ZAR1":          "#499894",
    "CCG10":         "#86BCB6",
    "TIR":           "#4E79A7",
    "TNP":           "#A0CBE8",
    "RPW8-NRG1":     "#B07AA1",
    "RPW8-ADR1":     "#D4A6C8",
    "Unclassified":  "#BAB0AC",
}

# ── NLR subclass annotation ────────────────────────────────────────────────
xl = pd.read_excel(
    os.path.join(nlr_d, "NLR_phylogeneti_classification.xlsx"),
    sheet_name="Sheet2", header=None, names=["protein_id", "raw_class"]
)
xl = xl[xl["protein_id"] != "ID"].copy()
xl["gene_id"] = xl["protein_id"].str.replace(r"-mRNA.*", "", regex=True)

CLASS_MAP = {
    "CC-clade2-NRC-helper":    "NRC-helper",
    "CC-clade2-NRC-sensor-01": "NRC-sensor-01",
    "CC-clade2-NRC-sensor-02": "NRC-sensor-02",
    "CC-clade2-ZAR1":          "ZAR1",
    "CC-clade2":               "CC-clade2",
    "CC-clade1":               "CC-clade1",
    "TIR":                     "TIR",
    "TNP":                     "TNP",
    "RPW8-NRG1":               "RPW8-NRG1",
    "RPW8-ADR1":               "RPW8-ADR1",
    "CCG10":                   "CCG10",
}
xl["class"] = xl["raw_class"].map(CLASS_MAP).fillna("Unclassified")
phylo = xl[["gene_id", "class"]].drop_duplicates("gene_id")

# ── Load NLR gene list ─────────────────────────────────────────────────────
nlr_genes = set(open(os.path.join(nlr_d, "nlr_gene_ids.txt")).read().split())
nlr_genes.discard("")

# ── Load all 9 results CSVs ────────────────────────────────────────────────
dfs = []
for f in sorted(os.listdir(res_d)):
    if not f.startswith("results_") or not f.endswith(".csv"):
        continue
    parts = f.replace(".csv", "").split("_")
    virus, dpi = parts[1], parts[2]
    d = pd.read_csv(os.path.join(res_d, f))
    d["virus"] = virus
    d["dpi"]   = dpi
    dfs.append(d)
all_res = pd.concat(dfs, ignore_index=True)

# ── Filter NLR DEGs ────────────────────────────────────────────────────────
nlr_degs = all_res[
    all_res["gene_id"].isin(nlr_genes) &
    all_res["padj"].notna() &
    (all_res["padj"] < 0.05) &
    (all_res["log2FoldChange"].abs() > 1)
].copy()
nlr_degs["abs_fc"] = nlr_degs["log2FoldChange"].abs()

idx = nlr_degs.groupby("gene_id")["abs_fc"].idxmax()
plot_df = nlr_degs.loc[idx].copy()
print(f"Unique NLR DEGs: {len(plot_df)}")

# ── Join subclass ──────────────────────────────────────────────────────────
plot_df = plot_df.merge(phylo, on="gene_id", how="left")
plot_df["class"] = plot_df["class"].fillna("Unclassified")

# ── Derived axes ──────────────────────────────────────────────────────────
plot_df["log10_baseMean"] = np.log10(plot_df["baseMean"])

# ── Determine Set A: genes in the green zone ──────────────────────────────
def in_green_zone(x, y):
    """Return True if the point (x, y) lies in the green (high-priority) zone."""
    threshold = max(np.log2(2000.0 / (10.0 ** x) + 1), 4.0)
    return bool((y > threshold) or (y < -threshold))

set_a_mask = plot_df.apply(
    lambda row: in_green_zone(row["log10_baseMean"], row["log2FoldChange"]), axis=1
)
set_a_genes = set(plot_df.loc[set_a_mask, "gene_id"])
print(f"Set A (green zone): {sorted(set_a_genes)}")

# ── Set B: specific genes ─────────────────────────────────────────────────
set_b_genes = {
    "NbT2T02g01323", "NbT2T02g01211", "NbT2T06g03156", "NbT2T15g02674",
    "NbT2T17g00315", "NbT2T05g01392", "NbT2T13g01305", "NbT2T13g01835",
}

# ── Union of Set A and Set B (only label genes present in plot_df) ─────────
label_genes = (set_a_genes | set_b_genes) & set(plot_df["gene_id"])
print(f"Total genes to label: {len(label_genes)}")
print(f"Labelled genes: {sorted(label_genes)}")

# ── Figure ─────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(11, 7.5))
fig.patch.set_facecolor("white")
ax = fig.add_axes([0.09, 0.10, 0.60, 0.76])
ax.set_facecolor("white")

# ── Compute x range for fills (before plotting so we can use data extent) ──
x_min = plot_df["log10_baseMean"].min()
x_max = plot_df["log10_baseMean"].max()
x_fill = np.linspace(x_min - 0.5, x_max + 0.5, 1000)

# Compute curve values on fill x array
y_curve_upper = np.log2(2000.0 / (10.0 ** x_fill) + 1)
y_curve_lower = -np.log2(2000.0 / (10.0 ** x_fill) + 1)

y_fc4 = np.full_like(x_fill, 4.0)

y_max_thresh = np.maximum(y_curve_upper, y_fc4)
y_min_thresh = np.minimum(y_curve_upper, y_fc4)

# y-axis limits for clipping fills
y_axis_max = 10.5
y_axis_min = -10.5

# ── Background fill zones (zorder=0, behind everything) ───────────────────
# Upper half (y > 0)
# Zone 1 — High priority: y > max(curve, 4)
ax.fill_between(x_fill,
                np.clip(y_max_thresh, y_axis_min, y_axis_max),
                y_axis_max,
                color="#D5F5E3", alpha=0.35, zorder=0, linewidth=0)
# Zone 2 — Moderate: min(curve, 4) < y < max(curve, 4)
ax.fill_between(x_fill,
                np.clip(y_min_thresh, y_axis_min, y_axis_max),
                np.clip(y_max_thresh, y_axis_min, y_axis_max),
                color="#FEF5E4", alpha=0.35, zorder=0, linewidth=0)
# Zone 3 — Low priority: 0 < y < min(curve, 4)
ax.fill_between(x_fill,
                0,
                np.clip(y_min_thresh, 0, y_axis_max),
                color="#F4F4F4", alpha=0.35, zorder=0, linewidth=0)

# Lower half (y < 0) — mirrored
y_max_thresh_lo = np.maximum(y_curve_lower, -y_fc4)
y_min_thresh_lo = np.minimum(y_curve_lower, -y_fc4)

# Zone 1 — High priority: y < min(curve_lower, -4)
ax.fill_between(x_fill,
                y_axis_min,
                np.clip(y_min_thresh_lo, y_axis_min, y_axis_max),
                color="#D5F5E3", alpha=0.35, zorder=0, linewidth=0)
# Zone 2 — Moderate: min(curve_lower, -4) < y < max(curve_lower, -4)
ax.fill_between(x_fill,
                np.clip(y_min_thresh_lo, y_axis_min, y_axis_max),
                np.clip(y_max_thresh_lo, y_axis_min, y_axis_max),
                color="#FEF5E4", alpha=0.35, zorder=0, linewidth=0)
# Zone 3 — Low priority: min(curve_lower, -4) < y < 0
ax.fill_between(x_fill,
                np.clip(y_max_thresh_lo, y_axis_min, 0),
                0,
                color="#F4F4F4", alpha=0.35, zorder=0, linewidth=0)

# Restore white background on top of fills
ax.set_facecolor("white")

ax.axhline(0, color="#CCCCCC", linewidth=0.8, linestyle="--", zorder=1)

# ── Horizontal dashed lines at log2FC = +4 and -4 ─────────────────────────
ax.axhline( 4, color="#888888", linewidth=0.8, linestyle="--", zorder=2)
ax.axhline(-4, color="#888888", linewidth=0.8, linestyle="--", zorder=2)

# ── Delta = 2000 absolute expression change curves ─────────────────────────
x_curve = np.linspace(x_min - 0.1, x_max + 0.1, 500)

y_upper = np.log2(2000.0 / (10.0 ** x_curve) + 1)
y_lower = -np.log2(2000.0 / (10.0 ** x_curve) + 1)

# Only plot where |y| > 1
mask_upper = y_upper > 1
mask_lower = y_lower < -1

if mask_upper.any():
    ax.plot(x_curve[mask_upper], y_upper[mask_upper],
            color="#444444", linewidth=0.9, linestyle="--", zorder=2)
if mask_lower.any():
    ax.plot(x_curve[mask_lower], y_lower[mask_lower],
            color="#444444", linewidth=0.9, linestyle="--", zorder=2)

# ── Scatter points (plotted after reference lines so they sit on top) ──────
present_classes = [c for c in CLASS_ORDER if c in plot_df["class"].values]

for cls in present_classes:
    sub = plot_df[plot_df["class"] == cls]
    ax.scatter(
        sub["log10_baseMean"],
        sub["log2FoldChange"],
        s=55,
        c=CLASS_COLS[cls],
        label=cls,
        alpha=0.88,
        edgecolors="none",
        zorder=3,
    )

ax.set_xlabel(r"$\log_{10}$(baseMean)", fontsize=10, color="#444444")
ax.set_ylabel(r"$\log_2$FC  (max |FC| contrast)", fontsize=10, color="#444444")

ax.yaxis.grid(True, color="#E8E8E8", linewidth=0.5, zorder=0)
ax.xaxis.grid(False)
ax.set_axisbelow(True)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_color("#AAAAAA")
ax.spines["bottom"].set_color("#AAAAAA")
ax.tick_params(colors="#555555", labelsize=9)

# ── Labels for horizontal dashed lines (right edge, just above each line) ──
x_right = ax.get_xlim()[1]
ax.text(x_right, 4,  "|log2FC| = 4",
        ha="right", va="bottom",
        fontsize=7.5, color="#888888", zorder=4,
        clip_on=False)
ax.text(x_right, -4, "|log2FC| = 4",
        ha="right", va="bottom",
        fontsize=7.5, color="#888888", zorder=4,
        clip_on=False)

# ── Labels for delta = 2000 curves ────────────────────────────────────────
if mask_upper.any():
    y_lim_top = ax.get_ylim()[1]
    within = mask_upper & (y_upper <= y_lim_top)
    if within.any():
        lx = x_curve[within][0]
        ly = y_upper[within][0]
    else:
        lx = x_curve[mask_upper][0]
        ly = y_lim_top
    ax.text(lx, ly, "Δ = 2000",
            ha="left", va="bottom",
            fontsize=7.5, color="#444444", zorder=4)
if mask_lower.any():
    y_lim_bot = ax.get_ylim()[0]
    within = mask_lower & (y_lower >= y_lim_bot)
    if within.any():
        lx = x_curve[within][0]
        ly = y_lower[within][0]
    else:
        lx = x_curve[mask_lower][0]
        ly = y_lim_bot
    ax.text(lx, ly, "Δ = 2000",
            ha="left", va="top",
            fontsize=7.5, color="#444444", zorder=4)

# ── Gene labels with adjustText ────────────────────────────────────────────
label_df = plot_df[plot_df["gene_id"].isin(label_genes)].copy()

texts = []
for _, row in label_df.iterrows():
    t = ax.text(
        row["log10_baseMean"],
        row["log2FoldChange"],
        row["gene_id"],
        fontsize=6.5,
        color="#1A1A1A",
        fontfamily="Arial",
        zorder=5,
    )
    texts.append(t)

adjust_text(
    texts,
    ax=ax,
    arrowprops=dict(arrowstyle="-", color="#888888", lw=0.5),
    expand_points=(1.4, 1.4),
    force_text=(0.6, 0.6),
)

# ── Title & subtitle ──────────────────────────────────────────────────────
fig.text(0.09, 0.970,
         "NLR DEG abundance vs fold change",
         ha="left", va="top",
         fontsize=13, fontweight="bold", color="#1A202C",
         fontfamily="Arial")

italic_prop = fm.FontProperties(fname="/users/fyp/fyp5/project/genome/Arial Italic.ttf", size=8.5)
fig.text(0.09, 0.935,
         "N. benthamiana",
         ha="left", va="top",
         color="#718096",
         fontproperties=italic_prop)

fig.text(0.09 + 0.076, 0.935,
         r"  |  98 NLR DEGs  |  padj < 0.05, |$\log_2$FC| > 1  |  max-FC contrast per gene",
         ha="left", va="top",
         fontsize=8.5, color="#718096")

# ── Class colour legend ────────────────────────────────────────────────────
class_handles = [
    Line2D([0], [0], marker="o", color="w",
           markerfacecolor=CLASS_COLS[c], markersize=7,
           label=c, markeredgecolor="white", markeredgewidth=0.3)
    for c in present_classes
]
ax.legend(
    handles=class_handles,
    title="NLR class",
    title_fontsize=8.5,
    fontsize=8.5,
    loc="upper left",
    bbox_to_anchor=(1.03, 1.0),
    bbox_transform=ax.transAxes,
    frameon=False,
    handletextpad=0.5,
    labelspacing=0.5,
    borderpad=0,
)

# ── Save without bbox tight so full figure is preserved ───────────────────
out_pdf = os.path.join(base, "figures/nlr_scatter_v5.pdf")
out_png = os.path.join(base, "figures/nlr_scatter_v5.png")
fig.savefig(out_pdf, dpi=300, facecolor="white")
os.system(
    f"gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 "
    f"-dGraphicsAlphaBits=4 -dTextAlphaBits=4 "
    f"-sOutputFile={out_png} {out_pdf}"
)
print(f"Saved: {out_pdf}")
print(f"Saved: {out_png}")

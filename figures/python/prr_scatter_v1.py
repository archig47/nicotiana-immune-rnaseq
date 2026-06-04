import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.lines import Line2D
import os
from adjustText import adjust_text

base  = "/users/fyp/fyp5/project/PRJNA945175/deseq2_clean"
res_d = os.path.join(base, "results")
nlr_d = os.path.join(base, "nlr_prr")

# Fonts
fm.fontManager.addfont("/users/fyp/fyp5/project/genome/Arial.ttf")
fm.fontManager.addfont("/users/fyp/fyp5/project/genome/Arial Bold.ttf")
fm.fontManager.addfont("/users/fyp/fyp5/project/genome/Arial Italic.ttf")
plt.rcParams["font.family"] = "Arial"
plt.rcParams["mathtext.fontset"] = "dejavusans"

# PRR ecto classification
ACCESSORY = {'EGF', 'PAN', 'SDOM', 'NEW', 'RCC1', 'TNFR', 'GP_PDE', 'UNIDENTIFIED'}

def get_ecto(rules):
    if not isinstance(rules, str): return 'Other'
    tokens = rules.split(',')[1:]
    primary = next((t for t in tokens if t not in ACCESSORY), None)
    if not primary: return 'Other'
    if primary == 'LRR':                                      return 'LRR'
    if primary == 'LysM':                                     return 'LysM'
    if primary == 'WAK':                                      return 'WAK'
    if primary in ('MAL', 'SPARK'):                           return 'Malectin/CrRLK'
    if primary in ('BLEC', 'LLEC', 'CLEC', 'GLEC', 'GNK2'): return 'Lectin'
    return 'Other'

ecto_cols = {
    'LRR':            '#5255A8',
    'Lectin':         '#6E8FD4',
    'WAK':            '#22AEAD',
    'Malectin/CrRLK': '#B8365A',
    'LysM':           '#E07848',
    'Other':          '#9A9A9A',
}
ecto_order = ['LRR', 'Lectin', 'WAK', 'Malectin/CrRLK', 'LysM', 'Other']

# Load PRR landscape
landscape = pd.read_excel(
    os.path.join(nlr_d, "nlr_prr_full_landscape.xlsx"),
    sheet_name="PRR_landscape"
)
landscape["ecto"] = landscape["matched_rules"].apply(get_ecto)
prr_genes = set(landscape["gene_id"].dropna().unique())
prr_ecto  = landscape[["gene_id", "ecto"]].drop_duplicates("gene_id").set_index("gene_id")["ecto"]
print(f"Total PRR genes in landscape: {len(prr_genes)}")

# Load all 9 results CSVs
dfs = []
for f in sorted(os.listdir(res_d)):
    if not f.startswith("results_") or not f.endswith(".csv"):
        continue
    d = pd.read_csv(os.path.join(res_d, f))
    stem  = f.replace(".csv", "")
    parts = stem.split("_")
    virus = parts[1]
    dpi   = parts[2]
    d["virus"] = virus
    d["dpi"]   = dpi
    dfs.append(d)
all_res = pd.concat(dfs, ignore_index=True)

# Filter PRR DEGs
prr_degs = all_res[
    all_res["gene_id"].isin(prr_genes) &
    all_res["padj"].notna() &
    (all_res["padj"] < 0.05) &
    (all_res["log2FoldChange"].abs() > 1)
].copy()
prr_degs["abs_fc"] = prr_degs["log2FoldChange"].abs()

idx     = prr_degs.groupby("gene_id")["abs_fc"].idxmax()
plot_df = prr_degs.loc[idx].copy()
n_degs  = len(plot_df)
print(f"Unique PRR DEGs: {n_degs}")

if n_degs < 500 or n_degs > 900:
    print(f"WARNING: Expected ~666 PRR DEGs but got {n_degs}. Check data before proceeding.")

# Join ecto class
plot_df["ecto"] = plot_df["gene_id"].map(prr_ecto).fillna("Other")

# Derived axes
plot_df["log10_baseMean"] = np.log10(plot_df["baseMean"])

# Clip extreme outliers for display (rare points with |log2FC| > 14)
Y_DISPLAY = 14.0
plot_df["log2FC_plot"] = plot_df["log2FoldChange"].clip(-Y_DISPLAY, Y_DISPLAY)
plot_df["clipped"] = (plot_df["log2FoldChange"].abs() > Y_DISPLAY)
n_clipped = plot_df["clipped"].sum()
if n_clipped > 0:
    print(f"Clipping {n_clipped} point(s) with |log2FC| > {Y_DISPLAY} for display")

# Green zone uses original (unclipped) values for correctness
def in_green_zone(x, y):
    threshold = max(np.log2(2000.0 / (10.0 ** x) + 1), 4.0)
    return bool((y > threshold) or (y < -threshold))

set_a_mask = plot_df.apply(
    lambda row: in_green_zone(row["log10_baseMean"], row["log2FoldChange"]), axis=1
)
set_a_genes = set(plot_df.loc[set_a_mask, "gene_id"])
print(f"Green zone genes: {len(set_a_genes)}")

# Label gene selection: if >25 in green zone, restrict to the most extreme
# First pass: |log2FC| > 7 OR log10(baseMean) > 3.0
# If still > 30, raise to |log2FC| > 8 OR log10(baseMean) > 3.2
if len(set_a_genes) > 25:
    green_df = plot_df[plot_df["gene_id"].isin(set_a_genes)]
    label_genes = set(green_df.loc[
        (green_df["abs_fc"] > 7) | (green_df["log10_baseMean"] > 3.0),
        "gene_id"
    ])
    if len(label_genes) > 30:
        label_genes = set(green_df.loc[
            (green_df["abs_fc"] > 8) | (green_df["log10_baseMean"] > 3.2),
            "gene_id"
        ])
    print(f"Green zone >25, restricting label set: {len(label_genes)} genes")
else:
    label_genes = set_a_genes
print(f"Total genes to label: {len(label_genes)}")

# Figure
fig = plt.figure(figsize=(10, 7))
fig.patch.set_facecolor("white")
ax = fig.add_axes([0.09, 0.10, 0.60, 0.76])
ax.set_facecolor("white")

# X range for fills
x_min = plot_df["log10_baseMean"].min()
x_max = plot_df["log10_baseMean"].max()
x_fill = np.linspace(x_min - 0.5, x_max + 0.5, 1000)

y_curve_upper = np.log2(2000.0 / (10.0 ** x_fill) + 1)
y_curve_lower = -np.log2(2000.0 / (10.0 ** x_fill) + 1)

y_fc4 = np.full_like(x_fill, 4.0)

y_max_thresh = np.maximum(y_curve_upper, y_fc4)
y_min_thresh = np.minimum(y_curve_upper, y_fc4)

y_axis_max = 10.5
y_axis_min = -10.5

# Background fill zones
# Upper half
ax.fill_between(x_fill,
                np.clip(y_max_thresh, y_axis_min, y_axis_max),
                y_axis_max,
                color="#D5F5E3", alpha=0.35, zorder=0, linewidth=0)
ax.fill_between(x_fill,
                np.clip(y_min_thresh, y_axis_min, y_axis_max),
                np.clip(y_max_thresh, y_axis_min, y_axis_max),
                color="#FEF5E4", alpha=0.35, zorder=0, linewidth=0)
ax.fill_between(x_fill,
                0,
                np.clip(y_min_thresh, 0, y_axis_max),
                color="#F4F4F4", alpha=0.35, zorder=0, linewidth=0)

# Lower half
y_max_thresh_lo = np.maximum(y_curve_lower, -y_fc4)
y_min_thresh_lo = np.minimum(y_curve_lower, -y_fc4)

ax.fill_between(x_fill,
                y_axis_min,
                np.clip(y_min_thresh_lo, y_axis_min, y_axis_max),
                color="#D5F5E3", alpha=0.35, zorder=0, linewidth=0)
ax.fill_between(x_fill,
                np.clip(y_min_thresh_lo, y_axis_min, y_axis_max),
                np.clip(y_max_thresh_lo, y_axis_min, y_axis_max),
                color="#FEF5E4", alpha=0.35, zorder=0, linewidth=0)
ax.fill_between(x_fill,
                np.clip(y_max_thresh_lo, y_axis_min, 0),
                0,
                color="#F4F4F4", alpha=0.35, zorder=0, linewidth=0)

ax.set_facecolor("white")
ax.set_ylim(y_axis_min, y_axis_max)

ax.axhline(0, color="#CCCCCC", linewidth=0.8, linestyle="--", zorder=1)
ax.axhline( 4, color="#888888", linewidth=0.8, linestyle="--", zorder=2)
ax.axhline(-4, color="#888888", linewidth=0.8, linestyle="--", zorder=2)

# Delta = 2000 curves
x_curve = np.linspace(x_min - 0.1, x_max + 0.1, 500)
y_upper = np.log2(2000.0 / (10.0 ** x_curve) + 1)
y_lower = -np.log2(2000.0 / (10.0 ** x_curve) + 1)

mask_upper = y_upper > 1
mask_lower = y_lower < -1

if mask_upper.any():
    ax.plot(x_curve[mask_upper], y_upper[mask_upper],
            color="#444444", linewidth=0.9, linestyle="--", zorder=2)
if mask_lower.any():
    ax.plot(x_curve[mask_lower], y_lower[mask_lower],
            color="#444444", linewidth=0.9, linestyle="--", zorder=2)

# Scatter points (use clipped y for plotting)
present_ectos = [e for e in ecto_order if e in plot_df["ecto"].values]

for ecto in present_ectos:
    sub = plot_df[plot_df["ecto"] == ecto]
    # Normal (non-clipped) points
    sub_norm = sub[~sub["clipped"]]
    if len(sub_norm):
        ax.scatter(
            sub_norm["log10_baseMean"],
            sub_norm["log2FC_plot"],
            s=20,
            c=ecto_cols[ecto],
            label=ecto,
            alpha=0.65,
            edgecolors="none",
            zorder=3,
        )
    # Clipped outlier points: plot as triangles at y_axis boundary
    sub_clip = sub[sub["clipped"]]
    if len(sub_clip):
        for _, row in sub_clip.iterrows():
            marker = "^" if row["log2FoldChange"] > 0 else "v"
            y_pos  = y_axis_max - 0.3 if row["log2FoldChange"] > 0 else y_axis_min + 0.3
            ax.scatter(
                row["log10_baseMean"], y_pos,
                s=25, c=ecto_cols[ecto],
                marker=marker, alpha=0.85, edgecolors="none", zorder=4
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

# Labels for dashed lines
x_right = ax.get_xlim()[1]
ax.text(x_right,  4, "|log2FC| = 4",
        ha="right", va="bottom",
        fontsize=7.5, color="#888888", zorder=4, clip_on=False)
ax.text(x_right, -4, "|log2FC| = 4",
        ha="right", va="bottom",
        fontsize=7.5, color="#888888", zorder=4, clip_on=False)

# Labels for delta = 2000 curves
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

# Gene labels with adjustText (use clipped y so labels stay in view)
label_df = plot_df[plot_df["gene_id"].isin(label_genes)].copy()

texts = []
for _, row in label_df.iterrows():
    y_label = float(np.clip(row["log2FoldChange"], y_axis_min, y_axis_max))
    t = ax.text(
        row["log10_baseMean"],
        y_label,
        row["gene_id"],
        fontsize=6.0,
        color="#1A1A1A",
        fontfamily="Arial",
        zorder=5,
    )
    texts.append(t)

if texts:
    adjust_text(
        texts,
        ax=ax,
        arrowprops=dict(arrowstyle="-", color="#888888", lw=0.5),
        expand_points=(1.4, 1.4),
        force_text=(0.6, 0.6),
    )

# Title & subtitle
fig.text(0.09, 0.970,
         "PRR DEG abundance vs fold change",
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
         r"  |  666 PRR DEGs  |  padj < 0.05, |$\log_2$FC| > 1  |  max-FC contrast per gene",
         ha="left", va="top",
         fontsize=8.5, color="#718096")

# Ecto colour legend
ecto_handles = [
    Line2D([0], [0], marker="o", color="w",
           markerfacecolor=ecto_cols[e], markersize=7,
           label=e, markeredgecolor="white", markeredgewidth=0.3)
    for e in present_ectos
]
ax.legend(
    handles=ecto_handles,
    title="PRR ectodomain",
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

# Save
out_pdf = os.path.join(base, "figures/prr_scatter_v1.pdf")
out_png = os.path.join(base, "figures/prr_scatter_v1.png")
fig.savefig(out_pdf, dpi=300, facecolor="white")
os.system(
    f"gs -dNOPAUSE -dBATCH -sDEVICE=png16m -r300 "
    f"-dGraphicsAlphaBits=4 -dTextAlphaBits=4 "
    f"-sOutputFile={out_png} {out_pdf}"
)
print(f"Saved: {out_pdf}")
print(f"Saved: {out_png}")

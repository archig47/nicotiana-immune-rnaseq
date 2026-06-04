#!/usr/bin/env python3
"""
Splicing analysis figures for NLR and PRR genes.
Figure 1: Summary heatmap — n genes with significant AS per condition × gene class
Figure 2: NLR gene-level ΔPSI dot plot for multi-condition AS genes
Figure 3: PRR ecto class distribution of AS events
"""
import os, sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from matplotlib import cm
import matplotlib.gridspec as gridspec

BASE    = "/users/fyp/fyp5/project/PRJNA945175/deseq2_clean"
RMATS   = "/users/fyp/fyp5/project/PRJNA945175/rmats_v2"
RES_DIR = os.path.join(BASE, "results")
NLR_IDS = os.path.join(BASE, "nlr_prr/nlr_gene_ids.txt")
PRR_XLS = os.path.join(BASE, "nlr_prr/nlr_prr_full_landscape.xlsx")
NLR_CLS = os.path.join(BASE, "nlr_prr/NLR_phylogeneti_classification.xlsx")
OUT_DIR = os.path.join(BASE, "figures")
FONT    = "/users/fyp/fyp5/project/genome/Arial.ttf"

os.makedirs(OUT_DIR, exist_ok=True)

# ── Font ──────────────────────────────────────────────────────────────────
from matplotlib import font_manager as fm
if os.path.exists(FONT):
    fm.fontManager.addfont(FONT)
    plt.rcParams["font.family"] = "Arial"
else:
    plt.rcParams["font.family"] = "sans-serif"

# ── Gene lists & annotation ────────────────────────────────────────────────
nlr_ids = set(open(NLR_IDS).read().strip().split())

prr_df  = pd.read_excel(PRR_XLS, sheet_name="PRR_landscape")
prr_ids = set(prr_df["gene_id"].dropna().astype(str))

nlr_cls_df = pd.read_excel(NLR_CLS)
nlr_subclass = dict(zip(
    nlr_cls_df["ID"].astype(str).str.replace(r'-mRNA.*$', '', regex=True),
    nlr_cls_df["NLR_class"].astype(str)
))

ACCESSORY = {"EGF","PAN","SDOM","NEW","RCC1","TNFR","GP_PDE","UNIDENTIFIED"}
def get_ecto(rules):
    if not isinstance(rules, str) or pd.isna(rules): return "Other"
    tokens = rules.split(",")[1:]
    primary = next((t for t in tokens if t not in ACCESSORY), None)
    if primary is None:         return "Other"
    if primary == "LRR":        return "LRR"
    if primary == "LysM":       return "LysM"
    if primary == "WAK":        return "WAK"
    if primary in ("MAL","SPARK"):           return "Malectin/CrRLK"
    if primary in ("BLEC","LLEC","CLEC","GLEC","GNK2"): return "Lectin"
    return "Other"
prr_ecto = {r["gene_id"]: get_ecto(r.get("matched_rules","")) for _,r in prr_df.iterrows()}

# ── Load significant events CSV ────────────────────────────────────────────
sig_df = pd.read_csv(os.path.join(OUT_DIR, "rmats_nlr_prr_sig_events.csv"))

COMPARISONS = [
    "HCRV_1dpi","HCRV_7dpi","HCRV_14dpi",
    "TSWV_1dpi","TSWV_7dpi","TSWV_14dpi",
    "TH_1dpi",  "TH_7dpi",  "TH_14dpi"
]
COMP_LABELS = {
    "HCRV_1dpi": "HCRV\n1 dpi", "HCRV_7dpi": "HCRV\n7 dpi", "HCRV_14dpi": "HCRV\n14 dpi",
    "TSWV_1dpi": "TSWV\n1 dpi", "TSWV_7dpi": "TSWV\n7 dpi", "TSWV_14dpi": "TSWV\n14 dpi",
    "TH_1dpi":   "TH\n1 dpi",   "TH_7dpi":   "TH\n7 dpi",   "TH_14dpi":   "TH\n14 dpi",
}

# ── Load DEG data for overlap annotation ──────────────────────────────────
all_res = []
for f in os.listdir(RES_DIR):
    if not f.startswith("results_"): continue
    parts = f.replace(".csv","").split("_")
    d = pd.read_csv(os.path.join(RES_DIR, f))
    d["virus"] = parts[1]; d["dpi"] = parts[2]
    all_res.append(d)
all_res = pd.concat(all_res, ignore_index=True)

def is_deg(gene_id, virus, dpi, df=all_res):
    row = df[(df["gene_id"]==gene_id) & (df["virus"]==virus) & (df["dpi"]==dpi)]
    if row.empty: return False
    r = row.iloc[0]
    return (not pd.isna(r.get("padj",np.nan))) and r.get("padj",1) < 0.05 and abs(r.get("log2FoldChange",0)) > 1

# ── NLR colours (from UpSet R script palette) ─────────────────────────────
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

PRR_COLS = {
    "LRR"           : "#5255A8",
    "Lectin"        : "#6E8FD4",
    "WAK"           : "#22AEAD",
    "Malectin/CrRLK": "#B8365A",
    "LysM"          : "#E07848",
    "Other"         : "#9A9A9A",
}

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 1: Summary heatmap — n unique genes with sig AS per condition
# ═══════════════════════════════════════════════════════════════════════════
fig1, axes = plt.subplots(1, 2, figsize=(9, 4.5), sharey=False)
fig1.patch.set_facecolor("white")

for ax, gene_class, title_tag in zip(axes, ["NLR","PRR"], ["NLR","PRR"]):
    sub = sig_df[sig_df["gene_class"] == gene_class]
    counts = sub.groupby("comparison")["gene_id"].nunique()
    counts = counts.reindex(COMPARISONS, fill_value=0)

    virus_groups = {
        "HCRV": ["HCRV_1dpi","HCRV_7dpi","HCRV_14dpi"],
        "TSWV": ["TSWV_1dpi","TSWV_7dpi","TSWV_14dpi"],
        "TH":   ["TH_1dpi",  "TH_7dpi",  "TH_14dpi"],
    }
    dpi_labels = ["1 dpi","7 dpi","14 dpi"]
    virus_order = ["HCRV","TSWV","TH"]

    matrix = np.array([[counts[f"{v}_{d}"] for d in ["1dpi","7dpi","14dpi"]] for v in virus_order])
    vmax = max(matrix.max(), 1)

    cmap = plt.cm.YlOrRd if gene_class == "PRR" else plt.cm.Blues
    im = ax.imshow(matrix, cmap=cmap, vmin=0, vmax=vmax, aspect="auto")

    ax.set_xticks(range(3)); ax.set_xticklabels(dpi_labels, fontsize=9)
    ax.set_yticks(range(3)); ax.set_yticklabels(virus_order, fontsize=10, fontweight="bold")
    ax.tick_params(length=0)

    for i in range(3):
        for j in range(3):
            n = int(matrix[i,j])
            ax.text(j, i, str(n), ha="center", va="center",
                    fontsize=13, fontweight="bold",
                    color="white" if n > vmax*0.6 else "#333333")

    cbar = fig1.colorbar(im, ax=ax, shrink=0.75, pad=0.04)
    cbar.set_label("n genes with\nsig AS event", fontsize=8)
    cbar.ax.tick_params(labelsize=8)

    ax.set_title(f"{title_tag} genes", fontsize=11, fontweight="bold", pad=8)
    for spine in ax.spines.values(): spine.set_visible(False)

fig1.suptitle(
    "Alternatively spliced NLR and PRR genes across viral infections\n"
    r"$\it{N. benthamiana}$ | FDR < 0.05, |ΔPSI| > 0.1",
    fontsize=11, y=1.02
)
fig1.tight_layout()
out1 = os.path.join(OUT_DIR, "rmats_summary_heatmap_v1.pdf")
fig1.savefig(out1, bbox_inches="tight", dpi=300)
print(f"Saved: {out1}")

import subprocess
png1 = out1.replace(".pdf",".png")
subprocess.run(["gs","-dNOPAUSE","-dBATCH","-sDEVICE=png16m","-r300",
                "-dGraphicsAlphaBits=4","-dTextAlphaBits=4",
                f"-sOutputFile={png1}", out1], check=True)
print(f"Saved: {png1}")

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 2: NLR gene ΔPSI dot plot — top NLR genes across comparisons
# ═══════════════════════════════════════════════════════════════════════════
nlr_sig = sig_df[sig_df["gene_class"] == "NLR"].copy()
gene_counts = nlr_sig.groupby("gene_id")["comparison"].nunique()
# Include genes with sig AS in ≥2 comparisons, + any with ΔPSI > 0.3 in any comparison
top_nlr = gene_counts[gene_counts >= 2].index.tolist()
high_dpsi = nlr_sig[nlr_sig["DPSI"].abs() > 0.3]["gene_id"].unique().tolist()
nlr_genes = sorted(set(top_nlr) | set(high_dpsi),
                   key=lambda g: -gene_counts.get(g, 0))

if nlr_genes:
    # Build matrix: rows = genes, cols = comparisons
    dpsi_mat = pd.DataFrame(index=nlr_genes, columns=COMPARISONS, dtype=float)
    fdr_mat  = pd.DataFrame(index=nlr_genes, columns=COMPARISONS, dtype=float)
    for gene in nlr_genes:
        for comp in COMPARISONS:
            rows = nlr_sig[(nlr_sig["gene_id"]==gene) & (nlr_sig["comparison"]==comp)]
            if rows.empty:
                # check if event exists but not significant
                # load raw data
                se_f = os.path.join(RMATS, comp, "SE.MATS.JC.txt")
                raw = pd.read_csv(se_f, sep="\t")
                raw["gene_id"] = raw["GeneID"].astype(str).str.strip('"').str.replace(r'^Nb(\d)', r'NbT2T\1', regex=True)
                raw_g = raw[raw["gene_id"]==gene]
                if not raw_g.empty:
                    dpsi_mat.loc[gene, comp] = pd.to_numeric(raw_g["IncLevelDifference"].iloc[0], errors="coerce")
                    fdr_mat.loc[gene, comp]  = pd.to_numeric(raw_g["FDR"].iloc[0], errors="coerce")
                else:
                    dpsi_mat.loc[gene, comp] = np.nan
                    fdr_mat.loc[gene, comp]  = np.nan
            else:
                # pick event with largest |ΔPSI|
                best = rows.reindex(rows["DPSI"].abs().sort_values(ascending=False).index).iloc[0]
                dpsi_mat.loc[gene, comp] = best["DPSI"]
                fdr_mat.loc[gene, comp]  = best["FDR"]

    # Sort genes by subclass then by n comparisons sig
    def gene_sort_key(g):
        sc = nlr_subclass.get(g, "Unknown")
        sc_order = list(NLR_COLS.keys())
        return (sc_order.index(sc) if sc in sc_order else 99, -gene_counts.get(g,0))
    nlr_genes_sorted = sorted(nlr_genes, key=gene_sort_key)

    n_genes = len(nlr_genes_sorted)
    fig2_h  = max(4.5, n_genes * 0.45 + 1.5)
    fig2, ax2 = plt.subplots(figsize=(7.5, fig2_h))
    fig2.patch.set_facecolor("white")

    x_positions = {c: i for i, c in enumerate(COMPARISONS)}
    y_positions = {g: i for i, g in enumerate(nlr_genes_sorted[::-1])}

    for gene in nlr_genes_sorted:
        sc  = nlr_subclass.get(gene, "Unknown")
        col = NLR_COLS.get(sc, "#BAB0AC")
        yp  = y_positions[gene]
        for comp in COMPARISONS:
            xp   = x_positions[comp]
            dpsi = dpsi_mat.loc[gene, comp]
            fdr  = fdr_mat.loc[gene, comp]
            if pd.isna(dpsi): continue
            sig  = (not pd.isna(fdr)) and fdr < 0.05 and abs(dpsi) > 0.1
            size = abs(dpsi) * 350
            marker = "o" if dpsi > 0 else "v"
            alpha  = 0.9 if sig else 0.25
            edge   = "black" if sig else "none"
            ew     = 0.6 if sig else 0
            ax2.scatter(xp, yp, s=size, c=col, marker=marker,
                        alpha=alpha, edgecolors=edge, linewidths=ew, zorder=3)

    # Gene labels
    ax2.set_yticks(range(n_genes))
    gene_labels = []
    for g in nlr_genes_sorted[::-1]:
        sc  = nlr_subclass.get(g, "Unknown")
        nc  = gene_counts.get(g, 0)
        gene_labels.append(f"{g}  [{sc}]")
    ax2.set_yticklabels(gene_labels, fontsize=8)

    # X axis
    ax2.set_xticks(range(9))
    ax2.set_xticklabels([COMP_LABELS[c] for c in COMPARISONS], fontsize=8)
    ax2.tick_params(axis="x", length=0)

    # Virus group separators
    for xsep in [2.5, 5.5]:
        ax2.axvline(xsep, color="#BBBBBB", lw=0.8, ls="--", zorder=1)

    # Virus group labels at top
    for label, centre in [("HCRV",1), ("TSWV",4), ("TH",7)]:
        ax2.text(centre, n_genes - 0.1, label, ha="center", va="bottom",
                 fontsize=9, fontweight="bold", transform=ax2.get_xaxis_transform())

    # Horizontal gene separators
    for y in range(n_genes):
        ax2.axhline(y, color="#EEEEEE", lw=0.4, zorder=0)

    # Legend
    legend_elements = []
    for sc, col in NLR_COLS.items():
        present = any(nlr_subclass.get(g,"Unknown") == sc for g in nlr_genes_sorted)
        if present:
            legend_elements.append(mpatches.Patch(facecolor=col, label=sc))
    # size legend
    for dpsi_ex, lbl in [(0.1,"0.1"), (0.2,"0.2"), (0.3,"0.3")]:
        legend_elements.append(
            plt.scatter([], [], s=dpsi_ex*350, c="grey", marker="o",
                        label=f"|ΔPSI|={lbl}")
        )

    ax2.legend(handles=legend_elements, loc="lower right", fontsize=7,
               frameon=True, framealpha=0.9, ncol=2,
               title="Subclass / |ΔPSI|", title_fontsize=8)

    ax2.set_xlim(-0.6, 8.6)
    ax2.set_ylim(-0.7, n_genes - 0.3)
    ax2.set_title(
        "Alternative splicing in NLR genes across viral infections\n"
        r"$\it{N. benthamiana}$ | filled = FDR < 0.05, |ΔPSI| > 0.1",
        fontsize=10, pad=20
    )
    for spine in ax2.spines.values(): spine.set_linewidth(0.5)
    ax2.spines["right"].set_visible(False)
    ax2.spines["top"].set_visible(False)
    ax2.grid(axis="x", color="#F0F0F0", lw=0.3)

    fig2.tight_layout()
    out2 = os.path.join(OUT_DIR, "rmats_nlr_dotplot_v1.pdf")
    fig2.savefig(out2, bbox_inches="tight", dpi=300)
    subprocess.run(["gs","-dNOPAUSE","-dBATCH","-sDEVICE=png16m","-r300",
                    "-dGraphicsAlphaBits=4","-dTextAlphaBits=4",
                    f"-sOutputFile={out2.replace('.pdf','.png')}", out2], check=True)
    print(f"Saved: {out2}")
    print(f"Saved: {out2.replace('.pdf','.png')}")
else:
    print("No NLR genes with sig AS in ≥2 comparisons")

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 3: PRR ecto class breakdown of sig AS events (bar plot)
# ═══════════════════════════════════════════════════════════════════════════
prr_sig = sig_df[sig_df["gene_class"] == "PRR"].copy()
prr_sig["ecto"] = prr_sig["gene_id"].map(prr_ecto).fillna("Other")

ecto_order = ["LRR","Lectin","WAK","Malectin/CrRLK","LysM","Other"]

# Unique genes per ecto per comparison
ecto_counts = (prr_sig.groupby(["comparison","ecto"])["gene_id"]
               .nunique().reset_index(name="n"))
ecto_counts["ecto"] = pd.Categorical(ecto_counts["ecto"], categories=ecto_order, ordered=True)

fig3, ax3 = plt.subplots(figsize=(9, 4))
fig3.patch.set_facecolor("white")

x = np.arange(9)
bottoms = np.zeros(9)
for ecto in ecto_order:
    vals = []
    for comp in COMPARISONS:
        row = ecto_counts[(ecto_counts["comparison"]==comp) & (ecto_counts["ecto"]==ecto)]
        vals.append(int(row["n"].iloc[0]) if not row.empty else 0)
    vals = np.array(vals)
    ax3.bar(x, vals, bottom=bottoms, color=PRR_COLS[ecto], label=ecto,
            width=0.65, edgecolor="white", linewidth=0.5)
    bottoms += vals

ax3.set_xticks(x)
ax3.set_xticklabels([COMP_LABELS[c] for c in COMPARISONS], fontsize=8.5)
ax3.tick_params(axis="x", length=0)
ax3.set_ylabel("Number of PRR genes with sig AS event", fontsize=9)
ax3.set_title(
    "PRR genes with significant alternative splicing across viral infections\n"
    r"$\it{N. benthamiana}$ | FDR < 0.05, |ΔPSI| > 0.1",
    fontsize=10
)

# Virus group separators
for xsep in [2.5, 5.5]:
    ax3.axvline(xsep, color="#BBBBBB", lw=0.8, ls="--", zorder=0)

for label, centre in [("HCRV",1), ("TSWV",4), ("TH",7)]:
    ax3.text(centre, ax3.get_ylim()[1]*1.01, label, ha="center", va="bottom",
             fontsize=9, fontweight="bold")

ax3.legend(title="PRR ectodomain", fontsize=8, title_fontsize=9,
           loc="upper right", frameon=True)
ax3.spines["right"].set_visible(False)
ax3.spines["top"].set_visible(False)
ax3.grid(axis="y", color="#EEEEEE", lw=0.4)

fig3.tight_layout()
out3 = os.path.join(OUT_DIR, "rmats_prr_bar_v1.pdf")
fig3.savefig(out3, bbox_inches="tight", dpi=300)
subprocess.run(["gs","-dNOPAUSE","-dBATCH","-sDEVICE=png16m","-r300",
                "-dGraphicsAlphaBits=4","-dTextAlphaBits=4",
                f"-sOutputFile={out3.replace('.pdf','.png')}", out3], check=True)
print(f"Saved: {out3}")
print(f"Saved: {out3.replace('.pdf','.png')}")

print("\nAll figures done.")

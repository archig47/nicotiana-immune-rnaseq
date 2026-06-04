#!/usr/bin/env python3
"""PRR AS bar chart v2 — fixed virus labels below x-axis."""
import os, subprocess
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

BASE    = "/users/fyp/fyp5/project/PRJNA945175/deseq2_clean"
OUT_DIR = os.path.join(BASE, "figures")
FONT    = "/users/fyp/fyp5/project/genome/Arial.ttf"

if os.path.exists(FONT):
    fm.fontManager.addfont(FONT)
    plt.rcParams["font.family"] = "Arial"

PRR_XLS = os.path.join(BASE, "nlr_prr/nlr_prr_full_landscape.xlsx")
prr_df  = pd.read_excel(PRR_XLS, sheet_name="PRR_landscape")

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

PRR_COLS = {
    "LRR"           : "#5255A8",
    "Lectin"        : "#6E8FD4",
    "WAK"           : "#22AEAD",
    "Malectin/CrRLK": "#B8365A",
    "LysM"          : "#E07848",
    "Other"         : "#9A9A9A",
}

COMPARISONS = [
    "HCRV_1dpi","HCRV_7dpi","HCRV_14dpi",
    "TSWV_1dpi","TSWV_7dpi","TSWV_14dpi",
    "TH_1dpi",  "TH_7dpi",  "TH_14dpi"
]
COMP_LABELS = ["1 dpi","7 dpi","14 dpi","1 dpi","7 dpi","14 dpi","1 dpi","7 dpi","14 dpi"]
ecto_order  = ["LRR","Lectin","WAK","Malectin/CrRLK","LysM","Other"]

sig_df  = pd.read_csv(os.path.join(OUT_DIR, "rmats_nlr_prr_sig_events.csv"))
prr_sig = sig_df[sig_df["gene_class"] == "PRR"].copy()
prr_sig["ecto"] = prr_sig["gene_id"].map(prr_ecto).fillna("Other")

ecto_counts = (prr_sig.groupby(["comparison","ecto"])["gene_id"]
               .nunique().reset_index(name="n"))

fig, ax = plt.subplots(figsize=(9.5, 4.8))
fig.patch.set_facecolor("white")

x        = np.arange(9)
bottoms  = np.zeros(9)
for ecto in ecto_order:
    vals = []
    for comp in COMPARISONS:
        row = ecto_counts[(ecto_counts["comparison"]==comp) & (ecto_counts["ecto"]==ecto)]
        vals.append(int(row["n"].iloc[0]) if not row.empty else 0)
    vals = np.array(vals, dtype=float)
    ax.bar(x, vals, bottom=bottoms, color=PRR_COLS[ecto], label=ecto,
           width=0.68, edgecolor="white", linewidth=0.5)
    bottoms += vals

# Total labels on bars
for i, tot in enumerate(bottoms):
    ax.text(i, tot + 0.3, str(int(tot)), ha="center", va="bottom",
            fontsize=8.5, fontweight="bold", color="#333333")

ax.set_xticks(x)
ax.set_xticklabels(COMP_LABELS, fontsize=9)
ax.tick_params(axis="x", length=0, pad=2)

# Virus group dividers
for xsep in [2.5, 5.5]:
    ax.axvline(xsep, color="#CCCCCC", lw=0.9, ls="--", zorder=0)

# Virus labels as x-axis secondary labels (below tick labels)
ax.annotate("HCRV", xy=(1, 0), xycoords=("data","axes fraction"),
            xytext=(0, -32), textcoords="offset points",
            ha="center", fontsize=10, fontweight="bold", color="#333333")
ax.annotate("TSWV", xy=(4, 0), xycoords=("data","axes fraction"),
            xytext=(0, -32), textcoords="offset points",
            ha="center", fontsize=10, fontweight="bold", color="#333333")
ax.annotate("TH", xy=(7, 0), xycoords=("data","axes fraction"),
            xytext=(0, -32), textcoords="offset points",
            ha="center", fontsize=10, fontweight="bold", color="#333333")

ax.set_ylabel("Number of PRR genes with\nsignificant AS event", fontsize=9)
ax.set_title(
    "PRR genes with significant alternative splicing across viral infections\n"
    r"$\it{N. benthamiana}$  |  FDR < 0.05, |" + u"Δ" + "PSI| > 0.1",
    fontsize=10, pad=8
)
ax.legend(title="PRR ectodomain", fontsize=8, title_fontsize=9,
          loc="upper left", frameon=True, framealpha=0.9)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.grid(axis="y", color="#EEEEEE", lw=0.4)
ax.set_xlim(-0.6, 8.6)

plt.subplots_adjust(left=0.10, right=0.97, top=0.88, bottom=0.21)

out3 = os.path.join(OUT_DIR, "rmats_prr_bar_v2.pdf")
fig.savefig(out3, dpi=300)
subprocess.run(["gs","-dNOPAUSE","-dBATCH","-sDEVICE=png16m","-r300",
                "-dGraphicsAlphaBits=4","-dTextAlphaBits=4",
                "-sOutputFile=" + out3.replace(".pdf",".png"), out3], check=True)
print("Saved " + out3)

import os, re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import font_manager as fm

ARIAL = '/users/fyp/fyp5/project/genome/Arial.ttf'
fm.fontManager.addfont(ARIAL)
plt.rcParams['font.family'] = 'Arial'

BASE     = '/users/fyp/fyp5/project/PRJNA945175/'
META     = BASE + 'deseq2_clean/sample_metadata_clean.csv'
STAR_DIR = BASE + 'aligned/'
FC_MAIN  = BASE + 'featurecounts/counts_v12_multimapper.txt.summary'
FC_081   = BASE + 'featurecounts/counts_081_multimapper.txt.summary'
MULTIQC  = BASE + 'multiqc_raw/multiqc_data/multiqc_general_stats.txt'
OUT      = BASE + 'figures/qc_'

# ── Colours ──────────────────────────────────────────────────────────────────
COND_COL = {'CK': '#9E9E9E', 'HCRV': '#2596be', 'TSWV': '#E07848', 'TH': '#6B2D8B'}

# ── Metadata ─────────────────────────────────────────────────────────────────
meta = pd.read_csv(META, index_col=0)
meta['rep'] = meta.groupby('group').cumcount() + 1
cond_ord = {'CK':0,'HCRV':1,'TSWV':2,'TH':3}
tp_ord   = {1:0,7:1,14:2}
meta['_co'] = meta['condition'].map(cond_ord)
meta['_to'] = meta['timepoint'].map(tp_ord)
meta = meta.sort_values(['_co','_to','rep'])
srrs   = meta.index.tolist()
labels = [f"{r['condition']} {r['timepoint']}dpi\nrep{r['rep']}" for _, r in meta.iterrows()]
colors = [COND_COL[r['condition']] for _, r in meta.iterrows()]

# ── Group separators ─────────────────────────────────────────────────────────
# We'll draw vertical lines between condition groups
groups = []
prev = None
for i, (_, r) in enumerate(meta.iterrows()):
    c = r['condition']
    if c != prev:
        groups.append((i, c))
        prev = c

# ── Raw read counts ───────────────────────────────────────────────────────────
mqc = pd.read_csv(MULTIQC, sep='\t', index_col=0)
r1  = mqc[mqc.index.str.endswith('_1')].copy()
r1.index = r1.index.str.replace('_1','')
reads_m = r1.reindex(srrs)['fastqc-total_sequences'].values   # already in millions

# ── STAR alignment ────────────────────────────────────────────────────────────
star = {}
for srr in srrs:
    log = os.path.join(STAR_DIR, srr, 'Log.final.out')
    txt = open(log).read()
    def g(pat): m = re.search(pat, txt); return float(m.group(1)) if m else 0.0
    uni  = g(r'Uniquely mapped reads % \|\s+([\d.]+)%')
    mul  = g(r'% of reads mapped to multiple loci \|\s+([\d.]+)%')
    tm   = g(r'% of reads mapped to too many loci \|\s+([\d.]+)%')
    ush  = g(r'% of reads unmapped: too short \|\s+([\d.]+)%')
    uot  = g(r'% of reads unmapped: other \|\s+([\d.]+)%')
    star[srr] = {'unique': uni, 'multi': mul+tm, 'unmapped': ush+uot}
star_df = pd.DataFrame(star).T.reindex(srrs)

# ── featureCounts ─────────────────────────────────────────────────────────────
fc = pd.read_csv(FC_MAIN, sep='\t', index_col=0)
fc.columns = [os.path.basename(os.path.dirname(c)) for c in fc.columns]
fc081 = pd.read_csv(FC_081, sep='\t', index_col=0)
fc081.columns = ['SRR23875081']
fc = pd.concat([fc, fc081], axis=1)
fc_tot  = fc.sum(axis=0)
fc_pct  = (fc.loc['Assigned'] / fc_tot * 100).reindex(srrs).values

n = len(srrs)
x = np.arange(n)
W = 0.6

def add_group_lines(ax, groups, n, ymax):
    for i, (start, cond) in enumerate(groups):
        end = groups[i+1][0] if i+1 < len(groups) else n
        cx  = (start + end - 1) / 2
        ax.text(cx, ymax*1.02, cond, ha='center', va='bottom', fontsize=9,
                fontweight='bold', color=COND_COL[cond])
        if start > 0:
            ax.axvline(start - 0.5, color='#cccccc', linewidth=0.8, zorder=1)

def style_ax(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('#aaaaaa')
    ax.spines['bottom'].set_color('#aaaaaa')
    ax.tick_params(axis='both', labelsize=7.5, color='#aaaaaa')
    ax.yaxis.grid(True, color='#eeeeee', linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)

# ━━━ Figure 1: Input reads ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
fig, ax = plt.subplots(figsize=(12, 3.8))
bars = ax.bar(x, reads_m, width=W, color=colors, edgecolor='none', zorder=2)
ax.axhline(30, color='#e74c3c', linewidth=1, linestyle='--', zorder=3, label='30M threshold')
ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=7)
ax.set_ylabel('Input reads (millions)', fontsize=10)
ax.set_xlim(-0.6, n-0.4)
ymax = reads_m.max() * 1.15
ax.set_ylim(0, ymax)
add_group_lines(ax, groups, n, ymax)
style_ax(ax)
legend_patches = [mpatches.Patch(color=v, label=k) for k,v in COND_COL.items()]
legend_patches.append(plt.Line2D([0],[0], color='#e74c3c', linestyle='--', label='30M threshold'))
ax.legend(handles=legend_patches, fontsize=8, frameon=False, bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0)
plt.tight_layout()
plt.savefig(OUT+'reads.png', dpi=250, bbox_inches='tight', facecolor='white')
plt.savefig(OUT+'reads.pdf', bbox_inches='tight', facecolor='white')
plt.close()
print("Figure 1 done")

# ━━━ Figure 2: STAR alignment ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
fig, ax = plt.subplots(figsize=(12, 3.8))
b1 = ax.bar(x, star_df['unique'],  width=W, color='#2596be', edgecolor='none', label='Uniquely mapped', zorder=2)
b2 = ax.bar(x, star_df['multi'],   width=W, bottom=star_df['unique'],
            color='#a8d5e8', edgecolor='none', label='Multi-mapped', zorder=2)
b3 = ax.bar(x, star_df['unmapped'], width=W, bottom=star_df['unique']+star_df['multi'],
            color='#e0e0e0', edgecolor='none', label='Unmapped', zorder=2)
ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=7)
ax.set_ylabel('Reads (%)', fontsize=10)
ax.set_xlim(-0.6, n-0.4)
ax.set_ylim(0, 110)
add_group_lines(ax, groups, n, 105)
style_ax(ax)
ax.legend(fontsize=8, frameon=False, bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0)
plt.tight_layout()
plt.savefig(OUT+'star.png', dpi=250, bbox_inches='tight', facecolor='white')
plt.savefig(OUT+'star.pdf', bbox_inches='tight', facecolor='white')
plt.close()
print("Figure 2 done")

# ━━━ Figure 3: featureCounts assignment ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
fig, ax = plt.subplots(figsize=(12, 3.8))
ax.bar(x, fc_pct, width=W, color=colors, edgecolor='none', zorder=2)
ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=7)
ax.set_ylabel('Reads assigned (%)', fontsize=10)
ax.set_xlim(-0.6, n-0.4)
ymax_fc = max(fc_pct) * 1.15
ax.set_ylim(0, ymax_fc)
add_group_lines(ax, groups, n, ymax_fc)
style_ax(ax)
legend_patches = [mpatches.Patch(color=v, label=k) for k,v in COND_COL.items()]
ax.legend(handles=legend_patches, fontsize=8, frameon=False, bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0)
plt.tight_layout()
plt.savefig(OUT+'fc.png', dpi=250, bbox_inches='tight', facecolor='white')
plt.savefig(OUT+'fc.pdf', bbox_inches='tight', facecolor='white')
plt.close()
print("Figure 3 done")
print("All done — saved to", OUT + "*.png/.pdf")

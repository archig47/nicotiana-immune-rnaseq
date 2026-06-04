import os, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib.cm as cm

font_manager.fontManager.addfont('/users/fyp/fyp5/project/genome/Arial.ttf')
plt.rcParams.update({'font.family':'Arial','font.size':11,'pdf.fonttype':42,'ps.fonttype':42})

BASE  = '/users/fyp/fyp5/project/PRJNA945175/promoters/homer_results_may2026_v2'
DRAFT = '/users/fyp/fyp5/project/PRJNA945175/figures/draft'

COMPS = [
    ('nlr_vs_nonNLR',       'NLRs\nvs genome'),
    ('prr_vs_nonPRR',       'PRRs\nvs genome'),
    ('degNLR_vs_nondegNLR', 'DEG NLRs\nvs non-DEG'),
    ('degPRR_vs_nondegPRR', 'DEG PRRs\nvs non-DEG'),
    ('degNLR_vs_degNonNLR', 'DEG NLRs\nvs DEG PRRs'),
    ('degPRR_vs_degNonPRR', 'DEG PRRs\nvs DEG NLRs'),
    ('upPRR_vs_downPRR',    'Up-PRRs\nvs Down-PRRs'),
    ('upNLR_vs_downNLR',    'Up-NLRs\nvs Down-NLRs'),
]

MOTIF_ROWS = [
    ('WRKY2',     ['WRKY2']),
    ('WRKY20',    ['WRKY20']),
    ('WRKY28/53', ['WRKY28','WRKY53']),
    ('WRKY42',    ['WRKY42']),
    ('WRKY71',    ['WRKY71']),
    ('NAC098',    ['NAC098']),
    ('MYB40',     ['MYB40']),
    ('MYB121',    ['MYB121']),
    ('MYB30',     ['MYB30']),
    ('MYB33',     ['MYB33']),
    ('DREB1/CBF', ['DREB1A','DREB1B','DREB1D','DREB1E']),
    ('DREB2F',    ['DREB2F']),
    ('RAP2-9',    ['RAP2-9']),
    ('ERF15',     ['ERF15']),
    ('TCP5',      ['TCP5']),
    ('TCP2',      ['TCP2']),
    ('DOF4.2',    ['DOF4.2']),
    ('PIF4',      ['PIF4']),
    ('BZIP18',    ['BZIP18']),
    ('HSFB2A',    ['HSFB2A']),
]
ROW_LABELS = [r[0] for r in MOTIF_ROWS]
N_Y, N_X = len(ROW_LABELS), len(COMPS)

TF_FAMILY = {
    'WRKY2':'WRKY','WRKY20':'WRKY','WRKY28/53':'WRKY','WRKY42':'WRKY','WRKY71':'WRKY',
    'NAC098':'NAC',
    'MYB40':'MYB','MYB121':'MYB','MYB30':'MYB','MYB33':'MYB',
    'DREB1/CBF':'AP2/ERF','DREB2F':'AP2/ERF','RAP2-9':'AP2/ERF','ERF15':'AP2/ERF',
    'TCP5':'TCP','TCP2':'TCP',
    'DOF4.2':'DOF','PIF4':'PIF','BZIP18':'bZIP','HSFB2A':'HSF',
}
FAM_COLOR = {
    'WRKY':'#1D3557','NAC':'#2A9D8F','MYB':'#E9C46A',
    'AP2/ERF':'#E76F51','TCP':'#F4A261',
    'DOF':'#8338EC','PIF':'#6C757D','bZIP':'#AAAAAA','HSF':'#C77DFF',
}

def extract_tf(name):
    n = name.replace('/JASPAR2024_plants','')
    m = re.match(r'^([A-Za-z0-9][A-Za-z0-9\.\-]*)', n)
    return m.group(1) if m else n

LOGP_THRESH = -3.0
logp_mat = np.full((N_Y, N_X), np.nan)
fe_mat   = np.full((N_Y, N_X), np.nan)

for xi, (comp_key, _) in enumerate(COMPS):
    fpath = f'{BASE}/{comp_key}/knownResults.txt'
    if not os.path.exists(fpath): continue
    df = pd.read_csv(fpath, sep='\t')
    df['TF']    = df['Motif Name'].apply(extract_tf)
    df['logP']  = df['Log P-value'].astype(float)
    df['tpct']  = df['% of Target Sequences with Motif'].str.replace('%','').astype(float)
    df['bgpct'] = df['% of Background Sequences with Motif'].str.replace('%','').astype(float)
    df['fe']    = np.log2((df['tpct']+0.1)/(df['bgpct']+0.1))
    for yi, (_, tf_list) in enumerate(MOTIF_ROWS):
        sub = df[df['TF'].isin(tf_list)]
        if sub.empty: continue
        best = sub.loc[sub['logP'].idxmin()]
        logp_mat[yi, xi] = best['logP']
        fe_mat[yi, xi]   = best['fe']

cmap = LinearSegmentedColormap.from_list('fe',['#F7F2FB','#C9B8E8','#9B72CF','#6B2D8B'])
norm = Normalize(vmin=0.5, vmax=4.5)
mapper = cm.ScalarMappable(norm=norm, cmap=cmap)
mapper.set_array([])

DOT_SIZE = 200

# ── Layout ────────────────────────────────────────────────────────────────────
FW, FH = 14.0, 11.0
LEFT = 2.65   # wider for larger y-axis text
RIGHT = 2.9
BOT = 0.75    # taller for larger x-axis text
TOP = 2.0
pw = FW - LEFT - RIGHT   # 8.45"
ph = FH - BOT - TOP      # 8.25"

fig = plt.figure(figsize=(FW, FH))
ax  = fig.add_axes([LEFT/FW, BOT/FH, pw/FW, ph/FH])

# Alternating family bands
fam_seq = [TF_FAMILY.get(lbl,'') for lbl in ROW_LABELS]
fam_changes = [0]
for i in range(1, N_Y):
    if fam_seq[i] != fam_seq[i-1]:
        fam_changes.append(i)
fam_changes.append(N_Y)

for k, start in enumerate(fam_changes[:-1]):
    end = fam_changes[k+1]
    col = '#F2F2F2' if k%2==0 else '#FFFFFF'
    ax.axhspan(N_Y-end-0.5, N_Y-start-0.5, facecolor=col, edgecolor='none', zorder=0)

for xi in range(N_X):
    ax.axvline(xi, color='#DDDDDD', lw=0.4, zorder=1)
for yi in range(N_Y):
    ax.axhline(yi, color='#DDDDDD', lw=0.4, zorder=1)

# Dots
for yi in range(N_Y):
    for xi in range(N_X):
        lp = logp_mat[yi, xi]
        fe = fe_mat[yi, xi]
        if np.isnan(lp) or lp >= LOGP_THRESH: continue
        fe_c = float(np.clip(fe if not np.isnan(fe) else 0.5, 0.5, 4.5))
        ax.scatter(xi, N_Y-1-yi, s=DOT_SIZE, c=[mapper.to_rgba(fe_c)],
                   edgecolors='#444444', linewidths=0.4, zorder=3, clip_on=False)

for xb in [1.5, 3.5, 5.5]:
    ax.axvline(xb, color='#888888', lw=1.0, ls='--', zorder=2)

ax.set_xlim(-0.5, N_X-0.5)
ax.set_ylim(-0.5, N_Y-0.5)
ax.set_xticks(range(N_X))
ax.set_xticklabels([c[1] for c in COMPS], fontsize=10, ha='center',
                   multialignment='center', linespacing=1.3)
ax.tick_params(axis='x', length=0, pad=9)
ax.set_yticks(range(N_Y))
ax.set_yticklabels(ROW_LABELS[::-1], fontsize=11, fontstyle='italic')
ax.tick_params(axis='y', length=0, pad=24)
for sp in ax.spines.values():
    sp.set_visible(False)

# ── TF family colour strip ────────────────────────────────────────────────────
STRIP_W = 0.18
strip_ax = fig.add_axes([(LEFT - STRIP_W - 0.08)/FW, BOT/FH, STRIP_W/FW, ph/FH])
strip_ax.set_xlim(0, 1); strip_ax.set_ylim(-0.5, N_Y-0.5)
strip_ax.axis('off')
for yi, lbl in enumerate(ROW_LABELS[::-1]):
    fam = TF_FAMILY.get(lbl,'')
    strip_ax.barh(yi, 1, height=0.82, left=0,
                  color=FAM_COLOR.get(fam,'#CCCCCC'), edgecolor='none')

# ── Group labels ──────────────────────────────────────────────────────────────
col_frac  = pw / FW / N_X
ax_left_f = LEFT / FW
ax_top_f  = (BOT + ph) / FH
groups = [
    ('Receptor identity',      0, 1),
    ('Virus responsiveness',   2, 3),
    ('Cross-class enrichment', 4, 5),
    ('Direction of change',    6, 7),
]
for label, x0, x1 in groups:
    xc  = ax_left_f + ((x0+x1)/2 + 0.5) * col_frac
    bly = ax_top_f + 0.020
    bx0 = ax_left_f + (x0 + 0.08) * col_frac
    bx1 = ax_left_f + (x1 + 0.92) * col_frac
    fig.add_artist(plt.Line2D([bx0,bx1],[bly,bly],
                               color='#555555', lw=1.1,
                               transform=fig.transFigure, clip_on=False))
    fig.text(xc, bly+0.010, label,
             ha='center', va='bottom', fontsize=11, color='#333333',
             transform=fig.transFigure)

# ── Colorbar ──────────────────────────────────────────────────────────────────
rp_left = (LEFT + pw + 0.38) / FW
cb_top  = (BOT + ph * 0.90) / FH
cb_bot  = (BOT + ph * 0.48) / FH
cbar_ax = fig.add_axes([rp_left, cb_bot, 0.20/FW, cb_top - cb_bot])
cb = plt.colorbar(mapper, cax=cbar_ax)
cb.set_label('log2 fold enrichment\n(% target / % background)', fontsize=9.5, labelpad=6)
cb.ax.tick_params(labelsize=9.5)
cb.set_ticks([1,2,3,4])

fig.text(rp_left + 0.10/FW, (BOT + ph*0.46)/FH,
         'circles: p < 0.001',
         ha='center', va='top', fontsize=9, color='#555555',
         transform=fig.transFigure)

# ── TF family legend ──────────────────────────────────────────────────────────
fam_order = ['WRKY','NAC','MYB','AP2/ERF','TCP','DOF','PIF','bZIP','HSF']
patches = [mpatches.Patch(facecolor=FAM_COLOR[f], label=f, edgecolor='none')
           for f in fam_order]
fig.legend(handles=patches, title='TF family', title_fontsize=10,
           loc='lower right',
           bbox_to_anchor=((FW-0.22)/FW, (BOT+0.18)/FH),
           bbox_transform=fig.transFigure,
           frameon=False, fontsize=9.5,
           handlelength=1.2, handleheight=0.9, handletextpad=0.5,
           borderpad=0, labelspacing=0.35)

# ── Title ─────────────────────────────────────────────────────────────────────
fig.text(0.5, 0.988,
         'Promoter motif enrichment across NLR and PRR immune receptor gene sets',
         ha='center', va='top', fontsize=13, fontweight='bold')
fig.text(0.5, 0.966,
         'JASPAR2024 plant TF motifs  |  2 kb upstream  |  HOMER hypergeometric test',
         ha='center', va='top', fontsize=10, color='#555555')

OUT = f'{DRAFT}/homer_motif_dotplot_v6.png'
plt.savefig(OUT, dpi=200, bbox_inches='tight', facecolor='white')
plt.savefig(OUT.replace('.png','.pdf'), bbox_inches='tight', facecolor='white')
print('Saved:', OUT)

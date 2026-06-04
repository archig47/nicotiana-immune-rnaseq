import re, os
import numpy as np
import pandas as pd
from Bio import SeqIO
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.stats import gaussian_kde

font_manager.fontManager.addfont('/users/fyp/fyp5/project/genome/Arial.ttf')
plt.rcParams.update({'font.family':'Arial','font.size':11,'pdf.fonttype':42,'ps.fonttype':42})

JASPAR_HOMER = '/tmp/JASPAR2024_plants_homer.motif'
FA_DIR = '/users/fyp/fyp5/project/PRJNA945175/promoters/timepoint_model_may2026_v2'
DRAFT  = '/users/fyp/fyp5/project/PRJNA945175/figures/draft'

FAM_COLOR = {
    'WRKY':'#1D3557','NAC':'#2A9D8F','MYB':'#E9C46A',
    'AP2/ERF':'#E76F51','TCP':'#F4A261',
    'DOF':'#8338EC','PIF':'#6C757D','bZIP':'#AAAAAA','HSF':'#C77DFF',
}

# (display_label, jaspar_key, family, fasta_file, gene_set_label)
TARGETS = [
    ('TCP2',      'TCP2',   'TCP',     'all_prr_promoters_2kb.fa', 'PRR'),
    ('WRKY42',    'WRKY42', 'WRKY',    'all_prr_promoters_2kb.fa', 'PRR'),
    ('MYB30',     'MYB30',  'MYB',     'deg_nlr_promoters_2kb.fa', 'DEG NLR'),
    ('DREB1/CBF', 'DREB1E', 'AP2/ERF', 'all_prr_promoters_2kb.fa', 'PRR'),
    ('DOF4.2',    'DOF4.2', 'DOF',     'all_prr_promoters_2kb.fa', 'PRR'),
]

def load_pwm(homer_file, tf_key):
    cur_name, cur_rows, capturing, thresh = None, [], False, 0
    with open(homer_file) as f:
        for line in f:
            if line.startswith('>'):
                if capturing and cur_rows:
                    return np.array(cur_rows), thresh
                parts = line.strip().split('\t')
                raw = parts[1].replace('/JASPAR2024_plants','')
                m = re.match(r'^([A-Za-z0-9][A-Za-z0-9\.\-]*)', raw)
                cur_name = m.group(1) if m else raw
                thresh = float(parts[2])
                capturing = (cur_name == tf_key)
                cur_rows = []
            elif capturing and line.strip():
                cur_rows.append([float(x) for x in line.strip().split('\t')])
    if capturing and cur_rows:
        return np.array(cur_rows), thresh
    return None, None

BASE_IDX = {b:i for b,i in zip('ACGTacgt',[0,1,2,3,0,1,2,3])}
SEQ_LEN = 2000

def scan_fasta(fa_path, pwm, thresh):
    lo = np.log2(np.clip(pwm, 1e-9, 1) / 0.25)
    L = len(lo)
    dists = []
    for rec in SeqIO.parse(fa_path, 'fasta'):
        seq = str(rec.seq).upper()
        for i in range(len(seq) - L + 1):
            s = 0.0; ok = True
            for j, b in enumerate(seq[i:i+L]):
                idx = BASE_IDX.get(b, -1)
                if idx < 0: ok = False; break
                s += lo[j, idx]
            if ok and s >= thresh:
                centre = i + L // 2
                dists.append(centre - SEQ_LEN)  # negative = upstream of TSS
    return dists

results = {}
for label, jkey, fam, fa_name, gs_label in TARGETS:
    pwm, thresh = load_pwm(JASPAR_HOMER, jkey)
    fa = f'{FA_DIR}/{fa_name}'
    dists = scan_fasta(fa, pwm, thresh)
    n_seq = sum(1 for _ in SeqIO.parse(fa, 'fasta'))
    print(f'{label}: {len(dists)} hits in {n_seq} {gs_label} sequences')
    results[label] = (fam, dists, gs_label)

# ── Plot ──────────────────────────────────────────────────────────────────────
FW, FH = 10.5, 6.5
fig, ax = plt.subplots(figsize=(FW, FH))
fig.subplots_adjust(left=0.10, right=0.70, top=0.87, bottom=0.13)

xs = np.linspace(-2000, 0, 600)
bw = 0.15

legend_lines = []
for label, (fam, dists, gs) in results.items():
    if len(dists) < 10:
        print(f'  skipping {label}: too few hits')
        continue
    kde = gaussian_kde(dists, bw_method=bw)
    dens = kde(xs)
    c = FAM_COLOR.get(fam, '#888888')
    ls = '--' if gs != 'PRR' else '-'
    line, = ax.plot(xs, dens, color=c, lw=2.2, linestyle=ls, zorder=3,
                    label=f'{label}  ({gs})')
    ax.fill_between(xs, dens, alpha=0.10, color=c, zorder=2)

# Reference lines and TSS
for xr in [-1500, -1000, -500]:
    ax.axvline(xr, color='#DDDDDD', lw=0.8, ls='--', zorder=1)
ax.axvline(0, color='#888888', lw=1.3, zorder=4)

ax.set_xlabel('Distance from TSS (bp)', fontsize=12)
ax.set_ylabel('Motif hit density', fontsize=12)
ax.set_xlim(-2000, 100)
ax.set_xticks([-2000,-1500,-1000,-500,0])
ax.set_xticklabels(['-2000','-1500','-1000','-500','0'], fontsize=10.5)
ax.tick_params(labelsize=10.5)
for sp in ['top','right']: ax.spines[sp].set_visible(False)

# TSS label — placed after ylim is determined
ymax = ax.get_ylim()[1]
ax.text(8, ymax*0.97, 'TSS', va='top', fontsize=10, color='#555555')

ax.legend(title='TF motif  (gene set)', title_fontsize=10, fontsize=10,
          loc='upper left', bbox_to_anchor=(1.02, 1.0),
          frameon=False, labelspacing=0.6,
          handlelength=1.5)

fig.text(0.40, 0.96,
         'Positional distribution of TF binding motif hits in promoter sequences',
         ha='center', va='top', fontsize=13, fontweight='bold')
fig.text(0.40, 0.935,
         '2 kb upstream regions  |  JASPAR2024 PWMs  |  KDE  ·  dashed = DEG NLR promoters',
         ha='center', va='top', fontsize=9.5, color='#555555')

OUT = f'{DRAFT}/homer_position_v3.png'
plt.savefig(OUT, dpi=200, bbox_inches='tight', facecolor='white')
plt.savefig(OUT.replace('.png','.pdf'), bbox_inches='tight', facecolor='white')
print('Saved:', OUT)

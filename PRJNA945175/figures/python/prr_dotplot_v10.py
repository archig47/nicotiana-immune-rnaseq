import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.font_manager as fm
fm.fontManager.addfont('/users/fyp/fyp5/project/genome/Arial.ttf')
fm.fontManager.addfont('/users/fyp/fyp5/project/genome/Arial Bold.ttf')
matplotlib.rcParams['font.family'] = 'Arial'
matplotlib.rcParams['font.sans-serif'] = ['Arial']
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import re, os, glob
from collections import Counter, defaultdict

base  = "/users/fyp/fyp5/project/PRJNA945175/deseq2_clean"
res_d = os.path.join(base, "results")
nlr_d = os.path.join(base, "nlr_prr")
out_v = 10

# ── 1. PRR annotation ────────────────────────────────────────────────────────
xl = pd.read_excel(os.path.join(nlr_d, "nlr_prr_full_landscape.xlsx"),
                   sheet_name="PRR_landscape")

ACCESSORY = {'EGF','PAN','SDOM','NEW','RCC1','TNFR','GP_PDE','UNIDENTIFIED'}

def get_ecto(rules):
    if not isinstance(rules, str): return 'Other'
    tokens = rules.split(',')[1:]
    primary = next((t for t in tokens if t not in ACCESSORY), None)
    if not primary: return 'Other'
    if primary == 'LRR':                                 return 'LRR'
    if primary == 'LysM':                                return 'LysM'
    if primary == 'WAK':                                 return 'WAK'
    if primary in ('MAL','SPARK'):                       return 'Malectin/CrRLK'
    if primary in ('BLEC','LLEC','CLEC','GLEC','GNK2'): return 'Lectin'
    return 'Other'

def get_backbone(rules):
    if not isinstance(rules, str): return 'Other'
    return rules.split(',')[0]

xl['ecto']     = xl['matched_rules'].apply(get_ecto)
xl['backbone'] = xl['matched_rules'].apply(get_backbone)
prr_genes = set(xl['gene_id'])
prr_anno  = xl.set_index('gene_id')[['ecto','backbone']]

# ── 2. Load results ──────────────────────────────────────────────────────────
frames = []
for f in glob.glob(os.path.join(res_d, "results_*.csv")):
    d = pd.read_csv(f)
    parts = os.path.splitext(os.path.basename(f))[0].split('_')
    d['contrast'] = f"{parts[1]}_{parts[2]}"
    frames.append(d)
all_res = pd.concat(frames, ignore_index=True)

prr_degs = all_res[
    all_res['gene_id'].isin(prr_genes) &
    all_res['padj'].notna() &
    (all_res['padj'] < 0.05) &
    (all_res['log2FoldChange'].abs() > 1)
].copy()

# ── 3. BLAST data (needed for receptor filter before gene selection) ──────────
blast_path = "/users/fyp/fyp5/project/genome/ath_blast/prr_vs_ath_swissprot.txt"
blast_df = pd.read_csv(blast_path, sep='\t', header=None,
                       names=['qseqid','sseqid','pident','length','evalue','bitscore','stitle'])
blast_df['gene_id'] = blast_df['qseqid'].str.replace('-mRNA', '', regex=False)
blast_best = (blast_df.sort_values(['evalue','bitscore'], ascending=[True, False])
                       .drop_duplicates('gene_id', keep='first')
                       .set_index('gene_id'))

# Descriptions that identify non-receptor annotation artifacts to skip
NON_RECEPTOR = re.compile(
    r'ring-h2|ubiquitin|nitrate transporter|plasmodesmata|'
    r'ep1-like glycoprotein|gpi-anchored|transcription factor|'
    r'cytochrome|ribosom|histone|aquaporin|'
    r'xyloglucan|endotransglucosylase|hydrolase|'
    r'extensin|dehydrin|lipid transfer',
    re.IGNORECASE
)

def is_receptor(gene_id):
    if gene_id not in blast_best.index:
        return True   # no BLAST hit → assume genuine PRR (in landscape by domain rules)
    desc = blast_best.loc[gene_id, 'stitle']
    return not bool(NON_RECEPTOR.search(str(desc)))

# ── 4. Gene selection: top-5 genuine receptors per ecto class ────────────────
ECTO_ORDER = ['LRR', 'Lectin', 'WAK', 'Malectin/CrRLK', 'LysM', 'Other']
max_lfc_all = prr_degs.groupby('gene_id')['log2FoldChange'].apply(lambda x: x.abs().max())
n_sig_all   = prr_degs.groupby('gene_id')['contrast'].nunique()
score_all   = max_lfc_all * np.log1p(n_sig_all)

top_genes = []
for ecto in ECTO_ORDER:
    ecto_genes  = prr_anno[prr_anno['ecto'] == ecto].index
    ecto_scores = score_all.reindex(ecto_genes).dropna().sort_values(ascending=False)
    n_picked = 0
    for gene_id in ecto_scores.index:
        if n_picked == 5: break
        if is_receptor(gene_id):
            top_genes.append(gene_id)
            n_picked += 1

print(f"Selected {len(top_genes)} genes across {len(ECTO_ORDER)} ecto classes")

# ── 5. Sort: ecto order (LRR first) then score descending ────────────────────
ecto_rank = {e: i for i, e in enumerate(ECTO_ORDER)}
anno_top = prr_anno.reindex(top_genes).fillna('Other').copy()
anno_top['score']     = score_all.reindex(top_genes)
anno_top['ecto_rank'] = anno_top['ecto'].map(ecto_rank).fillna(5).astype(int)
anno_top = anno_top.sort_values(['ecto_rank', 'score'], ascending=[True, False])
gene_order = anno_top.index.tolist()
gene_idx   = {g: i for i, g in enumerate(gene_order)}

# ── 6. Column order: TIMEPOINT-GROUPED ──────────────────────────────────────
col_order = [
    'HCRV_1dpi', 'TSWV_1dpi', 'TH_1dpi',
    'HCRV_7dpi', 'TSWV_7dpi', 'TH_7dpi',
    'HCRV_14dpi','TSWV_14dpi','TH_14dpi'
]

# ── 7. Significant observations only ────────────────────────────────────────
sig_data = prr_degs[
    prr_degs['gene_id'].isin(gene_order) &
    prr_degs['contrast'].isin(col_order)
].copy()
sig_data['direction'] = np.where(sig_data['log2FoldChange'] > 0, 'up', 'down')
sig_data['y'] = sig_data['gene_id'].map(gene_idx)
sig_data['x'] = sig_data['contrast'].map({c: i for i, c in enumerate(col_order)})

print(f"Genes: {len(gene_order)}, Dots: {len(sig_data)}, "
      f"Mean fill: {len(sig_data)/len(gene_order):.1f}/9")

# ── 8. Short gene labels + supplementary table ───────────────────────────────
def parse_gn(stitle):
    m = re.search(r'GN=(\S+)', stitle)
    if not m: return None
    gn = m.group(1)
    return None if re.match(r'^AT\dG\d{5,6}$', gn, re.IGNORECASE) else gn

def parse_atid(stitle):
    m = re.search(r'GN=(AT\dG\d{5,6})', stitle, re.IGNORECASE)
    return m.group(1) if m else None

def parse_full_desc(stitle):
    m = re.match(r'sp\|\S+\|\S+_ARATH\s+(.+?)\s+OS=', stitle)
    return m.group(1) if m else stitle

BB_SHORT = {'RLK': 'RLK', 'RLP': 'RLP', 'Secreted': 'Sec'}

def effective_backbone(gene_id):
    """Use annotation backbone; override Secreted→RLK when BLAST ortholog is a kinase."""
    bb = anno_top.loc[gene_id, 'backbone'] if gene_id in anno_top.index else 'RLK'
    if bb == 'Secreted' and gene_id in blast_best.index:
        if 'kinase' in blast_best.loc[gene_id, 'stitle'].lower():
            return 'RLK'
    return bb

# Pre-pass: count how many top genes share each named Arabidopsis gene
ath_counts = Counter()
for g in gene_order:
    if g in blast_best.index:
        gn = parse_gn(blast_best.loc[g, 'stitle'])
        if gn:
            ath_counts[gn] += 1

ath_seen      = Counter()
no_hit_counts = defaultdict(int)

def make_short_label(gene_id):
    bb  = effective_backbone(gene_id)
    bbs = BB_SHORT.get(bb, bb[:3])

    if gene_id in blast_best.index:
        stitle = blast_best.loc[gene_id, 'stitle']
        gn     = parse_gn(stitle)
        atid   = parse_atid(stitle)
        if gn:
            ath_seen[gn] += 1
            suffix = f".{ath_seen[gn]}" if ath_counts[gn] > 1 else ""
            return f"Nb{gn}{suffix} [{bbs}]"
        elif atid:
            return f"{atid} [{bbs}]"
    # No BLAST hit at all
    ecto = anno_top.loc[gene_id, 'ecto'] if gene_id in anno_top.index else 'Other'
    no_hit_counts[(ecto, bbs)] += 1
    n = no_hit_counts[(ecto, bbs)]
    return f"{ecto}-{bbs} #{n}"

gene_labels = [make_short_label(g) for g in gene_order]
print("\nGene labels:")
for g, lbl in zip(gene_order, gene_labels):
    print(f"  {g:25s}  {lbl}")

# Supplementary table
def blast_row(gene_id):
    if gene_id not in blast_best.index:
        return ('', '', '', '', '')
    row    = blast_best.loc[gene_id]
    stitle = row['stitle']
    return (
        row['sseqid'],
        parse_gn(stitle) or parse_atid(stitle) or '',
        f"{row['pident']:.1f}",
        f"{row['evalue']:.2e}",
        parse_full_desc(stitle)
    )

supp_rows = []
for g, lbl in zip(gene_order, gene_labels):
    ecto = anno_top.loc[g, 'ecto']
    bb   = anno_top.loc[g, 'backbone']
    acc, gn, pid, ev, desc = blast_row(g)
    supp_rows.append({
        'figure_label':   lbl,
        'gene_id':        g,
        'ecto_class':     ecto,
        'backbone':       bb,
        'ath_accession':  acc,
        'ath_gene':       gn,
        'pct_identity':   pid,
        'evalue':         ev,
        'full_description': desc
    })

supp_df = pd.DataFrame(supp_rows)
supp_path = os.path.join(base, f'figures/prr_dotplot_v{out_v}_gene_table.csv')
supp_df.to_csv(supp_path, index=False)
print(f"\nSupplementary table: {supp_path}")
print(supp_df[['figure_label','ath_gene','full_description']].to_string(index=False))

# ── 9. Figure layout ─────────────────────────────────────────────────────────
n_genes = len(gene_order)
n_cols  = 9

fig_w   = 11.0
fig_h   = n_genes * 0.38 + 4.0

fig = plt.figure(figsize=(fig_w, fig_h), facecolor='white')

top_frac  = 2.8 / fig_h
bot_frac  = 1.2 / fig_h
left_frac = 0.21
right_frac = 0.70

ax = fig.add_axes([left_frac, bot_frac,
                   right_frac - left_frac,
                   1 - top_frac - bot_frac])
ax.set_facecolor('white')

# ── 10. Background ecto bands ─────────────────────────────────────────────────
ecto_bg = {
    'LRR':            '#EFF5FF',
    'Lectin':         '#FFF6EE',
    'WAK':            '#FFEEEE',
    'Malectin/CrRLK': '#F6EEFF',
    'LysM':           '#EEFFF3',
    'Other':          '#F5F5F5'
}
ecto_cols = {
    'LRR':            '#4E79A7',
    'Lectin':         '#F28E2B',
    'WAK':            '#E15759',
    'Malectin/CrRLK': '#B07AA1',
    'LysM':           '#59A14F',
    'Other':          '#BAB0AC'
}

for ecto, grp in anno_top.groupby('ecto', sort=False):
    ys = [gene_idx[g] for g in grp.index]
    ax.axhspan(min(ys) - 0.5, max(ys) + 0.5,
               color=ecto_bg.get(ecto, '#F5F5F5'), alpha=0.6, zorder=0)

for y in range(n_genes):
    ax.axhline(y, color='#EEEEEE', linewidth=0.4, zorder=0)
for x in range(n_cols):
    ax.axvline(x, color='#EEEEEE', linewidth=0.4, zorder=0)
for x in [2.5, 5.5]:
    ax.axvline(x, color='#AAAAAA', linewidth=1.2, zorder=1)

# ── 11. Draw dots ─────────────────────────────────────────────────────────────
UP_COL   = '#6B2D8B'
DOWN_COL = '#0D7377'

def lfc_size(lfc_abs):
    return np.clip((lfc_abs - 1) / 7, 0, 1) ** 0.6 * 220 + 30

up_mask   = sig_data['direction'] == 'up'
down_mask = ~up_mask

ax.scatter(sig_data.loc[up_mask,   'x'], sig_data.loc[up_mask,   'y'],
           s=[lfc_size(v) for v in sig_data.loc[up_mask,   'log2FoldChange'].abs()],
           color=UP_COL,   alpha=0.85, zorder=3, linewidths=0)
ax.scatter(sig_data.loc[down_mask, 'x'], sig_data.loc[down_mask, 'y'],
           s=[lfc_size(v) for v in sig_data.loc[down_mask, 'log2FoldChange'].abs()],
           color=DOWN_COL, alpha=0.85, zorder=3, linewidths=0)

# ── 12. Axes ──────────────────────────────────────────────────────────────────
ax.set_xlim(-0.5, n_cols - 0.5)
ax.set_ylim(-0.5, n_genes - 0.5)
ax.invert_yaxis()

virus_xlabels = ['HCRV', 'TSWV', 'TH'] * 3
ax.set_xticks(range(n_cols))
ax.set_xticklabels(virus_xlabels, fontsize=9, color='#555555')
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
ax.tick_params(axis='x', which='both', length=0)

for label, xfrac in [('1 dpi', 1/9), ('7 dpi', 4/9), ('14 dpi', 7/9)]:
    ax.text(xfrac, 1.10, label,
            ha='center', va='bottom', fontsize=11, fontweight='bold',
            color='#1A202C', transform=ax.transAxes, clip_on=False)

for lo_frac, hi_frac in [(0/9, 3/9), (3/9, 6/9), (6/9, 9/9)]:
    ax.annotate('', xy=(hi_frac - 0.005, 1.07), xytext=(lo_frac + 0.005, 1.07),
                xycoords=ax.transAxes, textcoords=ax.transAxes,
                annotation_clip=False,
                arrowprops=dict(arrowstyle='-', color='#BBBBBB', lw=1.0))

ax.set_yticks(range(n_genes))
ax.set_yticklabels(gene_labels, fontsize=11, color='#1A202C')
ax.yaxis.tick_right()
ax.yaxis.set_label_position('right')
ax.tick_params(axis='y', which='both', length=0)

for spine in ax.spines.values():
    spine.set_visible(False)

# ── 13. Left ecto strip: one continuous rectangle per group ──────────────────
STRIP_X0 = -1.25
STRIP_W  = 0.50

for ecto, grp in anno_top.groupby('ecto', sort=False):
    ys    = [gene_idx[g] for g in grp.index]
    y_min = min(ys) - 0.5
    y_max = max(ys) + 0.5
    rect  = mpatches.Rectangle(
        (STRIP_X0, y_min), STRIP_W, y_max - y_min,
        linewidth=0, facecolor=ecto_cols.get(ecto, '#BAB0AC'),
        zorder=4, clip_on=False
    )
    ax.add_patch(rect)

for ecto, grp in anno_top.groupby('ecto', sort=False):
    ys    = [gene_idx[g] for g in grp.index]
    y_mid = np.mean(ys)
    ax.text(STRIP_X0 - 0.15, y_mid, ecto,
            ha='right', va='center', fontsize=9,
            fontweight='bold', color=ecto_cols.get(ecto, '#333333'),
            clip_on=False)

# ── 14. Legend — below axes ───────────────────────────────────────────────────
UP_PATCH   = mpatches.Patch(facecolor=UP_COL,   label='Upregulated')
DOWN_PATCH = mpatches.Patch(facecolor=DOWN_COL, label='Downregulated')

size_handles = []
for lfc_val, label in [(2, '|log2FC| = 2'), (4, '= 4'), (6, '= 6')]:
    ms = np.sqrt(lfc_size(lfc_val)) * 0.55
    size_handles.append(
        Line2D([0],[0], marker='o', color='w', markerfacecolor='#999999',
               markersize=ms, label=label)
    )

fig.legend(
    handles=[UP_PATCH, DOWN_PATCH],
    title='Direction', title_fontsize=8, fontsize=7.5,
    frameon=True, framealpha=0.95, edgecolor='#DDDDDD',
    loc='lower left',
    bbox_to_anchor=(left_frac, 0.01), bbox_transform=fig.transFigure,
    ncol=2, handlelength=0.8, handletextpad=0.4, borderpad=0.6, columnspacing=0.8
)
fig.legend(
    handles=size_handles,
    title='Magnitude', title_fontsize=8, fontsize=7.5,
    frameon=True, framealpha=0.95, edgecolor='#DDDDDD',
    loc='lower center',
    bbox_to_anchor=(0.50, 0.01), bbox_transform=fig.transFigure,
    ncol=3, handlelength=0.8, handletextpad=0.4, borderpad=0.6, columnspacing=0.8
)

# ── 15. Title & subtitle ──────────────────────────────────────────────────────
title_y    = 1 - (0.55 / fig_h)
subtitle_y = 1 - (1.15 / fig_h)

fig.text(0.02, title_y,
         'Top PRR DEGs by ecto class — response magnitude across 9 contrasts',
         fontsize=12, fontweight='bold', color='#1A202C', va='top')
fig.text(0.02, subtitle_y,
         'N. benthamiana  |  padj < 0.05, |log2FC| > 1  |  '
         'top 5 genuine receptors per ecto class by max|log2FC| × ln(1 + n sig. contrasts)  |  '
         'labels: NbGENE [backbone]; see supplementary table for full annotations',
         fontsize=8, color='#718096', va='top', style='italic')

# ── 16. Save ──────────────────────────────────────────────────────────────────
out_png = os.path.join(base, f'figures/prr_gene_dotplot_v{out_v}.png')
fig.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
print(f"\nSaved: {out_png}")
plt.close()

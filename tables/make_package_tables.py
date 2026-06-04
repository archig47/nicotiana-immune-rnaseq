#!/usr/bin/env python3
"""
Generate styled package tables for thesis appendix.
Table A.1: R packages  |  Table A.2: Python packages
Style: Arial 12pt, #2596be header, alternating #fce8ec / white rows.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import textwrap, os

# ── Colours ──────────────────────────────────────────────────────────────────
HEADER_BG  = '#2596be'
HEADER_FG  = '#ffffff'
ODD_ROW    = '#fce8ec'   # alternating pink (rows 1, 3, 5 …)
EVEN_ROW   = '#ffffff'
TEXT_FG    = '#000000'
GRID_COLOR = '#aaaaaa'
OUTER_CLR  = '#555555'

# ── Data ─────────────────────────────────────────────────────────────────────
HEADERS  = ['Package', 'Version', 'Purpose', 'Citation']
COL_W    = [1.80, 0.90, 2.90, 1.85]   # inches; total ≈ 7.45 in (fits portrait A4)

R_DATA = [
    ('DESeq2',       '1.50.2', 'Differential expression analysis and rlog normalisation',  'Love et al., 2014'),
    ('ggplot2',      '3.5.2',  'Data visualisation and figure generation',                 'Wickham, 2016'),
    ('pheatmap',     '1.0.13', 'Heatmap generation for NLR and PRR expression profiles',  'Kolde, 2019'),
    ('ggrepel',      '0.9.8',  'Non-overlapping gene labels on scatter plots',             'Slowikowski, 2024'),
    ('ggvenn',       '0.1.19', 'Venn diagram visualisation of DEG overlaps',               'Yan, 2023'),
    ('dplyr',        '1.1.4',  'Data manipulation and filtering',                          'Wickham et al., 2023'),
    ('tidyr',        '1.3.1',  'Data reshaping and pivoting',                              'Wickham et al., 2023'),
    ('readr',        '2.1.5',  'Reading delimited text files',                             'Wickham et al., 2024'),
    ('readxl',       '1.4.5',  'Reading Excel input files',                                'Wickham & Bryan, 2023'),
    ('patchwork',    '1.3.2',  'Combining multiple ggplot2 panels into single figures',    'Pedersen, 2024'),
    ('scales',       '1.4.0',  'Axis scaling and colour transformations for ggplot2',      'Wickham & Seidel, 2022'),
    ('RColorBrewer', '1.1.3',  'Colour palettes for heatmaps and figures',                 'Neuwirth, 2022'),
]

PY_DATA = [
    ('matplotlib',  '3.9.2', 'Figure generation and data visualisation',                  'Hunter, 2007'),
    ('pandas',      '2.2.2', 'Data manipulation and tabular analysis',                    'McKinney, 2010'),
    ('numpy',       '1.26.4','Numerical computation and array operations',                'Harris et al., 2020'),
    ('adjustText',  '1.3.0', 'Non-overlapping text label placement on scatter plots',     'Flyamer, 2023'),
    ('openpyxl',    '3.1.5', 'Reading and writing Excel files for HOMER result export',   'Gazoni & Clark, 2024'),
    ('scipy',       '1.13.1','Kernel density estimation for motif positional plots',      'Virtanen et al., 2020'),
    ('logomaker',   '0.8.7', 'TF binding site sequence logo generation',                  'Tareen & Kinney, 2020'),
]

# ── Drawing ───────────────────────────────────────────────────────────────────

def draw_table(headers, data, col_widths_in, output_path, dpi=250):
    FS_HDR   = 12      # header font size (pt)
    FS_BODY  = 12      # body font size (pt)
    CPR_IN   = 7.2     # approximate Arial 12pt chars per inch
    LH_IN    = 0.245   # line height in inches at 12pt
    PAD_IN   = 0.16    # vertical padding (top+bottom) per cell
    MIN_RH   = 0.44    # minimum row height
    HDR_H    = 0.46    # header row height
    LEFT_PAD = 0.13    # text left margin inside each cell

    # Pre-wrap cell text
    wrapped, row_heights = [], []
    for row in data:
        wr, max_lines = [], 1
        for cell, cw in zip(row, col_widths_in):
            nchars = max(10, int(cw * CPR_IN))
            wt     = textwrap.fill(str(cell), width=nchars)
            nl     = wt.count('\n') + 1
            max_lines = max(max_lines, nl)
            wr.append(wt)
        wrapped.append(wr)
        row_heights.append(max(MIN_RH, max_lines * LH_IN + PAD_IN))

    total_w = sum(col_widths_in)
    total_h = HDR_H + sum(row_heights)

    fig, ax = plt.subplots(figsize=(total_w, total_h))
    ax.set_xlim(0, total_w)
    ax.set_ylim(0, total_h)
    ax.axis('off')

    # ── Header ──
    y, x = total_h - HDR_H, 0
    for header, cw in zip(headers, col_widths_in):
        ax.add_patch(mpatches.Rectangle((x, y), cw, HDR_H,
                     linewidth=0, facecolor=HEADER_BG, zorder=2))
        ax.text(x + LEFT_PAD, y + HDR_H / 2, header,
                color=HEADER_FG, fontsize=FS_HDR, fontweight='bold',
                fontfamily='Arial', va='center', ha='left', zorder=3)
        x += cw

    # ── Data rows ──
    y_cur = total_h - HDR_H
    for ri, (wr_row, rh) in enumerate(zip(wrapped, row_heights)):
        y_cur -= rh
        bg = ODD_ROW if ri % 2 == 0 else EVEN_ROW   # row 1 = pink, row 2 = white …
        x = 0
        for cell_text, cw in zip(wr_row, col_widths_in):
            ax.add_patch(mpatches.Rectangle((x, y_cur), cw, rh,
                         linewidth=0, facecolor=bg, zorder=2))
            ax.text(x + LEFT_PAD, y_cur + rh / 2, cell_text,
                    color=TEXT_FG, fontsize=FS_BODY, fontfamily='Arial',
                    va='center', ha='left', zorder=3,
                    multialignment='left', linespacing=1.3)
            x += cw

    # ── Horizontal grid lines ──
    y_cur = total_h - HDR_H
    for rh in row_heights:
        ax.axhline(y_cur, color=GRID_COLOR, linewidth=0.5, zorder=4)
        y_cur -= rh

    # ── Vertical grid lines (internal only) ──
    x = 0
    for cw in col_widths_in[:-1]:
        x += cw
        ax.axvline(x, color=GRID_COLOR, linewidth=0.5, zorder=4)

    # ── Outer border ──
    ax.add_patch(mpatches.Rectangle((0, 0), total_w, total_h,
                 linewidth=1.2, edgecolor=OUTER_CLR, facecolor='none', zorder=5))

    plt.savefig(output_path, dpi=dpi, bbox_inches='tight',
                facecolor='white', pad_inches=0.06)
    plt.close()
    print(f"Saved: {output_path}")


if __name__ == '__main__':
    out = os.path.expanduser('~/Downloads')
    draw_table(HEADERS, R_DATA,  COL_W, os.path.join(out, 'table_A1_r_packages.png'))
    draw_table(HEADERS, PY_DATA, COL_W, os.path.join(out, 'table_A2_python_packages.png'))
    print("Done.")

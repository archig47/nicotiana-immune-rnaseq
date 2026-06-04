#!/usr/bin/env python3
"""
HOMER comparison design table for thesis appendix (Table A.4).
Style: Arial 11pt, #2596be header, alternating #fce8ec / white rows,
light-blue group column with bold labels shown once per group.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import textwrap, os

HEADER_BG  = '#2596be'
HEADER_FG  = '#ffffff'
ODD_ROW    = '#fce8ec'
EVEN_ROW   = '#ffffff'
GROUP_BG   = '#deeef5'
TEXT_FG    = '#000000'
GRID_COLOR = '#aaaaaa'
OUTER_CLR  = '#555555'

HEADERS = ['Group', 'Foreground gene set (n)', 'Background gene set (n)', 'Biological question']
COL_W   = [1.70, 1.80, 1.80, 2.85]   # total ≈ 8.15 in

# Raw data — plain strings, no manual \n
DATA = [
    ('Receptor identity',
     'All NLR promoters (n = 284)',
     'Non-DEG expressed genome (n = 15,738)',
     'TF motifs constitutively enriched in NLR genes relative to the expressed genome'),
    ('Receptor identity',
     'All PRR promoters (n = 1,252)',
     'Non-DEG expressed genome (n = 15,738)',
     'TF motifs constitutively enriched in PRR genes relative to the expressed genome'),
    ('Virus responsiveness',
     'DEG NLR promoters (n = 98)',
     'Non-DEG expressed NLRs (n = 88)',
     'TF motifs distinguishing virus-responsive NLRs from non-responsive NLR genes'),
    ('Virus responsiveness',
     'DEG PRR promoters (n = 666)',
     'Non-DEG expressed PRRs (n = 328)',
     'TF motifs distinguishing virus-responsive PRRs from non-responsive PRR genes'),
    ('Cross-class enrichment',
     'DEG NLR promoters (n = 98)',
     'DEG PRR promoters (n = 666)',
     'TF motifs enriched in virus-responsive NLRs relative to virus-responsive PRRs'),
    ('Cross-class enrichment',
     'DEG PRR promoters (n = 666)',
     'DEG NLR promoters (n = 98)',
     'TF motifs enriched in virus-responsive PRRs relative to virus-responsive NLRs (mirror comparison)'),
    ('Direction of change',
     'Upregulated NLR DEGs (n = 55)',
     'Downregulated NLR DEGs (n = 50)',
     'TF motifs preferentially associated with NLR upregulation versus downregulation'),
    ('Direction of change',
     'Upregulated PRR DEGs (n = 463)',
     'Downregulated PRR DEGs (n = 357)',
     'TF motifs preferentially associated with PRR upregulation versus downregulation'),
]

def draw_table(headers, data, col_widths_in, output_path, dpi=250):
    FS_HDR   = 11
    FS_BODY  = 11
    CPR_IN   = 7.5     # Arial 11pt chars per inch
    LH_IN    = 0.230
    PAD_IN   = 0.18
    MIN_RH   = 0.44
    HDR_H    = 0.44
    LEFT_PAD = 0.12

    # Pre-wrap every cell and compute row heights
    wrapped_rows = []
    row_heights  = []
    for row in data:
        wr, max_lines = [], 1
        for cell, cw in zip(row, col_widths_in):
            nchars = max(8, int(cw * CPR_IN))
            wt     = textwrap.fill(str(cell), width=nchars)
            nl     = wt.count('\n') + 1
            max_lines = max(max_lines, nl)
            wr.append(wt)
        wrapped_rows.append(wr)
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
    prev_group  = None
    group_ri    = 0
    y_cur       = total_h - HDR_H

    for ri, (wr_row, rh) in enumerate(zip(wrapped_rows, row_heights)):
        y_cur -= rh
        group = data[ri][0]
        if group != prev_group:
            group_ri = 0
            prev_group = group
        else:
            group_ri += 1

        bg = ODD_ROW if ri % 2 == 0 else EVEN_ROW
        x = 0

        for ci, (cell_text, cw) in enumerate(zip(wr_row, col_widths_in)):
            cell_bg = GROUP_BG if ci == 0 else bg
            ax.add_patch(mpatches.Rectangle((x, y_cur), cw, rh,
                         linewidth=0, facecolor=cell_bg, zorder=2))

            # Group column: only show text on first row of each group
            if ci == 0:
                if group_ri == 0:
                    ax.text(x + LEFT_PAD, y_cur + rh / 2, cell_text,
                            color=TEXT_FG, fontsize=FS_BODY, fontfamily='Arial',
                            fontweight='bold', va='center', ha='left', zorder=3,
                            multialignment='left', linespacing=1.3)
            else:
                ax.text(x + LEFT_PAD, y_cur + rh / 2, cell_text,
                        color=TEXT_FG, fontsize=FS_BODY, fontfamily='Arial',
                        va='center', ha='left', zorder=3,
                        multialignment='left', linespacing=1.3)
            x += cw

    # ── Horizontal row lines ──
    y_cur = total_h - HDR_H
    for rh in row_heights:
        ax.axhline(y_cur, color=GRID_COLOR, linewidth=0.5, zorder=4)
        y_cur -= rh

    # ── Thicker separator between groups ──
    y_cur = total_h - HDR_H
    prev_group = None
    for row, rh in zip(data, row_heights):
        g = row[0]
        if prev_group is not None and g != prev_group:
            ax.axhline(y_cur, color='#888888', linewidth=1.1, zorder=5)
        prev_group = g
        y_cur -= rh

    # ── Vertical column lines ──
    x = 0
    for cw in col_widths_in[:-1]:
        x += cw
        ax.axvline(x, color=GRID_COLOR, linewidth=0.5, zorder=4)

    # ── Outer border ──
    ax.add_patch(mpatches.Rectangle((0, 0), total_w, total_h,
                 linewidth=1.2, edgecolor=OUTER_CLR, facecolor='none', zorder=6))

    plt.savefig(output_path, dpi=dpi, bbox_inches='tight',
                facecolor='white', pad_inches=0.06)
    plt.close()
    print(f"Saved: {output_path}")


if __name__ == '__main__':
    import matplotlib.transforms
    out = os.path.expanduser('~/Downloads')
    draw_table(HEADERS, DATA, COL_W,
               os.path.join(out, 'table_A4_homer_design.png'))
    print("Done.")

#!/usr/bin/env python3
"""
Generate tool comparison/justification table for thesis appendix.
Style: Arial, #2596be header, alternating #fce8ec / white rows.
Note: wide table — recommended for landscape orientation in appendix.
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
TEXT_FG    = '#000000'
GRID_COLOR = '#aaaaaa'
OUTER_CLR  = '#555555'

HEADERS = ['Tool considered', 'Purpose', 'Tool selected', 'Rationale for selection']
COL_W   = [1.35, 1.85, 1.40, 3.50]   # total ≈ 8.1 in — use landscape in Word

DATA = [
    (
        'FastQC',
        'Per-sample sequencing quality control',
        'Falco v1.2.5',
        'Falco produces concordant results across all QC modules with a ~7× speed advantage '
        '(9 s vs ~2 min per file). FastQC was run in parallel on a representative sample to '
        'confirm concordance; the one disagreement (duplication FAIL vs WARN on R1) reflects '
        'a known threshold difference between tools, not a data discrepancy.',
    ),
    (
        'Trimmomatic',
        'Adapter trimming and quality filtering',
        'Trim Galore v0.6.11',
        'Trimmomatic requires manual specification of adapter sequences, introducing user error '
        'risk. It is Java-based and comparatively slow. Trim Galore auto-detects adapters, '
        'wraps Cutadapt with integrated paired-end handling, and separates QC and trimming '
        'outputs to support independent verification.',
    ),
    (
        'fastp',
        'Adapter trimming with integrated QC reporting',
        'Trim Galore v0.6.11',
        "fastp's primary advantages — polyG tail correction and overlap-based adapter trimming "
        '— are relevant for NovaSeq/NextSeq instruments but not for this dataset (HiSeq PE150 '
        'TruSeq chemistry). Trim Galore was preferred for auditability; output quality is '
        'expected to be equivalent for this library type.',
    ),
    (
        'Cutadapt',
        'Adapter trimming',
        'Trim Galore v0.6.11',
        'Trim Galore wraps Cutadapt internally and adds automatic adapter detection and '
        'paired-end handling. Using Cutadapt directly would require manual configuration of '
        'these features with no performance benefit.',
    ),
    (
        'Bowtie2',
        'Short-read alignment',
        'STAR v2.7.11b',
        'Bowtie2 is not splice-aware and cannot correctly map reads spanning exon-exon '
        'junctions — a fundamental requirement for RNA-seq analysis of eukaryotic transcriptomes.',
    ),
    (
        'HISAT2',
        'Splice-aware read alignment',
        'STAR v2.7.11b',
        'STAR provides superior sensitivity at splice junctions and in repetitive genomic '
        'regions, both critical for this dataset. NLR immune receptor genes cluster in '
        'repeat-rich regions of the N. benthamiana genome. Server RAM (313 GB) comfortably '
        'accommodated STAR\'s higher memory requirements (~50–60 GB for index generation).',
    ),
    (
        'Salmon',
        'Alignment-free transcript quantification',
        'featureCounts v2.1.1',
        'featureCounts operates on STAR-generated BAM files and supports the fractional '
        'multimapper counting strategy (-M --fraction) required for the highly repetitive '
        'N. benthamiana genome (78.87% repeat content). '
        'Salmon\'s pseudoalignment approach does not support this multimapper handling strategy.',
    ),
    (
        'edgeR',
        'Count-based differential expression analysis',
        'DESeq2 v1.50.2',
        'DESeq2\'s apeglm shrinkage estimator improves stability of low-magnitude fold-change '
        'estimates, particularly important given the batch-condition confound in this dataset '
        'where conservative estimates reduce spurious DEG calls. Both tools are appropriate '
        'for count-based RNA-seq; DESeq2 was preferred for this reason.',
    ),
    (
        'UMI-tools',
        'UMI extraction and deduplication',
        'Not applied',
        'The TruSeq stranded library protocol does not incorporate unique molecular identifiers '
        '(UMIs), confirmed from SRA metadata. Applying UMI-based processing to non-UMI data '
        'would silently corrupt count data by misinterpreting adapter sequences as barcodes.',
    ),
    (
        'Picard MarkDuplicates',
        'PCR duplicate removal',
        'Not applied',
        'Duplicate removal is not appropriate for standard bulk RNA-seq without UMIs. '
        'PCR duplicates cannot be distinguished from reads arising from highly expressed '
        'transcripts; deduplication would remove genuine biological signal and '
        'systematically bias differential expression results.',
    ),
    (
        'AME (MEME Suite)',
        'Known motif enrichment in a set of sequences',
        'HOMER v5.1',
        'AME uses rank-based statistics (Wilcoxon or Fisher) that are not designed for '
        'discrete pairwise foreground/background comparisons. HOMER\'s hypergeometric test '
        'directly quantifies enrichment of a motif in a foreground set relative to a '
        'defined background, with fold enrichment output — the appropriate design for '
        'comparing receptor class and DEG gene sets.',
    ),
    (
        'STREME (MEME Suite)',
        'De novo motif discovery in a set of sequences',
        'HOMER v5.1',
        'STREME discovers novel motifs from scratch and is not designed for systematic '
        'enrichment testing of known binding sites. This analysis used the JASPAR2024 CORE '
        'plants library (805 experimentally characterised TF binding matrices); known-motif '
        'enrichment with HOMER -mknown was the appropriate strategy.',
    ),
    (
        'PlantPAN',
        'Web-based plant promoter and cis-regulatory element analysis',
        'HOMER v5.1',
        'PlantPAN is a browser-based tool and cannot scale to genome-wide pairwise '
        'comparisons across multiple gene sets. Its cis-regulatory database is smaller '
        'and less comprehensively curated than JASPAR2024 CORE plants, and it offers '
        'no support for the custom foreground/background comparison design used here.',
    ),
]


def draw_table(headers, data, col_widths_in, output_path, dpi=250):
    FS_HDR   = 11
    FS_BODY  = 10
    CPR_IN   = 8.2     # Arial 10pt chars per inch
    LH_IN    = 0.205
    PAD_IN   = 0.14
    MIN_RH   = 0.38
    HDR_H    = 0.42
    LEFT_PAD = 0.11

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

    y, x = total_h - HDR_H, 0
    for header, cw in zip(headers, col_widths_in):
        ax.add_patch(mpatches.Rectangle((x, y), cw, HDR_H,
                     linewidth=0, facecolor=HEADER_BG, zorder=2))
        ax.text(x + LEFT_PAD, y + HDR_H / 2, header,
                color=HEADER_FG, fontsize=FS_HDR, fontweight='bold',
                fontfamily='Arial', va='center', ha='left', zorder=3)
        x += cw

    y_cur = total_h - HDR_H
    for ri, (wr_row, rh) in enumerate(zip(wrapped, row_heights)):
        y_cur -= rh
        bg = ODD_ROW if ri % 2 == 0 else EVEN_ROW
        x = 0
        for cell_text, cw in zip(wr_row, col_widths_in):
            ax.add_patch(mpatches.Rectangle((x, y_cur), cw, rh,
                         linewidth=0, facecolor=bg, zorder=2))
            ax.text(x + LEFT_PAD, y_cur + rh / 2, cell_text,
                    color=TEXT_FG, fontsize=FS_BODY, fontfamily='Arial',
                    va='center', ha='left', zorder=3,
                    multialignment='left', linespacing=1.3)
            x += cw

    y_cur = total_h - HDR_H
    for rh in row_heights:
        ax.axhline(y_cur, color=GRID_COLOR, linewidth=0.5, zorder=4)
        y_cur -= rh

    x = 0
    for cw in col_widths_in[:-1]:
        x += cw
        ax.axvline(x, color=GRID_COLOR, linewidth=0.5, zorder=4)

    ax.add_patch(mpatches.Rectangle((0, 0), total_w, total_h,
                 linewidth=1.2, edgecolor=OUTER_CLR, facecolor='none', zorder=5))

    plt.savefig(output_path, dpi=dpi, bbox_inches='tight',
                facecolor='white', pad_inches=0.06)
    plt.close()
    print(f"Saved: {output_path}")


if __name__ == '__main__':
    out = os.path.expanduser('~/Downloads')
    draw_table(HEADERS, DATA, COL_W,
               os.path.join(out, 'table_A3_tool_comparison.png'))
    print("Done.")

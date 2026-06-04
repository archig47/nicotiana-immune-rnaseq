import pandas as pd
from pathlib import Path

configfile: "config.yaml"

# ── Load samples ───────────────────────────────────────────────────────────────
samples_df = pd.read_csv(config["samples_file"])
samples_df = samples_df.set_index("sample_id", drop=False)
SAMPLES = samples_df["sample_id"].tolist()

def get_r1(wildcards):
    return samples_df.loc[wildcards.sample, "fastq_r1"]

def get_r2(wildcards):
    return samples_df.loc[wildcards.sample, "fastq_r2"]

# ── STAR index: use existing or build ─────────────────────────────────────────
BUILD_INDEX = not bool(config["genome"]["star_index"])
STAR_INDEX  = "results/star_index" if BUILD_INDEX else config["genome"]["star_index"]

# ── Final outputs ──────────────────────────────────────────────────────────────
rule all:
    input:
        "results/qc/raw/multiqc_report.html",
        "results/qc/trimmed/multiqc_report_trimmed.html",
        "results/counts/counts_matrix.tsv",
        "results/deseq2/DEGs_full.csv",
        "results/deseq2/DEGs_significant.csv",
        "results/deseq2/plots/pca_plot.pdf",
        "results/deseq2/plots/volcano_plot.pdf",
        "results/deseq2/plots/ma_plot.pdf",
        "results/deseq2/plots/heatmap_top50.pdf",
        "results/deseq2/deseq2_summary.txt"

# ── (Optional) Build STAR index ───────────────────────────────────────────────
if BUILD_INDEX:
    rule build_star_index:
        input:
            fasta = config["genome"]["fasta"],
            gtf   = config["genome"]["gtf"]
        output:
            directory("results/star_index")
        params:
            sa_index_nbases = config["star"]["genome_sa_index_nbases"],
            overhang        = config["star"]["sjdb_overhang"]
        threads: 16
        log: "logs/star_index.log"
        shell:
            """
            mkdir -p {output}
            STAR --runMode genomeGenerate \
                 --genomeDir {output} \
                 --genomeFastaFiles {input.fasta} \
                 --sjdbGTFfile {input.gtf} \
                 --genomeSAindexNbases {params.sa_index_nbases} \
                 --sjdbOverhang {params.overhang} \
                 --runThreadN {threads} \
                 > {log} 2>&1
            """

# ── Raw FastQC ─────────────────────────────────────────────────────────────────
rule raw_fastqc:
    input:
        r1 = get_r1,
        r2 = get_r2
    output:
        touch("results/qc/raw/.done_{sample}")
    log: "logs/fastqc/raw_{sample}.log"
    shell:
        """
        mkdir -p results/qc/raw
        fastqc {input.r1} {input.r2} -o results/qc/raw -t 4 > {log} 2>&1
        """

rule raw_multiqc:
    input:
        expand("results/qc/raw/.done_{sample}", sample=SAMPLES)
    output:
        "results/qc/raw/multiqc_report.html"
    log: "logs/multiqc_raw.log"
    shell:
        "multiqc results/qc/raw/ -o results/qc/raw --filename multiqc_report.html -f > {log} 2>&1"

# ── Trim Galore ────────────────────────────────────────────────────────────────
rule trim_galore:
    input:
        r1 = get_r1,
        r2 = get_r2
    output:
        r1 = "results/trimmed/{sample}_val_1.fq.gz",
        r2 = "results/trimmed/{sample}_val_2.fq.gz"
    params:
        quality    = config["trim_galore"]["quality"],
        length     = config["trim_galore"]["length"],
        stringency = config["trim_galore"]["stringency"]
    log: "logs/trim_galore/{sample}.log"
    shell:
        """
        trim_galore \
            --quality {params.quality} \
            --length {params.length} \
            --stringency {params.stringency} \
            --paired \
            --basename {wildcards.sample} \
            -o results/trimmed \
            {input.r1} {input.r2} \
            > {log} 2>&1
        """

# ── Post-trim FastQC ───────────────────────────────────────────────────────────
rule post_trim_fastqc:
    input:
        r1 = "results/trimmed/{sample}_val_1.fq.gz",
        r2 = "results/trimmed/{sample}_val_2.fq.gz"
    output:
        touch("results/qc/trimmed/.done_{sample}")
    log: "logs/fastqc/trimmed_{sample}.log"
    shell:
        """
        mkdir -p results/qc/trimmed
        fastqc {input.r1} {input.r2} -o results/qc/trimmed -t 4 > {log} 2>&1
        """

rule post_trim_multiqc:
    input:
        expand("results/qc/trimmed/.done_{sample}", sample=SAMPLES)
    output:
        "results/qc/trimmed/multiqc_report_trimmed.html"
    log: "logs/multiqc_trimmed.log"
    shell:
        """
        multiqc results/trimmed/ results/qc/trimmed/ \
            -o results/qc/trimmed \
            --filename multiqc_report_trimmed.html \
            -f > {log} 2>&1
        """

# ── STAR Alignment ─────────────────────────────────────────────────────────────
rule star_align:
    input:
        r1    = "results/trimmed/{sample}_val_1.fq.gz",
        r2    = "results/trimmed/{sample}_val_2.fq.gz",
        index = STAR_INDEX
    output:
        bam       = "results/alignment/{sample}/Aligned.sortedByCoord.out.bam",
        log_final = "results/alignment/{sample}/Log.final.out",
        counts    = "results/alignment/{sample}/ReadsPerGene.out.tab"
    params:
        multimapNmax = config["star"]["multimapNmax"],
        overhang     = config["star"]["sjdb_overhang"],
        gtf          = config["genome"]["gtf"]
    threads: 16
    log: "logs/star/{sample}.log"
    shell:
        """
        mkdir -p results/alignment/{wildcards.sample}
        STAR \
            --genomeDir {input.index} \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat \
            --outFileNamePrefix results/alignment/{wildcards.sample}/ \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS NM \
            --outFilterMultimapNmax {params.multimapNmax} \
            --sjdbGTFfile {params.gtf} \
            --sjdbOverhang {params.overhang} \
            --quantMode GeneCounts \
            --runThreadN {threads} \
            > {log} 2>&1
        """

rule index_bam:
    input:  "results/alignment/{sample}/Aligned.sortedByCoord.out.bam"
    output: "results/alignment/{sample}/Aligned.sortedByCoord.out.bam.bai"
    shell:  "samtools index {input}"

# ── featureCounts ──────────────────────────────────────────────────────────────
rule featurecounts:
    input:
        bams    = expand("results/alignment/{sample}/Aligned.sortedByCoord.out.bam", sample=SAMPLES),
        bai     = expand("results/alignment/{sample}/Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES),
        gtf     = config["genome"]["gtf"]
    output:
        counts  = "results/counts/counts_matrix.tsv",
        summary = "results/counts/counts_matrix.tsv.summary"
    params:
        strandedness = config["featurecounts"]["strandedness"],
        feature_type = config["featurecounts"]["feature_type"],
        attribute    = config["featurecounts"]["attribute"]
    threads: 8
    log: "logs/featurecounts.log"
    shell:
        """
        mkdir -p results/counts
        featureCounts \
            -T {threads} \
            -t {params.feature_type} \
            -g {params.attribute} \
            -s {params.strandedness} \
            -p -B -C \
            -a {input.gtf} \
            -o {output.counts} \
            {input.bams} \
            > {log} 2>&1
        """

# ── DESeq2 ─────────────────────────────────────────────────────────────────────
rule deseq2:
    input:
        counts  = "results/counts/counts_matrix.tsv",
        samples = config["samples_file"]
    output:
        full        = "results/deseq2/DEGs_full.csv",
        significant = "results/deseq2/DEGs_significant.csv",
        pca         = "results/deseq2/plots/pca_plot.pdf",
        volcano     = "results/deseq2/plots/volcano_plot.pdf",
        ma          = "results/deseq2/plots/ma_plot.pdf",
        heatmap     = "results/deseq2/plots/heatmap_top50.pdf",
        summary     = "results/deseq2/deseq2_summary.txt"
    log: "logs/deseq2.log"
    shell:
        """
        mkdir -p results/deseq2/plots
        Rscript scripts/deseq2_analysis.R \
            {input.counts} \
            {input.samples} \
            config.yaml \
            > {log} 2>&1
        """

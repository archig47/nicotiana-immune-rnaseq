SAMPLES = ["SRR6676954", "SRR6676955", "SRR6676956", "SRR6676957", "SRR6676958", "SRR6676959"]
MOCK = ["SRR6676954", "SRR6676955", "SRR6676956"]
PST = ["SRR6676957", "SRR6676958", "SRR6676959"]

rule all:
    input:
        "results/deseq2/volcano_plot.pdf",
        "results/deseq2/pca_plot.pdf"

# Trimming
rule trim_galore:
    input:
        r1="data/reads/{sample}_1.fastq",
        r2="data/reads/{sample}_2.fastq"
    output:
        r1="results/trimmed/{sample}_1_val_1.fq.gz",
        r2="results/trimmed/{sample}_2_val_2.fq.gz"
    shell:
        """
        source /users/fyp/fyp5/miniconda3/etc/profile.d/conda.sh
        conda activate nlr_pipeline
        trim_galore --quality 20 --length 30 --stringency 3 --paired -o results/trimmed {input.r1} {input.r2}
        gzip results/trimmed/{wildcards.sample}_1_val_1.fq
        gzip results/trimmed/{wildcards.sample}_2_val_2.fq
        """

# DESeq2 Analysis (using trimmed files as input for now)
rule deseq2:
    input:
        expand("results/trimmed/{sample}_1_val_1.fq.gz", sample=SAMPLES)
    output:
        pca="results/deseq2/pca_plot.pdf",
        volcano="results/deseq2/volcano_plot.pdf",
        degs="results/deseq2/DEGs_significant.csv"
    shell:
        """
        mkdir -p results/deseq2
        source /users/fyp/fyp5/miniconda3/etc/profile.d/conda.sh
        conda activate nlr_pipeline
        Rscript scripts/deseq2_analysis.R
        """

# Snakemake workflow for RNA-seq analysis
# Full version will be generated

rule all:
    input:
        "results/deseq2/volcano_plot.pdf"

rule download_test_data:
    output:
        expand("data/reads/SRR{srr}.fastq", srr=["6676954", "6676955", "6676956", "6676957", "6676958", "6676959"])
    shell:
        "bash test_data.sh"

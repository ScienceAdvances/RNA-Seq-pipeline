rule fastqc:
    input:
        get_fastq
    output:
        html = "{O}/QC/{S}{U}_fastqc.html",
        zip = "{O}/QC/{S}{U}_fastqc.zip"
    params:
        extra = "--quiet"
    log:
        "{O}/logs/fastqc/{S}{U}.log"
    wrapper:
        config["warpper_mirror"]+"bio/fastqc"


rule samtools_flagstat:
    input:
        "{O}/Count/{S}_sortedByCoord.bam"
    output:
        temp("{O}/QC/{S}.bam.flagstat")
    log:
        "{O}/logs/samtools_flagstat/{S}.log"
    wrapper:
        config["warpper_mirror"]+"bio/samtools/flagstat"

rule align_multiqc:
    input:
        expand("{{O}}/QC/{s}.bam.flagstat", s=SAMPLES)
    output:
        "{O}/Report/align_multiqc.html"
    log:
        "{O}/logs/align_multiqc.log"
    wrapper:
        config["warpper_mirror"]+"bio/multiqc"

rule fastqc_multiqc:
    input:
        expand("{{O}}/QC/{s}{u}_fastqc.zip", s=SAMPLES,u=UNITS)
    output:
        "{O}/Report/fastqc_multiqc.html"
    log:
        "{O}/logs/fastp/fastqc_multiqc.log"
    wrapper:
        config["warpper_mirror"]+"bio/multiqc"

rule qc_plot:
    input:
        "{O}/Count/Count.csv.gz"
    output:
        report("{O}/Report/qc_plot.pdf", caption="../Report/qc_plot.rst", category="Quanlity Control")
    params:
        n_gene=config["qc_plot"].get("n_genes"),
        ellipses=config["qc_plot"].get("ellipses"),
        ellipse_size=config["qc_plot"].get("ellipse_size"),
        group=samples.Group,
        sample_order=samples.Sample,
    conda:
        "../envs/qc_plot.yaml"
    script:
        "../scripts/qc_plot.R"
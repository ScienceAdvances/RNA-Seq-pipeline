rule star_align:
    input:
        unpack(get_trimmed_fastq),
        idx=config['index'],
    output:
        aln= "{O}/Count/{S}_sortedByCoord.bam",
        log="{O}/logs/Count/{S}.log",
        sj="{O}/Count/{S}_splice_junctions.tsv",
        reads_per_gene="{O}/Count/{S}_ReadsPerGene.tsv",
    params:
        extra="--outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --twopassMode Basic"
    threads: 16
    log:
        "{O}/logs/STAR/{S}.log"
    script: "../scripts/star_align.py"
    # wrapper:
    #     config["warpper_mirror"]+"bio/star/align"

rule merge_count:
    input:
        expand("{{O}}/Count/{_}_ReadsPerGene.tsv", _=samples.Sample)
    output:
        temp("{O}/Count/star_count.csv.gz")
    conda:
        "../envs/gene.yaml"
    params:
        samples.Sample
    script:
        "../scripts/merge_count.py"

rule normalise:
    input:
        counts="{O}/Count/star_count.csv.gz",
        info=config['gene_info']
    output:
        report("{O}/Count/Count.csv.gz", caption="../report/normalise.rst",category="STRA-Count")
    conda:
        "../envs/gene.yaml"
    script:
        "../scripts/normalise.py"
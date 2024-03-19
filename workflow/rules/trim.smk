if config["is_pair_end"]:
    rule fastp_pe:
        input:
            sample=get_fastq
        output:
            trimmed=[temp("{O}/Clean/{S}{U}_1.fastq.gz"), temp("{O}/Clean/{S}{U}_2.fastq.gz")],
            html="{O}/QC/{S}{U}.fastp.html",
            json="{O}/QC/{S}{U}.fastp.json",
        log:
            "{O}/logs/fastp/{S}{U}.log"
        threads: 16
        wrapper:
            config["warpper_mirror"]+"bio/fastp"
else:
    rule fastp_se:
        input:
            sample=get_fastq
        output:
            trimmed=temp("{O}/results/Clean/{S}{U}.fastq.gz"),
            html="{O}/report/{S}{U}.fastp.html",
            json="{O}/report/{S}{U}.fastp.json",
        log:
            "{O}/logs/fastp/{S}{U}.log"
        threads: 16
        wrapper:
            config["warpper_mirror"]+"bio/fastp"
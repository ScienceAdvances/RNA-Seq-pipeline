rule dea:
    input:
        "{O}/Count/Count.csv.gz"
    output:
        dea_res=report(expand("{{O}}/DEA/DEA_res_{_}.csv.gz",_=get_contrast()), caption="../report/dea.rst", category="DEA"),
        volcano=report(expand("{{O}}/DEA/Volcano_plot_{_}.pdf",_=get_contrast()), caption="../report/dea.rst", category="DEA"),
        ma=report(expand("{{O}}/DEA/MA_plot_{_}.pdf",_=get_contrast()), caption="../report/dea.rst", category="DEA"),
        norm_count="{O}/DEA/star_deseq_norm_count.csv.gz",
    params:
        contrast=config["compare"],
        lfc=config["lfc"],
        pvalue=config["pvalue"],
        samples=samples,
        use_adjust_p=config["use_adjust_p"]
    conda:
        "../envs/gene.yaml"
    threads: 8
    script:
        "../scripts/dea.py"
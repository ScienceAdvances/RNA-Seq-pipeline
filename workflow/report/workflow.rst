Bulk RNA-seq quantification and quality control workflow
Reads were mapped onto {{ snakemake.config["genome"]["species"] }} build {{ snakemake.config["genome"]["build"] }} with `STAR`_, and adapters were removed with `fastp`_.

For nomalisztion, `gtftools`_ was used to calculate gene_length and `bioninfokit`_ was used to give TPM, FPKM and CPM results.

For quality control, PCA plot, dendrogram plot and heatmap were used to show differences among samples or groups.

{% if snakemake.config["dea"]["active"] %}
`PyDESeq2`_ was used to do differential gene anlysis.
{% endif %}

.. _fastp: https://github.com/OpenGene/fastp
.. _STAR: https://github.com/alexdobin/STAR
.. _gtftools: https://github.com/RacconC/gtftools
.. _bioninfokit: https://github.com/reneshbedre/bioinfokit
.. _PyDESeq2: https://github.com/owkin/PyDESeq2
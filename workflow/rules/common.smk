import numpy as np
import pandas as pd
from snakemake.utils import validate
from snakemake.utils import min_version
from pathlib import Path

min_version("7.25.0")

report: "../report/workflow.rst"

container: "mambaorg/micromamba:1.5.7"

#=====================================================
# validate config.yaml file and samples.csv file
#=====================================================

configfile: "config/config.yaml"

# validate(config, schema="../schemas/config.schema.yaml")
OUTDIR = config['outdir']
samples = pd.read_csv(config["samples"], dtype=str,sep='\t',header=0).fillna(value="")
if not "Unit" in samples.columns:
    samples.loc[:,"Unit"]=""
samples["Unit"] = [f"_{x}" if x else x for x in samples.Unit]
samples.set_index(keys=["Sample", "Unit"], drop=False,inplace=True)

samples.index = samples.index.set_levels(
	[i.astype(str) for i in samples.index.levels]
)  # enforce str in index
# if units are not none, add a _ prefix
fastqs = config["fastqs"]
if config["is_pair_end"]:
	fq1=[f"{fastqs}/{x}{y}_1.fastq.gz" for x,y in zip(samples.Sample,samples.Unit)]
	fq2=[f"{fastqs}/{x}{y}_2.fastq.gz" for x,y in zip(samples.Sample,samples.Unit)]
	samples.insert(loc=0,column="fq2",value=fq2)
	samples.insert(loc=0,column="fq1",value=fq1)
else:
	fq1=[f"{fastqs}/{x}{y}.fastq.gz" for x,y in zip(samples.Sample,samples.Unit)]
	samples.insert(loc=0,column="fq1",value=fq1)

validate(samples, schema="../schemas/samples.schema.yaml")
# validate(samples, schema="workflow/schemas/samples.schema.yaml")
SAMPLES=samples.index.get_level_values(0)
UNITS=samples.index.get_level_values(1)
#=====================================================
# Wildcard constraints
#=====================================================

wildcard_constraints:
	S="|".join(SAMPLES),
	U="|".join(UNITS),
	O=config['outdir']

#=====================================================
# Helper functions
#=====================================================

def get_fastq(wildcards):
	"""Get fastq files of given sample and unit."""
	fastqs = samples.loc[(wildcards.S, wildcards.U), ]
	if config["is_pair_end"]:
		return [fastqs.fq1, fastqs.fq2]
	return [fastqs.fq1]

def get_trimmed_fastq(wildcards):
	"""Get trimmed reads of given sample and unit."""
	fastqs = samples.loc[wildcards.S, :]
	if config["is_pair_end"]:
		# paired-end sample
		fq1=expand("{{O}}/Clean/{{S}}{_}_1.fastq.gz",_=UNITS)
		fq2=expand("{{O}}/Clean/{{S}}{_}_2.fastq.gz",_=UNITS)
		return {"fq1":fq1,"fq2":fq2}
	# single end sample
	return {"fq1":expand("{{O}}/trimmed/{{S}}{_}.fastq.gz",_=UNITS)}

def get_contrast():
	itertuples = pd.read_csv(config['compare'],sep='\t').itertuples(index=False)
	return [f"{x}_VS_{y}" for x,y in itertuples]
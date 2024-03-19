import numpy as np
import pandas as pd
from bioinfokit import analys
import functools

merge = functools.partial(pd.merge,left_index=True, right_index=True)

counts=pd.read_csv(snakemake.input.counts,index_col=0)
info=pd.read_csv(snakemake.input.info,index_col=0)
counts = merge(info,counts)

info = counts.loc[:,["Chr","Start","Stop","Strand","Symbol","GeneType","Length"]]

df = counts.drop(columns=["Chr","Start","Stop","Strand","Symbol","GeneType"])

nm = analys.norm()
nm.tpm(df=df, gl='Length')
nm.rpkm(df=df, gl='Length')
nm.cpm(df=df.drop(columns=["Length"]))

rpkm=nm.rpkm_norm
pkm = "fpkm" if snakemake.config["fastqs"].get("pe") else "rpkm"
rpkm.columns=[pkm+x for x in rpkm.columns]

tpm=nm.tpm_norm
tpm.columns=["tpm_"+x for x in tpm.columns]

cpm=nm.cpm_norm
cpm.columns=["cpm_"+x for x in cpm.columns]

df.drop(columns=["Length"],inplace=True)
df.columns=["count_"+x for x in df.columns]

star=functools.reduce(merge,[info,df,tpm,rpkm,cpm])

star.to_csv(snakemake.output[0])
import pandas as pd
import functools

reader = functools.partial(pd.read_csv,sep='\t',skiprows=4,header=None)
count_list=[reader(x,usecols=[1],dtype="int") for x in snakemake.input]

gene_id =  reader(snakemake.input[0],usecols=[0],index_col=0,dtype="str")

df = pd.concat(count_list,axis=1)
df.columns = snakemake.params[0]
df.insert(loc=0,value=gene_id.index.values,column="GeneID")

df.to_csv(snakemake.output[0],index=False)
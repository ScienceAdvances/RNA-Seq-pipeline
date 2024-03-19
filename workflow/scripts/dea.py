import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.preprocessing import deseq2_norm

from bioinfokit import visuz,analys
from bioinfokit.analys import norm, get_data
import matplotlib.pylab as plt

# snakemake input
n_jobs=snakemake.threads
samples=snakemake.params.samples
contrast=snakemake.params.contrast

# params
use_adjust_p=snakemake.params.use_adjust_p
lfc_thr = (snakemake.params.lfc,snakemake.params.lfc)
pv_thr = (snakemake.params.pvalue,snakemake.params.pvalue)


vs = pd.read_csv(contrast,sep='\t')
vss=[[x,y] for x,y in vs.itertuples(index=False) ]

volcano_path=[x.split(".pdf")[0] for x in snakemake.output.volcano]
ma_path=[x.split(".pdf")[0] for x in snakemake.output.ma]
dea_path=snakemake.output.dea_res
norm_count_path=snakemake.output.get("norm_count","")
is_equal=len(volcano_path)==len(ma_path)==len(vss)==len(dea_path)
assert is_equal,"number of output files must be equal"

#=====================================================
# prepare counts
#=====================================================
counts = pd.read_csv(snakemake.input[0],
            usecols=["Symbol"]+["count_"+ x for x in samples.Sample.unique()],
            index_col="Symbol",
            header=0
            )
counts.columns = [x.split('_')[1] for x in counts.columns]

def unique_exprs(frame, reductions=np.median):
    """
    基因去重复
    frame: row_index = genenames, columns = samples
    """
    frame.index.name=None
    frame.loc[:,'Ref'] = frame.apply(reductions, axis=1)
    frame.sort_values(by='Ref', ascending=False, inplace=True)
    frame.loc[:,'Ref'] = frame.index
    frame.drop_duplicates(subset='Ref', inplace=True)
    frame.drop(columns='Ref', inplace=True)
    return frame
counts=unique_exprs(counts)
counts=counts.transpose()

clinical_df=samples.drop_duplicates(subset="Sample").loc[:,["Sample","Group"]].set_index(keys="Sample",drop=True)

genes_to_keep = counts.columns[counts.sum(axis=0) >= 10]
counts = counts[genes_to_keep]

#=====================================================
# deseq2 norm counts
#=====================================================
if norm_count_path:
	deseq2_counts,size_factors = deseq2_norm(counts.values)
	deseq2_counts=pd.DataFrame(data=deseq2_counts,index=counts.index,columns = counts.columns)
	deseq2_counts.to_csv(norm_count_path)

#=====================================================
# run deseq
#=====================================================
def run_dea(vs,dea_res_path,volcano_path,ma_path):
	dds = DeseqDataSet(
		counts=counts,
		clinical=clinical_df,
		design_factors="Group",
		refit_cooks=True,
		n_cpus=n_jobs,
		ref_level=["Group", vs[1]]
	)
	dds.deseq2()

	stat_res = DeseqStats(dds, n_cpus=n_jobs,contrast=['Group', vs[0], vs[1]])
	stat_res.summary()
	res=stat_res.results_df
	res.insert(loc=0,value=res.index,column="Gene")
	res.to_csv(dea_res_path,index=False)

	#=====================================================
	# volcano plot
	#=====================================================
	if use_adjust_p:
		pv = "paj"
	else:
		pv = "pvalue"
	visuz.GeneExpression.volcano(
		df=res
		,lfc='log2FoldChange'
		,lfc_thr=lfc_thr
		,pv=pv
		,pv_thr=pv_thr
		,sign_line=True
		,gstyle=2
		,show=False
		,plotlegend=True
		,color=("red","grey","green")
		,legendpos='upper right'
		,legendanchor=(1.46,1)
		,figtype='pdf'
		,figname=volcano_path
		)

	#=====================================================
	# MA plot
	#=====================================================

	res.loc[:,"baseMean"]=np.log2(res.baseMean+1)
	visuz.GeneExpression.ma(
		df=res
		,lfc='log2FoldChange'
		,basemean="baseMean"
		,pv=pv
		,gstyle=2
		,show=False
		,plotlegend=True
		,color=("red","grey","green")
		,figtype='pdf'
		,figname=ma_path
		)

for a,b,c,d in zip(vss,dea_path,volcano_path,ma_path):
	run_dea(a,b,c,d)
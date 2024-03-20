import argparse
import os
from typing import Optional
from pathlib import Path

def main():
    parse = argparse.ArgumentParser(description="RNA-seq")
    parse.add_argument("--fastqs", help="dir of fastq files", type=str, required=True)
    parse.add_argument(
        "--species",
        help="species ABBR. one of hsa, mmu, mtu",
        choices=["hsa", "mmu", "mtu"],
        default="hsa",
        required=True,
    )
    parse.add_argument("--is_pair_end", help="out dir path",default=True, choices=[True, False], type=bool, required=True)
    parse.add_argument("--genome", help="genome dir path",default="/home/data/wd/victor", type=str, required=True)
    parse.add_argument("--pipeline", help="pipeline dir path",default="/home/data/wd/victor/RNA-seq-pipeline", type=str, required=True)
    parse.add_argument("--samples", help="sample.tsv file path", type=str, required=True)
    parse.add_argument("--outdir", help="out_dir", type=str, required=True)
    parse.add_argument("--compare", help="compare.tsv file path", type=str, required=False)
    parse.add_argument("--pvalue", help="pvalue cutoff", type=float, default=0.05, required=False)
    parse.add_argument("--use_adjust_p", help="use adjust pvalue or not", type=bool, default=False, required=False)
    parse.add_argument("--lfc", help="log2FoldChange cutoff", type=float, default=0.585, required=False)

    parse.add_argument("--n_jobs", help="no. of cpu to use", type=int, default=16, required=False)
    parse.add_argument("--rerun_incomplete", help="rerun incomplete pipline",action="store_true", required=False)
    parse.add_argument("--dry_run", help="fake run",action="store_true", required=False)

    args = parse.parse_args()
    n_jobs = args.n_jobs
    is_pair_end = args.is_pair_end
    rerun_incomplete = "--rerun-incomplete" if args.rerun_incomplete else ""
    dry_run = "--dry-run" if args.dry_run else ""
    fastqs = str(Path(args.fastqs).absolute())
    samples = str(Path(args.samples).absolute())
    outdir = str(Path(args.outdir).absolute())
    compare = str(Path(args.compare).absolute()) if args.compare else ""
    genome = str(Path(args.genome).absolute())
    pipeline = str(Path(args.pipeline).absolute())
    lfc = args.lfc
    pvalue = args.pvalue
    use_adjust_p = args.use_adjust_p
    species = args.species
    index=f"{genome}/index/{species}"
    gene_info=f"{genome}/{species}.gene_info.tsv"
    os.chdir(pipeline)
    print(f" ============ 输出文件路径为 {outdir} ============ ")
    
    os.system(f"snakemake --use-conda -c {n_jobs} {rerun_incomplete} {dry_run} \
        --config fastqs={fastqs} \
        species={species} is_pair_end={is_pair_end} index={index} samples={samples} \
        gene_info={gene_info} compare={compare} outdir={outdir} \
        use_adjust_p={use_adjust_p} lfc={lfc} pvalue={pvalue}"
        )

if __name__ == '__main__':
    main()
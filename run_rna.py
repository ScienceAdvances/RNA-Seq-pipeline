import argparse
import os
from typing import Optional
from pathlib import Path

def main():
    parse = argparse.ArgumentParser(description="RNA-seq")
    parse.add_argument("-f", "--fastqs", help="dir of fastq files", type=str, required=True)
    parse.add_argument(
        "-s",
        "--species",
        help="species ABBR. one of hsa, mmu, mtu",
        choices=["hsa", "mmu", "mtu"],
        default="hsa",
        required=True,
    )
    parse.add_argument("-p", "--is_pair_end", help="out_dir",default=True, choices=[True, False], type=bool, required=True)
    parse.add_argument("-m", "--samples", help="sample.tsv file", type=str, required=True)
    parse.add_argument("-o", "--outdir", help="out_dir", type=str, required=True)
    parse.add_argument("-c", "--compare", help="compare.tsv file", type=str, required=False)
    parse.add_argument("-v", "--pvalue", help="pvalue", type=float, default=0.05, required=False)
    parse.add_argument("-u", "--use_adjust_p", help="use_adjust_p", type=bool, default=False, required=False)
    parse.add_argument("-l", "--lfc", help="lfc", type=float, default=0.585, required=False)

    parse.add_argument("-n", "--n_jobs", help="n cpu", type=int, default=16, required=False)
    parse.add_argument("-r", "--rerun_incomplete", help="sample.tsv file",action="store_true", required=False)
    parse.add_argument("-d", "--dry_run", help="dry_run",action="store_true", required=False)

    args = parse.parse_args()
    n_jobs = args.n_jobs
    is_pair_end = args.is_pair_end
    rerun_incomplete = "--rerun-incomplete" if args.rerun_incomplete else ""
    dry_run = "--dry-run" if args.dry_run else ""
    fastqs = str(Path(args.fastqs).absolute())
    samples = str(Path(args.samples).absolute())
    outdir = str(Path(args.outdir).absolute())
    compare = str(Path(args.compare).absolute())
    lfc = args.lfc
    pvalue = args.pvalue
    use_adjust_p = args.use_adjust_p
    species = args.species
    index=f"/home/data/wd/victor/index/{species}"
    gene_info=f"/home/data/wd/victor/index/{species}.gene_info.tsv"
    os.chdir("/home/data/wd/victor/RNA-seq-pipeline")
    print(f" ============ 输出文件路径为 {outdir} ============ ")
    
    os.system(f"snakemake --use-conda -c {n_jobs} {rerun_incomplete} {dry_run} \
        --config fastqs={fastqs} \
        species={species} is_pair_end={is_pair_end} index={index} samples={samples} \
        gene_info={gene_info} compare={compare} outdir={outdir} \
        use_adjust_p={use_adjust_p} lfc={lfc} pvalue={pvalue}"
        )

if __name__ == '__main__':
    main()
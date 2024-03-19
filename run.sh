python /home/data/wd/victor/RNA-seq-pipeline/run_rna.py \
    --fastqs fastq \
    --species mtu \
    --is_pair_end True \
    --samples sample.tsv \
    --outdir ./ \
    --compare compare.tsv \
    --pvalue 0.05 \
    --use_adjust_p True \
    --lfc 0.585 \
    --n_jobs 8 \
    --rerun_incomplete
    # --dry_run
    

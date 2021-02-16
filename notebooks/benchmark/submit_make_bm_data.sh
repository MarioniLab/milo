#!/bin/bash

pop_file=$1

outdir=/nfs/team205/ed6/bin/milo_benchmarking/outfiles

## Activate conda env
conda activate milo_bm

## Run 
cat $pop_file | \
while read pop
    do
    for pop_enr in $(seq 0.6 0.1 0.9)
        do
        for seed in $(seq 43 1 45)
            do
            echo "Rscript --pop_enrichment $pop_enr ./make_bm_data.R /nfs/team205/ed6/data/milo_benchmark/embryo_data_bm.RDS $seed $pop" | bsub -o ${outdir}/milo_make_bm_data_${seed}.out -e ${outdir}/milo_make_bm_data_${seed}.err -G team283 -R"select[mem>3500] rusage[mem=3500]" -M3500
            done
        done
    done

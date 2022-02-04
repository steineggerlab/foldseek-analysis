#!/bin/bash
#SBATCH --job-name=scop-clesw
#SBATCH --ntasks=1
#SBATCH --array=[1-175]%1
#SBATCH -t 0
#SBATCH -c 1
#SBATCH -w super003

BATCH=$1

#conda activate ssw

## store file name to a query array
LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p /storage/stephkim/prj-foldseek/organized_data/benchmark/data/batch_submission/core$BATCH.list)
id=$LINE

## timing begins
echo starttime: $BATCH
date

#/data3/stephanie/scope-benchmark/ssw-orig/src/ssw_test -o 22 -e 5 -a cle-opt.mat -p -f 50 /data3/stephanie/scope-benchmark/scope-orig/cle.fin.fasta $id 2> /dev/null | mawk '/^target/{target=$2} /^query/{query=$2} /^optimal_alignment_score/{score=$2; print query,target,score}'  |sort -k1,1 -k3,3nr > /data3/stephanie/prj-foldseek/SCOP-benchmark/CLE-SW/result/${1##*/}".m8";	
/home/stephkim/lib/ssw/src/ssw_test -o 22 -e 5 -a /storage/stephkim/prj-foldseek/organized_data/benchmark/data/clesw/cle-opt.mat -p -f 50 -c /storage/stephkim/prj-foldseek/organized_data/benchmark/data/clesw/cle.fin.fasta /storage/stephkim/prj-foldseek/organized_data/benchmark/data/clesw/cle/$id.cle > /storage/stephkim/prj-foldseek/organized_data/benchmark/results/CLE-SW/$id.m8 

## timing ends
echo endtime: $BATCH
date

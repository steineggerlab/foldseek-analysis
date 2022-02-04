#!/bin/bash
#SBATCH --job-name=scop-CE
#SBATCH --ntasks=1
#SBATCH -t 0
#SBATCH -c 1
#SBATCH -w super003

BATCH=$1

## start timing
echo start time $BATCH
date

## beginning of code
java -jar /storage/stephkim/prj-foldseek/organized_data/benchmark/scripts/CEalign.jar /storage/stephkim/prj-foldseek/organized_data/benchmark/data/timebenchmark/batch_submission/core$BATCH.list /mnt/scratch/scop-pdb-ent/ /storage/stephkim/prj-foldseek/organized_data/benchmark/results/CE/timebenchmark/ > /storage/stephkim/prj-foldseek/organized_data/benchmark/results/CE/timebenchmark/srun$BATCH.out 

## end timing
echo end time $BATCH
date

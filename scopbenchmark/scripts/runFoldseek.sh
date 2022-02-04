#!/bin/bash
#SBATCH --job-name=scop-foldseek
#SBATCH --ntasks=1
#SBATCH -t 0
#SBATCH -c 64
#SBATCH -w super003

## timing the job run
echo start timing
date

## foldseek convert2db
#/data3/stephanie/scope-benchmark/foldseek/build/src/foldseek convert2db /data3/stephanie/scope-benchmark/pdb/ /data3/stephanie/prj-foldseek/SCOP-benchmark/foldseek/db --threads 2

## foldseek prefilter
#/data2/martin/tmp/foldseek/build/src/foldseek prefilter /storage/databases/alphafold2/alphafolddb_ss /storage/databases/alphafold2/alphafolddb_ss /data3/stephanie/prj-foldseek/alphafoldDB/alphafolddb_pref --threads 128

## foldseek align
#/data2/martin/tmp/foldseek/build/src/foldseek align /storage/databases/alphafold2/alphafolddb_ss /storage/databases/alphafold2/alphafolddb_ss /data3/stephanie/prj-foldseek/alphafoldDB/alphafolddb_pref /data3/stephanie/prj-foldseek/alphafoldDB/alphafolddb_aln -a --threads 128



## foldseek easy-search (combination of all convert2db + prefilter + align)
/storage/databases/alphafold2/foldseek/build/src/foldseek easy-search /mnt/scratch/scop-pdb/ /mnt/scratch/scop-pdb/ /storage/stephkim/prj-foldseek/organized_data/benchmark/results/foldseek/swaln /storage/stephkim/prj-foldseek/organized_data/benchmark/results/foldseek/tmp/ --threads 64 -s 9 --max-seqs 2000 -e 10000 --alignment-mode 0

## end time
echo end timing
date


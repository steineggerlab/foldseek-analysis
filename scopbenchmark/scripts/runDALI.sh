#!/bin/bash
#SBATCH --job-name=scop-DALI
#SBATCH --ntasks=1
#SBATCH --array=[1-33]%1
#SBATCH -t 0
#SBATCH -c 1
#SBATCH -w super003

BATCH=$1

## beginning of the code

## store file name to a query array
LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p /storage/stephkim/prj-foldseek/organized_data/benchmark/data/dali/batch_submission/core$BATCH.list)
id=$LINE

## store file name to a target array
while read line ; do
        targetArray+=($line)
done < /storage/stephkim/prj-foldseek/organized_data/benchmark/data/dali/scop.full.dali.list

## nested for loop to run FASTalign
#cd /data3/stephanie/prj-foldseek/SCOP-benchmark/DALI/timebenchmark/

## formatting input files (DALI/import.pl)
#for i in "${array[@]}"; do
#	pdbid=$(echo $i | cut -c2-5)
#	/data3/stephanie/external-lib/DaliLite.v5/bin/import.pl --pdbfile /data3/stephanie/prj-foldseek/SCOP-benchmark/DALI/newpdb/$i.ent --pdbid $pdbid --dat /data3/stephanie/prj-foldseek/SCOP-benchmark/DALI/DAT --clean
#done

## start timing
echo start timing $BATCH
date

## beginning of code
/home/stephkim/lib/DaliLite.v5/bin/dali.pl --cd1 $id --db /storage/stephkim/prj-foldseek/organized_data/benchmark/data/dali/scop.full.dali.list --TITLE systematic --dat1 /mnt/scratch/scop-DAT/DAT --dat2 /mnt/scratch/scop-DAT/DAT --outfmt "summary" --clean

## end timing
echo end timing $BATCH
date



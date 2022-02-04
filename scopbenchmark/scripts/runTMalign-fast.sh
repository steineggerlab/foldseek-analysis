#!/bin/bash
#SBATCH --job-name=scop-tmfast
#SBATCH --ntasks=1
#SBATCH --array=[1]
#SBATCH -t 0
#SBATCH -c 1
#SBATCH -w super003

BATCH=$1


## store file name to a query array
LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p /storage/stephkim/prj-foldseek/organized_data/benchmark/data/batch_submission/core$BATCH.list)
id=$LINE

## store file name to a target array
while read line ; do
        array+=($line)
done < /storage/stephkim/prj-foldseek/organized_data/benchmark/data/scop.file.full.name.list

## start timing
echo start timing: $BATCH
date

## nested for loop to run FASTalign
for i in  "${array[@]}"; do
	/home/stephkim/lib/TM-align/TMalign /storage/stephkim/prj-foldseek/organized_data/benchmark/data/scop-pdb/$id /storage/stephkim/prj-foldseek/organized_data/benchmark/data/scop-pdb/$i -fast >> /storage/stephkim/prj-foldseek/organized_data/benchmark/results/tmalign-fast/all-against-all/output2/$id.out
done

## end timing
echo end timing: $BATCH
date

#rm *_tmaln.out

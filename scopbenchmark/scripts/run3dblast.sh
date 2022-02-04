#!/bin/bash
#SBATCH --job-name=scop-3dblast
#SBATCH --ntasks=1
#SBATCH --array=[1-11]%1
#SBATCH -t 0
#SBATCH -c 1
#SBATCH -w super003

BATCH=$1

#./3d-blast -p <query protein> <chain id> -d <SADB file> -o <output file>

## store file name to a query array
LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p /storage/stephkim/prj-foldseek/organized_data/benchmark/data/batch_submission/core$BATCH.list)
id=$LINE


## timing begins
echo starttime: $BATCH
date

## formatDB (already converted to SA from local computer - dssp)
#/home/stephkim/lib/blast-2.2.26/bin/formatdb -i $database -o T -p T

## blastall
/home/stephkim/lib/blast-2.2.26/bin/blastall -p blastp -M BLOSUM62 -d /mnt/scratch/scop-3dblast/3dblastSA -i /mnt/scratch/scop-fasta/$id.fasta -o /storage/stephkim/prj-foldseek/organized_data/benchmark/results/3dblast/$BATCH.$id.out -v 50 -b 50 -e 10 -F F -G 8 -E 2
#/home/stephkim/lib/Standalone_3d-blast_Linux/blast/bin/blastall -p blastp -M BLOSUM62 -d /mnt/scratch/scop-3dblast/3dblastSA -i /mnt/scratch/scop-pdb/$id -o /storage/stephkim/prj-foldseek/organized_data/benchmark/results/3dblast/$BATCH.$id.out -v 50 -b 50 -e 10 -F F -G 8 -E 2
## timing ends
echo endtime: $BATCH
date

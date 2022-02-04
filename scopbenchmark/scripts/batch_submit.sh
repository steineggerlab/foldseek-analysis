#!/bin/bash

## run batch (in various directories)
#for i in {1..64}; do
#	if [ -z "$(ls -A /data3/stephanie/prj-foldseek/alphafoldDB/DALI-segment/result/tmp_$i)" ]; then
#		cd /data3/stephanie/prj-foldseek/alphafoldDB/DALI-segment/result/tmp_$i
#		sbatch /data3/stephanie/prj-foldseek/alphafoldDB/DALI-segment/runDALI.sh $i
#	fi
#done


## run batch (single directory)
for i in {224..224}; do
	cd /storage/stephkim/prj-foldseek/organized_data/benchmark/results/tmalign/all-against-all/slurmoutput
	#mkdir core$i
	#cd /storage/stephkim/prj-foldseek/organized_data/benchmark/results/dali/timebenchmark/core$i/
	#rm *
	sbatch /storage/stephkim/prj-foldseek/organized_data/benchmark/scripts/runTMalign.sh $i
done
	

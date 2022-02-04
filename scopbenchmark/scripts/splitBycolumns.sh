#!/bin/bash

echo TMalign
## timing the job run
date

## beginning of the code
## store file name to an array
while read line ; do
        array+=($line)
done < /data3/stephanie/prj-foldseek/SCOP-benchmark/SCOPpdb.fileName.list

for i in "${array[@]}"; do
	#echo awk '$1 ~ /^"'"$i"'"/' /data3/stephanie/prj-foldseek/SCOP-benchmark/TMalign/tmaln-all.sort.out
	awk -v var=$i '$1 ~var {print $0}' /storage/stephkim/prj-foldseek/organized_data/benchmark/results/tmalign/TMraw.scop.complete.sort.fix.out >> /storage/stephkim/prj-foldseek/organized_data/benchmark/results/tmalign/fullresult/$i.m8
done

## end time
date



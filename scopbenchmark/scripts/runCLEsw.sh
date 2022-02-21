#!/bin/bash

## query file array
while read line ; do
	array+=($line)
done < ./cle.file.list

## timing starts
date

## running SSW
for id in "${array[@]}"; do
	./ssw_test -o 100 -e 10 -a ./cle.shen.mat -p -c ./cle.fin.fasta ./cle/$id.cle > ./$id.m8 
done

## timing ends
date

## extracting alignment result from raw file
for id in "${array[@]}"; do
        grep -a -n "^target_name:\|^query_name:\|^optimal_alignment_score:" ./$id.m8 | awk '{print $1,$2}' | awk -F ":" '{print $1,$3}' | sort -nk1,1 | awk '{print $2}' | xargs -n 3 | awk '{print $2,$1,$3}' | sort -nrk3,3 > ./$id.clean.m8;
done

## concat individual CLEsw alignment file
cat ./*.m8 | awk '{print $1, $2, $3}' | sort -k1,1 -k3,3nr > ../alignResults/rawoutput/cleswaln

## generate ROCX file
./bench.awk ../data/scop_lookup.fix.tsv <(cat ../alignResults/rawoutput/cleswaln) > ../alignResults/rocx/clesw.rocx

## calculate auc
 awk '{ famsum+=$3; supfamsum+=$4; foldsum+=$5}END{print famsum/NR,supfamsum/NR,foldsum/NR}' ../alignresults/rocx/clesw.rocx

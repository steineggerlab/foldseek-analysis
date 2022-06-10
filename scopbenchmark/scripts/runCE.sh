#!/bin/bash

## start timing
date

## beginning of code
java -jar ./CEalign.jar ./core$BATCH.list ./scop-pdb/ ./output/ > ./output/srun.out 

## end timing
date

## extracting CE Z-score, rmsd, tmscore
for f in ./output/*.result; do
        awk '{print$1 "\t" $2 "\t" $8 "\t" $9 "\t" $14}' $f | sed 's|:|\t|g' | awk '{print $1 "\t" $2 "\t" $4 "\t" $6 "\t" $8}' | sort -nrk3 > ./$f.extracted.out
done

## concat all extracted result file
cat ./*.extracted.out | awk '{print $1, $2, $3}' | sort -k1,1 -k3,3nr > ../alignResults/rawoutput/cealn

## generate ROCX file
./bench.awk ../data/scop_lookup.fix.tsv <(cat ../alignResults/rawoutput/cealn) > ../alignResults/rocx/ce.rocx

## calculate auc
 awk '{ famsum+=$3; supfamsum+=$4; foldsum+=$5}END{print famsum/NR,supfamsum/NR,foldsum/NR}' ../alignResults/rocx/ce.rocx

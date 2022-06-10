#!/bin/bash

## start timing
date 

## foldseek easy-search (combination of all convert2db + prefilter + align)
./foldseek/build/src/foldseek easy-search ./benchmark/data/scop-pdb/ ./benchmark/data/scop-pdb/ ./benchmark/alignResults/rawoutput/foldseekaln ./benchmark/alignResults/tmp/ --threads 64 -s 9.5 --max-seqs 2000 -e 10

## end timing
date

## generate ROCX file
./bench.awk ../data/scop_lookup.fix.tsv <(cat ../alignResults/rawoutput/foldseekaln) > ../alignResults/rocx/foldseek.rocx

## calculate auc
 awk '{ famsum+=$3; supfamsum+=$4; foldsum+=$5}END{print famsum/NR,supfamsum/NR,foldsum/NR}' ../alignResults/rocx/foldseek.rocx



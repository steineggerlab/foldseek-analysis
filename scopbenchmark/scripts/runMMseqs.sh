#!/bin/bash

## start timing
date 

## mmseqs search module
mmseqs search ./scop-fasta ./scop-fasta ./benchmark/alignResults/rawoutput/mmseqs.m8 tmp -a --threads 64 -s 7.5 -e 10000 --max-seqs 2000

## end timing
date


## generate ROCX file
./bench.awk ../data/scop_lookup.fix.tsv <(cat ../alignResults/rawoutput/mmseqsaln) > ../alignResults/rocx/mmseqs.rocx

## calculate auc
 awk '{ famsum+=$3; supfamsum+=$4; foldsum+=$5}END{print famsum/NR,supfamsum/NR,foldsum/NR}' ../alignResults/rocx/mmseqs.rocx

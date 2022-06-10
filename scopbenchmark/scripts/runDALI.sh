#!/bin/bash

## store file name to a target array
while read line ; do
        array+=($line)
done < ./scop.full.dali.list

## converting pdb files to dat files  (DALI/import.pl)
for i in "${array[@]}"; do
#	pdbid=$(echo $i | cut -c2-5)
	./import.pl --pdbfile ./scop-pdb/$i --pdbid $i --dat ./DAT --clean
done

## start timing
date

## run DALI all-versus-all
for i in "${array[@]}"; do
	./dali.pl --cd1 $i --db ./scop.full.dali.list --TITLE systematic --dat1 ./DAT --dat2 ./DAT --outfmt "summary" --clean
done

## end timing
date

## extract DALI Z-score
for f in ./*.txt; do
	awk -v query="$f" '(NF == 7 && $1 != "ident") {print query, $2, $3, $4, $5, $6, $7}' $f | sed 's/-//g' | sed 's/.txt//g' >> ./$f.out
done

## concat individual 3dblast alignment file
cat ./*.out | awk '{print $1, $2, $3}' | sort -k1,1 -k3,3nr > ../alignResults/rawoutput/dalialn

## generate ROCX file
./bench.awk ../data/scop_lookup.fix.tsv <(cat ../alignResults/rawoutput/dalialn) > ../alignResults/rocx/dali.rocx

## calculate auc
 awk '{ famsum+=$3; supfamsum+=$4; foldsum+=$5}END{print famsum/NR,supfamsum/NR,foldsum/NR}' ../alignResults/rocx/dali.rocx
  

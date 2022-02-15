#!/bin/bash

while read line ; do
        array+=($line)
done < ../data/scop.file.full.name.list


##### 3DBLAST #####
## start timing (3dblast) 
date

## blastall all-versus-all
for id "${array[@]}"; do
	./blastall -p blastp -M BLOSUM62 -d ../data/3dblastSA -i ../data/3dblast-fasta/$id.fasta -o ./$id.out -v 12000 -b 12000 -e 10 -F F -G 8 -E 2
done

## end timing (3dblast)
date

## extract 3dblast value from raw
for file in *.out; do
        grep -n -a "^./scop-pdb/" ./$file | sed 's|./scop-pdb/| |g' | sed 's|:||g' | sort -nk1,1 | awk -v var="$file" '{print var,$2,$5,$4}' | sed 's|.out||g' > ./$file;
done


##### 3dblast + SSW #####
## start timing (3dblast + ssw)
date

## 3dblast aligned with Smith-Waterman (ssw) 
for id "${array[@]}"; do
	./ssw_test -o 8 -e 2 -a ./BLOSUM62 -p -c ../data/3dblast.fin.fasta ../data/3dblast-fasta/$id.fasta > ./$id.m8
done

## end timing (3dblast + ssw)
date

## extracting alignment result from raw file
for file in "${array[@]}"; do
        grep -a -n "^target_name:\|^query_name:\|^optimal_alignment_score:" ./$file.m8 | awk '{print $1,$2}' | awk -F ":" '{print $1,$3}' | sort -nk1,1 | awk '{print $2}' | xargs -n 3 | awk '{print $2,$1,$3}' | sort -nrk3,3 > ./$file.m8;
done

## concat individual 3dblast alignment file
cat ./*.m8 | awk '{print $1, $2, $3}' | sort -k1,1 -k3,3nr > ../alignResults/rawoutput/3dblastswaln

## generate ROCX file
./bench.awk ../data/scop_lookup.fix.tsv <(cat ../alignResults/rawoutput/3dblastswaln) > ../alignResults/rocx/3dblastsw.rocx

## calculate auc
awk '{ famsum+=$3; supfamsum+=$4; foldsum+=$5}END{print famsum/NR,supfamsum/NR,foldsum/NR}' ../alignResults/rocx/3dblastsw.rocx

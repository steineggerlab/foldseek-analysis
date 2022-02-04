#!/bin/bash
#./bench.awk scop_lookup.tsv <(cat 3di.v3.result/*.m8) > 3di_v3.adjopt.rocx
## sort raw score file
awk '{print $1, $2, $3, $4}' /storage/stephkim/prj-foldseek/organized_data/benchmark/results/tmalign-fast/all-against-all/totalRaw.out | sort -k1,1 -k3,3nr > /storage/stephkim/prj-foldseek/organized_data/benchmark/results/tmalign-fast/all-against-all/TMfast.raw.scop.merged.rocx.sort.out

## run bench.awk script for a rocx file
./bench.awk scop_lookup.fix.tsv <(cat /storage/stephkim/prj-foldseek/organized_data/benchmark/results/tmalign-fast/all-against-all/TMfast.raw.scop.merged.rocx.sort.out) > /storage/stephkim/prj-foldseek/organized_data/benchmark/results/tmalign-fast/all-against-all/TMfast.v2.rocx

#!/bin/bash

# Pacth SSW first to handle larger alphabet sizes:
patch -u tmp/ssw/src/ssw.c -i ssw.patch
(cd tmp/ssw/src && make && cp ssw_test ../../ssw_test)

for k in {4..40..4}
do
  RUN=$RUN ./learnAlphabet.sh $k 20 "data/pdbs_train.txt" "data/pdbs_val.txt" "tmp/kmodels/$k"
done


cat data/pdbs_train.txt data/pdbs_val.txt | sort | uniq > tmp/pdbs_all.txt

for k in {4..40..4}
do
  echo -n "$k " >> tmp/koptimization.log
  INVST=$(cat tmp/kmodels/$k/invalid_state.txt)
  RUN=$RUN ./run-benchmark.sh \
  "tmp/kmodels/$k/encoder.pt" \
  "tmp/kmodels/$k/states.txt" \
  "tmp/kmodels/$k/sub_score.mat" \
  tmp/pdbs_all.txt data/scop_lookup.tsv 270 0 2 $INVST \
  | tee -a tmp/koptimization.log

  cp tmp/result.rocx "tmp/kmodels/$k/result.rocx"
done



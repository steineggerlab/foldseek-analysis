
mkdir -p tmp

# Split PDBs by fold into 4 partitions
python3 split_folds.py

# Build 4x train/test sets
rm tmp/pdbs_train?.txt
rm tmp/pdbs_val?.txt

for k in {0..3} # each fold
do
  for s in {0..3} # each spit
  do
    if [ $k = $s ]; then
      awk '/^d/{print $1}' "tmp/fold_split$s.txt" >> "tmp/pdbs_val$k.txt"
    else
      awk '/^d/{print $1}' "tmp/fold_split$s.txt" >> "tmp/pdbs_train$k.txt"
    fi
  done
done

# Train splitmodels
for k in {0..3}
do
  RUN=$RUN ./learnAlphabet.sh 20 \
    "tmp/pdbs_train$k.txt" "tmp/pdbs_train$k.txt" "tmp/splitmodels/sp$k"
done

# Compare splitmodels against foldseek
for k in {0..3}
do
  echo "Fold $k" >> tmp/crossval.log
  echo -n "Ref: " >> tmp/crossval.log
  RUN=$RUN ./run-benchmark.sh \
    data/foldseek_v1/encoder.pt \
    data/foldseek_v1/states.txt \
    data/foldseek_v1/sub_score.mat \
    "tmp/pdbs_val$k.txt" data/scop_lookup.tsv 270 0 2 D \
    | tee -a tmp/crossval.log
  cp tmp/result.rocx "tmp/crossval_ref$k.rocx"

  echo -n "Splitmodel: " >> tmp/crossval.log
  INVST=$(cat tmp/splitmodels/sp$k/invalid_state.txt)
  RUN=$RUN ./run-benchmark.sh \
    "tmp/splitmodels/sp$k/encoder.pt"\
    "tmp/splitmodels/sp$k/states.txt"\
    "tmp/splitmodels/sp$k/sub_score.mat"\
    "tmp/pdbs_val$k.txt" data/scop_lookup.tsv 270 0 2 $INVST \
    | tee -a tmp/crossval.log
  cp tmp/result.rocx "tmp/crossval_splitmodel$k.rocx"
done


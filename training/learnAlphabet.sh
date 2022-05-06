#!/bin/bash -x

if [ "$#" -ne 5 ]; then
    echo "Illegal number of parameters"
    exit
fi

K=$1
TRIES=$2
PDBS_TRAIN=$3
PDBS_VAL=$4

OUTPUT_DIR=$5 # save alphabet here

mkdir -p tmp
mkdir -p $OUTPUT_DIR

# Fetch PDBs
if [ ! -d tmp/pdb ]; then
    curl https://wwwuser.gwdg.de/~compbiol/foldseek/scp40pdb.tar.gz | tar -xz -C tmp
fi

# Compile ssw_test
if [ ! -f tmp/ssw_test ]; then
    git clone --depth 1 https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library tmp/ssw
    (cd tmp/ssw/src && make)
    cp tmp/ssw/src/ssw_test tmp/ssw_test
fi

# Filter alignments for training
awk 'FNR==NR {pdbs[$1]=1; next}
     ($1 in pdbs) && ($2 in pdbs) {print $1,$2,$10}' \
         $PDBS_TRAIN data/tmaln-06.out > tmp/pairfile_train.out

# Benchmark virtual center positions
#for D in {1,1.5,2,2.5,3}; do
#    for THETA in {45..360..45}; do
#        for TAU in {0..180..45}; do
#            echo -n $THETA $TAU $D >> "$OUTPUT_DIR/log.txt"
#            # create alphabet and benchmark
#        done
#    done
#done

THETA=270
TAU=0
D=2

# Create training data
./create_vqvae_training_data.py \
    tmp/pdb tmp/pairfile_train.out $THETA $TAU $D tmp/vaevq_training_data.npy

for ((seed=0;seed<$TRIES;seed++))
do
    echo -n "$seed " >> "$OUTPUT_DIR/log.txt"

    ./train_vqvae.py $seed tmp/vaevq_training_data.npy tmp $K \
    | awk '/opt_loss=/{printf "%s ", $2}' >> "$OUTPUT_DIR/log.txt"
    
    $RUN \
      ./encode_pdbs.py tmp/encoder.pt tmp/states.txt \
      --pdb_dir tmp/pdb --virt $THETA $TAU $D \
      < $PDBS_TRAIN > tmp/seqs.csv
    
    ./create_submat.py tmp/pairfile_train.out tmp/seqs.csv --mat tmp/sub_score.mat
    
    ./run-benchmark.sh tmp/encoder.pt tmp/states.txt tmp/sub_score.mat \
      $PDBS_VAL data/scop_lookup.tsv $THETA $TAU $D X >> "$OUTPUT_DIR/log.txt"

done

# Find best seed
# TMalign.rocx => 0.928162 0.662063 0.275436
SEED=$(awk '{print $1, ($3/0.928162 + $4/0.662063 + $5/0.275436) / 3}' "$OUTPUT_DIR/log.txt" \
        | sort -rk 2,2 | head -n 1 | awk '{print $1}')

# Create final alphabet
./train_vqvae.py "$SEED" tmp/vaevq_training_data.npy $OUTPUT_DIR $K

# Create final submat
cp $PDBS_TRAIN tmp/pdbs_submat.txt
awk 'FNR==NR {pdbs[$1]=1; next}
     ($1 in pdbs) && ($2 in pdbs) {print $1,$2,$10}' \
         tmp/pdbs_submat.txt data/tmaln-06.out > tmp/pairfile_submat.out

$RUN ./encode_pdbs.py $OUTPUT_DIR/encoder.pt $OUTPUT_DIR/states.txt \
  --pdb_dir tmp/pdb --virt $THETA $TAU $D \
    < tmp/pdbs_submat.txt > tmp/seqs.csv

./create_submat.py tmp/pairfile_submat.out tmp/seqs.csv \
  --mat tmp/sub_score.mat --merge_state X \
  | tee tmp/create_submat.log

awk '/^assign_invalid_states_to/{printf "%s", $3}' tmp/create_submat.log > "$OUTPUT_DIR/invalid_state.txt"

# Add X to submat TODO: adapt to k
awk 'NR==1 {printf "%s   X\n", $0}
NR!=1 {printf "%s   0\n", $0}
END{print "X   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0"}' \
 tmp/sub_score.mat > "$OUTPUT_DIR/sub_score.mat"

# Merge X state in seq. file
#awk 'FNR==NR && /^assign_invalid_states_to/{newxstate=$3;nextfile}
#     FNR!=NR {gsub(/X/,newxstate,$2);print $1,$2}' \
#         tmp/create_submat.log tmp/seqs.csv > tmp/seqs_no_x.csv


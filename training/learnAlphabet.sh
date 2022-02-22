#!/bin/bash -x

mkdir -p tmp

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
         data/pdbs_train.txt data/tmaln-06.out > tmp/pairfile_train.out

# Filter scop_lookup file for validation
awk 'FNR==NR {pdbs[$1]=1; next}
     ($1 in pdbs) {print $0}' \
        data/pdbs_val.txt data/scop_lookup.tsv > tmp/scop_lookup_val.tsv

# Benchmark virtual center positions
#for D in {1,1.5,2,2.5,3}; do
#    for ALPHA in {45..360..45}; do
#        for BETA in {0..180..45}; do
#            echo -n $ALPHA $BETA $D >> tmp/log.txt
#            # create alphabet and benchmark
#        done
#    done
#done

ALPHA=270
BETA=0
D=2

# Create training data
./create_vqvae_training_data.py \
    tmp/pdb tmp/pairfile_train.out $ALPHA $BETA $D tmp/vaevq_training_data.npy

rm tmp/log.txt
for seed in {1..100}
do
    echo -n "$seed " >> tmp/log.txt

    ./train_vqvae.py $seed tmp/vaevq_training_data.npy tmp
    
    ./encode_pdbs.py tmp/encoder.pt tmp/states.txt --pdb_dir tmp/pdb --virt $ALPHA $BETA $D \
        < data/pdbs_train.txt > tmp/seqs.csv
    
    ./create_submat.py tmp/pairfile_train.out tmp/seqs.csv --mat tmp/sub_score.mat
    
    # Encode validation PDBs
    ./encode_pdbs.py tmp/encoder.pt tmp/states.txt --pdb_dir tmp/pdb --virt $ALPHA $BETA $D \
        < data/pdbs_val.txt > tmp/seqs_val.csv
    
    # Prepare for benchmark
    mkdir -p tmp/splits
    mkdir -p tmp/alignments
    awk '{print ">" $1} {print $2}' < tmp/seqs_val.csv > tmp/target.fasta
    split -n 30 -d tmp/target.fasta tmp/splits/split_ --additional-suffix=.fasta
    
    ./run-smithwaterman.sh 8 2
    
    ../scopbenchmark/scripts/bench.awk tmp/scop_lookup_val.tsv \
        <(cat tmp/alignments/*.m8) > tmp/result.rocx
    
    # Calculate AUC
    awk '{famsum+=$3; supfamsum+=$4; foldsum+=$5} END{print famsum/NR,supfamsum/NR,foldsum/NR}' \
        tmp/result.rocx >> tmp/log.txt
    
done

# Find best seed
# TMalign.rocx => 0.928162 0.662063 0.275436
SEED=$(awk '{print $1, ($2/0.928162 + $3/0.662063 + $4/0.275436) / 3}' tmp/log.txt \
        | sort -rk 2,2 | head -n 1 | awk '{print $1}')

# Create final alphabet
./train_vqvae.py "$SEED" tmp/vaevq_training_data.npy tmp
find tmp/pdb -type f -name 'd*' -printf "%f\n" > tmp/pdbs.txt
awk '{print $1,$2,$10}' data/tmaln-06.out > tmp/pairfile.out
./encode_pdbs.py tmp/encoder.pt tmp/states.txt --pdb_dir tmp/pdb --virt $ALPHA $BETA $D \
    < tmp/pdbs.txt > tmp/seqs.csv
./create_submat.py tmp/pairfile.out tmp/seqs.csv --mat tmp/sub_score.mat --merge_state X \
    | tee tmp/create_submat.log
# Merge X state in seq. file
awk 'FNR==NR && /^assign_invalid_states_to/{newxstate=$3;nextfile}
     FNR!=NR {gsub(/X/,newxstate,$2);print $1,$2}' \
         tmp/create_submat.log tmp/seqs.csv > tmp/seqs_no_x.csv


#!/bin/bash -x

ENCODER=$1
STATES=$2
SUBMAT=$3

PDBS=$4 # list of SIDs
SCOPLOOKUP=$5

THETA=$6 # virtc
TAU=$7
D=$8

INVALIDSTATES=$9

# Filter scop_lookup.tsv
awk 'FNR==NR {pdbs[$1]=1; next}
     ($1 in pdbs) {print $0}' \
        $PDBS $SCOPLOOKUP > tmp/scop_lookup_filtered.tsv

# Encode validation PDBs
$RUN ./encode_pdbs.py $ENCODER $STATES --pdb_dir tmp/pdb \
  --virt $THETA $TAU $D \
  --invalid-state $INVALIDSTATES \
    < $PDBS > tmp/seqs.csv

cp $SUBMAT tmp/sub_score.mat

# Prepare for benchmark
mkdir -p tmp/splits
mkdir -p tmp/alignments
awk '{print ">" $1} {print $2}' < tmp/seqs.csv > tmp/target.fasta
split -n 30 -d tmp/target.fasta tmp/splits/split_ --additional-suffix=.fasta

./run-smithwaterman.sh 8 2  # sub_score.mat, target.fasta

./roc1.awk tmp/scop_lookup_filtered.tsv \
    <(cat tmp/alignments/*.m8) > tmp/result.rocx

# Calculate AUC
awk '{famsum+=$3; supfamsum+=$4; foldsum+=$5} END{print famsum/NR,supfamsum/NR,foldsum/NR}' \
    tmp/result.rocx

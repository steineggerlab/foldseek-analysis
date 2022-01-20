#!/bin/bash

# Find optimal virtual center
#for D in {2..4}
#do
#    for ALPHA in {45..360..45}
#    do
#        for BETA in {0..180..45}
#        do
#            echo $ALPHA $BETA $D
#            ./create_vqvae_training_data.py data/pdbs $ALPHA $BETA $D
#            ./train_vqvae.py
#            # run benchmark ...
#        done
#    done
#done

PDBDIR=$1

./create_vqvae_training_data.py "$PDBDIR" 270 0 2 

SEED=42 ./train_vqvae.py

echo 'd12asa_' \
    | ./encode_pdbs.py encoder.pt states.txt --pdb_dir "$PDBDIR" --virt 270 0 2

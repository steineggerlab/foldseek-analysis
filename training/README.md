
# Learn 3Di states

```
    ./learnAlphabet.sh
```

Creates 3Di states (encoder.pt, states.txt) and a substitution matrix (sub_score.mat).
Tests 100 seeds for weight initialization (training an alphabet + benchmarking) and uses the best one.

Note: downloads PDBs and Smith-Waterman aligner.

## Requirements
- python3: torch, biopython

## Included data
- tmaln-06.out: TMalign alignments all-against-all in SCOPe 2.07, where TMscore >= 0.6
- pdbs_train.txt: SCOPe SIDs of trainings set
- pdbs_val.txt: SCOPe SIDs of validation set
- scop_lookup.tsv: SCOPe classification per domain


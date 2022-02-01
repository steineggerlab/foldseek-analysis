# E-Value related scripts
---
##  Fitting the e-value parameter for mu and lambda
The jupyter notebook shows the workflow how to fit mu/lambda and how to train a NN from the query composition.
The fitevd.cpp is mostly code from HMMer calibrate and is licensed as GPLv3.

## Parsing and plotting E-Value data
To parse the E-Values from Foldsee, MMseqs2, or 3D-BLAST for either SCOP data or
shuffled data into a different format for plotting use get_evalues_3dblast.py,
get_FP_shuffle.py, get_TP_FP_SCOP.py.
To then use that data to adjust the bin sizes of a frequency histrogram, use
prep_plot_evalues.py.

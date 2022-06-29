# Foldseek (foldseek Version: fdc503ea7fb6c34991146ccf9bdec07125f1b6ba)
# make target db
foldseek createdb spici_target foldseek_spici/target
# call foldseek
foldseek easy-search spici_query  foldseek_spici/target foldseek_multidomain_results tmp --threads 64 --format-output query,target,qstart,tstart,cigar -s 9.5 -e 1000 --max-seqs 10000

# MMseqs2 Version: 13.45111
# call mmseqs2
mmseqs search foldseek_spici/query foldseek_spici/target mmseqs_multidomain_results tmp --threads 64 -a -s 7.5 -e 1000 --max-seqs 10000
# convertalis
mmseqs convertalis foldseek_spici/query foldseek_spici/target mmseqs_spici_results mmseqs_spici_results.m8 --threads 1 --format-output query,target,qstart,tstart,cigar

# TODO STEPH fill in commands and versions
# TM-align

# DALi

# CLE-SW

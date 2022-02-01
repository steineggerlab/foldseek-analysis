#! /usr/bin/env python

'''
This script reads in the results of shuffled SCOP analyses and outputs a list of
query ID, target ID, E-Value

Usage:
python3 get_FP_shuffle.py SHUFFLED_RESULTS OUTPUT_FILE

'''

import sys
import os


result_file = sys.argv[1]
output_file = sys.argv[2]

fp = []

with open(result_file) as input:
    for line in input:
        result_line = line.split('\t')
        id_q = result_line[0].split('_')[0].split('.ent')[0]
        id_t = result_line[1].split('_')[0].split('.ent')[0]
        evalue = result_line[10].strip() # for mmseqs2 shuffle, change for other shuffled data

        fp.append(id_q + ',' + id_t+ ',' +  evalue + ',0\n')

file = open(output_file, 'w')
for evalue in fp:
    file.write(evalue)
file.close()

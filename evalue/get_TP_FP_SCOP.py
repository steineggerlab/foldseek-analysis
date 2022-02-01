#! /usr/bin/env python

'''
This script requires the result file for the SCOP analysis, the SCOP
fold/superfamily/familiy lookup file and the TM-results file, it then determins
true/false positives

Usage:
python3 get_TP_FP_SCOP.py RESULT_FILE SCOP_LOOKUP_FILE TMSCORE_FILE OUTPUT_FILE

'''

import sys
import os

result_file = sys.argv[1]
lookup = sys.argv[2] #'scop_lookup.fix.tsv'
tmalign = sys.argv[3] #'TMraw.scop.complete.both.out'
output_file = sys.argv[4]

class_lookup = {}

# get scop lookup info
with open(lookup) as input:
    for line in input:
        input_line = line.split('\t')
        id = input_line[0].split('_')[0].split('.')[0]
        cat = input_line[1].split('.')[1].strip(), input_line[1].split('.')[2].strip(), input_line[1].split('.')[3].strip(), input_line[1].split('.')[0].strip()
        class_lookup[id] = cat

# read in tm-align results file
tmscore = {}
with open(tmalign) as input:
    for line in input:
        input_line = line.split('\t')
        id_q = input_line[0].split('_')[0].split('.')[0]
        id_t = input_line[1].split('_')[0].split('.')[0]
        score = max(float(input_line[2]),float(input_line[3]))
        tmscore[id_q + '-' + id_t] = score

tpfp = []
exclude = 0

with open(result_file) as input:
    for line in input:
        # make sure the correct delim is used
        result_line = line.split('\t')
        id_q = result_line[0].split('_')[0].split('.')[0] # query id
        id_t = result_line[1].split('_')[0].split('.')[0] # target id
        # make sure the correct columns are used
        # evalue is 11th column for foldseek
        evalue = result_line[10].strip()
        # exclude sequences who are beta-propellers and certain
        # parallel beta-sheets a FP
        if( (class_lookup[id_q][3] == 'b' and int(class_lookup[id_q][0]) <  71 and int(class_lookup[id_q][0]) > 65)
        or (class_lookup[id_q][3] == 'c' and ((int(class_lookup[id_q][0]) <  6 and int(class_lookup[id_q][0]) > 1)
        or (int(class_lookup[id_q][0]) == 27) or (int(class_lookup[id_q][0]) == 28)  or (int(class_lookup[id_q][0]) == 30)  or (int(class_lookup[id_q][0]) == 31)))
        ):
            exclude += 1
        # consider anything with a tmscore >= 0.3 a TP
        elif((class_lookup[id_q][0] == class_lookup[id_t][0]) or (tmscore[id_q + '-' + id_t] >= 0.3)):
            tpfp.append(id_q + ',' + id_t+ ',' +  evalue + ',1\n')
        else:
            tpfp.append(id_q + ',' + id_t+ ',' +  evalue + ',0\n')

file = open(output_file, 'w')
for evalue in tpfp:
    file.write(evalue)
file.close()

#! /usr/bin/env python

'''
This script reads in every file in the folder (sys.argv[1]) that contains the
3D-BLAST results and extracts the query and target ID and the e-Value and
output a single file (sys.argv[2]) that contains only that information.

usage:
python3 get_evalues_3dblast.py 3D-BLAST_RESULTS_FOLDER OUTPUT_FILE

'''

import sys
import os

results_dir =  sys.argv[1]
final_results_file = sys.argv[2]

results = []

arr = os.listdir(results_dir)
for file in arr:
    with open(results_dir + file) as input:
        check1 = 0
        check2 = 0
        id_q = file.split('.')[1]
        for line in input:
            if 'Sequences producing significant alignments:' in line:
                check1 = 1
            elif check1 == 1 and check2 == 1 and len(line.strip()) == 0:
                check1 = 0
            elif check1 == 1 and len(line.strip()) > 0:
                check2 = 1
                result_line = line.split('   ')
                id_t = result_line[0].split(' ')[0] .split('/')[-1].strip()
                evalue = result_line[-1].strip()
                if evalue[0] == 'e' or evalue[0] == 'E':
                    evalue = '1' + str(evalue)
                results.append(id_q + ',' + id_t + ',' + evalue + '\n')

file = open(final_results_file, 'w')

for evalue in results:
    file.write(evalue)
file.close()

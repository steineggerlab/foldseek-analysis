import argparse

import numpy as np
import re
import os
import sys

input = sys.argv[1]
tp_threshold = sys.argv[2]
hits = int(sys.argv[3])

query = ''
counter = 0
sens = 0
selec = 0
total_sens = 0
total_selec = 0
total_query = 0
total_ali = 0

with open(input, 'r') as f:
    for line in f:
        if line[0] == 'q':
            info_list = line.split(',')
            for info in info_list:
                lala = info.split('=')
                if lala[0] == 'query ':
                    query_id = lala[1]
                elif lala[0] == 'qLength ':
                    query_length = int(lala[1])
                else:
                    ali_length = int(lala[1])
            if query ==  '':
                query = query_id
            elif query != query_id:
                total_sens += sens
                total_selec += selec
                counter = 0
                sens = 0
                selec = 0
                query = query_id
        else:
            if counter == 0:
                total_query = total_query + query_length
            if query == query_id and counter < hits:
                lddt = line.strip().split(',')
                count = 0
                counter += 1
                ave_lddt = 0
                total_ali += ali_length
                for value in lddt:
                    if value == '':
                        continue
                    if float(value) >= float(tp_threshold):
                        count += 1
                if query_length > 0:
                    sens += count
                if ali_length > 0:
                    selec += count


total_sens += sens
total_selec += selec

total_sens = total_sens/total_query
total_selec = total_selec/total_ali

print(total_sens/hits, ",", total_selec)
# print("sens:", total_sens/hits)
# print("selec:", total_selec)

import argparse

import numpy as np
import re
import os
import sys

input = sys.argv[1]

line_count = 0
with open(input, 'r') as f:
    for line in f:
        line_count = line_count + 1
        ali = []
        bla = line.strip().split(',')
        for lddt in bla:
            if len(lddt) > 0:
                ali.append(float(lddt))
        if len(ali) > 0:
            average = sum(ali)/len(ali)
            print(average)

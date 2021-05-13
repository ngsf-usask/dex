#!/usr/bin/env python3

import sys
import re

orig = sys.argv[1] # original gtf file from ensembl

x=1
with open(orig, 'r') as infile:
    for line in infile:
        if line[0] == '#':
            print(line.strip())
        elif '\texon\t' not in line:
            print(line.strip())
        else:
            newl = re.sub('ENS...E...........',f'{x}',line.strip())
            print(newl)
            x+=1

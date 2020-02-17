#!/usr/bin/env python3

import pickle
import argparse
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description =
                                 '- Pickle to CSV  -',
                                 formatter_class=RawTextHelpFormatter)
parser.add_argument('-f', required= True, help='Input file')
ns = parser.parse_args()

file=ns.f
with open(file,'rb') as f:
    data = pickle.load(f)


for key in data:
    line=key
    #if 'rna' not in data[key]: continue
    #if 'GC%' not in data[key]: continue
    for field in data[key]:
        line+=","+str(data[key][field])
    print (line)
    

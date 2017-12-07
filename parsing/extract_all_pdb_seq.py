#!/usr/bin/env python

from parse_pdb import *
import sys
import operator
import numpy as np
from collections import defaultdict
import gzip

try:
    pdbfile = gzip.open(sys.argv[1], 'rb')
except:
    pdbfile = open(sys.argv[1], 'r')
code = get_acc(pdbfile)
if (code == 'XXXX' or code == '' or code == 'xxxx' or code == '    '):
    code=sys.argv[1]
pdbfile.close()
try:
    pdbfile = gzip.open(sys.argv[1], 'rb')
except:
    pdbfile = open(sys.argv[1], 'r')
name = sys.argv[1].replace(".pdb","_")
try:
    chain = sys.argv[2]
except:
    chains = get_all_chains(pdbfile)
else:
    chains = chain
print ":"+code+":"
print pdbfile
print chains
pdbfile.close()
#print chains
for chain in chains:
    print (code, chain)
    try:
        pdbfile = gzip.open(sys.argv[1], 'r')
    except:
        pdbfile = open(sys.argv[1], 'r')
    outname=code+chain+".fa"
    print outname
    outfile = open(outname,'w')
    outfile.write("> " + code + "_" + chain + "\n")
    outfile.write(get_atom_seq(pdbfile, chain) + "\n")
    pdbfile.close()
    outfile.close()
#pdbfile = open(sys.argv[1], 'r')
#print get_coordinates(pdbfile)
#pdbfile.close()

                                

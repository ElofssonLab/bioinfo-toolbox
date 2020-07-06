#!/usr/bin/env python

from parse_pdb_resid import *
import sys
import operator
import numpy as np
from collections import defaultdict
import gzip

filename=sys.argv[1]
if filename.endswith(".gz"):
    pdbfile = gzip.open(sys.argv[1], 'rb')
else:
    pdbfile = open(sys.argv[1], 'r')

code = get_acc(pdbfile)
if (code == 'XXXX' or code == '' or code == 'xxxx' or code == '    '):
    code=sys.argv[1]
pdbfile.close()
if filename.endswith(".gz"):
    pdbfile = gzip.open(sys.argv[1], 'rb')
else:
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
    if filename.endswith(".gz"):
        pdbfile = gzip.open(sys.argv[1], 'r')
    else:
        pdbfile = open(sys.argv[1], 'r')
    outname=name+chain+".fasta"
    print outname
    outfile = open(outname,'w')
    outfile.write("> " + code + "_" + chain + "\n")
    outfile.write(get_atom_seq(pdbfile, chain) + "\n")
    pdbfile.close()
    outfile.close()
#pdbfile = open(sys.argv[1], 'r')
#print get_coordinates(pdbfile)
#pdbfile.close()

                                

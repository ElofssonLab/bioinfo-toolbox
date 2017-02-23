#!/usr/bin/env python

from parse_pdb import *
import sys
import operator
import numpy as np
from collections import defaultdict



pdbfile = open(sys.argv[1], 'r')
#chain = sys.argv[2]
code = get_acc(pdbfile)
if (code == 'XXXX' or code == '' or code == 'xxxx' or code == '    '):
    code=sys.argv[1]
#print ":"+code+":"
pdbfile.close()
pdbfile = open(sys.argv[1], 'r')
name = sys.argv[1].replace(".pdb","_")
chains = get_all_chains(pdbfile)
pdbfile.close()
#print chains
for chain in chains:
    pdbfile = open(sys.argv[1], 'r')
    coord=read_chain(pdbfile, chain)
    pdbfile.close()
    outname=name+chain+".pdb"
    outfile = open(outname,'w')
    write(coord,outfile)

    #pdbfile = open(sys.argv[1], 'r')
#print get_coordinates(pdbfile)
#pdbfile.close()

                                

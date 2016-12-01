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
for chainA in chains:
    for chainB in chains:
        if chainA >= chainB:
            continue
        pdbfile = open(sys.argv[1], 'r')
        coordA=read_chain(pdbfile, chainA)
        pdbfile.close()
        pdbfile = open(sys.argv[1], 'r')
        coordB=read_chain(pdbfile, chainB)
        pdbfile.close()
        coord=[]
        coord.append(coordA[0])
        coord.append(coordA[1]+coordB[1])
        coord.append(coordA[2])
        outname=name+chainA+"-"+chainB+".pdb"
        outfile = open(outname,'w')
        write(coord,outfile)

    #pdbfile = open(sys.argv[1], 'r')
#print get_coordinates(pdbfile)
#pdbfile.close()

                                

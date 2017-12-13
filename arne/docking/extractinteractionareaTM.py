#!/usr/bin/env python
import sys
import operator
import numpy as np
import os
from collections import defaultdict




def parse_atm_record(line):

    record = defaultdict()
    record['name'] = line[0:6].strip()
    record['atm_no'] = int(line[6:11])
    record['atm_name'] = line[12:16].strip()
    record['res_name'] = line[17:20].strip()
    record['chain'] = line[21]
    record['res_no'] = int(line[22:26])
    record['insert'] = line[26].strip()
    record['x'] = float(line[30:38])
    record['y'] = float(line[38:46])
    record['z'] = float(line[46:54])
    record['occ'] = float(line[54:60])
    record['B'] = float(line[60:66])
    
    return record

def get_area(bvaluefile):
    file = open(bvaluefile, 'r')
    TMarea = 0.
    nonTMarea =0. 
    TMboundary=15.
    for line in file:
        if line.startswith('ATOM'):
            atm_record = parse_atm_record(line)
            if (atm_record['z'] > TMboundary or -1*atm_record['z'] > TMboundary):
                nonTMarea=nonTMarea+atm_record['B']
            else:
                TMarea=TMarea+atm_record['B']
    file.close()
    return (TMarea,nonTMarea)
        
if __name__ == '__main__':

    extension=".pdb.pdb-Bvalue.pdb"
    chaindir = "tm-chains/"
    pairdir = "tm-pairs/"
    minarea = 500.
    
    for pair in os.listdir(pairdir):
        if os.path.isfile(pairdir+pair) and pair.endswith(extension) and pair.find('*'):
            chainA=pair[0:4]+"_"+pair[5]
            chainB=pair[0:4]+"_"+pair[7]
            code=pair[0:4]
            codeA=pair[5]
            codeB=pair[7]
            (TMarea,nonTMarea)= get_area(pairdir+pair)
            (TMareaA,nonTMareaA)= get_area(chaindir+chainA+extension)
            (TMareaB,nonTMareaB)= get_area(chaindir+chainB+extension)
            TMdiff=float(TMareaA)+float(TMareaB)-float(TMarea)
            nonTMdiff=float(nonTMareaA)+float(nonTMareaB)-float(nonTMarea)
            if (TMdiff+nonTMdiff)>minarea:
                print chainA,chainB,TMarea,TMareaA,TMareaB,TMdiff,nonTMarea,nonTMareaA,nonTMareaB,nonTMdiff
            # Three inputs - com


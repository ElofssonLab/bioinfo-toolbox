#!/usr/bin/env python
import sys
import operator
import numpy as np
import os
from collections import defaultdict


def get_area(sasafile):
    file = open(sasafile, 'r')
    area = 0.
    for line in file:
        if line.startswith('Total'):
            (name,colon,area)=line.split()
            file.close()
            return area
            
if __name__ == '__main__':

    chaindir = "tm-chains/"
    pairdir = "tm-pairs/"
    minarea = 500.
    
    for pair in os.listdir(pairdir):
        if os.path.isfile(pairdir+pair) and pair.endswith('.sasa') and pair.find('*'):
            chainA=pair[0:4]+"_"+pair[5]
            chainB=pair[0:4]+"_"+pair[7]
            code=pair[0:4]
            codeA=pair[5]
            codeB=pair[7]
            area= get_area(pairdir+pair)
            areaA= get_area(chaindir+chainA+".sasa")
            areaB= get_area(chaindir+chainB+".sasa")
            diff=float(areaA)+float(areaB)-float(area)
            if diff>minarea:
                print chainA,chainB,area,areaA,areaB,diff
            # Three inputs - com


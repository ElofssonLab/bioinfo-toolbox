#!/usr/bin/env python
import sys
import operator
import numpy as np
import os
from collections import defaultdict


if __name__ == '__main__':

    listfile = open(sys.argv[1], 'r')
    cdhitfile = open(sys.argv[2], 'r')
    chains={}
    for line in cdhitfile:
        if line.startswith('>'):
            splitted=line.split()
            chains[splitted[1]]=1
            
    print chains
    lists={}       
    for line in listfile:
            splitted=line.split(' ',3)
            chainA=splitted[0]
            chainB=splitted[1]
            if (chainA in chains and chainB in chains):
                print line
            elif (chainA in chains and ( not chainA in lists )):
                print line
                lists[chainA]=1
            elif (chainB in chains and ( not chainB in lists )):
                print line
                lists[chainB]=1
#            else:
                # print "SKIP",line
                

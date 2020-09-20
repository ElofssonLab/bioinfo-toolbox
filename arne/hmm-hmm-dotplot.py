#!/usr/bin/env python3

import sys
import re
import matplotlib.pyplot as plt
import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import pandas as pd 

# Should calculate KL divergence
#from scipy.stats import entropy
#

# quick funcation fromthe stacjexchange
def KL(a, b):
    a = np.asarray(a, dtype=np.float)
    b = np.asarray(b, dtype=np.float)
    return np.sum(np.where(a != 0, a * np.log(a / b), 0))


# Gives a hmmerfile, returns array of
def parse_hmm(file):
    array=[]
    for line in open(file, 'r'):
        if re.match(r'^\s+\d+\s+\d+\.\d+',line):
            #print(line)
            vector=line.split()
            pos=int(vector[0])
            values=[]
            for i in vector[1:21]:
                values+=[np.exp(-1*float(i))]
            #print (pos,values)
            array+=[values]
    #print (array)
    nparray=np.array(array)
    return nparray
    
p = argparse.ArgumentParser(description = '- plotting trRosetta maps-',
                            formatter_class=RawTextHelpFormatter)
p.add_argument('-fileA','--inputA','-i', required= True, help='Input File A')
p.add_argument('-fileB','--inputB','-j', required= True, help='Input file B')
p.add_argument('-out','--output','-o', required= True, help='Output CSV file')
p.add_argument('-png','--plot','-p', required= False, help='Outplot plot file (optional)')

ns = p.parse_args()

#handleA = open(ns.inputA, 'r')
#handleB = open(ns.inputB, 'r')
#handleout = open(fileA, 'w')

a=parse_hmm(ns.inputA)
#print (a)
b=parse_hmm(ns.inputB)

res = np.zeros((len(a), len(b)))
for i in range(len(a)):
    for j in range(len(b)):
        #print (i,j,a[i],b[j],KL(a[i],b[j]))
        res[i,j]=np.log(KL(a[i],b[j]))

pd.DataFrame(res).to_csv(ns.output)
if (ns.plot):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #fig, (ax1, ax2) = plt.subplots(ncols=2)
    #ax2=ax.twin()
    cax = ax.matshow(res, cmap="hot")
    #print (res)
    ax.set(title="Dotplot between two HMMs")
    fig.colorbar(cax)
    fig.savefig(ns.plot)

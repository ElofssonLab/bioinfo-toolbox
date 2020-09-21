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


def hhscore(a,b,f):
    sum=0
    for i in np.arange(len(a)):
        #print (i,a[i],b[i],f[i])
        sum+=a[i]*b[i]/f[i]
    return np.log(sum)
    #return sum
# quick funcation fromthe stacjexchange
def KL(a, b):
    a = np.asarray(a, dtype=np.float)
    b = np.asarray(b, dtype=np.float)
    return np.sum(np.where(a != 0, a * np.log(a / b), 0))

def shannon_entropy(A, mode="auto", verbose=False):
    """
    https://stackoverflow.com/questions/42683287/python-numpy-shannon-entropy-array
    """
    A = np.asarray(A)

    # Determine distribution type
    if mode == "auto":
        condition = np.all(A.astype(float) == A.astype(int))
        if condition:
            mode = "discrete"
        else:
            mode = "continuous"
    if verbose:
        print(mode, file=sys.stderr)
    # Compute shannon entropy
    pA = A / A.sum()
    # Remove zeros
    pA = pA[np.nonzero(pA)[0]]
    if mode == "continuous":
        return -np.sum(pA*np.log2(A))  
    if mode == "discrete":
        return -np.sum(pA*np.log2(pA))   

def mutual_information(x,y, mode="auto", normalized=False):
    """
    I(X, Y) = H(X) + H(Y) - H(X,Y)
    https://stackoverflow.com/questions/20491028/optimal-way-to-compute-pairwise-mutual-information-using-numpy
    """
    x = np.asarray(x)
    y = np.asarray(y)
    # Determine distribution type
    if mode == "auto":
        condition_1 = np.all(x.astype(float) == x.astype(int))
        condition_2 = np.all(y.astype(float) == y.astype(int))
        if all([condition_1, condition_2]):
            mode = "discrete"
        else:
            mode = "continuous"

    H_x = shannon_entropy(x, mode=mode)
    H_y = shannon_entropy(y, mode=mode)
    H_xy = shannon_entropy(np.concatenate([x,y]), mode=mode)

    # Mutual Information
    I_xy = H_x + H_y - H_xy
    if normalized:
        return I_xy/np.sqrt(H_x*H_y)
    else:
        return  I_xy

    
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
p.add_argument('-log','--log','-l', required= False, help='use log',action='store_true',default=False)
p.add_argument('-norm','--normalize','-n', required= False, help='Normalize',action='store_true',default=False)
p.add_argument('-verb','--verbose','-v', required= False, help='verbose output',action='store_true',default=False)

p.add_argument('-win','--window','-w', required= False, default=0,help='Windows size for diagonal')
p.add_argument('-type','--type','-t', required= False, default="HH",help='Type (valid choises are KL, CC, MI, Shannon,HH')

ns = p.parse_args()

#handleA = open(ns.inputA, 'r')
#handleB = open(ns.inputB, 'r')
#handleout = open(fileA, 'w')

a=parse_hmm(ns.inputA)
#print (a)
b=parse_hmm(ns.inputB)

res = np.zeros((len(a), len(b)))
tiny=1.e-10


# AA freq A C D E  F G H I  K L M N  P Q R S  T V W Y
f=[0.074,0.033,0.059,0.058, 0.040,0.074,0.029,0.038,
   0.072,0.076,0.018,0.044,0.050, 0.037,0.042,0.081,
   0.062,0.068,0.013,0.033]


#   Alanine 7.4 %
#   Arginine 4.2 %
#   Asparagine 4.4 %
#   Aspartic Acid 5.9 %
#   Cysteine 3.3 %
#   Glutamic Acid 5.8 %
#   Glutamine 3.7 %
#   Glycine 7.4 %
#   Histidine 2.9 %
#   Isoleucine 3.8 %
#   Leucine 7.6 %
#   Lysine 7.2 %
#   Methionine 1.8 %
#   Phenylalanine 4.0 %
#   Proline 5.0 %
#   Serine 8.1 %
#   Threonine 6.2 %
#   Tryptophan 1.3 %
#   Tyrosine 3.3 %
#   Valine 6.8 %



#print (i,j,a[i],b[j],KL(a[i],b[j]))
#res[i,j]=np.log(tiny+KL(a[i],b[j]))
if (ns.type=="KL"):
    for i in range(len(a)):
        for j in range(len(b)):
            res[i,j]=KL(a[i],b[j])
elif (ns.type=="Shannon"):
    for i in range(len(a)):
        for j in range(len(b)):
            res[i,j]=shannon_entropy(a[i],b[j])
elif (ns.type=="MI"):
    for i in range(len(a)):
        for j in range(len(b)):
            res[i,j]=mutual_information(a[i],b[j])
elif(ns.type=="CC"):
    for i in range(len(a)):
        for j in range(len(b)):
            res[i,j]=np.corrcoef(a[i],b[j])[0,1]
elif(ns.type=="HH"):
    for i in range(len(a)):
        for j in range(len(b)):
            res[i,j]=hhscore(a[i],b[j],f)
            
else:sys.exit()
   
if ns.verbose:
    print (res)
# Normalise along diagonal (length W)
newres = np.zeros((len(a), len(b)))
if int(ns.window)>0:
    w=ns.window
    for i in range(len(a)):
        for j in range(len(b)):
            sum=0
            num=0
            for k in np.arange(-1*int(w),int(w)):
                if ((i+k)>=0)and ((j+k)>=0) and ((i+k)<len(a))and ((j+k)<len(b)):
                    sum+=res[i+k,j+k]
                    num+=1
            if (ns.log):
                newres[i,j] = np.log((tiny+sum)/num)
            else:
                newres[i,j] = sum/num
else:
    if (ns.log):
        newres=np.log(res)
    else:
        newres=res

if ns.verbose:
    print (newres)
if (ns.normalize):
    tempres=np.zeros((len(a), len(b)))
    row_sums = newres.sum(axis=1)
    col_sums = newres.sum(axis=0)
    for i in range(len(a)):
        for j in range(len(b)):
            tempres[i,j] = 2*newres[i,j] / (row_sums[i]/len(a)+col_sums[j]/len(b))
            if ns.verbose:
                print (i,j,newres[i,j],row_sums[i]/len(a),col_sums[j]/len(b))
    newres=tempres

    
if ns.verbose:
    print (newres)
pd.DataFrame(newres).to_csv(ns.output)
if (ns.plot):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #fig, (ax1, ax2) = plt.subplots(ncols=2)
    #ax2=ax.twin()
    cax = ax.matshow(newres, cmap="hot")
    #print (res)
    ax.set(title="Dotplot between two HMMs")
    ax.set_xlabel(ns.inputA)
    ax.set_ylabel(ns.inputB)
    fig.colorbar(cax)
    fig.savefig(ns.plot)
   
   

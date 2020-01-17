#!/usr/bin/env python3
import os
import h5py
import pickle
import random
import argparse
import tensorflow
import numpy as np
import pandas as pd
import seaborn as sb
import datetime as day
import matplotlib.pyplot as plt

if __name__ == '__main__':

    with open('outpred_rna.pickle','rb') as f:
        rnadata = pickle.load(f)

    with open('outpred_pro.pickle','rb') as f:
        prodata = pickle.load(f)

    with open('mobidata_K.pickle','rb') as f:
        mobidata = pickle.load(f)

    flist = open('formatted_list','r')
    disxgcR = {'A':{'dis1':[], 'gc':[]}, 'B':{'dis1':[], 'gc':[]}, 'E':{'dis1':[], 'gc':[]}}
    disxgcI = {'A':{'dis1':[], 'gc':[]}, 'B':{'dis1':[], 'gc':[]}, 'E':{'dis1':[], 'gc':[]}}
    for code in flist:
        fasta = open('FASTA/'+code.rstrip(),'r')
        iudis = open('disorder_IUpred/'+code.rstrip(),'r')
         
        fastaseq = ''
        for fline in fasta: 
            if fline.startswith('>'): continue
            fastaseq += fline.rstrip()

        iudisseq = ''
        for fline in iudis:
            if fline.startswith('#'): continue
            if float(fline.rstrip().split('\t')[2]) >= 0.4: iudisseq += 'D'
            else: iudisseq += 'S'

        pos = 0
        nounknown = 0
        disR = 0
        disI = 0
        first = mobidata[code.rstrip()]['disorder'][0][0]
        for el in mobidata[code.rstrip()]['disorder']:
            for n in range(el[0], el[1]+1):
                if el[2] == 'D': disR += 1
                if iudisseq[pos] == 'D': disI += 1
                pos += 1
                if el[2] != 'X' and el[2] != 'C': nounknown += 1

        if nounknown == 0: continue
        if mobidata[code.rstrip()]['kingdom'] not in ['A', 'B', 'E']: continue
        if mobidata[code.rstrip()]['GC%'] < 60 and mobidata[code.rstrip()]['GC%'] > 40: continue
        disxgcR[mobidata[code.rstrip()]['kingdom']]['gc'].append(mobidata[code.rstrip()]['GC%'])
        disxgcR[mobidata[code.rstrip()]['kingdom']]['dis1'].append((float(disR)/nounknown)*100)
        disxgcI[mobidata[code.rstrip()]['kingdom']]['gc'].append(mobidata[code.rstrip()]['GC%'])
        disxgcI[mobidata[code.rstrip()]['kingdom']]['dis1'].append((float(disI)/nounknown)*100)

    plot = 0
    fig, axes = plt.subplots(2, 4, sharex=True, sharey=True)
    for dic in [disxgcR, rnadata, prodata, disxgcI]:
        print (plot)
        print (len(dic['A']['gc']), len(dic['A']['dis1'])) 
        print (len(dic['B']['gc']), len(dic['B']['dis1'])) 
        print (len(dic['E']['gc']), len(dic['E']['dis1']))
        dfA = pd.DataFrame(dic['A'])
        dfB = pd.DataFrame(dic['B'])
        dfE = pd.DataFrame(dic['E'])
        
        sb.regplot(dfA['gc'], dfA['dis1'], color="blue", scatter_kws={'s':1}, ax=axes[0, plot], label='Arc')
        sb.regplot(dfB['gc'], dfB['dis1'], color="red", scatter_kws={'s':1}, ax=axes[0, plot], label='Bac')
        sb.regplot(dfE['gc'], dfE['dis1'], color="green", scatter_kws={'s':1}, ax=axes[0, plot], label='Euk')
        plt.ylim(0, 100)
        plt.xlim(0, 100)
        plt.legend()

        sb.kdeplot(dfA['gc'], dfA['dis1'], n_levels=6, ax=axes[1, plot], cmap="Blues")
        sb.kdeplot(dfB['gc'], dfB['dis1'], n_levels=6, ax=axes[1, plot], cmap="Reds")
        sb.kdeplot(dfE['gc'], dfE['dis1'], n_levels=6, ax=axes[1, plot], cmap="Greens")
        plt.ylim(0, 100)
        plt.xlim(0, 100)

        plot += 1 

    plt.show()


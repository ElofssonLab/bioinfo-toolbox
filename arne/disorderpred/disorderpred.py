#!/usr/bin/env python3 
import os
import h5py
import pickle
import random
import re
import argparse
import tensorflow as tf
import numpy as np
import pandas as pd
import seaborn as sb
import datetime as day
import nnlab as nl
from keras import backend as K
import matplotlib.pyplot as plt
from argparse import RawTextHelpFormatter
from keras.models import load_model, Model
import csv


#!/usr/bin/env python3
import pickle
import h5py
import numpy as np

codons = [
'ATA', 'ATC', 'ATT', 'ATG', 'ACA', 'ACC', 'ACG', 'ACT', 
'AAC', 'AAT', 'AAA', 'AAG', 'AGC', 'AGT', 'AGA', 'AGG',                  
'CTA', 'CTC', 'CTG', 'CTT', 'CCA', 'CCC', 'CCG', 'CCT', 
'CAC', 'CAT', 'CAA', 'CAG', 'CGA', 'CGC', 'CGG', 'CGT', 
'GTA', 'GTC', 'GTG', 'GTT', 'GCA', 'GCC', 'GCG', 'GCT', 
'GAC', 'GAT', 'GAA', 'GAG', 'GGA', 'GGC', 'GGG', 'GGT', 
'TCA', 'TCC', 'TCG', 'TCT', 'TTC', 'TTT', 'TTA', 'TTG', 
'TAC', 'TAT', 'TGC', 'TGT', 'TGG']

nuc_encode = {
'A':[1,0,0,0],
'T':[0,1,0,0],
'C':[0,0,1,0],
'G':[0,0,0,1]}

res_encode = {
'A':[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
'R':[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
'N':[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
'D':[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
'C':[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
'Q':[0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
'E':[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
'G':[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
'H':[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
'I':[0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
'L':[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
'K':[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
'M':[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
'F':[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
'P':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
'S':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
'T':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
'W':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
'Y':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
'V':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1]}

label = {'S':0, 'X':0.5, 'C':0.5, 'D':1}

##### codon wise one-hot encoding table
pos = 0
cod_encode = {}
for codon in codons: 
    cod_encode[codon] = []
    for n in range(61): 
        if n == pos: cod_encode[codon].append(1)
        else: cod_encode[codon].append(0)
    pos += 1
#####


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description =
    '- Neural network build manager -',
    formatter_class=RawTextHelpFormatter)
    #parser.add_argument('-t', required= True, help='path to training file list')
    #parser.add_argument('-d', required= True, help='path to data folder')
    parser.add_argument('-f', required= True, help='feature kind (pro, rna)')
    parser.add_argument('-m', required= True, help='model')
    parser.add_argument('-gc', required= False,  help=' GC', action='store_true')
    parser.add_argument('-seq', required= False,  help='sequence file (fasta format)', action='store_true')
    parser.add_argument('-kingdom', required= False,  help=' Kingdom', action='store_true')
    ns = parser.parse_args()



    
    rocstep=0.01
    cutoff =0.4 # Prediction cutoff.
    iupredcutoff =0.4 # Prediction cutoff.
    seed = 42
    tiny=1.e-10
    os.environ['PYTHONHASHSEED'] = '0'
    np.random.seed(seed)
    random.seed(seed)
    config = tf.ConfigProto()
    config.gpu_options.allow_growth = True
    tf.keras.backend.set_session(tf.Session(config=config))
    session_conf = tf.ConfigProto(intra_op_parallelism_threads=1, inter_op_parallelism_threads=1)
    tf.set_random_seed(seed)
    sess = tf.Session(graph=tf.get_default_graph(), config=session_conf)
    K.set_session(sess)

    if ns.f == 'pro': feat_len=20
    elif ns.f == 'rna': feat_len=61
    else: sys.exit('Unknown ')
    if ns.gc : feat_len+=1
    if ns.kingdom : feat_len+=3


    ##### Dataset #####
    print ('Loading data ...')
    data = h5py.File(ns.d+'formatted_data_GC_kingdom.h5py','r')
    with open(ns.d+'mobidata_K.pickle','rb') as f:
        gc = pickle.load(f)
    #test_list = []
    #for line in open(ns.t, 'r'): test_list.append(line.rstrip())

    model =ns.m


    pred = {}
    i=0
    
    for m in models:
    nt (m,models[m])
    
    _list = []
    line in open(models[m], 'r'): test_list.append(line.rstrip())
   
    l = load_model(m)

    t_cm = {}
    tive=np.zeros(101)
    tive=np.zeros(101)
    sitive=np.zeros(101)
    gative=np.zeros(101)
    
     thr in np.arange(0.,1.,rocstep): test_cm[thr] = test_cm.get(thr, {'PP':{'TP':0,'FP':0},'PN':{'TN':0,'FN':0}})

    # Prediction #####
    d = {'Name':[], 'kingdom':[], 'gc':[], 'TP':[], 'FP':[], 'FN':[], 'TN':[], 'Pred':[], 'Diso':[]}
    d in iupred
    is = {}
    
     code in test_list:
     iudis = open('disorder_IUpred/'+code.rstrip(),'r')
     iudis.append(code)
     for fline in iudis:
         if fline.startswith('#'): continue
         iudis[code].append(float(fline.rstrip().split('\t')[2]))
         #if float(fline.rstrip().split('\t')[2]) >= iupredcutoff: iudisseq += 'D'
         #else: iudisseq += 'S'
    
    protein in test_list:
    iudis = open('disorder_IUpred/'+protein.rstrip(),'r')
    iupred=[]
    for fline in iudis:
        if fline.startswith('#'): continue
        #if float(fline.rstrip().split('\t')[2]) >= iupredcutoff:
        #    iupred += 1
        #else:
        #    iupred += 0
        iupred.append( float(fline.rstrip().split('\t')[2]))
    iudis.close()
    #i+=1
    #if i>10: continue
    if ns.gc :
        if (ns.kingdom ):
            sample=np.concatenate((data[protein]['kingdom'],data[protein]['GC'],data[protein][ns.f]),axis=1)
        else:
            sample=np.concatenate((data[protein]['GC'],data[protein][ns.f]),axis=1)
    else:
        if (ns.kingdom ):
            sample=np.concatenate((data[protein]['kingdom'],data[protein]['GC'],data[protein][ns.f]),axis=1)
        else:
            sample = np.array(data[protein][ns.f], dtype=np.float64)
    X = sample[:,:-1].reshape(1, len(sample), len(sample[0])-1)
    Y = sample[:,-1]
    #print (protein,X,Y)
    prediction = model.predict_on_batch(X)

    ##### Confusion Matrix #####  (This is way too slow to implement)..
    #for thr in test_cm:
    #    for pos in range(len(prediction[0])):
    #        if Y[pos] == 0.5: continue
    #        if prediction[0][pos][0] >= thr:
    #            if Y[pos] == 1: test_cm[thr]['PP']['TP'] += 1
    #            else: test_cm[thr]['PP']['FP'] += 1
    #        else:
    #            if Y[pos] == 1: test_cm[thr]['PN']['FN'] += 1
    #            else: test_cm[thr]['PN']['TN'] += 1
    

    for pos in range(len(prediction[0])):
        if Y[pos] == 0.5: continue
        thr=int(tiny+prediction[0][pos][0]/rocstep)
        iuthr=int(tiny+iupred[pos]/rocstep)
        if Y[pos] == 1:
            positive[thr]+=1
            iupositive[iuthr]+=1
        else:
            negative[thr]+=1
            iunegative[iuthr]+=1

            ##### Data selection #####
    if gc[protein]['kingdom'] not in ['A', 'B', 'E']: continue
    if (not gc[protein]['GC%']): continue 
    
    #acc1 = 0
    TP=0
    FP=0
    FN=0
    TN=0
    nounknown = 0
    iudiso=0
    for pos in range(len(prediction[0])):
        if Y[pos] == 0.5: continue
        #if prediction[0][pos][0] >= cutoff: acc1 += 1
        if prediction[0][pos][0] >= cutoff:
            if Y[pos] == 1: TP += 1
            else: FP += 1
        else:
            if Y[pos] == 1: FN += 1
            else: TN += 1
        nounknown += 1
        if iupred[pos] > iupredcutoff:
            iudiso+=1
    if nounknown == 0: continue
    else:
        pred[protein]=[]
        pred[protein].append(protein)
        pred[protein].append(gc[protein]['kingdom'])
        pred[protein].append(gc[protein]['GC%'])
        pred[protein].append(TP/nounknown)
        pred[protein].append(FP/nounknown)
        pred[protein].append(FN/nounknown)
        pred[protein].append(TN/nounknown)
        TPR=TP/(tiny+TP+FN)
        FPR=FP/(tiny+FP+TN)
        Spec=TN/(tiny+FP+TN)
        PPV=TP/(tiny+FP+TP)
        F1=2*TP/(tiny+2*TP+FP+FN)
        MCC=((TP*TN)-FP*FN)/np.sqrt(tiny+(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
        pred[protein].append(TPR)
        pred[protein].append(FPR)
        pred[protein].append(Spec)
        pred[protein].append(PPV)
        pred[protein].append(F1)
        pred[protein].append(MCC)
        pred[protein].append((TP+FP)/nounknown)
        pred[protein].append((TP+FN)/nounknown)
        pred[protein].append(iudiso/nounknown)


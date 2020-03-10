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

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description =
    '- Neural network build manager -',
    formatter_class=RawTextHelpFormatter)
    parser.add_argument('-t', required= True, help='path to training file list')
    parser.add_argument('-d', required= True, help='path to data folder')
    parser.add_argument('-f', required= True, help='feature kind (pro, rna)')
    parser.add_argument('-m', required= True, help='model')
    parser.add_argument('-gc', required= False,  help=' GC', action='store_true')
    ns = parser.parse_args()
    if (ns.gc): ns.gc='GC'
    else: ns.gc='noGC'
    cutoff =0.4 # Prediction cutoff.
    iupredcutoff =0.4 # Prediction cutoff.
    seed = 42
    tiny=1.e-20
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
    if ns.gc == 'GC': feat_len+=1

    
    ##### Dataset #####
    print ('Loading data ...')
    data = h5py.File(ns.d+'formatted_data_GC.h5py','r')
    with open(ns.d+'mobidata_K.pickle','rb') as f:
        gc = pickle.load(f)
    test_list = []
    for line in open(ns.t, 'r'): test_list.append(line.rstrip())

    model = load_model(ns.m)

    test_cm = {}
    for thr in nl.thrlist: test_cm[thr] = test_cm.get(thr, {'PP':{'TP':0,'FP':0},'PN':{'TN':0,'FN':0}})

# This is not complete as we have to check if we have the data for that residue as well
    ##### Prediction #####
    #pred = {'Name':[], 'kingdom':[], 'gc':[], 'TP':[], 'FP':[], 'FN':[], 'TN':[], 'Pred':[], 'Diso':[]}
    pred = {}
    i=0

    for protein in test_list:
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
        if ns.gc == 'GC':
            sample=np.concatenate((data[protein]['GC'],data[protein][ns.f]),axis=1)
        else:
            sample = np.array(data[protein][ns.f], dtype=np.float64)
        X = sample[:,:-1].reshape(1, len(sample), len(sample[0])-1)
        Y = sample[:,-1]
        print (X,Y)
        prediction = model.predict_on_batch(X)

    ##### Confusion Matrix #####
        for thr in test_cm:
            for pos in range(len(prediction[0])):
                if Y[pos] == 0.5: continue
                if prediction[0][pos][0] >= thr:
                    if Y[pos] == 1: test_cm[thr]['PP']['TP'] += 1
                    else: test_cm[thr]['PP']['FP'] += 1
                else:
                    if Y[pos] == 1: test_cm[thr]['PN']['FN'] += 1
                    else: test_cm[thr]['PN']['TN'] += 1
                    
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
        protlength=len(prediction[0])
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
            pred[protein].append(protlength)

    model=re.sub(r'.*\/','',ns.m)

    with open('predictions/outpred_'+model+'.pickle','wb') as f:
        pickle.dump(pred, f)
    #print (pred)
    keys=[['Name', 'kingdom', 'gc', 'TP', 'FP', 'FN', 'TN','TPR','FPR','Spec','PPV','F1','MCC', 'Pred', 'Diso','IUPRED']]
    
    with open('predictions/outpred_'+model+'.csv','w',newline="") as f:
        w = csv.writer(f)
        w.writerows(keys)
        for prot in pred:
            print (prot)
            print (pred[prot])
            w.writerows([pred[prot]])
    #print (disxgc)
    #print (test_cm.keys())
    #print (test_cm.values())
    #print (test_cm)
    rockeys=[['Thr','TP','FP','TN','TPR','FPR','Spec','PPV','F1','MCC']]
    with open('predictions/outpred_'+model+'.roc','w',newline="") as f:
        w = csv.writer(f)
        w.writerows(rockeys)
        for thr in test_cm.keys():
            line=[]
            line.append(thr)
            for k in ["TP","FP"]:
                line.append(test_cm[thr]['PP'][k])
            for k in ["TN","FN"]:
                line.append(test_cm[thr]['PN'][k])

            TP=test_cm[thr]['PP']['TP']
            FP=test_cm[thr]['PP']['FP']
            TN=test_cm[thr]['PN']['TN']
            FN=test_cm[thr]['PN']['FN']
            TPR=TP/(tiny+TP+FN)
            FPR=FP/(tiny+FP+TN)
            Spec=TN/(tiny+FP+TN)
            PPV=TP/(tiny+FP+TP)
            F1=2*TP/(tiny+2*TP+FP+FN)
            MCC=((TP*TN)-FP*FN)/np.sqrt(tiny+(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
            line.append(TPR)
            line.append(FPR)
            line.append(Spec)
            line.append(PPV)
            line.append(F1)
            line.append(MCC)
            #print (line)    
            w.writerows([line])
    print (test_cm)

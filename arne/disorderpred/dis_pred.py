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
    formatter_class=RawTextHelpFormatter)
    #parser.add_argument('-t', required= True, help='path to training file list')
    #parser.add_argument('-d', required= True, help='path to data folder')
    parser.add_argument('-f', required= True, help='feature kind (pro, rna)')
    #parser.add_argument('-m', required= True, help='model')
    #parser.add_argument('-final', required= False,  help=' Use Final Model', action='store_true')
    parser.add_argument('-gc', required= False,  help=' GC', action='store_true')
    parser.add_argument('-gcgenomic', required= False,  help=' GCgenomic', action='store_true')
    parser.add_argument('-kingdom', required= False,  help=' Kingdom', action='store_true')
    ns = parser.parse_args()
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
    if ns.gc: feat_len+=1
    if ns.gcgenomic:  feat_len+=1

    
    ##### Dataset #####
    #print ('Loading data ...')
    #data = h5py.File(ns.d+'formatted_data_GC.h5py','r')
    #with open(ns.d+'mobidata_K.pickle','rb') as f:
    #    gc = pickle.load(f)
    #test_list = []
    #for line in open(ns.t, 'r'): test_list.append(line.rstrip())

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
        #print (protein,X,Y)
        prediction = model.predict_on_batch(X)


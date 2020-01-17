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
import nnlab as nl
from keras import backend as K
import matplotlib.pyplot as plt
from argparse import RawTextHelpFormatter
from keras.models import load_model, Model

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description =
    '- Neural network build manager -',
    formatter_class=RawTextHelpFormatter)
    parser.add_argument('-t', required= True, help='path to training file list')
    parser.add_argument('-d', required= True, help='path to data folder')
    parser.add_argument('-f', required= True, help='feature kind (pro, rna)')
    parser.add_argument('-m', required= True, help='model')

    ns = parser.parse_args()

    seed = 42
    os.environ['PYTHONHASHSEED'] = '0'
    np.random.seed(seed)
    random.seed(seed)
    config = tensorflow.ConfigProto()
    config.gpu_options.allow_growth = True
    tensorflow.keras.backend.set_session(tensorflow.Session(config=config))
    session_conf = tensorflow.ConfigProto(intra_op_parallelism_threads=1, inter_op_parallelism_threads=1)
    tensorflow.set_random_seed(seed)
    sess = tensorflow.Session(graph=tensorflow.get_default_graph(), config=session_conf)
    K.set_session(sess)

    if ns.f == 'pro': feat_len=20
    else: feat_len=61

    ##### Dataset #####
    print ('Loading data ...')
    data = h5py.File(ns.d+'formatted_data.h5py','r')
    with open(ns.d+'mobidata_K.pickle','rb') as f:
        gc = pickle.load(f)
    test_list = []
    for line in open(ns.t, 'r'): test_list.append(line.rstrip())

    model = load_model(ns.m)

    test_cm = {}
    for thr in nl.thrlist: test_cm[thr] = test_cm.get(thr, {'PP':{'TP':0,'FP':0},'PN':{'TN':0,'FN':0}})

    ##### Prediction #####
    disxgc = {'A':{'dis1':[], 'gc':[]}, 'B':{'dis1':[], 'gc':[]}, 'E':{'dis1':[], 'gc':[]}}
    for protein in test_list:
        sample = np.array(data[protein][ns.f], dtype=np.float64)
        X = sample[:,:-1].reshape(1, len(sample), len(sample[0])-1)
        Y = sample[:,-1]
        #print (protein,X,Y)
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
        if gc[protein]['GC%'] < 60 and gc[protein]['GC%'] > 40: continue 
        
        acc1 = 0
        nounknown = 0
        for pos in range(len(prediction[0])):
            if Y[pos] == 0.5: continue
            if prediction[0][pos][0] >= 0.4: acc1 += 1
            nounknown += 1

        if nounknown == 0: continue
        else: disxgc[gc[protein]['kingdom']]['dis1'].append((acc1/nounknown)*100)
        disxgc[gc[protein]['kingdom']]['gc'].append(gc[protein]['GC%'])

    with open('outpred_'+ns.f+'.pickle','wb') as f:
        pickle.dump(disxgc, f)


#!/usr/bin/env python3
import os
import h5py
import re
import pickle
import random
import argparse
import numpy as np
import datetime as day
import nnlab as nl
from argparse import RawTextHelpFormatter
import tensorflow as tf
from keras import losses
from keras.layers import *
from keras import backend as K
from keras.regularizers import l2
from keras.optimizers import Adam
from keras.models import load_model, Model
from keras.engine.topology import Input, Layer
from keras.layers.advanced_activations import ELU
from multiprocessing.pool import ThreadPool

def model_set(feat_len, units, act, reg, drp, bn=False):
    inlayer = Input(shape=(None, feat_len))
    a = nl.bi_LSTM(units, reg, drp, True, bn, inlayer)
    a = nl.bi_LSTM(units, reg, drp, True, bn, a)
    out = Dense(1, activation='sigmoid', use_bias='True')(a)
    model = Model(inputs=inlayer, outputs=out)
    adam = Adam(lr=0, beta_1=0.9, beta_2=0.999, epsilon=None, decay=0, amsgrad=False)
    model.summary()
    model.compile(loss='binary_crossentropy', optimizer = adam, sample_weight_mode='temporal')
    return model

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description =
    '- Neural network build manager -',
    formatter_class=RawTextHelpFormatter)
    parser.add_argument('-t', required= True, help='path to training file list')
    parser.add_argument('-v', required= True, help='path to validation file list')
    parser.add_argument('-d', required= True, help='path to data folder')
    parser.add_argument('-f', required= True, help='feature kind (pro, rna)')
    parser.add_argument('-gc', required= False,  help=' GC', action='store_true')
    parser.add_argument('-gcgenomic', required= False,  help=' GCgenomic', action='store_true')
    parser.add_argument('-kingdom', required= False,  help=' Kingdom', action='store_true')

    parser.add_argument('-ep', required= False, default= '100', help='epoch number')
    parser.add_argument('-bs', required= False, default= '10', help='mini-batch size')
    parser.add_argument('-lr', required= False, default= '0.001', help='learning rate')

    parser.add_argument('-id', required= False, default= '1', help='train model id')
    ns = parser.parse_args()
    type=''
    if (ns.gc):
        type+='_GC'
    if(ns.gcgenomic):
        type+='_GCgenomic'

    if (ns.kingdom): ns.kingdom='KINGDOM_'
    else: ns.kingdom=''
    tiny=1.e-10
    testset=re.sub(r'.*\/','',ns.t)
    seed = 42
    thr=0.4
    seed = random.randint(1,99)
    os.environ['PYTHONHASHSEED'] = '0'
    np.random.seed(seed)
    random.seed(seed)
    config = tf.ConfigProto()
    config.gpu_options.allow_growth = True
    config.log_device_placement=True
    tf.keras.backend.set_session(tf.Session(config=config))
    session_conf = tf.ConfigProto(intra_op_parallelism_threads=1, inter_op_parallelism_threads=1)
    session_conf.log_device_placement=True
    tf.set_random_seed(seed)
    sess = tf.Session(graph=tf.get_default_graph(), config=session_conf)
    K.set_session(sess)

    hrun_id = str(ns.ep)+"-"+testset+"-"+str(ns.bs)+'-'+str(ns.lr)+'_'+str(ns.id)+'_'+str(ns.f)+str(type)+"_"+str(ns.kingdom)+str(seed)

    epochs = int(ns.ep)
    batch = int(ns.bs)
    LR = float(ns.lr)
    act = 'relu'
    reg = 0.000000001
    drp = 0.3
    if ns.f == 'pro': feat_len=20
    elif ns.f == 'rna': feat_len=61
    else: sys.exit('Unknown ')
    if (ns.gc): feat_len+=1
    if (ns.gcgenomic): feat_len+=1
    if ns.kingdom == 'KINGDOM_': feat_len+=3

    
    ##### Dataset #####
    print ('Loading data ...')
    data = h5py.File(ns.d,'r')
    train_list = []
    for line in open(ns.t, 'r'): train_list.append(line.rstrip())
    val_list = []
    for line in open(ns.v, 'r'): val_list.append(line.rstrip())

    model = model_set(feat_len, 64, act, reg, drp)

    ##### Epoch cycle #####
    print ('Start training ...')
    mb = []
    best = 0
    K.set_value(model.optimizer.lr, LR)
    for n in range(1, epochs+1):

        score = 0
        epoch_losses = []
        print('\nEpoch '+str(n)+'/'+str(epochs)+' --- Base LR >>> '+str(K.get_value(model.optimizer.lr)))
        run_id = str(n)+'-'+str(ns.bs)+'-'+str(ns.lr)+'_'+str(ns.id)
        ##### Batch formatting #####
        for protein in train_list:
            # This is actually too slow...
            mb.append(np.array(data[protein][ns.f], dtype=np.float64))
                
            # Z is not used - skip
            if len(mb) > batch-1:
                mbX = []
                mbY = []
                mbZ = []
                timesteps = 0
                for sample in mb: timesteps = max(len(sample), timesteps)
                for sample in mb:
                    while len(sample) < timesteps: 
                        sample = np.append(sample, np.zeros(len(sample[0])).reshape(1,len(sample[0])), axis=0)
                        sample[-1][-1] = 0.5
                    X = sample[:,:-1]
                    mbX.append(X)
                    Y = sample[:,-1]
                    mbY.append(Y)
                    #Z = []
                    #for value in Y: 
                    #    if value == 1: Z.append(10)
                    #    if value == 0: Z.append(1)
                    #    if value == 0.5: Z.append(0)
                    #mbZ.append(Z)
                mbX = np.array(mbX, dtype=np.float64).reshape(batch, timesteps, len(X[0]))
                mbY = np.array(mbY, dtype=np.float64).reshape(batch, timesteps, 1)
                #mbZ = np.array(mbZ, dtype=np.float64).reshape(batch, timesteps, 1)
                ##### Model train #####
                loss = model.train_on_batch(mbX, mbY)
                epoch_losses.append(loss)
                mb = []

                
        ##### Prediction/evaluation over validation set #####
        print ('Epoch '+str(n)+' complete! Evaluation...')
        TP=0
        FP=0
        TN=0
        FN=0
        for protein in val_list:
            sample = np.array(data[protein][ns.f], dtype=np.float64)
            X = sample[:,:-1].reshape(1, len(sample), len(sample[0])-1)
            Y = sample[:,-1]
            prediction = model.predict_on_batch(X)
            for pos in range(len(prediction[0])):
                if Y[pos] == 0.5: continue
                if prediction[0][pos][0] >= thr: 
                    if Y[pos] == 1:
                        TP+=1
                    else:
                        FP+=1
                else:
                    if Y[pos] == 1:
                        FN+=1
                    else:
                        TN+=1

        #TPR=TP/(tiny+TP+FN)
        #FPR=FP/(tiny+FP+TN)
        F1=2*TP/(tiny+2*TP+FP+FN)
        #score = np.trapz(TPR,FPR)
        score = F1

        ##### Model save, stats display #####
        acc = 0
        for el in epoch_losses: acc += el

        outfile = open('logs/history_'+hrun_id+".log",'a')
        if score >= best:
            model.save('models/model_'+hrun_id+".ann")
            print ('Loss: '+str(acc/len(epoch_losses))+' -  Val score: '+str(score)+' >>> Partial saved!')
            outfile.write('Epoch_'+str(n)+' Loss: '+str(acc/len(epoch_losses))+' - Val score: '+str(score)+' >>> Partial saved!\n')
            best = score
        else:
            print ('Loss: '+str(acc/len(epoch_losses))+' -  Val score: '+str(score))
            outfile.write('Epoch_'+str(n)+' Loss: '+str(acc/len(epoch_losses))+' - Val score: '+str(score)+'\n')
            K.set_value(model.optimizer.lr, K.get_value(model.optimizer.lr)*0.9)
        outfile.close()
    model.save('models/model_final_'+hrun_id+".ann")
    print ('Loss: '+str(acc/len(epoch_losses))+' -  Val score: '+str(score))
    outfile.write('Epoch_'+str(n)+' Loss: '+str(acc/len(epoch_losses))+' - Val score: '+str(score)+'\n')
    outfile.close()


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
    #parser.add_argument('-t', required= True, help='path to training file list')
    parser.add_argument('-d', required= True, help='path to data folder')
    parser.add_argument('-f', required= True, help='feature kind (pro, rna)')
    #parser.add_argument('-m', required= True, help='model')
    parser.add_argument('-final', required= False,  help=' Use Final Model', action='store_true')
    parser.add_argument('-gc', required= False,  help=' GC', action='store_true')
    parser.add_argument('-gcgenomic', required= False,  help=' GCgenomic', action='store_true')
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
    if (ns.f =='pro'):
        field='_pro'
    elif (ns.f =='rna'):
        field='_rna'
    if (ns.gc):
        field+="_GC"
    if(ns.gcgenomic):
        field+="_GCgenomic"
    if (ns.kingdom):
        field+="_kingdom"

    print (field)
    data = h5py.File(ns.d+'formatted_data'+field+".h5py",'r')
    with open(ns.d+'mobidata_K.pickle','rb') as f:
        gc = pickle.load(f)
    #test_list = []
    #for line in open(ns.t, 'r'): test_list.append(line.rstrip())
    anndir="models/"
    if (ns.f =='pro'):
        if (ns.gcgenomic):
            #models={"model_100-test1234.sorted.GCgenomic-10-0.001_1_pro_GCgenomic_83.ann":"DIScv/test5.sorted.GCgenomic",
            #        "model_100-test1235.sorted.GCgenomic-10-0.001_1_pro_GCgenomic_95.ann":"DIScv/test4.sorted.GCgenomic",
            #        "model_100-test1245.sorted.GCgenomic-10-0.001_1_pro_GCgenomic_25.ann":"DIScv/test3.sorted.GCgenomic",
            #        "model_100-test1345.sorted.GCgenomic-10-0.001_1_pro_GCgenomic_80.ann":"DIScv/test2.sorted.GCgenomic",
            #        "model_100-test2345.sorted.GCgenomic-10-0.001_1_pro_GCgenomic_6.ann" :"DIScv/test1.sorted.GCgenomic"
            #}
            models={"model_100-test1234.sorted.GCgenomic-10-0.001_1_pro_GCgenomic_95.ann":"DIScv/test5.sorted.GCgenomic",
                    "model_100-test1235.sorted.GCgenomic-10-0.001_1_pro_GCgenomic_4.ann":"DIScv/test4.sorted.GCgenomic",
                    "model_100-test1245.sorted.GCgenomic-10-0.001_1_pro_GCgenomic_34.ann":"DIScv/test3.sorted.GCgenomic",
                    "model_100-test1345.sorted.GCgenomic-10-0.001_1_pro_GCgenomic_18.ann":"DIScv/test2.sorted.GCgenomic",
                    "model_100-test2345.sorted.GCgenomic-10-0.001_1_pro_GCgenomic_43.ann" :"DIScv/test1.sorted.GCgenomic"
            }
            if (ns.final):
                models={"model_final_100-test1234.sorted.GCgenomic-10-0.001_1_pro_GCgenomic_95.ann":"DIScv/test5.sorted.GCgenomic",
                        "model_final_100-test1235.sorted.GCgenomic-10-0.001_1_pro_GCgenomic_4.ann":"DIScv/test4.sorted.GCgenomic",
                        "model_final_100-test1245.sorted.GCgenomic-10-0.001_1_pro_GCgenomic_34.ann":"DIScv/test3.sorted.GCgenomic",
                        "model_final_100-test1345.sorted.GCgenomic-10-0.001_1_pro_GCgenomic_18.ann":"DIScv/test2.sorted.GCgenomic",
                        "model_final_100-test2345.sorted.GCgenomic-10-0.001_1_pro_GCgenomic_43.ann" :"DIScv/test1.sorted.GCgenomic"
                }
            if (ns.gc):
                models={"model_100-test1234.sorted.GCgenomic-10-0.001_1_pro_GC_GCgenomic_62.ann":"DIScv/test5.sorted.GCgenomic",
                        "model_100-test1235.sorted.GCgenomic-10-0.001_1_pro_GC_GCgenomic_75.ann":"DIScv/test4.sorted.GCgenomic",
                        "model_100-test1245.sorted.GCgenomic-10-0.001_1_pro_GC_GCgenomic_63.ann":"DIScv/test3.sorted.GCgenomic",
                        "model_100-test1345.sorted.GCgenomic-10-0.001_1_pro_GC_GCgenomic_53.ann":"DIScv/test2.sorted.GCgenomic",
                        "model_100-test2345.sorted.GCgenomic-10-0.001_1_pro_GC_GCgenomic_88.ann" :"DIScv/test1.sorted.GCgenomic"
                }
        elif (ns.gc):
            models={"model_100-test1234.sorted.GCgenomic-10-0.001_1_pro_GC_16.ann":"DIScv/test5.sorted.GCgenomic",
                    "model_100-test1235.sorted.GCgenomic-10-0.001_1_pro_GC_21.ann":"DIScv/test4.sorted.GCgenomic",
                    "model_100-test1245.sorted.GCgenomic-10-0.001_1_pro_GC_86.ann":"DIScv/test3.sorted.GCgenomic",
                    "model_100-test1345.sorted.GCgenomic-10-0.001_1_pro_GC_68.ann":"DIScv/test2.sorted.GCgenomic",
                    "model_100-test2345.sorted.GCgenomic-10-0.001_1_pro_GC_93.ann":"DIScv/test1.sorted.GCgenomic"
            }
            if (ns.final):
                models={"model_final_100-test1234.sorted.GCgenomic-10-0.001_1_pro_GC_16.ann":"DIScv/test5.sorted.GCgenomic",
                        "model_final_100-test1235.sorted.GCgenomic-10-0.001_1_pro_GC_21.ann":"DIScv/test4.sorted.GCgenomic",
                        "model_final_100-test1245.sorted.GCgenomic-10-0.001_1_pro_GC_86.ann":"DIScv/test3.sorted.GCgenomic",
                        "model_final_100-test1345.sorted.GCgenomic-10-0.001_1_pro_GC_68.ann":"DIScv/test2.sorted.GCgenomic",
                        "model_final_100-test2345.sorted.GCgenomic-10-0.001_1_pro_GC_93.ann":"DIScv/test1.sorted.GCgenomic"
                }
        else:
            models={"model_100-test1234.sorted.GCgenomic-10-0.001_1_pro_54.ann":"DIScv/test5.sorted.GCgenomic",
                    "model_100-test1235.sorted.GCgenomic-10-0.001_1_pro_50.ann":"DIScv/test4.sorted.GCgenomic",
                    "model_100-test1245.sorted.GCgenomic-10-0.001_1_pro_11.ann":"DIScv/test3.sorted.GCgenomic",
                    "model_100-test1345.sorted.GCgenomic-10-0.001_1_pro_49.ann":"DIScv/test2.sorted.GCgenomic",
                    "model_100-test2345.sorted.GCgenomic-10-0.001_1_pro_26.ann":"DIScv/test1.sorted.GCgenomic" 
            }
    elif(ns.f=='rna'):
        if (ns.gcgenomic):
            models={
            }
        elif (ns.gc):
            models={
            }
        else:
            models={
            }


    pred = {}
    i=0
    
    for m in models:
        #print (m,models[m])
        
        test_list = []
        for line in open(models[m], 'r'): test_list.append(line.rstrip())
    
        model = load_model(anndir+m)

        #test_cm = {}
        positive=np.zeros(101)
        negative=np.zeros(101)
        iupositive=np.zeros(101)
        iunegative=np.zeros(101)
        
        #for thr in np.arange(0.,1.,rocstep): test_cm[thr] = test_cm.get(thr, {'PP':{'TP':0,'FP':0},'PN':{'TN':0,'FN':0}})

        ##### Prediction #####
        #pred = {'Name':[], 'kingdom':[], 'gc':[], 'TP':[], 'FP':[], 'FN':[], 'TN':[], 'Pred':[], 'Diso':[]}
        #read in iupred
        #iudis = {}
        #
        #for code in test_list:
        #    iudis = open('disorder_IUpred/'+code.rstrip(),'r')
        #    iudis.append(code)
        #    for fline in iudis:
        #        if fline.startswith('#'): continue
        #        iudis[code].append(float(fline.rstrip().split('\t')[2]))
        #        #if float(fline.rstrip().split('\t')[2]) >= iupredcutoff: iudisseq += 'D'
        #        #else: iudisseq += 'S'
    
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
                pred[protein].append(gc[protein]['GC%'])  # Replace with GC genomic?
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
    #mod=re.sub(r'.*\/','',m)
    #set=re.sub(r'.*\/','',ns.t)

        
    with open('predictions/outpred'+field+'.pickle','wb') as f:
        pickle.dump(pred, f)
    #print (pred)
    keys=[['Name', 'kingdom', 'gc', 'TP', 'FP', 'FN', 'TN','TPR','FPR','Spec','PPV','F1','MCC', 'Pred', 'Diso','IUPRED','Length']]
    
    with open('predictions/outpred'+field+'.csv','w',newline="") as f:
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

    #print(positive,negative)
    
    rockeys=[['Thr','TP','FP','FN','TN','TPR','FPR','Spec','PPV','F1','MCC',
              'iuTP','iuFP','iuFN','iuTN','iuTPR','iuFPR','iuSpec','iuPPV','iuF1','iuMCC',]]
    with open('predictions/outpred'+field+'.roc','w',newline="") as f:
        w = csv.writer(f)
        w.writerows(rockeys)
        TP=np.sum(positive)
        FP=np.sum(negative)
        TN=0
        FN=0
        iuTP=np.sum(positive)
        iuFP=np.sum(negative)
        iuTN=0
        iuFN=0
        for thr in np.arange(0.,1.+tiny,rocstep):
            line=[]
            line.append(thr)
            index=int(tiny+thr/rocstep)
            TP-=positive[index]
            FP-=negative[index]
            TN+=negative[index]
            FN+=positive[index]
            line.append(TP)            
            line.append(FP)            
            line.append(FN)            
            line.append(TN)            
            TPR=TP/(tiny+TP+FN)
            FPR=FP/(tiny+FP+TN)
            Spec=TN/(tiny+FP+TN)
            PPV=TP/(tiny+FP+TP)
            F1=2*TP/(tiny+2*TP+FP+FN)
            #print (thr,index,positive[index],negative[index],TP,FP,FN,TN,TP+FP,np.sum(positive),
            #       np.sum(positive[0:index]),TN+FN,np.sum(negative),np.sum(negative[0:index]))
            MCC=((TP*TN)-FP*FN)/np.sqrt(tiny+(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
            line.append(TPR)
            line.append(FPR)
            line.append(Spec)
            line.append(PPV)
            line.append(F1)
            line.append(MCC)


            iuTP-=iupositive[index]
            iuFP-=iunegative[index]
            iuTN+=iunegative[index]
            iuFN+=iupositive[index]
            line.append(iuTP)            
            line.append(iuFP)            
            line.append(iuFN)            
            line.append(iuTN)            
            iuTPR=iuTP/(tiny+iuTP+iuFN)
            iuFPR=iuFP/(tiny+iuFP+iuTN)
            iuSpec=iuTN/(tiny+iuFP+iuTN)
            iuPPV=iuTP/(tiny+iuFP+iuTP)
            iuF1=2*iuTP/(tiny+2*iuTP+iuFP+iuFN)
            #print (thr,index,positive[index],negative[index],TP,FP,FN,TN,TP+FP,np.sum(positive),
            #       np.sum(positive[0:index]),TN+FN,np.sum(negative),np.sum(negative[0:index]))
            iuMCC=((iuTP*iuTN)-iuFP*iuFN)/np.sqrt(tiny+(iuTP+iuFP)*(iuTP+iuFN)*(iuTN+iuFP)*(iuTN+iuFN))
            line.append(iuTPR)
            line.append(iuFPR)
            line.append(iuSpec)
            line.append(iuPPV)
            line.append(iuF1)
            line.append(iuMCC)

            #print (line)    
            w.writerows([line])
    #print (test_cm)

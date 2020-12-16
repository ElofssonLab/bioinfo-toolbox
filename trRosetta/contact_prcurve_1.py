#!/usr/bin/env python3:q
import os
import sys
import math
import argparse
import numpy as np
from Bio.PDB import *
import matplotlib.pyplot as plt
from sklearn.metrics import auc
from scipy.stats import pearsonr
import multiprocessing as mp
import seaborn as sb
import matplotlib.pyplot as plt
import pandas as pd


def get_sep(seq1):

    three2one = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D',
                 'CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
                 'HIS':'H','ILE':'I','LEU':'L','LYS':'K',
                 'MET':'M','PHE':'F','PRO':'P','SER':'S',
                 'THR':'T','TRP':'W','TYR':'Y','VAL':'V',
                 'MSE':'M'}
    seq = ''
    prv = ''
    chain = ''
    for line in open(seq1):
        if line.startswith('ATOM'):
            if chain == '': chain = line[21]
            elif chain != line[21]: break 
            
            if line[22:27].strip() != prv:
                seq += three2one[line[17:20]]
                prv = line[22:27].strip()
        if line.startswith('TER'): break

    return seq

def prob_to_dist(distpred):
    #nres = distpred.shape[0]

    bin_step = 0.5
    smallest_bin = 2
    nbins = distpred.shape[-1] - 1
    bins = np.array([(smallest_bin+bin_step/2)+bin_step*i for i in range(nbins)])

    raw_predictions = distpred[:,:,1:nbins+1]
    binary_predictions = distpred[:,:,0]

    contact_prob = np.sum(raw_predictions, axis=-1)

    norm_predictions = np.divide(raw_predictions, contact_prob[:,:,None])
    values = np.multiply(norm_predictions, bins[None, None, :])
    mean_values = np.sum(values, axis=-1)

    return mean_values

def interface_contacts(code1, code2):

    p = PDBParser(QUIET=True)
    str1 = p.get_structure('', code1)
    str2 = p.get_structure('', code2)
    str1 = Selection.unfold_entities(str1[0], 'R')
    str2 = Selection.unfold_entities(str2[0], 'R')
    bins = [2+(0.5*n) for n in range(1,37)]

    cmap = np.zeros((len(str1), len(str2), 1))
    dcmap = np.zeros((len(str1), len(str2), 37))

    for pos1, res1 in enumerate(str1):
        CA1 = CB1 = ''
        for atom in res1:
            if atom.get_id() == 'CB': CB1 = atom
            elif atom.get_id() == 'CA': CA1 = atom
        if CB1 != '': C1 = CB1
        else: C1 = CA1

        for pos2, res2 in enumerate(str2):
            CA2 = CB2 = ''
            for atom in res2:
                if atom.get_id() == 'CB': CB2 = atom
                elif atom.get_id() == 'CA': CA2 = atom
            if CB2 != '': C2 = CB2
            else: C2 = CA2

            d = C1-C2

            cmap[pos1, pos2] = d

            for pos, thr in enumerate(bins):
                if d <= thr: 
                    dcmap[pos1,pos2,pos] = 1.0  
                    break
            if np.all(dcmap[pos1,pos2]==0.0): dcmap[pos1,pos2,0] = 1.0

    return cmap, dcmap

def mergesort_pred(pred):
    linear = [ [x,y,np.max(pred[x,y,1:])] for x in range(pred.shape[0]) for y in range(pred.shape[1]) ]
    
    elem = 1
    while len(linear) > elem:
        joblist = []
        for idx in range(0, len(linear)+elem*2, elem*2):
            ida = idx+elem
            idb = idx+elem*2
            if len(linear) >= idb:
                a = linear[idx:ida] 
                b = linear[ida:idb]
            elif len(linear) >= ida:
                a = linear[idx:ida]
                b = linear[ida:]
            elif len(linear) >= idx:
                a = linear[idx:]
                b = []
            else: continue
            joblist.append([a, b])

        pool = mp.Pool(processes=6)
        results = pool.map(merge, joblist)
        pool.close()
        pool.join()
    
        linear = [ el for res in results for el in res ]
        elem *= 2
   
    return linear

def merge(job):
    l = []
    l1 = job[0]
    l2 = job[1]
    p1 = p2 = 0
    while len(l) < len(l1)+len(l2):
        if p1 == len(l1): l.extend(l2[p2:])
        elif p2 == len(l2): l.extend(l1[p1:])
        elif l1[p1][2] >= l2[p2][2]: 
            l.append(l1[p1])
            p1 += 1
        else:
            l.append(l2[p2])
            p2 += 1
    return l

def prob_auc(pred, real, sortp, tol, bins):

    TPRs = []
    PPVs = []
    TP = TN = FP = FN = 0
    for x, y, p in sortp:
        realbin = np.argmax(real[x,y,1:])
        predbin = np.argmax(pred[x,y,1:])
        if realbin > bins and predbin > bins: continue
        tolrange = [el for el in range(realbin-tol, realbin+tol+1) if el > 0 and el <= 36]
        if predbin in tolrange: FN += 1
        else: TN += 1

    prev = ''
    for x, y, p in sortp:
        realbin = np.argmax(real[x,y,1:])
        predbin = np.argmax(pred[x,y,1:])
        if realbin > bins and predbin > bins: continue
        tolrange = [el for el in range(realbin-tol, realbin+tol+1) if el > 0 and el <= 36]

        if predbin in tolrange:
            FN -= 1
            TP += 1
        else:
            TN -= 1
            FP += 1

        if p != prev:
            if (TP != 0 or FN != 0) and (TP != 0 or FP != 0): 
                #print (TP, FP, FN, TN, TP/(TP+FN), TP/(TP+FP))
                TPRs += [TP/(TP+FN)]
                PPVs += [TP/(TP+FP)]
        prev = p

    if (TP != 0 or FN != 0) and (TP != 0 or FP != 0): 
        TPRs += [TP/(TP+FN)]
        PPVs += [TP/(TP+FP)]

    df = {'ppv':PPVs, 'tpr':TPRs}
    df = pd.DataFrame(df)
    sb.lineplot(x='tpr', y='ppv', data=df)
    plt.show()
    
    if len(TPRs) < 2 or len(PPVs) < 2: return 0.0
    else: return auc(TPRs, PPVs)


if __name__ == "__main__":

    p = argparse.ArgumentParser(description = '- Plot PPV stats for a Cmap or a list of Cmaps')
    p.add_argument('-c', required= True, help='trRosetta file with .npz format (folder path only if code list is specified)')
    p.add_argument('-s1', required= True, help='structure file 1 with .pdb format (folder path only if code list is specified)')
    p.add_argument('-s2', required= True, help='structure file 2 with .pdb format (folder path only if code list is specified)')
    p.add_argument('-t', required= False, default=2, help='bin tolerance to define True Positive (max value bin in pred +/- tol\
                                                           must match max value bin in real)')
    p.add_argument('-d', required= False, default=20, help='positions in real or pred cmap need to have maxval\
                                                            in a bin <= than this to be counted in AUC calculations')
    ns = p.parse_args()

    str1 = ns.s1
    str2 = ns.s2
    npz = ns.c
    tol = int(ns.t)
    dis = int(ns.d)

    #Extract the interface contacts from the real structures
    real_cmap_dist, real_cmap = interface_contacts(str1, str2)
    sep = len(get_sep(str1)) 

    #Extract the predicted inter-protein contacts
    with np.load(npz) as npz_file: pred_cmap = npz_file['dist'][:sep, sep:, :]

    assert real_cmap.shape == pred_cmap.shape,\
           'Cmaps doesn\'t match, real - {} vs pred - {}'.format(real_cmap.shape, pred_cmap.shape)

    sorted_list = mergesort_pred(pred_cmap)

    auc = prob_auc(pred_cmap, real_cmap, sorted_list, tol, dis)

    print ('#Maximum value allowed bins: first {} (within {} Armstrong)'.format(dis, 2+(dis*0.5)))
    print ('#Bin tolerance: {} (allowing error of {} Armstrong)'.format(tol, (tol*0.5)+0.5))
    print (auc)

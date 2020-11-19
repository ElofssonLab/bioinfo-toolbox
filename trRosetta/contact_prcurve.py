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

    cmap = np.zeros((len(str1), len(str2), 37))

    for pos1, res1 in enumerate(str1):
        CA1 = CB1 = ''
        for atom in res1:
            if atom.get_id() == 'CB': CB1 = atom
            elif atom.get_id() == 'CA': CA1 = atom
        if CB1 != '': C1 = CB1
        elif CA1!= '' :  C1 = CA1
        else: continue

        for pos2, res2 in enumerate(str2):
            CA2 = CB2 = ''
            for atom in res2:
                if atom.get_id() == 'CB': CB2 = atom
                elif atom.get_id() == 'CA': CA2 = atom
            if CB2 != '': C2 = CB2
            elif CA2 != '': C2 = CA2
            else: continue
            #print (pos2,res2,C1,":",C2)

            d = C1-C2
            for pos, thr in enumerate(bins):
                if d <= thr: 
                    cmap[pos1,pos2,pos] = 1.0  
                    break
            if np.all(cmap[pos1,pos2]==0.0): 
                cmap[pos1,pos2,0] = 1.0
                cmap[pos1,pos2,-1] = 1.0
    return cmap

def sort_pred(pred):
    sorted_list = []
    for px, x in enumerate(pred):
        for py, y in enumerate(x):
            prob = np.max(y[1:])
            if sorted_list == []:
                sorted_list.append([px,py,prob])
                continue
            for p, el in enumerate(sorted_list):
                if prob > el[2]: 
                    sorted_list = sorted_list[:p]+[[px,py,prob]]+sorted_list[p:]
                    break
            if [px,py,prob] not in sorted_list: sorted_list.append([px,py,prob])

    return sorted_list

def bin_auc(pred, real, sortp, tol, int_bins):

    TPRs = []
    PPVs = []
    TP = TN = FP = FN = 0
    for x, y, p in sortp:
        realbin = np.argmax(real[x,y,1:])
        predbin = np.argmax(pred[x,y,1:])
        if realbin > int_bins and predbin > int_bins: continue
        tolrange = [el for el in range(realbin-tol, realbin+tol+1) if el > 0 and el <= 35]
        if predbin in tolrange: FN += 1
        else: TN += 1
    
    prev = 1.1
    for x, y, p in sortp:
        realbin = np.argmax(real[x,y,1:])
        predbin = np.argmax(pred[x,y,1:])
        if realbin > int_bins and predbin > int_bins: continue
        tolrange = [el for el in range(realbin-tol, realbin+tol+1) if el > 0 and el <= 35]

        if predbin in tolrange:
            FN -= 1
            TP += 1
        else:
            FP += 1
            TN -= 1

        if p != prev:
            if (TP != 0 or FN != 0) and (TP != 0 or FP != 0): 
                TPRs += [TP/(TP+FN)]
                PPVs += [TP/(TP+FP)]
        prev = p

    if (TP != 0 or FN != 0) and (TP != 0 or FP != 0): 
        TPRs += [TP/(TP+FN)]
        PPVs += [TP/(TP+FP)]

    #print (PPVs, TPRs)
    df = {'ppv':PPVs, 'tpr':TPRs}
    df = pd.DataFrame(df)
    #print (df)
    sb.lineplot(x='tpr', y='ppv', data=df)
    plt.show()
    
    if len(TPRs) < 2 or len(PPVs) < 2: return 0.0
    else: return auc(TPRs, PPVs)


if __name__ == "__main__":

    p = argparse.ArgumentParser(description = '- Plot PPV stats for a Cmap or a list of Cmaps')
    p.add_argument('-c', required= True, help='trRosetta file with .npz format (folder path only if code list is specified)')
    p.add_argument('-s1', required= True, help='structure file 1 with .pdb format (folder path only if code list is specified)')
    p.add_argument('-s2', required= True, help='structure file 2 with .pdb format (folder path only if code list is specified)')
    ns = p.parse_args()

    aucs12 = aucs8 = {}
    tolrange = range(2, 8, 2)
    for t in tolrange: aucs12[t] = 0
    for t in tolrange: aucs8[t] = 0

    header = '#\t\t\t'
    for t in tolrange: header += 'Int.12A-Tol.{}\t'.format((t*0.5)+0.5)
    for t in tolrange: header += 'Int.8A-Tol.{}\t'.format((t*0.5)+0.5)
    print (header)

    str1 = ns.s1
    str2 = ns.s2
    npz = ns.c

    #Extract the interface contacts from the real structures
    real_cmap = interface_contacts(str1, str2)
    sep = len(get_sep(str1)) 

    #Extract the predicted inter-protein contacts
    with np.load(npz) as npz_file: pred_cmap = npz_file['dist'][:sep, sep:, :]
    
    #assert real_cmap.shape == pred_cmap.shape,\
    #       'Cmaps doesn\'t match, real - {} vs pred - {}'.format(real_cmap.shape, pred_cmap.shape)
    if real_cmap.shape != pred_cmap.shape:
           print ('Cmaps doesn\'t match, real - {} vs pred - {}'.format(real_cmap.shape, pred_cmap.shape))
    else:
        sorted_list = sort_pred(pred_cmap)

        codeauc8 = []
        codeauc12 = []
        for t in tolrange: 
            aucval12 = bin_auc(pred_cmap, real_cmap, sorted_list, t, 35)
            aucs12[t] += aucval12
            codeauc12.append(aucval12)
        for t in tolrange: 
            aucval8 = bin_auc(pred_cmap, real_cmap, sorted_list, t, 12)
            aucs8[t] += aucval8
            codeauc8.append(aucval8)

        values = ns.c.split('/')[-1].rstrip('.npz')+'\t'
        for aucval in codeauc12: values += '{}\t\t'.format(round(aucval, 3))
        for aucval in codeauc8: values += '{}\t\t'.format(round(aucval, 3))
        print (values)


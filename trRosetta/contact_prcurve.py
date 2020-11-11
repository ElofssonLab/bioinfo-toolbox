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

    no_contact_prob = distpred[:,:,0]
    contact_prob = np.expand_dims(np.sum(raw_predictions, axis=-1), axis=-1)
    contact_split = np.subtract(contact_prob, no_contact_prob)

    norm_predictions = np.divide(raw_predictions, contact_prob)
    values = np.multiply(norm_predictions, bins[None, None, :])
    mean_values = np.expand_dims(np.sum(values, axis=-1), axis=-1)
    mean_values[np.where(contact_split<=0)] == 21

    data = np.concatenate((contact_prob, mean_values), axis=-1)

    return data

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
        else: C1 = CA1

        for pos2, res2 in enumerate(str2):
            CA2 = CB2 = ''
            for atom in res2:
                if atom.get_id() == 'CB': CB2 = atom
                elif atom.get_id() == 'CA': CA2 = atom
            if CB2 != '': C2 = CB2
            else: C2 = CA2

            d = C1-C2
            for pos, thr in enumerate(bins):
                if d <= thr: 
                    cmap[pos1][pos2][pos] = 1.0  
                    break
            if np.all(cmap[pos1][pos2]==0.0): 
                cmap[pos1][pos2][0] = 1.0

    return cmap

def sort_pred(pred):
    sorted_list = []
    for px, x in enumerate(pred):
        for py, y in enumerate(x):
            prob = np.max(y)
            if sorted_list == []:
                sorted_list.append([px,py,prob])
                continue
            for p, el in enumerate(sorted_list):
                if prob > el[2]: sorted_list = sorted_list[:p]+[[px,py,prob]]+sorted_list[p:]
                break
            if [px,py,prob] not in sorted_list: sorted_list.append([px,py,prob])

    return sorted_list

def bin_auc(real, pred, sortp, tol):

    TPRs = []
    PPVs = []
    TP = TN = FP = FN = 0
    for x, y, p in sortp:
        realbin = np.argmax(real[x][y])
        if max(real[x][y][realbin], real[x][y][0]) != real[x][y][0]: 
            tolrange = [el for el in range(realbin-tol, realbin+tol+1) if el > 0]
        else: tolrange = [0]+[el for el in range(36, 36-(tol+1), -1)]

        if np.argmax(pred[x][y]) in tolrange: FN += 1
        else: TN += 1
   
    prev = 1.1
    for x, y, p in sortp:
        realbin = np.argmax(real[x][y])
        if max(real[x][y][realbin], real[x][y][0]) != real[x][y][0]:
            tolrange = [el for el in range(realbin-tol, realbin+tol+1) if el > 0]
        else: tolrange = [0]+[el for el in range(36, 36-(tol+1), -1)]

        if np.argmax(pred[x][y]) in tolrange:
            FN -= 1
            TP += 1
        else:
            FP += 1
            TN -= 1

        if p != prev:
            if TP != 0 or FN != 0: TPRs += [TP/(TP+FN)]
            else: TPRs += [0.0]
            if TP != 0 or FP != 0: PPVs += [TP/(TP+FP)]
            else: PPVs += [0.0]

        prev = p

    if TP != 0 or FN != 0: TPRs += [TP/(TP+FN)]
    else: TPRs += [0.0]
    if TP != 0 or FP != 0: PPVs += [TP/(TP+FP)]
    else: PPVs += [0.0]

    return auc(TPRs, PPVs)


if __name__ == "__main__":

    p = argparse.ArgumentParser(description = '- Plot PPV stats for a Cmap or a list of Cmaps')
    p.add_argument('-l', required= False, default=None, help='codes list (codes must be identical for all files)')
    p.add_argument('-c', required= True, help='trRosetta file with .npz format (folder path only if code list is specified)')
    p.add_argument('-s1', required= True, help='structure file 1 with .pdb format (folder path only if code list is specified)')
    p.add_argument('-s2', required= True, help='structure file 2 with .pdb format (folder path only if code list is specified)')
    ns = p.parse_args()

    if ns.l is None: runlist = ['']
    else: runlist = [code.rstrip() for code in open(ns.l)]
    
    aucs = {}
    for t in range(0, 6, 2): aucs[t] = 0

    for code in runlist:
        if code == '':
            str1 = ns.s1
            str2 = ns.s2
            npz = ns.c
        else: 
            str1 = ns.s1+code+'.pdb'
            str2 = ns.s2+code+'.pdb'
            npz = ns.c+code+'.npz'

        #Extract the interface contacts from the real structures
        real_cmap = interface_contacts(str1, str2)
        sep = len(get_sep(str1)) 

        #Extract the predicted inter-protein contacts
        with np.load(npz) as npz_file: pred_cmap = npz_file['dist'][:sep, sep:, :]
        #pdmap = npz_to_casp(pred_cmap)
    
        if real_cmap.shape != pred_cmap.shape:
            print ('Cmaps doesn\'t match, real - {} vs pred - {}'.format(real_cmap.shape, pred_cmap.shape))
            continue

        sorted_list = sort_pred(pred_cmap)
        for t in range(0, 6, 2): aucs[t] += bin_auc(pred_cmap, real_cmap, sorted_list, t)
        
    for t in range(0, 6, 2): print ('AUC (tolerance of {}A): {}'.format(t*0.5, aucs[t]/len(runlist)))

    

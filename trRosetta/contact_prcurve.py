#!/usr/bin/env python3:q
import os
import sys
import math
import numpy as np
from Bio.PDB import *
import matplotlib.pyplot as plt
from sklearn.metrics import auc

if len(sys.argv) != 4:
    print("Usage: {} code_list npz_path npz_suffix str_path str1_suffix str2_suffix".format(os.path.basename(__file__)))
    sys.exit()

three2one = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D',
             'CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
             'HIS':'H','ILE':'I','LEU':'L','LYS':'K',
             'MET':'M','PHE':'F','PRO':'P','SER':'S',
             'THR':'T','TRP':'W','TYR':'Y','VAL':'V',
             'MSE':'M'}

def pdb_scan(pdb):
    """
    Extract fasta sequence from a pdb formatted file
    """
    prv = ''
    seq = ''
    for line in pdb:
        if line.startswith('ATOM'):
            if line[22:27].strip() != prv:
                seq += three2one[line[17:20]]
                prv = line[22:27].strip()

    return seq


def npz_to_casp(raw_data, first_seq, Pthr):
    """
    Extract inter-protein distances from trRosetta npz-file.
    Returns a dictionary where keys are tuples containing the 
    index of the two positions in contact, while the corresponding 
    value is a list containing calculated CB-CB distance and sum 
    of all bins probabilities for that residue prediction;
    """
    pred_contacts = {}

    # Length of the protein
    nres = raw_data.shape[0]
    print ("#Map size: "+ str(nres))

    # Number of bins, minus one for the first contact binary bin
    nbins = raw_data.shape[-1] - 1

    # Generate up the midpoint of each distance bin to be able to
    # calculate correct distances. Measure in Ångström
    bins = np.array([(smallest_bin+bin_step/2)+bin_step*i for i in range(nbins)])

    # Pull out the predictions, first array element is 
    # probability of non-contact, only pull out as many bins as we want
    # [L, L, nbins] -> [L, L, wanted_bins]
    raw_predictions = raw_data[:,:,1:nbins+1]
    binary_predictions = raw_data[:,:,0]

    # What is the total probability of contact?
    contact_probabilities = np.sum(raw_predictions, axis=-1)

    # Normalize all probabilities by dividing by the sum of the raw_predictions
    # across the last dimension, [L, L, nbins] -> [L, L, nbins]
    norm_predictions = np.divide(raw_predictions, np.sum(raw_predictions, axis=-1)[:,:,None])

    # Multiply normalized predictions with the bin (sizes) to get the correct value distribution
    values = np.multiply(norm_predictions, bins[None, None, :])

    # Make some assertions to confirm array operations har occured correctly
    # Has the normalization worked?
    assert np.array_equal(np.divide(raw_predictions[0, 5, :],np.sum(raw_predictions[0,5,:])), norm_predictions[0,5])
    # Has the multiplication with the bins worked?
    assert np.array_equal(np.multiply(np.divide(raw_predictions[0, 5, :],np.sum(raw_predictions[0, 5, :])), bins),
                          values[0,5, :])

    # Calculate mean value and errors
    mean_values = np.sum(values, axis=-1)
    variances = np.sum(
                    np.multiply(
                        np.square(
                            np.subtract(bins[None,None,:], mean_values[:,:,None])),
                            norm_predictions), axis=-1)
    sd = np.sqrt(variances)

    for i in range(len(first_seq)-1):
        for j in range(len(first_seq)+20, nres-1):
            if contact_probabilities[i][j] < Pthr: continue
            pred_contacts[(i+1, j+1)] = [mean_values[i,j], contact_probabilities[i,j]]

    return pred_contacts


def interface_extract(str1, str2, idx2_shift, Ithr):
    """
    Extract a list of CB-CB contacts between the two structures
    numbered like in the npz array. Ithr is the minimal CB-CB
    distance to define a contact. Returns a dictionary where keys 
    are tuples containing the index of the two positions in contact, 
    while the corresponding value is the CB-CB distance;
    """
    pos1 = 1
    true_contacts = {}
    for residue1 in str1:
        alpha1 = ''
        beta1 = ''
        atom_list1 = []
        for atom1 in residue1:
            if atom1.get_id() == 'CB': beta1 = atom1
            if atom1.get_id() == 'CA': alpha1 = atom1

        if beta1 != '': atom1=beta1
        elif alpha1 != '': atom1=alpha1
        else: continue

        pos2 = idx2_shift
        for residue2 in str2:
            alpha2 = ''
            beta2 = ''
            for atom2 in residue2:
                if atom2.get_id() == 'CB': beta2 = atom2 
                if atom2.get_id() == 'CA': alpha2 = atom2 
            
            if beta2 != '': atom2=beta2
            elif alpha2 != '': atom2=alpha2
            else: continue

            betabeta_dist = atom1-atom2
            if betabeta_dist < Ithr: true_contacts[(pos1, pos2)] = betabeta_dist
            pos2 += 1
        pos1 += 1

    return true_contacts

def find_distance(atom1, atom2):
    """
    calculate distance between two atoms
    """ 
    seta=atom1.coord
    setb=atom2.coord
    dbb = math.sqrt((seta[0]-setb[0])**2+(seta[1]-setb[1])**2+(seta[2]-setb[2])**2)

    return dbb


if __name__ == "__main__":

    code = sys.argv[1].split('/')[-1][:4]

    Ithr = 12
    Pthr = 0.5
    bin_step = 0.5
    smallest_bin = 2

    npz_path = sys.argv[1]
    str1_path = sys.argv[2]
    str2_path = sys.argv[3]
    thresholds = np.arange(9, 19, 0.1)

    #Extract the interface contacts from the real structures
    try: 
        with open(str1_path) as pdb: first_seq = pdb_scan(pdb)
    except: 
        print ("#" + str1_path + 'not found!')
        sys.exit()
    try: 
        with open(str2_path) as pdb: second_seq = pdb_scan(pdb)
    except:
        print ("#" + str2_path + 'not found!')
        sys.exit()

    p = PDBParser(QUIET=True)
    str1 = p.get_structure('', str1_path)
    res_list1 = Selection.unfold_entities(str1, 'R')
    str2 = p.get_structure('', str2_path)
    res_list2 = Selection.unfold_entities(str2, 'R')

    real = interface_extract(res_list1, res_list2, len(first_seq)+20, Ithr)
    print ("#Structure pair size: "+ str(len(first_seq)+len(second_seq)+20))

    #Extract the predicted inter-protein contacts
    try:
        with np.load(npz_path) as npz_file: raw_data = npz_file['dist']
    except:
        print("#Error reading npz file: {}".format(input_file), file=sys.stderr)
        sys.exit()

    pred = npz_to_casp(raw_data, first_seq, Pthr)

    #Compare predicted vs real contacts over different thresholds
    rocdata = {}
    for t in thresholds: rocdata[t] = {'TP':0, 'FP':0, 'FN':0}
    positive_number = 0

    for p in pred: 
        if p in list(real.keys()):
            for t in thresholds:
                if pred[p][0] <= t: rocdata[t]['TP'] += 1
                else: rocdata[t]['FN'] += 1
        else:
            for t in thresholds:
                if pred[p][0] <= t: rocdata[t]['FP'] += 1

    #Compute ppv and tpr for each threshold skipping useless points
    PPV = []
    TPR = []
    for t in thresholds:
        ppv = 0.0
        tpr = 0.0
        if rocdata[t]['TP'] != 0 or rocdata[t]['FP'] != 0: 
            ppv = rocdata[t]['TP']/(rocdata[t]['TP']+rocdata[t]['FP'])
        if rocdata[t]['TP'] != 0 or rocdata[t]['FN'] != 0:
            tpr = rocdata[t]['TP']/(rocdata[t]['TP']+rocdata[t]['FN'])

        if ppv == 0.0 and tpr == 0.0: continue 
        PPV.append(ppv)
        TPR.append(tpr)

    #Compute AUC, plot Precision-Recall curve and save in a figure
    values = 'n'
    if len(TPR) != 0: 
        for val in TPR:
            if val>0: values = 'y'
    if values=='n': 
        print (code.rstrip()+" - Pthr:"+str(Pthr)+" - Ithr:"+str(Ithr)+" - AUC:", 0.0)
        sys.exit()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(TPR, PPV)
    plt.ylabel("Precision")
    plt.ylim(0,1)
    plt.xlabel("Recall")
    plt.xlim(0,1)
    AUC = round(auc(TPR, PPV), 2)
    print (code.rstrip()+" - Pthr:"+str(Pthr)+" - Ithr:"+str(Ithr)+" - AUC:", AUC)

    plt.title(code.rstrip()+" - AUC:"+str(AUC)+" - Pthr:"+str(Pthr)+" - Ithr:"+str(Ithr))
    fig.savefig('pr-curve_'+code.rstrip()+'.png')



#!/usr/bin/env python3
import os
import sys
import math
import numpy as np
from Bio.PDB import *


if len(sys.argv) != 4:
    #print("Usage: {} code_list npz_path npz_suffix str_path str1_suffix str2_suffix".format(os.path.basename(__file__)))

    print("Usage: {}  npz_path str_path2 str_path2".format(os.path.basename(__file__)))

    sys.exit()

three2one = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D',
             'CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
             'HIS':'H','ILE':'I','LEU':'L','LYS':'K',
             'MET':'M','PHE':'F','PRO':'P','SER':'S',
             'THR':'T','TRP':'W','TYR':'Y','VAL':'V',
             'MSE':'M'}

#npz_path = sys.argv[2]
#str1_path = sys.argv[4]
#str2_path = sys.argv[4]
smallest_bin = 2
bin_step = 0.5

def pdb_scan(pdb):

    prv = ''
    seq = ''
    for line in pdb:
        if line.startswith('ATOM'):
            if line[22:27].strip() != prv:
                seq += three2one[line[17:20]]
                prv = line[22:27].strip()

    return seq


def npz_to_casp(raw_data, first_seq, Pthr=0.5):
    """
    Extract inter-protein distances from trRosetta npz-file .
    """
    pred_contacts = {}

    # Length of the protein
    nres = raw_data.shape[0]
    print (nres)

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
    max_probability = np.amax(raw_predictions, axis=-1)
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
            if contact_probabilities[i][j] > Pthr:
                pred_contacts[(i+1, j+1)] = [mean_values[i,j], contact_probabilities[i,j], max_probability[i,j]]

    return pred_contacts


def interface_extract(str1, str2, idx2_shift):

    pos1 = 1
    true_contacts = {}
    for residue1 in str1:
        atom_list1 = []
        for atom1 in residue1:
            atom_list1.append(atom1)

        ns = NeighborSearch(atom_list1)

        pos2 = idx2_shift
        for residue2 in str2:
            for atom2 in residue2:
                neighbours = ns.search(atom2.coord, 8, level='A')
                if len(list(neighbours)) != 0:
                    min_dist = find_distance(residue1, residue2)
                    true_contacts[(pos1, pos2)] = min_dist
                    break
            pos2 += 1
        pos1 += 1

    return true_contacts

def find_distance(res1, res2):

    min_dist = 9999
    for atom1 in res1:
        for atom2 in res2:
            seta=atom1.coord
            setb=atom2.coord
            dab = math.sqrt((seta[0]-setb[0])**2+(seta[1]-setb[1])**2+(seta[2]-setb[2])**2)
            min_dist=min(dab,min_dist)

    return min_dist


if __name__ == "__main__":

    data = {}

    #for code in open(sys.argv[1]):
        
    #npz_path = sys.argv[2]+line.rstrip()+sys.argv[3]
    #str1_path = sys.argv[4]+line.rstrip()+sys.argv[5]
    #str2_path = sys.argv[4]++line.rstrip()+sys.argv[6]
    npz_path = sys.argv[1]
    str1_path = sys.argv[2]
    str2_path = sys.argv[3]
    smallest_bin = 2
    bin_step = 0.5

    try: 
        with open(str1_path) as pdb: first_seq = pdb_scan(pdb)
    except: 
        print (str1_path + 'not found!')
        sys.exit()
    try: 
        with open(str2_path) as pdb: second_seq = pdb_scan(pdb)
    except:
        print (str2_path + 'not found!')
        sys.exit()

    try:
        with np.load(npz_path) as npz_file: raw_data = npz_file['dist']
    except:
        print("Error reading npz file: {}".format(input_file), file=sys.stderr)
        sys.exit()

    pred = npz_to_casp(raw_data, first_seq)

    p = PDBParser()
    str1 = p.get_structure('', str1_path)
    str2 = p.get_structure('', str2_path) 
    print ("HEJ",str1_path,str2_path)
    real = interface_extract(str1[0]['A'], str2[0]['A'], len(first_seq)+20)

    #data[code] = [real, pred]
    print (real,pred)
    #P = set(real.keys())

    #Pr = set(pred.keys())
    #TP = list(Pr.intersection(P))
    
    #for key in s1: 
    #    if key in real: print (key, real[key], pred[key])
    #    else: print (key, "WRONG", pred[key])


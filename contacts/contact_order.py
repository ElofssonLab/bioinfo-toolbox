import sys
import argparse
from math import *

import numpy as np

from os.path import expanduser
home = expanduser("~")
sys.path.append(home + '/bioinfo-toolbox')

from parsing import parse_contacts
from parsing import parse_pdb


def get_cb_contacts(cb_lst):

    seqlen = len(cb_lst)
    dist_mat = np.zeros((seqlen, seqlen), np.float)
    dist_mat.fill(float('inf'))
    
    for i, cb1 in enumerate(cb_lst):
        for j, cb2 in enumerate(cb_lst):
            diff_vec = cb1 - cb2
            dist_mat[i,j] = np.sqrt(np.sum(diff_vec * diff_vec))
    return dist_mat


def get_co_pdb(pdb_filename, chain, cb_cutoff=8):

    cb_lst = parse_pdb.get_cb_coordinates(open(pdb_filename, 'r'), chain)
    dist_mat = get_cb_contacts(cb_lst)
    ref_contact_map = dist_mat < cb_cutoff
    atom_seq = parse_pdb.get_atom_seq(open(pdb_filename, 'r'), chain)
    np.fill_diagonal(ref_contact_map, 0)

    co = 0
    L = len(atom_seq)
    N = sum((ref_contact_map != 0).sum(0)/2)

    for (i,j), is_contact in np.ndenumerate(ref_contact_map):
        if i < j and is_contact:
            S_ij = j - i
            co += S_ij
    co = float(co)/float(N*L)
    print co

    return co

    
if __name__ == "__main__":

    p = argparse.ArgumentParser(description='Calculate contact order for given PDB structure or contact map.')
    p.add_argument('pdb')#, required=True)
    p.add_argument('chain')#, required=True)

    args = vars(p.parse_args(sys.argv[1:]))
    
    """
    c_filename = args['contact_file']

    # guessing separator of constraint file
    line = open(c_filename,'r').readline()
    if len(line.split(',')) != 1:
        sep = ','
    elif len(line.split(' ')) != 1:
        sep = ' '
    else:
        sep = '\t'
    """
    co = get_co_pdb(args['pdb'], args['chain'])



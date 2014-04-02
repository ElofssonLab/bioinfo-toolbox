import sys
import argparse
from math import *

# on UPPMAX only
sys.path.append('/sw/apps/bioinfo/biopython/1.59/tintin/lib/python')

import Bio.PDB
from Bio import pairwise2

import numpy as np

from os.path import expanduser
home = expanduser("~")
sys.path.append(home + '/bioinfo-toolbox')

from parsing import parse_contacts
from parsing import parse_fasta
from parsing import parse_pdb


def get_cb_contacts(gapped_cb_lst):

    seqlen = len(gapped_cb_lst)
    dist_mat = np.zeros((seqlen, seqlen), np.float)
    dist_mat.fill(float('inf'))
    
    for i, cb1 in enumerate(gapped_cb_lst):
        if cb1 == '-':
            continue
        for j, cb2 in enumerate(gapped_cb_lst):
            if cb2 == '-':
                continue
            diff_vec = cb1 - cb2
            dist_mat[i,j] = np.sqrt(np.sum(diff_vec * diff_vec))
    return dist_mat


def get_ppv_helper(contacts_x, contacts_y, ref_contact_map, ref_len, factor, atom_seq_ali=[]):

    num_c = min(len(contacts_x), ceil(ref_len * factor))
    TP = 0.0
    FP = 0.0
    PPV = 0.0
    for i in range(num_c):
        c_x = contacts_x[i]
        c_y = contacts_y[i]
        if atom_seq_ali:
            if atom_seq_ali[c_x] == '-':
                continue
            if atom_seq_ali[c_y] == '-':
                continue
        if ref_contact_map[c_x, c_y] > 0:
            TP += 1.0
        else:
            FP += 1.0

    if TP > 0 and FP > 0:
        TP = TP / (ref_len*factor)
        FP = FP / (ref_len*factor)
        PPV = TP / (TP + FP)

    return (PPV, TP, FP)


def get_ppv(fasta_filename, c_filename, pdb_filename, factor=1.0, chain='', sep=' ', outfilename='', noalign=False):  
    
    acc = fasta_filename.split('.')[-2][-5:-1]

    ### get sequence
    seq = parse_fasta.read_fasta(open(fasta_filename, 'r')).values()[0][0]
    ref_len = len(seq)

    ### get top "factor" * "ref_len" predicted contacts
    contacts = parse_contacts.parse(open(c_filename, 'r'), sep)

    contacts_x = []
    contacts_y = []
    scores = []
    contact_dict = {}

    count = 0
    for i in range(len(contacts)):
        score = contacts[i][0]
        c_x = contacts[i][1] - 1
        c_y = contacts[i][2] - 1

        pos_diff = abs(c_x - c_y)
        too_close = pos_diff < 20

        if not too_close:
            contacts_x.append(c_x)
            contacts_y.append(c_y)
            scores.append(score)
            count += 1
           
        if count >= ref_len * factor:
            break
 
    cb_lst = parse_pdb.get_cb_coordinates(open(pdb_filename, 'r'), chain)

    if noalign:
        dist_mat = get_cb_contacts(cb_lst)
        cb_cutoff = 8
        ref_contact_map = dist_mat < cb_cutoff
        PPV, TP, FP = get_ppv_helper(contacts_x, contacts_y, ref_contact_map, ref_len, factor)
    else:
        atom_seq = parse_pdb.get_atom_seq(open(pdb_filename, 'r'), chain)
                
        align = pairwise2.align.globalms(atom_seq, seq, 2, -1, -0.5, -0.1)

        atom_seq_ali = align[-1][0]
        seq_ali = align[-1][1]

        j = 0
        gapped_cb_lst = []

        for i in xrange(len(atom_seq_ali)):
            if atom_seq_ali[i] == '-':
                gapped_cb_lst.append('-')
            elif seq_ali[i] == '-':
                j += 1
                continue
            else:
                gapped_cb_lst.append(cb_lst[j])
                j += 1

        dist_mat = get_cb_contacts(gapped_cb_lst)
        cb_cutoff = 8
        ref_contact_map = dist_mat < cb_cutoff
   
        PPV, TP, FP = get_ppv_helper(contacts_x, contacts_y, ref_contact_map, ref_len, factor, atom_seq_ali=atom_seq_ali)

    print '%s %s %s %s' % (pdb_filename, PPV, TP, FP)
    return (pdb_filename, PPV, TP, FP)
  

    
if __name__ == "__main__":

    p = argparse.ArgumentParser(description='Plot protein residue contact maps.')
    p.add_argument('fasta_file')#, required=True)
    p.add_argument('contact_file')#, required=True)
    p.add_argument('pdb')
    p.add_argument('-o', '--outfile', default='')
    p.add_argument('-f', '--factor', default=1.0, type=float)
    p.add_argument('--chain', default='')
    p.add_argument('--noalign', action='store_true')

    args = vars(p.parse_args(sys.argv[1:]))

    fasta_filename = args['fasta_file']
    c_filename = args['contact_file']

    # guessing separator of constraint file
    line = open(c_filename,'r').readline()
    if len(line.split(',')) != 1:
        sep = ','
    elif len(line.split(' ')) != 1:
        sep = ' '
    else:
        sep = '\t'
    
    get_ppv(args['fasta_file'], args['contact_file'], args['pdb'], args['factor'], chain=args['chain'], sep=sep, outfilename=args['outfile'], noalign=args['noalign'])


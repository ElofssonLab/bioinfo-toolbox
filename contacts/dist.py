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


def get_min_dist(res1, res2):
    
    min_dist = float('inf')

    for atm1 in res1:
        for atm2 in res2:
            diff_vec = atm1 - atm2
            dist = np.sqrt(np.sum(diff_vec * diff_vec))
            if dist < min_dist:
                min_dist = dist

    return min_dist


def get_dist_mat_heavy(gapped_res_lst):

    seqlen = len(gapped_res_lst)
    dist_mat = np.zeros((seqlen, seqlen), np.float)
    dist_mat.fill(float('inf'))
    
    for i, res1 in enumerate(gapped_res_lst):
        if res1 == '-':
            continue
        for j, res2 in enumerate(gapped_res_lst):
            if res2 == '-':
                continue
            dist_mat[i,j] = get_min_dist(res1[1], res2[1])
    return dist_mat



def get_dist_mat(gapped_lst):

    seqlen = len(gapped_lst)
    dist_mat = np.zeros((seqlen, seqlen), np.float)
    dist_mat.fill(float('inf'))
    
    for i, cb1 in enumerate(gapped_lst):
        if cb1 == '-':
            continue
        for j, cb2 in enumerate(gapped_lst):
            if cb2 == '-':
                continue
            diff_vec = cb1 - cb2
            dist_mat[i,j] = np.sqrt(np.sum(diff_vec * diff_vec))
    return dist_mat


def get_dist_helper(contacts_x, contacts_y, dist_mat, atom_seq_ali=[]):

    contacts_dist = []
    num_c = len(contacts_x)
    for i in range(num_c):
        c_x = contacts_x[i]
        c_y = contacts_y[i]
        dist = dist_mat[c_x, c_y]
        if atom_seq_ali:
            if atom_seq_ali[c_x] == '-' or atom_seq_ali[c_y] == '-':
                contacts_dist.append(-1)
                continue
        contacts_dist.append(dist)

    return contacts_dist


def get_dist(fasta_filename, c_filename, pdb_filename, chain='', 
sep='', outfilename='', noalign=False, dist_type='CB'):  
    
    acc = fasta_filename.split('.')[-2][-5:-1]

    ### get sequence
    seq = parse_fasta.read_fasta(open(fasta_filename, 'r')).values()[0][0]
    ref_len = len(seq)

    ### get top "factor" * "ref_len" predicted contacts
    contacts = parse_contacts.parse(open(c_filename, 'r'), sep,
            min_dist=5)

    contacts_x = []
    contacts_y = []
    scores = []

    count = 0
    for i in range(len(contacts)):
        score = contacts[i][0]
        c_x = contacts[i][1] - 1
        c_y = contacts[i][2] - 1

        contacts_x.append(c_x)
        contacts_y.append(c_y)
        scores.append(score)
        count += 1
           
    res_lst = parse_pdb.get_coordinates(open(pdb_filename, 'r'), chain)
    cb_lst = parse_pdb.get_cb_coordinates(open(pdb_filename, 'r'), chain)
    ca_lst = parse_pdb.get_ca_coordinates(open(pdb_filename, 'r'), chain)

    if noalign:
        if dist_type == 'CB':
            dist_mat = get_dist_mat(cb_lst)
        elif dist_type == 'CA':
            dist_mat = get_dist_mat(ca_lst)
        else:
            dist_mat = get_dist_mat_heavy(res_lst)

        contacts_dist = get_dist_helper(contacts_x, contacts_y, dist_mat)

    else:
        atom_seq = parse_pdb.get_atom_seq(open(pdb_filename, 'r'), chain)
                
        align = pairwise2.align.globalms(atom_seq, seq, 2, -1, -0.5, -0.1)
        atom_seq_ali = align[-1][0]
        seq_ali = align[-1][1]
        j = 0
        gapped_res_lst = []
        gapped_cb_lst = []
        gapped_ca_lst = []

        for i in xrange(len(atom_seq_ali)):
            if atom_seq_ali[i] == '-':
                gapped_res_lst.append('-')
                gapped_cb_lst.append('-')
                gapped_ca_lst.append('-')
            elif seq_ali[i] == '-':
                j += 1
                continue
            else:
                gapped_res_lst.append(res_lst[j])
                gapped_cb_lst.append(cb_lst[j])
                gapped_ca_lst.append(ca_lst[j])
                j += 1
    
        assert(len(gapped_ca_lst) == len(gapped_cb_lst) == len(gapped_res_lst))

        if dist_type == 'CB':
            dist_mat = get_dist_mat(gapped_cb_lst)
        elif dist_type == 'CA':
            dist_mat = get_dist_mat(gapped_ca_lst)
        else:
            dist_mat = get_dist_mat_heavy(gapped_res_lst)

        contacts_dist = get_dist_helper(contacts_x, contacts_y, dist_mat, atom_seq_ali=atom_seq_ali)

    assert(len(contacts_dist) == len(contacts_x) == len(contacts_y) == len(scores))

    num_c = len(contacts_dist)

    if outfilename:
        with open(outfilename, 'w') as outfile:
            for i in xrange(num_c):
                outfile.write('%s %s %f %f\n' % (contacts_x[i],
                    contacts_y[i], scores[i], contacts_dist[i]))

    return (contacts_x, contacts_y, scores, contacts_dist)
  
    

if __name__ == "__main__":

    p = argparse.ArgumentParser(description='Plot protein residue contact maps.')
    p.add_argument('fasta_file')#, required=True)
    p.add_argument('contact_file')#, required=True)
    p.add_argument('pdb')
    p.add_argument('-o', '--outfile', default='')
    p.add_argument('-d', '--disttype', default='CB')
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
    
    (contacts_x, contacts_y, scores, contacts_dist) = get_dist(
            args['fasta_file'], args['contact_file'], args['pdb'],
            chain=args['chain'], sep=sep, outfilename=args['outfile'],
            dist_type=args['disttype'], noalign=args['noalign'])

    


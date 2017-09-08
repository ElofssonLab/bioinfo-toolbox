import sys
import argparse
from math import *

# on UPPMAX only
#sys.path.append('/sw/apps/bioinfo/biopython/1.59/tintin/lib/python')

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


def print_distances(contacts_x, contacts_y, scores, dist_mat, atom_seq_ali=[], outfile=""):

    num_c = len(contacts_x)
    outstr = ""
    for i in range(num_c):
        c_x = contacts_x[i]
        c_y = contacts_y[i]
        if atom_seq_ali:
            if atom_seq_ali[c_x] == '-':
                continue
            if atom_seq_ali[c_y] == '-':
                continue
        if outfile:
            outstr += "%s %s %s %s\n" % (c_x, c_y, scores[i], dist_mat[c_x, c_y])
        else:
            print c_x+1, c_y+1, scores[i], dist_mat[c_x, c_y]
    if outfile:
        with open(outfile, 'w') as outf:
            outf.write(outstr)


def get_ppv_helper(contacts_x, contacts_y, ref_contact_map, ref_len, factor, atom_seq_ali=[]):

    num_c = len(contacts_x)
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
            TP += 1.0 / num_c
        else:
            FP += 1.0 / num_c

    if TP > 0:
        PPV = TP / (TP + FP)

    return (PPV, TP, FP)


def get_ppv(fasta_filename, c_filename, pdb_filename, factor=1.0,
        min_score=-1.0, chain='', sep=' ', outfilename='', name='', noalign=False, min_dist=5, print_dist=False):  
    
    acc = fasta_filename.split('.')[-2][-5:-1]

    ### get sequence
    seq = parse_fasta.read_fasta(open(fasta_filename, 'r')).values()[0][0]
    ref_len = len(seq)

    ### get top ranked predicted contacts
    contacts = parse_contacts.parse(open(c_filename, 'r'), sep, min_dist=min_dist)

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
        too_close = pos_diff < min_dist

        if not too_close:
            contacts_x.append(c_x)
            contacts_y.append(c_y)
            scores.append(score)
            count += 1
           
        if min_score == -1.0 and count >= ref_len * factor:
            break
        if score < min_score:
            break
    
    assert(len(contacts_x) == len(contacts_y) == len(scores))

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
        if print_dist:
            print_distances(contacts_x, contacts_y, scores, dist_mat, atom_seq_ali=atom_seq_ali, outfile=outfilename)
        cb_cutoff = 8
        ref_contact_map = dist_mat < cb_cutoff
   
        PPV, TP, FP = get_ppv_helper(contacts_x, contacts_y, ref_contact_map, ref_len, factor, atom_seq_ali=atom_seq_ali)
    if name:
        print '%s %s %s %s' % (name, PPV, TP, FP)
    else:
        print '%s %s %s %s %s' % (fasta_filename, c_filename, PPV, TP, FP)
    return (pdb_filename, PPV, TP, FP)
  
    
def get_ppv_hbond(fasta_filename, c_filename, hbond_filename,
        factor=1.0, min_score=-1.0, sep=' ', outfilename=''):  

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
        too_close = pos_diff < 5

        if not too_close:
            contacts_x.append(c_x)
            contacts_y.append(c_y)
            scores.append(score)
            count += 1
           
        if min_score == -1.0 and count >= ref_len * factor:
            break
        if score < min_score:
            break
   
    ref_contact_map = np.zeros((ref_len, ref_len))

    hbonds_raw = open(hbond_filename).readlines()
    hbonds = [line.strip().split(' ') for line in hbonds_raw] #map(split(' '), map(strip, hbonds_raw))

    for h in hbonds:
        i = int(h[0]) - 1
        j = int(h[1]) - 1
        val = float(h[2])
        ref_contact_map[i,j] = -val
        ref_contact_map[j,i] = -val
    
    PPV, TP, FP = get_ppv_helper(contacts_x, contacts_y, ref_contact_map, ref_len, factor)

    print '%s %s %s %s' % (hbond_filename, PPV, TP, FP)
    return (hbond_filename, PPV, TP, FP)



if __name__ == "__main__":

    p = argparse.ArgumentParser(description='Plot protein residue contact maps.')
    p.add_argument('fasta_file')#, required=True)
    p.add_argument('contact_file')#, required=True)
    p.add_argument('pdb')
    p.add_argument('-o', '--outfile', default='')
    p.add_argument('-f', '--factor', default=1.0, type=float)
    p.add_argument('-s', '--score', default=-1.0, type=float)
    p.add_argument('--chain', default='')
    p.add_argument('--noalign', action='store_true')
    p.add_argument('--name', default='')
    p.add_argument('--min_dist', default=5, type=int)
    p.add_argument('--print_dist', action='store_true')

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
    
    #if len(open(args['pdb']).readline().split(' ')) != 3:
    if True:
        get_ppv(args['fasta_file'], args['contact_file'], args['pdb'],
                args['factor'], chain=args['chain'], sep=sep,
                outfilename=args['outfile'], noalign=args['noalign'],
                min_score=args['score'], name=args['name'], min_dist=args['min_dist'], print_dist=args['print_dist'])
    else:
        get_ppv_hbond(args['fasta_file'], args['contact_file'],
                args['pdb'], args['factor'], sep=sep,
                outfilename=args['outfile'], min_score=args['score'])


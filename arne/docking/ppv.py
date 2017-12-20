#!/usr/bin/env python2


import sys
import argparse
from math import *
import Bio.PDB
from Bio import pairwise2

from os.path import expanduser
home = expanduser("~")
sys.path.append(home + '/git/bioinfo-toolbox/')
from parsing import parse_contacts
from parsing import parse_fasta
from parsing import parse_pdb

from  ppvlibrary import * 


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
    bfactor = parse_pdb.get_area(open(pdb_filename, 'r'), chain)
    surf = parse_pdb.get_dist_to_surface(open(pdb_filename, 'r'), chain)

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
        gapped_cb_lst = []

        ali_lst =[]
        j = 0
        k = 0
        for i in xrange(len(atom_seq_ali)):
            #print i,j,k,seq_ali[i],atom_seq_ali[i]
            if atom_seq_ali[i] == '-':
                gapped_cb_lst.append(['-'])
                ali_lst.append(-9999)
                k += 1
            elif seq_ali[i] == '-':
                j += 1
                continue
            else:
                ali_lst.append(j)
                gapped_cb_lst.append(cb_lst[j])
                k += 1
                j += 1

        dist_mat = get_cb_contacts(gapped_cb_lst)
        area = parse_pdb.get_area(open(pdb_filename, 'r'), chain)
        surf = parse_pdb.get_dist_to_surface(open(pdb_filename, 'r'), chain)
        if print_dist:
            print_distances(contacts_x, contacts_y, scores, dist_mat,
                                                area, surf, ref_len,ref_len,
                                                seq, fasta_filename, fasta_filename, ali_lst=ali_lst, atom_seq=atom_seq,
                                                outfile=outfilename)
        cb_cutoff = 8
        ref_contact_map = dist_mat < cb_cutoff
   
        PPV, TP, FP = get_ppv_helper(contacts_x, contacts_y, ref_contact_map, ref_len, factor, atom_seq_ali=atom_seq_ali)
    if name:
        print '%s %s %s %s' % (name, PPV, TP, FP)
    else:
        print '%s %s %s %s %s' % (fasta_filename, c_filename, PPV, TP, FP)
    return (pdb_filename, PPV, TP, FP)
  
    


if __name__ == "__main__":

    p = argparse.ArgumentParser(description='Plot protein residue contact maps.')
    p.add_argument('fasta_file')#, required=True)
    p.add_argument('pdb')
    p.add_argument('contact_file')#, required=True)
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


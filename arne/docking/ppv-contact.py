#!/usr/bin/env python2

import sys
import argparse
from math import *
import Bio.PDB
from Bio import pairwise2
import numpy as np


from os.path import expanduser
home = expanduser("~")
sys.path.append(home + '/git/bioinfo-toolbox/')

from parsing import parse_contacts
from parsing import parse_fasta
from parsing import parse_pdb


    
def get_cb_contacts(gapped_cb_lst):

    seqlen = len(gapped_cb_lst)
    dist_mat = np.zeros((seqlen, seqlen), np.float)
    dist_mat.fill(float('inf'))
    for i, cb1 in enumerate(gapped_cb_lst):
        if cb1[0] ==  '-':
            continue
        for j, cb2 in enumerate(gapped_cb_lst):
            if cb2[0] == '-':
                continue
            diff_vec = cb1 - cb2
            dist_mat[i,j] = np.sqrt(np.sum(diff_vec * diff_vec))
    return dist_mat


def print_distances(contacts_x, contacts_y, scores, dist_mat, area, lenA,lenB,atom_seq_ali=[], outfile=""):
    num_c = len(contacts_x)
    outstr = ""
    domain = ""
    for i in range(num_c):
        c_x = contacts_x[i]
        c_y = contacts_y[i]
        if atom_seq_ali:
            if atom_seq_ali[c_x] == '-':
                continue
            if atom_seq_ali[c_y] == '-':
                continue
        if c_x < lenA and c_y < lenA:
            domain="A"
        elif c_x > lenA and c_y > lenA:
            domain="B"
        else:
            domain="I"
        areaX=area[c_x][1]
        areaY=area[c_y][1]
        outstr += "%s %s %s %.2f %.2f %s %s\n" % (domain,c_x, c_y,areaX,areaY, scores[i], dist_mat[c_x, c_y])
    if outfile:
        with open(outfile, 'w') as outf:
            outf.write(outstr)
    else:
        print outstr

def get_ppv_helper(contacts_x, contacts_y, ref_contact_map, ref_len, factor, atom_seq_ali=[]):
    num_c = int(ref_len*factor)
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


def get_ppv_helper_interface(contacts_x, contacts_y, ref_contact_map, area, chainlenA, chainlenB,numcontacts,cutoff,atom_seq_ali=[]):
    num_c = numcontacts
    count=0
    TP = 0.0
    FP = 0.0
    PPV = 0.0
    countE=0
    TPe = 0.0
    FPe = 0.0
    PPVe = 0.0
    i=0
    for i in range(len(contacts_x)):
        c_x = contacts_x[i]
        c_y = contacts_y[i]
        if atom_seq_ali:
            if atom_seq_ali[c_x] == '-':
                continue
            if atom_seq_ali[c_y] == '-':
                continue
        
        if ( c_x < chainlenA and c_y > chainlenA) or ( c_x > chainlenA and c_y < chainlenA):
            #print "I: ",i,c_x,c_y,ref_contact_map[c_x, c_y]
            if (count < numcontacts):
                count += 1
                if ref_contact_map[c_x, c_y] > 0:
                    TP += 1.0 
                else:
                    FP += 1.0 
            if (countE < numcontacts):
                if area[c_x][1]>cutoff and area[c_y][1]>cutoff:
                    countE += 1
                    if ref_contact_map[c_x, c_y] > 0:
                        TPe += 1.0
                    else:
                        FPe += 1.0 

    if TP > 0:
        PPV = TP / (TP + FP)
    if TPe > 0:
        PPVe = TPe / (TPe + FPe)
    return (PPV, TP, FP,PPVe, TPe, FPe)


    
def get_ppv(fasta_filenameA, c_filename, pdb_filenameA,
            fasta_filenameB, pdb_filenameB,factor=1.0, min_score=-1.0,
            chainA='',chainB='', sep=' ', outfilename='', name='', noalign=False,
            min_dist=5, interfacelen=10, print_dist=False,cutoff=0.25):  
    

    ### get sequence
    seqA = parse_fasta.read_fasta(open(fasta_filenameA, 'r')).values()[0][0]
    seqB = parse_fasta.read_fasta(open(fasta_filenameB, 'r')).values()[0][0]
    seq=seqA+seqB

    ref_lenA = len(seqA)
    ref_lenB = len(seqB)
    ref_len = len(seq)

    ### get top ranked predicted contacts
    contacts = parse_contacts.parse(open(c_filename, 'r'), sep, min_dist=min_dist)

    contacts_x = []
    contacts_y = []
    scores = []
    contactsA_x = []
    contactsA_y = []
    scoresA = []
    contactsB_x = []
    contactsB_y = []
    scoresB = []
    contactsI_x = []
    contactsI_y = []
    scoresI = []
    contact_dict = {}

    count = 0
    countA = 0
    countB = 0
    countI = 0
    for i in range(len(contacts)):
        score = contacts[i][0]
        c_x = contacts[i][1] - 1
        c_y = contacts[i][2] - 1
        #print i,c_x,c_y,score

        pos_diff = abs(c_x - c_y)
        too_close = pos_diff < min_dist

        if not too_close:
            # Check if it is protein A or B or interchain contact
            contacts_x.append(c_x)
            contacts_y.append(c_y)
            scores.append(score)
            count += 1
            if (c_x < ref_lenA and c_y < ref_lenA):
                contactsA_x.append(c_x)
                contactsA_y.append(c_y)
                scoresA.append(score)
                countA += 1
            elif (c_x >= ref_lenA and c_y >= ref_lenA):
                contactsB_x.append(c_x-ref_lenA)
                contactsB_y.append(c_y-ref_lenA)
                scoresB.append(score)
                countB += 1
            else:
                contactsI_x.append(c_x)
                contactsI_y.append(c_y)
                scoresI.append(score)
                countI += 1
           
                #        if min_score == -1.0 and count >= ref_len * factor:
                #            break
                #        if score < min_score:
                #            break
                
    assert(len(contactsA_x) == len(contactsA_y) == len(scoresA))
    assert(len(contactsB_x) == len(contactsB_y) == len(scoresB))
    assert(len(contactsI_x) == len(contactsI_y) == len(scoresI))

    cb_lstA = parse_pdb.get_cb_coordinates(open(pdb_filenameA, 'r'), chainA)
    cb_lstB = parse_pdb.get_cb_coordinates(open(pdb_filenameB, 'r'), chainB)
    cb_lst=cb_lstA+cb_lstB
    bfactorA = parse_pdb.get_area(open(pdb_filenameA, 'r'), chainA)
    bfactorB = parse_pdb.get_area(open(pdb_filenameB, 'r'), chainA)
    bfactor = bfactorA+bfactorB
    #print cb_lst,noalign
    if noalign:
        dist_mat = get_cb_contacts(cb_lst)
        dist_matA = get_cb_contacts(cb_lstA)
        dist_matB = get_cb_contacts(cb_lstB)
        #PPV, TP, FP = get_ppv_helper(contacts_x, contacts_y, ref_contact_map, ref_len, factor)
    else:
        atom_seqA = parse_pdb.get_atom_seq(open(pdb_filenameA, 'r'), chainA)
        atom_seqB = parse_pdb.get_atom_seq(open(pdb_filenameB, 'r'), chainB)
        atom_seq = atom_seqA + atom_seqB
        align = pairwise2.align.globalms(atom_seq, seq, 2, -1, -0.5, -0.1)
        alignA = pairwise2.align.globalms(atom_seqA, seqA, 2, -1, -0.5, -0.1)
        alignB = pairwise2.align.globalms(atom_seqB, seqB, 2, -1, -0.5, -0.1)
        atom_seq_ali = align[-1][0]
        seq_ali = align[-1][1]
        atom_seq_aliA = alignA[-1][0]
        seq_aliA = alignA[-1][1]
        atom_seq_aliB = alignB[-1][0]
        seq_aliB = alignB[-1][1]
        j = 0
        gapped_cb_lst = []
        gapped_cb_lstA = []
        gapped_cb_lstB = []

        for i in xrange(len(atom_seq_ali)):
            if atom_seq_ali[i] == '-':
                gapped_cb_lst.append(['-'])
            elif seq_ali[i] == '-':
                j += 1
                continue
            else:
                gapped_cb_lst.append(cb_lst[j])
                j += 1
        j = 0
        for i in xrange(len(atom_seq_aliA)):
            if atom_seq_aliA[i] == '-':
                gapped_cb_lstA.append(['-'])
            elif seq_aliA[i] == '-':
                j += 1
                continue
            else:
                gapped_cb_lstA.append(cb_lst[j])
                j += 1
        j = 0
        for i in xrange(len(atom_seq_aliB)):
            if atom_seq_aliB[i] == '-':
                gapped_cb_lstB.append(['-'])
            elif seq_aliB[i] == '-':
                j += 1
                continue
            else:
                gapped_cb_lstB.append(cb_lst[j])
                j += 1

        dist_mat = get_cb_contacts(gapped_cb_lst)
        dist_matA = get_cb_contacts(gapped_cb_lstA)
        dist_matB = get_cb_contacts(gapped_cb_lstB)
    cb_cutoff = 8
    ref_contact_map = dist_mat < cb_cutoff
    ref_contact_mapA = dist_matA < cb_cutoff
    ref_contact_mapB = dist_matB < cb_cutoff

    if print_dist:
        print_distances(contacts_x, contacts_y, scores, dist_mat, bfactor, ref_lenA,ref_lenB,atom_seq_ali=atom_seq_ali, outfile=outfilename)

    PPV, TP, FP = get_ppv_helper(contacts_x, contacts_y, ref_contact_map, ref_len, factor, atom_seq_ali=atom_seq_ali)
    PPVa, TPa, FPa = get_ppv_helper(contactsA_x, contactsA_y, ref_contact_mapA, interfacelen, factor, atom_seq_ali=atom_seq_aliA)
    PPVb, TPb, FPb = get_ppv_helper(contactsB_x, contactsB_y, ref_contact_mapB, interfacelen, factor, atom_seq_ali=atom_seq_aliB)
    PPVi, TPi, FPi, PPViE, TPiE, FPiE = get_ppv_helper_interface(contacts_x, contacts_y, ref_contact_map, bfactor, ref_lenA,ref_lenB, interfacelen, cutoff, atom_seq_ali=atom_seq_ali)
    #for i in range(10):
    #    print "I: ",i,contactsI_x[i],contactsI_y[i],scoresI[i],dist_mat[contactsI_x[i]][contactsI_y[i]],ref_contact_map[contactsI_x[i]][contactsI_y[i]]
    #    print "A: ",i,contactsA_x[i],contactsA_y[i],scoresA[i],dist_mat[contactsA_x[i]][contactsA_y[i]],ref_contact_map[contactsA_x[i]][contactsA_y[i]]
    #    print "B: ",i,contactsB_x[i],contactsB_y[i],scoresB[i],dist_mat[contactsB_x[i]][contactsB_y[i]],ref_contact_map[contactsB_x[i]][contactsB_y[i]]

    if name:
        print '%s %s %s %s' % (name, PPVa, TPa, FPa)
        print '%s %s %s %s' % (name, PPVb, TPb, FPb)
        print '%s %s %s %s' % ("BOTH", PPV, TP, FP)
        print '%s %s %s %s' % ("Interface", PPVi, TPi, FPi)
        print '%s %s %s %s' % ("Interface Exposed", PPViE, TPiE, FPiE)
    else:
        print '%s %s %s %s %s' % (fasta_filenameA, c_filename, PPVa, TPa, FPa)
        print '%s %s %s %s %s' % (fasta_filenameB, c_filename, PPVb, TPb, FPb)
        print '%s %s %s %s %s' % ("BOTH", c_filename, PPV, TP, FP)
        print '%s %s %s %s %s' % ("Interface", c_filename, PPVi, TPi, FPi)
        print '%s %s %s %s %s' % ("Interface Exposed", c_filename, PPViE, TPiE, FPiE)
    return (pdb_filenameA, PPV, TP, FP)
  
    


if __name__ == "__main__":

    p = argparse.ArgumentParser(description='Plot protein residue contact maps.')
    p.add_argument('fasta_fileA')#, required=True)
    p.add_argument('fasta_fileB')#, required=True)
    p.add_argument('pdbA')
    p.add_argument('pdbB')
    p.add_argument('contact_file')#, required=True)
    p.add_argument('-o', '--outfile', default='')
    p.add_argument('-f', '--factor', default=1.0, type=float)
    p.add_argument('-s', '--score', default=-1.0, type=float)
    p.add_argument('--chainA', default='')
    p.add_argument('--chainB', default='')
    p.add_argument('--noalign', action='store_true')
    p.add_argument('--name', default='')
    p.add_argument('--min_dist', default=5, type=int)
    p.add_argument('--interfacelen', default=10, type=int)
    p.add_argument('--cutoff', default=0.25, type=float)
    p.add_argument('--print_dist', action='store_true')

    args = vars(p.parse_args(sys.argv[1:]))

    fasta_filenameA = args['fasta_fileA']
    fasta_filenameB = args['fasta_fileB']
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
    get_ppv(args['fasta_fileA'], args['contact_file'], args['pdbA'],
            args['fasta_fileB'], args['pdbB'], args['factor'],
            chainA=args['chainA'], chainB=args['chainB'], sep=sep,
            outfilename=args['outfile'], noalign=args['noalign'],
            min_score=args['score'], name=args['name'],
            min_dist=args['min_dist'],
            interfacelen=args['interfacelen'],print_dist=args['print_dist'],cutoff=args['cutoff'])


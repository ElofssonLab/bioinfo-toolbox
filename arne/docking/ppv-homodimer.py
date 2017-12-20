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

    
def get_ppv(fasta_filenameA, c_filename, pdb_filenameA,
            fasta_filenameB, pdb_filenameB,factor=1.0, min_score=-1.0,
            chainA='',chainB='', sep=' ', outfilename='', name='', noalign=False,
            min_dist=5, interfacelen=10, print_dist=False,cutoff=0.25):  
    

    ### get sequence
    seqA = parse_fasta.read_fasta(open(fasta_filenameA, 'r')).values()[0][0]
    seqB = parse_fasta.read_fasta(open(fasta_filenameB, 'r')).values()[0][0]
    seq=seqA+seqA # Actually the contact map sequence is just two copies of seqA

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
            # The contacts only covers 
            contacts_x.append(c_x)
            contacts_y.append(c_y)
            scores.append(score)
            #contacts_x.append(c_x+ref_lenA)
            #contacts_y.append(c_y+ref_lenA)
            #scores.append(score)
            contactsA_x.append(c_x)
            contactsA_y.append(c_y)
            scoresA.append(score)
            contactsB_x.append(c_x)
            contactsB_y.append(c_y)
            scoresB.append(score)

           
                #        if min_score == -1.0 and count >= ref_len * factor:
                #            break
                #        if score < min_score:
                #            break
                
    assert(len(contacts_x) == len(contacts_y) == len(scores))
    assert(len(contactsA_x) == len(contactsA_y) == len(scoresA))
    assert(len(contactsB_x) == len(contactsB_y) == len(scoresB))
    assert(len(contactsI_x) == len(contactsI_y) == len(scoresI))

    cb_lstA = parse_pdb.get_cb_coordinates(open(pdb_filenameA, 'r'), chainA)
    cb_lstB = parse_pdb.get_cb_coordinates(open(pdb_filenameB, 'r'), chainB)
    cb_lst=cb_lstA+cb_lstB
    bfactorA = parse_pdb.get_area(open(pdb_filenameA, 'r'), chainA)
    bfactorB = parse_pdb.get_area(open(pdb_filenameB, 'r'), chainB)
    bfactor = bfactorA+bfactorB
    surfA = parse_pdb.get_dist_to_surface(open(pdb_filenameA, 'r'), chainA)
    surfB = parse_pdb.get_dist_to_surface(open(pdb_filenameB, 'r'), chainB)
    surf = surfA+surfB
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
        alignB = pairwise2.align.globalms(atom_seqB, seqA, 2, -1, -0.5, -0.1) # Align to seq A
        atom_seq_ali = align[-1][0]
        seq_ali = align[-1][1]
        atom_seq_aliA = alignA[-1][0]
        seq_aliA = alignA[-1][1]
        atom_seq_aliB = alignB[-1][0]
        seq_aliB = alignB[-1][1]
        gapped_cb_lst = []
        gapped_cb_lstA = []
        gapped_cb_lstB = []
        ali_lst =[]
        ali_lstA =[]
        ali_lstB =[]
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
        j = 0
        k = 0
        for i in xrange(len(atom_seq_aliA)):
            if atom_seq_aliA[i] == '-':
                gapped_cb_lstA.append(['-'])
                ali_lstA.append(-9999)
                k += 1
            elif seq_aliA[i] == '-':
                j += 1
                continue
            else:
                ali_lstA.append(j)
                gapped_cb_lstA.append(cb_lstA[j])
                k += 1
                j += 1
        j = 0
        k = 0
        for i in xrange(len(atom_seq_aliB)):
            #print "B",i,j,k,seq_aliB[i],atom_seq_aliB[i]
            if atom_seq_aliB[i] == '-':
                gapped_cb_lstB.append(['-'])
                ali_lstB.append(-9999)
                k += 1
            elif seq_aliB[i] == '-':
                j += 1
                continue
            else:
                ali_lstB.append(j)
                gapped_cb_lstB.append(cb_lstB[j])
                k += 1
                j += 1


        #print len(gapped_cb_lst),len(gapped_cb_lstA),len(gapped_cb_lstB)
        dist_mat = get_cb_contacts(gapped_cb_lst)
        dist_matA = get_cb_contacts(gapped_cb_lstA)
        dist_matB = get_cb_contacts(gapped_cb_lstB)
    cb_cutoff = 8
    #ref_contact_map = dist_mat < cb_cutoff
    # This routine adds all interface and B chain contacts
    contacts_x,contacts_y,scores = get_interface_contacts(contacts_x, contacts_y, scores, dist_mat, ref_lenA,
                                                          factor, cb_cutoff+4,atom_seq_ali=atom_seq_ali)
    ref_contact_map = dist_mat < cb_cutoff  
    ref_contact_mapA = dist_matA < cb_cutoff
    ref_contact_mapB = dist_matB < cb_cutoff
    # Here we need to append
    if print_dist:
        print_distances(contacts_x, contacts_y, scores, dist_mat,
                        bfactor, surf, ref_lenA,ref_lenB,
                        seq, fasta_filenameA, fasta_filenameB, ali_lst=ali_lst, atom_seq=atom_seq,
                        outfile=outfilename)

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
    print 'PPV %s %s %s %s %s %s' % (c_filename, PPV, PPVa, PPVb, PPVi, PPViE)
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


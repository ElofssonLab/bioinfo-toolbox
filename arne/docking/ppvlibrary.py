import sys
import re
from math import *
import Bio.PDB
import numpy as np

from os.path import expanduser
home = expanduser("~")
sys.path.append(home + '/git/bioinfo-toolbox/')



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


def print_distances(contacts_x, contacts_y, scores, dist_mat, area, dist, lenA,lenB,seq,nameA,nameB,atom_seq_ali=[], outfile=""):
    num_c = len(contacts_x)
    outstr = ""
    domain = ""
    codeA=nameA
    codeB=nameB
    
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
        distX=dist[c_x]
        distY=dist[c_y]
        outstr += "%s %s %s %s %s %.2f %.2f %.2f %.2f %s %s %s %s \n" % (domain,c_x,seq[c_x], c_y,seq[c_y],
                                                                  areaX,areaY,distX,distY, scores[i],
                                                                         dist_mat[c_x, c_y],codeA,codeB)
    if outfile:
        with open(outfile, 'w') as outf:
            outf.write(outstr)
    else:
        print outstr

def get_interface_contacts(contacts_x, contacts_y, scores, dist_mat, ref_len, factor, cutoff, atom_seq_ali=[]):
    for i in range(len(contacts_x)):
        c_x = contacts_x[i]
        c_y = contacts_y[i]
        if atom_seq_ali:
            if atom_seq_ali[c_x] == '-':
                continue
            if atom_seq_ali[c_y] == '-':
                continue
        if dist_mat[c_x, c_y] > cutoff and c_x < ref_len:
            contacts_x.append(c_x)
            contacts_y.append(c_y+ref_len)
            scores.append(scores[i])
    return (contacts_x,contacts_y,scores)


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


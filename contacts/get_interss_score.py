from math import *
import sys
sys.path.append('/bubo/sw/apps/bioinfo/biopython/1.59/tintin/lib/python')

import Bio.PDB
from Bio import pairwise2

import numpy as np
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform

#import rpy2

import matplotlib
matplotlib.use('Agg')
from matplotlib import pylab
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.cbook as cbook

import parse_contacts
import parse_psipred
import parse_fasta
#import parse_jones
import parse_pdb

# check: two residues are in the same helix iff
#        they are in a helix and all residues inbetween are in a helix
def in_same_helix(res1, res2, ss):
    result = True
    first_res = min(res1, res2) + 3
    sec_res = max(res1, res2) - 3
    for i in range(first_res - 1, sec_res):
        if ss[i] != 'H':
            result = False
            break
    return result


def in_dom(i, pos_lst):
    #offset = max(1, pos_lst[0][0] + 10)
    pos_i = 1
    for pos in pos_lst:
        #if i >= (pos[0] - offset - 1) and i < (pos[1] - offset):
        if i >= (pos[0] - 1) and i < pos[1]:
            return pos_i
        pos_i += 1
    pos_i = 0
    return pos_i


def get_min_dist(res1, res2):
    
    min_dist = float('inf')

    for atm1 in res1:
        for atm2 in res2:
            #print '%s - %s' % (atm1, atm2)
            diff_vec = atm1 - atm2
            dist = np.sqrt(np.sum(diff_vec * diff_vec))
            if dist < min_dist:
                min_dist = dist

    return min_dist


def get_heavy_contacts(gapped_res_lst):

    seqlen = len(gapped_res_lst)
    #print seqlen
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


def get_cb_contacts(gapped_cb_lst):

    seqlen = len(gapped_cb_lst)
    #print seqlen
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


def get_dom_seq(acc, ref_len, query_file):
    
    dom_seq = ''
    dom_lst = []

    for line in query_file:
        line_arr = line.split('\t')
        curr_acc = line_arr[0]

        if curr_acc != acc:
            continue
        
        pos_lst = eval(line_arr[1])
        offset = pos_lst[0][0] + 10
        for pos in pos_lst:
            dom_lst.append((max(1, pos[0] - offset), pos[1] - offset))

        for i in xrange(ref_len):
            dom_i = in_dom(i, pos_lst)
            dom_seq += str(dom_i)

        return dom_seq, dom_lst


def get_ss_pos(ss_seq):
    
    ss_pos_lst = []
    prev_ss = ss_seq[0]
    start = 1
    end = len(ss_seq)
    in_ss_flag = False

    for i, ss in enumerate(ss_seq):
        if ss != 'C' and ss != prev_ss: # change in ss => new element
            if in_ss_flag: # no coil between secstruct elements
                end = i
                ss_pos_lst.append((start, end))
            start = i + 1
            in_ss_flag = True
        elif in_ss_flag and ss != prev_ss: # change + already in ss => element ends
            end = i 
            ss_pos_lst.append((start, end))
            in_ss_flag = False
        prev_ss = ss

    if len(ss_pos_lst) == 0:
        ss_pos_lst.append(start, end)
    return ss_pos_lst


def plot_map(acc, fasta_filename, contact_filename, psipred_filename, sep=',', pdb_filename='', chain='A'):  
   
    pdb_flag = pdb_filename.strip() != ''

    rep_len = 1000
    psipred_filename = '%s.horiz' % '.'.join(fasta_filename.split('.')[:-1])
    ss = parse_psipred.horizontal(open(psipred_filename, 'r'))
    ss_lst = get_ss_pos(ss)
    print ss
    print ss_lst
    seq = parse_fasta.read_fasta(open(fasta_filename, 'r')).values()[0][0]
    #ss = seq
    ref_len = len(seq)

    if pdb_flag:
        #pdb_code = pdb_filename.split('/')[-1].split('.')[0]
        res_lst = parse_pdb.get_coordinates(open(pdb_filename, 'r'), chain)
        #cb_lst = parse_pdb.get_cb_coordinates(open(pdb_filename, 'r'), chain)
        atom_seq = parse_pdb.get_atom_seq(open(pdb_filename, 'r'), chain)

        #print seq
        #print atom_seq
        align = pairwise2.align.globalms(atom_seq, seq, 2, -1, -0.5, -0.1)
        #print align[-1]
        atom_seq_ali = align[-1][0]
        seq_ali = align[-1][1]
        print atom_seq_ali
        print seq_ali

        """
        # get 1D-domain assignments as char sequence
        # 'D' = in domain / 'N' = not in domain
        acc = fasta_filename.split('/')[-1].split('.')[0]
        dom_seq, dom_lst = get_dom_seq(acc, ref_len, open('/bubo/home/h9/mircomic/glob/2013-05-29_human_repeats/query.txt', 'r'))
        dom_seq_lst = map(int, list(dom_seq))
        print dom_seq
        """

        #dist_mat = calc_dist_matrix_heavy(ref_chain, ref_chain)
        j = 0
        gapped_res_lst = []
        #gapped_cb_lst = []
        for i in xrange(len(atom_seq_ali)):
            if atom_seq_ali[i] == '-':
                gapped_res_lst.append('-')
                #gapped_cb_lst.append('-')
            elif seq_ali[i] == '-':
                j += 1
                continue
            else:
                gapped_res_lst.append(res_lst[j])
                #gapped_cb_lst.append(cb_lst[j])
                j += 1
        #print gapped_res_lst
        print len(gapped_res_lst)
        print len(res_lst)
        #print len(gapped_cb_lst)
        #print len(cb_lst)
        print len(atom_seq)
        print len(seq)
        dist_mat = get_heavy_contacts(gapped_res_lst)
        #dist_mat = get_cb_contacts(gapped_cb_lst)
        heavy_cutoff = 5
        cb_cutoff = 8

        #ref_contact_map = (dist_mat < 8) & (dist_mat > 4)
        ref_contact_map = dist_mat < heavy_cutoff
        ref_contacts = np.where(dist_mat < heavy_cutoff)
        #ref_contact_map = dist_mat < cb_cutoff
        #ref_contacts = np.where(dist_mat < cb_cutoff)

        ref_contacts_x = ref_contacts[0]
        ref_contacts_y = ref_contacts[1]
        
    """
    tmp_x = []
    tmp_y = []

    for i in range(len(ref_contacts_x)):
        x = ref_contacts_x[i]
        y = ref_contacts_y[i]
        if y > x:
            tmp_x.append(x)
            tmp_y.append(y)

    ref_contacts_x = tmp_x
    ref_contacts_y = tmp_y
    """
    #print dist_mat
    #print ref_contact_map
    #print 'ref_contacts=' + str(ref_contacts[1])

    contacts = parse_contacts.parse(open(contact_filename, 'r'), sep)
    #contacts_cut = contacts[0:ref_len]

    contacts_x = []
    contacts_y = []
    scores = []
    contact_dict = {}

    count = 0
    #for i in range(ref_len * 1):
    for i in range(len(contacts)):
        score = contacts[i][0]
        c_x = contacts[i][1] - 1
        c_y = contacts[i][2] - 1

        pos_diff = abs(c_x - c_y)
        too_close = pos_diff < 5
        too_far = pos_diff > rep_len * 1.5

        if not too_close:
            contacts_x.append(c_x)
            contacts_y.append(c_y)
            scores.append(score)
            count += 1
           
        if count >= ref_len * 1.0:
            break

    if pdb_flag:
        PPVs = []
        inter_PPVs = []
        intra_PPVs = []
        TPs = []
        inter_TPs = []
        intra_TPs = []
        FPs = []
        inter_FPs = []
        intra_FPs = []

        for num_c in range(min(len(contacts_x), ref_len * 1))[1:]:
            TP = 0.0
            intra_TP = 0.0
            inter_TP = 0.0
            FP = 0.0
            intra_FP = 0.0
            inter_FP = 0.0
            for i in range(num_c):
                c_x = contacts_x[i]
                c_y = contacts_y[i]
                if atom_seq_ali[c_x] == '-':
                    continue
                if atom_seq_ali[c_y] == '-':
                    continue
                c_x_dom = in_dom(c_x, dom_lst)
                c_y_dom = in_dom(c_y, dom_lst)
                #print c_x_dom
                #print c_y_dom
                if ref_contact_map[c_x, c_y] > 0:
                    TP += 1.0
                    if c_x_dom != c_y_dom and c_x_dom != 0 and c_y_dom != 0:
                        inter_TP += 1.0
                    if c_x_dom == c_y_dom and c_x_dom != 0 and c_y_dom != 0:
                        intra_TP += 1.0
                else:
                    FP += 1.0
                    if c_x_dom != c_y_dom and c_x_dom != 0 and c_y_dom != 0:
                        inter_FP += 1.0
                    if c_x_dom == c_y_dom and c_x_dom != 0 and c_y_dom != 0:
                        intra_FP += 1.0

            #print '%s, %s, %s, %s, %s, %s' % (TP, FP, inter_TP, inter_FP, intra_TP, intra_FP)
            if TP > 0 and FP > 0:
                PPVs.append(TP / (TP + FP))
                if inter_TP > 0 or inter_FP > 0:
                    inter_PPVs.append(inter_TP / (inter_TP + inter_FP))
                if intra_TP > 0 or intra_FP > 0:
                    intra_PPVs.append(intra_TP / (intra_TP + intra_FP))
                TPs.append(TP / ref_len)
                FPs.append(FP / ref_len)

        print len(PPVs)
        if len(PPVs) > 0:
            print PPVs[-1]
        else:
            PPVs.append(0.0)
            print PPVs[-1]
        
        if len(inter_PPVs) > 0:
            print inter_PPVs[-1]
        else:
            inter_PPVs.append(0.0)
            print inter_PPVs[-1]
        
        if len(intra_PPVs) > 0:
            print intra_PPVs[-1]
        else:
            intra_PPVs.append(0.0)
            print intra_PPVs
        
        #print TPs[-1]
        #print FPs[-1]
        
        PPV_file = open('/bubo/home/h9/mircomic/glob/2013-05-29_human_repeats/pconsc_predictions/PPV.txt', 'a')
        PPV_file.write('%s\t%s\n' % (acc, PPVs[-1]))
        PPV_file.close()
        intra_PPV_file = open('/bubo/home/h9/mircomic/glob/2013-05-29_human_repeats/pconsc_predictions/intra_PPV.txt', 'a')
        intra_PPV_file.write('%s\t%s\n' % (acc, intra_PPVs[-1]))
        intra_PPV_file.close()
        inter_PPV_file = open('/bubo/home/h9/mircomic/glob/2013-05-29_human_repeats/pconsc_predictions/inter_PPV.txt', 'a')
        inter_PPV_file.write('%s\t%s\n' % (acc, inter_PPVs[-1]))
        inter_PPV_file.close()


    ### get pairwise distances from the alignment
    #alinum = parse_jones.get_numeric('%s.jones' % '.'.join(contact_filename.split('.')[:-1]))
    #alidist = squareform(pdist(alinum.T, 'euclidean'))
    #alidistlog = np.log(alidist)

    """
    dot_matrix = dotter.calc_dot_matrix(seq)
    tmp_dot_matrix = dot_matrix
    for (i,j), score in np.ndenumerate(dot_matrix):
        if i > j:
            tmp_dot_matrix[i, j] = 0.0
    dot_matrix = tmp_dot_matrix
    """

    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    for i in range(len(ss)):
        if ss[i] == 'H':
            plt.plot(i, i, 'o', c='#8B0043', mec="#8B0043", markersize=2)
            #plt.plot(i, i, 'o', c='#999999', mec="#444444")
        if ss[i] == 'E':
            plt.plot(i, i, 'D', c='#0080AD', mec="#0080AD", markersize=2)
            #plt.plot(i, i, 'D', c='#999999', mec="#444444", markersize=6)
        if ss[i] == 'C':
            continue
            #plt.plot(i, i, 'D', c='#999999', mec='#999999', markersize=3)
            #plt.plot(i, i, 'D', c='#999999', mec='#999999', markersize=3)

    """
    for dom in dom_lst:
        start = dom[0] - 1
        end = dom[1] - 1
        dom_len = end - start
        print '%d - %d' % (start, end)
        ax.add_patch(plt.Rectangle((start, start), dom_len, dom_len, facecolor='#EEEEEE', edgecolor='black', lw=0.5, zorder=0))
        for dom2 in dom_lst:
            if dom2 == dom:
                continue
            start2 = dom2[0] - 1
            end2 = dom2[1] - 1
            dom_len2 = end2 - start2
            ax.add_patch(plt.Rectangle((start, start2), dom_len, dom_len2, facecolor='#EEEEEE', edgecolor='black', lw=0.5, zorder=0, alpha=0.3))
    """

    ss_col_dict = {'H': '#8B0043', 'E': '#0080AD', 'HE': '#3F467D'}
    for ss_elem in ss_lst:
        start = ss_elem[0] - 1
        end = ss_elem[1] - 1
        elem_len = end - start
        curr_ss = ss[start]
        #print '%d - %d' % (start, end)
        ax.add_patch(plt.Rectangle((start, start), elem_len, elem_len, facecolor=ss_col_dict[curr_ss], edgecolor='black', lw=0.5, zorder=0, alpha=0.5))
        for ss_elem2 in ss_lst:
            if ss_elem2 == ss_elem:
                continue
            start2 = ss_elem2[0] - 1
            end2 = ss_elem2[1] - 1
            elem_len2 = end2 - start2
            curr_ss2 = ss[start2]
            if curr_ss != curr_ss2: # interaction between sheet and helix
                curr_ss2 = 'HE'
            ax.add_patch(plt.Rectangle((start, start2), elem_len, elem_len2, facecolor=ss_col_dict[curr_ss2], edgecolor='black', lw=0.5, zorder=0, alpha=0.2))

    #ax.imshow(dot_matrix, origin='lower', cmap=cm.binary)
    if pdb_flag:
        ax.scatter(ref_contacts_x, ref_contacts_y, marker='o', c='#CCCCCC', lw=0)
        #ax.scatter(range(ref_len), range(ref_len), marker='d', c=dom_seq_lst, lw=0, edgecolor=dom_seq_lst, cmap=cm.spectral_r)
        #plt.plot(ref_contacts_x, ref_contacts_y, 'o', c='#CCCCCC', mec='#CCCCCC')
        fig.suptitle('%s\nPPV = %.2f   intra-PPV = %.2f   inter-PPV = %.2f' % (contact_filename, PPVs[-1], intra_PPVs[-1], inter_PPVs[-1]))
    else: 
        fig.suptitle('%s - %s' % (acc, contact_filename))
    
    sc = ax.scatter(contacts_x[::-1], contacts_y[::-1], marker='o', c=scores[::-1], s=4, alpha=0.75, cmap=cm.jet, linewidths=0.5)
    #sc = ax.scatter(contacts_nf_y[::-1], contacts_nf_x[::-1], marker='o', c=scores_nf[::-1], s=8, alpha=0.75, cmap=cm.jet, linewidths=0.5)

    plt.gca().set_xlim([0,ref_len])
    plt.gca().set_ylim([0,ref_len])
    plt.colorbar(sc)
    #cbar = plt.colorbar(ax, ticks=[min(scores), max(scores)])
    #cbar.ax.set_yticklabels(['Low', 'High'])
    #cbar.set_label(r'Contact Score')

    pp = PdfPages('%s_ContactMap.pdf' % contact_filename)
    pp.savefig(fig)
    pp.close()

    #outfile = open('%s.contacts' % '.'.join(contact_filename.split('.')[0:-1]),'w')
    #for i in range(len(scores)):
    #    outfile.write('%s,%s,%s\n' % (int(contacts_x[i] + 1), int(contacts_y[i] + 1), scores[i]))

    if pdb_flag:
        fig.clf()
        ax2 = fig.add_subplot(111)
        ax2.plot(PPVs)
        pp = PdfPages('%s_PPVs.pdf' % contact_filename)
        pp.savefig(fig)
        pp.close() 


if __name__ == "__main__":

    if len(sys.argv) < 3:
        sys.stderr.write('Usage: python plot_contact_map.py <fasta_filename> <contact_filename> <psipred_filename> <pdb_filename>\n')

    acc = sys.argv[1]
    fasta_filename = sys.argv[2]
    contact_filename = sys.argv[3]
    psipred_filename = sys.argv[4]

    # guessing separator of constraint file
    line = open(contact_filename,'r').readline()
    if len(line.split(',')) != 1:
        sep = ','
    elif len(line.split(' ')) != 1:
        sep = ' '
    else:
        sep = '\t'

    if len(sys.argv) == 6:
        pdb_filename = sys.argv[5]
        chain = sys.argv[6]
        plot_map(acc, fasta_filename, contact_filename, psipred_filename, sep, pdb_filename, chain)
    else:
        #rep_len = int(sys.argv[3])
        #print "sep=" + sep
        plot_map(acc, fasta_filename, contact_filename, psipred_filename, sep)


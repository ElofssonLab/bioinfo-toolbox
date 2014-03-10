import sys
import argparse
from math import *

import Bio.PDB
from Bio import pairwise2

import numpy as np

import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.cm as cm

sys.path.append('/home/x_mirmi/bioinfo-toolbox')
sys.path.append('/home/mircomic/bioinfo-toolbox')

from parsing import parse_contacts
from parsing import parse_psipred
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


def get_heavy_contacts(gapped_res_lst):

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



def get_ppvs(contacts_x, contacts_y, ref_contact_map, atom_seq_ali, ref_len, factor):

    PPVs = []
    TPs = []
    FPs = []

    for num_c in range(min(len(contacts_x), ref_len * factor))[1:]:
        TP = 0.0
        FP = 0.0
        for i in range(num_c):
            c_x = contacts_x[i]
            c_y = contacts_y[i]
            if atom_seq_ali[c_x] == '-':
                continue
            if atom_seq_ali[c_y] == '-':
                continue
            if ref_contact_map[c_x, c_y] > 0:
                TP += 1.0
            else:
                FP += 1.0

        if TP > 0 and FP > 0:
            PPVs.append(TP / (TP + FP))
            TPs.append(TP / ref_len)
            FPs.append(FP / ref_len)

    if len(PPVs) == 0:
        PPVs.append(0.0)

    return PPVs


def get_tp_colors(contacts_x, contacts_y, ref_contact_map, atom_seq_ali):

    tp_colors = []

    for i in range(len(contacts_x)):
        c_x = contacts_x[i]
        c_y = contacts_y[i]
        if atom_seq_ali[c_x] == '-':
            #tp_colors.append('green')
            tp_colors.append('red')
            continue
        if atom_seq_ali[c_y] == '-':
            #tp_colors.append('green')
            tp_colors.append('red')
            continue
        if ref_contact_map[c_x, c_y] > 0:
            tp_colors.append('blue')
        else:
            tp_colors.append('red')

    return tp_colors
 

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


def plot_map(fasta_filename, c_filename, factor, c2_filename='', psipred_filename='', pdb_filename='', is_heavy=False, chain='', sep=',', outfilename=''):  
    #acc = c_filename.split('.')[0]
    acc = fasta_filename.split('.')[0][:4]

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
           
        if count >= ref_len * factor:
            break
 

    ### start plotting
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ### plot secondary structure on the diagonal if given
    if psipred_filename:
        ss = parse_psipred.horizontal(open(psipred_filename, 'r'))
        ss_lst = get_ss_pos(ss)
        for i in range(len(ss)):
            if ss[i] == 'H':
                plt.plot(i, i, 'o', c='#8B0043', mec="#8B0043", markersize=2)
            if ss[i] == 'E':
                plt.plot(i, i, 'D', c='#0080AD', mec="#0080AD", markersize=2)
            if ss[i] == 'C':
                continue

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



    ### plot reference contacts in the background if given
    if pdb_filename:
        res_lst = parse_pdb.get_coordinates(open(pdb_filename, 'r'), chain)
        cb_lst = parse_pdb.get_cb_coordinates(open(pdb_filename, 'r'), chain)
        atom_seq = parse_pdb.get_atom_seq(open(pdb_filename, 'r'), chain)
                
        align = pairwise2.align.globalms(atom_seq, seq, 2, -1, -0.5, -0.1)

        atom_seq_ali = align[-1][0]
        seq_ali = align[-1][1]

        j = 0
        gapped_res_lst = []
        gapped_cb_lst = []

        for i in xrange(len(atom_seq_ali)):
            if atom_seq_ali[i] == '-':
                gapped_res_lst.append('-')
                gapped_cb_lst.append('-')
            elif seq_ali[i] == '-':
                j += 1
                continue
            else:
                gapped_res_lst.append(res_lst[j])
                gapped_cb_lst.append(cb_lst[j])
                j += 1

        if is_heavy:
            dist_mat = get_heavy_contacts(gapped_res_lst)
            heavy_cutoff = 5
            ref_contact_map = dist_mat < heavy_cutoff
            ref_contacts = np.where(dist_mat < heavy_cutoff)
        else:
            dist_mat = get_cb_contacts(gapped_cb_lst)
            cb_cutoff = 8
            ref_contact_map = dist_mat < cb_cutoff
            ref_contacts = np.where(dist_mat < cb_cutoff)
        
        ref_contacts_x = ref_contacts[0]
        ref_contacts_y = ref_contacts[1]
       
        PPVs = get_ppvs(contacts_x, contacts_y, ref_contact_map, atom_seq_ali, ref_len, factor)
        tp_colors = get_tp_colors(contacts_x, contacts_y, ref_contact_map, atom_seq_ali)
   
        print '%s\t%s' % (acc, PPVs[-1])
      
        ax.scatter(ref_contacts_x, ref_contacts_y, marker='o', c='#CCCCCC', lw=0, edgecolor='#CCCCCC')


    ### plot predicted contacts from second contact map if given
    if c2_filename:
        contacts2 = parse_contacts.parse(open(c2_filename, 'r'), sep)
        contacts2_x = []
        contacts2_y = []
        scores2 = []
        contact_dict2 = {}

        count = 0

        for i in range(len(contacts2)):
            score = contacts2[i][0]
            c_x = contacts2[i][1] - 1
            c_y = contacts2[i][2] - 1

            pos_diff = abs(c_x - c_y)
            too_close = pos_diff < 5

            if not too_close:
                contacts2_x.append(c_x)
                contacts2_y.append(c_y)
                scores2.append(score)
                count += 1
               
            if count >= ref_len * factor:
                break

        ### use TP/FP color coding if reference contacts given
        if pdb_filename:
            PPVs2 = get_ppvs(contacts2_x, contacts2_y, ref_contact_map, atom_seq_ali, ref_len, factor)
            tp2_colors = get_tp_colors(contacts2_x, contacts2_y, ref_contact_map, atom_seq_ali)
            print '%s\t%s' % (acc, PPVs2[-1])
            fig.suptitle('%s\nPPV (upper left) = %.2f | PPV (lower right) = %.2f' % (acc, PPVs[-1], PPVs2[-1]))
            sc = ax.scatter(contacts2_y[::-1], contacts2_x[::-1], marker='o', c=tp2_colors[::-1], s=6, alpha=0.75, linewidths=0.0)
            sc = ax.scatter(contacts_x[::-1], contacts_y[::-1], marker='o', c=tp_colors[::-1], s=6, alpha=0.75, linewidths=0.0)
        else:
            sc = ax.scatter(contacts2_y[::-1], contacts2_x[::-1], marker='o', c='#D70909', edgecolor='#D70909', s=4, linewidths=0.5)
            sc = ax.scatter(contacts_x[::-1], contacts_y[::-1], marker='o', c='#004F9D', edgecolor='#004F9D', s=4, linewidths=0.5)


    ### plot predicted contacts from first contact map on both triangles
    ### if no second contact map given
    else:
        if pdb_filename:
            fig.suptitle('%s\nPPV = %.2f' % (acc, PPVs[-1]))
            sc = ax.scatter(contacts_x[::-1], contacts_y[::-1], marker='o', c=tp_colors[::-1], s=6, alpha=0.75, linewidths=0.0)
            sc = ax.scatter(contacts_y[::-1], contacts_x[::-1], marker='o', c=tp_colors[::-1], s=6, alpha=0.75, linewidths=0.0)
        else:
            sc = ax.scatter(contacts_x[::-1], contacts_y[::-1], marker='o', c=scores[::-1], s=4, alpha=0.75, cmap=cm.jet, linewidths=0.1)
            sc = ax.scatter(contacts_y[::-1], contacts_x[::-1], marker='o', c=scores[::-1], s=4, alpha=0.75, cmap=cm.jet, linewidths=0.1)
            plt.colorbar(sc)

    plt.gca().set_xlim([0,ref_len])
    plt.gca().set_ylim([0,ref_len])

    if outfilename:
        if outfilename.endswith('.pdf'):
            pp = PdfPages(outfilename)
            pp.savefig(fig)
            pp.close()
        else:
            pp = PdfPages('%s.pdf' % outfilename)
            pp.savefig(fig)
            pp.close()
    else:
        pp = PdfPages('%s_ContactMap.pdf' % c_filename)
        pp.savefig(fig)
        pp.close()


    
if __name__ == "__main__":

    p = argparse.ArgumentParser(description='Plot protein residue contact maps.')
    p.add_argument('fasta_file')#, required=True)
    p.add_argument('contact_file')#, required=True)
    p.add_argument('-o', '--outfile', default='')
    p.add_argument('-f', '--factor', default=2.0, type=float)
    p.add_argument('--c2', default='')
    p.add_argument('--psipred_horiz', default='')
    p.add_argument('--pdb', default='')
    p.add_argument('--heavy', action='store_true')
    p.add_argument('--chain', default='')

    args = vars(p.parse_args(sys.argv[1:]))

    fasta_filename = args['fasta_file']
    c_filename = args['contact_file']
    psipred_filename = args['psipred_horiz']

    # guessing separator of constraint file
    line = open(c_filename,'r').readline()
    if len(line.split(',')) != 1:
        sep = ','
    elif len(line.split(' ')) != 1:
        sep = ' '
    else:
        sep = '\t'
    
    plot_map(args['fasta_file'], args['contact_file'], args['factor'], c2_filename=args['c2'], psipred_filename=args['psipred_horiz'], pdb_filename=args['pdb'], is_heavy=args['heavy'], chain=args['chain'], sep=sep, outfilename=args['outfile'])

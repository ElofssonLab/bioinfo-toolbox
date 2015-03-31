import sys, os, re, string
import argparse
from math import *

# on UPPMAX only
sys.path.append('/sw/apps/bioinfo/biopython/1.59/tintin/lib/python')

import Bio.PDB
from Bio import pairwise2

import numpy as np

import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

from os.path import expanduser
home = expanduser("~")
sys.path.append(home + '/bioinfo-toolbox')

from parsing import parse_contacts
from parsing import parse_psipred
from parsing import parse_fasta
from parsing import parse_pdb


def s_score(d, d0):
    return 1/(1+pow(d/d0, 2))

def s_score_vec(d, d0):
    f = np.vectorize(s_score, otypes=[np.float])
    return f(d, d0)


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
    
    #offset = 0
    #first_i = gapped_cb_lst[0].keys()[0]
    #if first_i < 0:
    #    offset = abs(first_i)

    for i, cb1 in enumerate(gapped_cb_lst):
        if cb1 == '-':
            continue
        for j, cb2 in enumerate(gapped_cb_lst):
            if cb2 == '-':
                continue
            diff_vec = cb1 - cb2
            #dist_mat[i+offset,j+offset] = np.sqrt(np.sum(diff_vec * diff_vec))
            dist_mat[i,j] = np.sqrt(np.sum(diff_vec * diff_vec))
    return dist_mat



def get_ppvs(contacts_x, contacts_y, ref_contact_map, atom_seq_ali, ref_len, factor):

    PPVs = []
    TPs = []
    FPs = []

    for num_c in range(min(len(contacts_x), ceil(ref_len * factor)) + 1)[1:]:
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
                TP += 1.0 / (ref_len*factor)
            else:
                FP += 1.0 / (ref_len*factor)

        if TP > 0.0:
            PPVs.append(TP / (TP + FP))
            TPs.append(TP)
            FPs.append(FP)


    if len(PPVs) == 0:
        PPVs.append(0.0)
    if len(TPs) == 0:
        TPs.append(0.0)
    if len(FPs) == 0:
        FPs.append(0.0)

    return PPVs, TPs, FPs


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
 

def get_seqlen(filename):
    alifile = open(filename, 'r')
    l = alifile.readline()
    if l.startswith('>'):
        L = len(alifile.readline().strip())
    else:
        L = len(l.strip())
    alifile.close()
    return L


def get_ali_coverage(filename):
    L = get_seqlen(filename)
    N = 0
    alifile = open(filename, 'r')
    coverage = dict.fromkeys(range(L), 0)
    for line in alifile:
        # ignore headers
        if line.startswith('>'):
            continue
        # remove possible inserts
        line = line.translate(None,string.ascii_lowercase)
        pos_lst = [m.start() for m in re.finditer('-', line)]
        for pos in pos_lst:
            coverage[pos] += 1
        N += 1
    alifile.close()
    coverage_lst = [1-(coverage[i]/float(N)) for i in range(L)]
    return coverage_lst


def plot_map(fasta_filename, c_filename, factor, c2_filename='', psipred_horiz_fname='', psipred_vert_fname='', pdb_filename='', is_heavy=False, chain='', sep=',', outfilename='', ali_filename=''):  
   
    #acc = c_filename.split('.')[0]
    #acc = fasta_filename.split('.')[0][:4]
    acc = '.'.join(os.path.basename(fasta_filename).split('.')[:-1])

    ### get sequence
    seq = parse_fasta.read_fasta(open(fasta_filename, 'r')).values()[0][0]
    ref_len = len(seq)

    ### get top "factor" * "ref_len" predicted contacts
    contacts = parse_contacts.parse(open(c_filename, 'r'), sep)
    contacts_np = parse_contacts.get_numpy_cmap(contacts)

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
    fig = plt.figure(figsize=(8, 8), dpi=96, facecolor='w')
    #ax = fig.add_subplot(111)
    ax = plt.subplot2grid((8,8), (1, 1), colspan=7, rowspan=7)
    ax.tick_params(labelleft='off')
    ax.set_xlim([-3,ref_len])
    ax.set_ylim([-3,ref_len])


    ### plot secondary structure on the diagonal if given
    if psipred_horiz_fname or psipred_vert_fname:
        if psipred_horiz_fname:
            ss = parse_psipred.horizontal(open(psipred_horiz_fname, 'r'))
        else:
            ss = parse_psipred.vertical(open(psipred_vert_fname, 'r'))

        assert len(ss) == ref_len
 
        for i in range(len(ss)):
            if ss[i] == 'H':
                ax.plot(-1.5, i, 'o', c='#8B0043', mec="#8B0043", markersize=2)
                ax.plot(i, -1.5, 'o', c='#8B0043', mec="#8B0043", markersize=2)
            if ss[i] == 'E':
                ax.plot(-1.5, i, 'D', c='#0080AD', mec="#0080AD", markersize=2)
                ax.plot(i, -1.5, 'D', c='#0080AD', mec="#0080AD", markersize=2)
            if ss[i] == 'C':
                continue
        ax2.axis('off')

    
    ### plot alignment coverage if alignemnt given
    if ali_filename:
        coverage_lst = get_ali_coverage(ali_filename)
        ax2 = plt.subplot2grid((8,8), (1,0), rowspan=7, sharey=ax)
        ax2.plot([0]+coverage_lst+[0], [0]+range(ref_len)+[ref_len-1], 'k', lw=0)
        ax2.axvline(x=0.25, lw=0.5, c='black', ls=':')
        ax2.axvline(x=0.5, lw=0.5, c='black', ls=':')
        ax2.axvline(x=0.75, lw=0.5, c='black', ls=':')
        ax2.fill([0]+coverage_lst+[0], [0]+range(ref_len)+[ref_len-1], facecolor='gray', lw=0, alpha=0.5)
        ax2.set_ylim([-3,ref_len])
        ax2.set_xticks([0, 1])
        ax2.invert_xaxis()
        #ax2.spines['top'].set_visible(False)
        #ax2.spines['left'].set_visible(False)
        #ax.get_xaxis().tick_bottom()
        #ax.get_yaxis().tick_right()
        ax2.grid()

        ax3 = plt.subplot2grid((8,8), (0,1), colspan=7, sharex=ax)
        ax3.plot([0]+range(ref_len)+[ref_len-1], [0]+coverage_lst+[0], 'k', lw=0)
        ax3.axhline(y=0.25, lw=0.5, c='black', ls=':')
        ax3.axhline(y=0.5, lw=0.5, c='black', ls=':')
        ax3.axhline(y=0.75, lw=0.5, c='black', ls=':')
        ax3.fill([0]+range(ref_len)+[ref_len-1], [0]+coverage_lst+[0], facecolor='gray', lw=0, alpha=0.5)
        #ax3.xaxis.tick_top()
        ax3.set_xlim([-3,ref_len])
        ax3.set_yticks([0, 1])
        ax3.tick_params(labelbottom='off')
        #ax3.spines['top'].set_visible(False)
        #ax3.spines['right'].set_visible(False)
        #ax.get_xaxis().tick_top()
        #ax.get_yaxis().tick_left()
        ax3.grid()


    ### plot reference contacts in the background if given
    if pdb_filename:
        res_lst = parse_pdb.get_coordinates(open(pdb_filename, 'r'), chain)
        cb_lst = parse_pdb.get_cb_coordinates(open(pdb_filename, 'r'), chain)
        atom_seq = parse_pdb.get_atom_seq(open(pdb_filename, 'r'), chain)
                
        align = pairwise2.align.globalms(atom_seq, seq, 2, -1, -0.5, -0.1)

        atom_seq_ali = align[-1][0]
        seq_ali = align[-1][1]

        #print atom_seq_ali
        #print seq_ali

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
            #ref_contacts = np.where(np.ma.array(dist_mat, mask=np.tri(dist_mat.shape[0]), fill_value=float("inf")) < cb_cutoff)
        
        ref_contacts_x = ref_contacts[0]
        ref_contacts_y = ref_contacts[1]

        PPVs, TPs, FPs = get_ppvs(contacts_x, contacts_y, ref_contact_map, atom_seq_ali, ref_len, factor)
        tp_colors = get_tp_colors(contacts_x, contacts_y, ref_contact_map, atom_seq_ali)
   
        print '%s %s %s %s' % (pdb_filename, PPVs[-1], TPs[-1], FPs[-1])
      
        cmap = cm.get_cmap("binary")
        cmap.set_bad([1,1,1,0])
        dist_mat_masked = np.ma.array(dist_mat, mask=np.tri(dist_mat.shape[0], k=-1))
        sc = ax.imshow(s_score_vec(dist_mat_masked, 5), cmap=cmap, interpolation='none')
        
        ref_contacts_diag_x = []
        ref_contacts_diag_y = []
        for i in range(len(ref_contacts_x)):
            x_i = ref_contacts_x[i]
            y_i = ref_contacts_y[i]
            if not dist_mat_masked.mask[x_i, y_i]:
                ref_contacts_diag_x.append(x_i)
                ref_contacts_diag_y.append(y_i)
       
        ax.scatter(ref_contacts_diag_x, ref_contacts_diag_y, marker='o', c='#CCCCCC', lw=0, edgecolor='#CCCCCC')


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
            PPVs2, TPs2, FPs2 = get_ppvs(contacts2_x, contacts2_y, ref_contact_map, atom_seq_ali, ref_len, factor)
            tp2_colors = get_tp_colors(contacts2_x, contacts2_y, ref_contact_map, atom_seq_ali)
            print '%s %s %s %s' % (pdb_filename, PPVs2[-1], TPs2[-1], FPs2[-1])
            fig.suptitle('%s\nPPV (upper left) = %.2f | PPV (lower right) = %.2f' % (acc, PPVs[-1], PPVs2[-1]))
            sc = ax.scatter(contacts2_y[::-1], contacts2_x[::-1], marker='o', c=tp2_colors[::-1], s=6, alpha=0.75, lw=0)
            sc = ax.scatter(contacts_x[::-1], contacts_y[::-1], marker='o', c=tp_colors[::-1], s=6, alpha=0.75, lw=0)
        else:
            sc = ax.scatter(contacts2_y[::-1], contacts2_x[::-1], marker='o', c='#D70909', edgecolor='#D70909', s=4, linewidths=0.5)
            sc = ax.scatter(contacts_x[::-1], contacts_y[::-1], marker='o', c='#004F9D', edgecolor='#004F9D', s=4, linewidths=0.5)


    ### plot predicted contacts from first contact map on both triangles
    ### if no second contact map given
    else:
        if pdb_filename:
            pdb_acc = parse_pdb.get_acc(open(pdb_filename))
            if pdb_acc:
                if chain:
                    fig.suptitle('%s (PDB: %s, chain %s)\nPPV = %.2f' % (acc, pdb_acc, chain, PPVs[-1]))
                else:
                    fig.suptitle('%s (PDB: %s)\nPPV = %.2f' % (acc, pdb_acc, PPVs[-1]))
            else:
                fig.suptitle('%s\nPPV = %.2f' % (acc, PPVs[-1]))
            cmap = cm.get_cmap("hot_r")
            cmap.set_bad([1,1,1,0])
            contacts_np_masked = np.ma.array(contacts_np, mask=np.tri(contacts_np.shape[0], k=-1))
            #sc = ax.imshow(contacts_np_masked.T, cmap=cmap)
            sc = ax.scatter(contacts_x[::-1], contacts_y[::-1], marker='o', c=tp_colors[::-1], s=6, alpha=0.75, linewidths=0.0)
            #sc = ax.scatter(contacts_y[::-1], contacts_x[::-1], marker='o', c=tp_colors[::-1], s=6, alpha=0.75, linewidths=0.0)
        else:
            #if c_filename.startswith('data'):
            #    acc = c_filename.split('/')[1]
            #else:
            #    acc = c_filename.split('/')[-1]
            fig.suptitle('%s' % acc)
            #sc = ax.imshow(contacts_np + contacts_np.T, cmap=cm.hot_r)
            sc = ax.imshow(contacts_np, cmap=cm.hot_r, vmin=0.2, vmax=1.0, interpolation='none')
            #divider1 = make_axes_locatable(ax)
            #cax1 = divider1.append_axes("right", size="2%", pad=0.05)
            #plt.colorbar(sc, cax=cax1)
            #plt.colorbar(sc, ax=ax)
            sc = ax.scatter(contacts_x[::-1], contacts_y[::-1],
                    marker='o', c="black", s=4, alpha=0.75,
                    linewidths=0.1, edgecolors='none')
            #sc = ax.scatter(contacts_y[::-1], contacts_x[::-1], marker='o', c=scores[::-1], s=4, alpha=0.75, cmap=cm.hot_r, linewidths=0.1, edgecolors='none')

    #plt.gca().set_xlim([0,ref_len])
    #plt.gca().set_ylim([0,ref_len])

    ax.grid()
    #ax.invert_yaxis()

    if outfilename:
        if outfilename.endswith('.pdf'):
            pp = PdfPages(outfilename)
            pp.savefig(fig)
            pp.close()
        elif outfilename.endswith('.png'):
            plt.savefig(outfilename)
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
    p.add_argument('--psipred_vert', default='')
    p.add_argument('--pdb', default='')
    p.add_argument('--heavy', action='store_true')
    p.add_argument('--chain', default='')
    p.add_argument('--alignment', default='')

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
    
    plot_map(args['fasta_file'], args['contact_file'], args['factor'], c2_filename=args['c2'], psipred_horiz_fname=args['psipred_horiz'], psipred_vert_fname=args['psipred_vert'], pdb_filename=args['pdb'], is_heavy=args['heavy'], chain=args['chain'], sep=sep, outfilename=args['outfile'], ali_filename=args['alignment'])


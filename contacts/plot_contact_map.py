#!/usr/bin/env python
import sys, os, re, string
import argparse
from math import *

# on UPPMAX only
sys.path.append('/sw/apps/bioinfo/biopython/1.59/tintin/lib/python')

from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

import numpy as np

import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.cm as cm
#from mpl_toolkits.axes_grid1 import make_axes_locatable

from os.path import expanduser
home = expanduser("~")
sys.path.append(home + '/bioinfo-toolbox/parsing')
sys.path.append(home + '/git/bioinfo-toolbox/parsing')

import parse_contacts
import parse_psipred
import parse_fasta
import parse_pdb
import parse_hhblits_hhr
import parse_a3m


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



def get_ppvs(contacts_x, contacts_y, ref_contact_map, ref_len, factor, atom_seq_ali=[]):

    PPVs = []
    TPs = []
    FPs = []

    for num_c in range(min(len(contacts_x), int(ceil(ref_len * factor))) + 1)[1:]:
        TP = 0.0
        FP = 0.0
        for i in range(num_c):
            c_x = contacts_x[i]
            c_y = contacts_y[i]
            if atom_seq_ali:
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


def get_tp_colors(contacts_x, contacts_y, ref_contact_map, atom_seq_ali=[]):

    tp_colors = []

    for i in range(len(contacts_x)):
        c_x = contacts_x[i]
        c_y = contacts_y[i]
        if atom_seq_ali:
            if atom_seq_ali[c_x] == '-':
                #tp_colors.append('green')
                tp_colors.append('red')
                continue
            if atom_seq_ali[c_y] == '-':
                #tp_colors.append('green')
                tp_colors.append('red')
                continue
        if ref_contact_map[c_x, c_y] > 0:
            tp_colors.append('#00ff00')
        else:
            tp_colors.append('red')

    return tp_colors
 

def get_colors(contacts_np, ref_contact_map=[], th=0.5, binary=False):

    N = contacts_np.shape[0]
    img = np.ones((N,N,4))

    for i in xrange(N):
        for j in xrange(N):
            if j < (i+5):
                continue
            sc = contacts_np[i,j]
            if len(ref_contact_map) > 0:
                #print N, ref_contact_map.shape[0]
                assert N == ref_contact_map.shape[0]
                # FN
                if sc <= th and ref_contact_map[i,j] < 8:
                    #img[i,j] = [0.5,0.5,1,1]
                    img[i,j] = [0.5,0.5,0.5,1]
                    img[j,i] = [0.5,0.5,0.5,1]
                # TP
                elif sc > th and ref_contact_map[i,j] < 8:
                    img[i,j] = [0,1,0,1]
                    img[j,i] = [0,1,0,1]
                    #img[j,i] = [1-sc,1-sc,1-sc,1]
                # FP
                #elif contacts_np[i,j] > th and ref_contact_map[i,j] >= 8:
                elif sc > th and ref_contact_map[i,j] >= 12:
                    img[i,j] = [1,0,0,1]
                    img[j,i] = [1,0,0,1]
                    #img[j,i] = [1-sc,1-sc,1-sc,1]
                # grey zone between 8 and 12 Angstroem
                elif sc > th and (ref_contact_map[i,j] < 12 or ref_contact_map[i,j] >= 8):
                    if binary:
                        img[i,j] = [1,0,0,1]
                        img[j,i] = [1,0,0,1]
                    else:
                        val = (ref_contact_map[i,j] - 8)/(12 - 8)
                        img[i,j] = [0.5+val/2,1-val/2,0,1]
                        img[j,i] = [0.5+val/2,1-val/2,0,1]
                        #img[j,i] = [1-sc,1-sc,1-sc,1]
            else:
                if sc > th:
                    img[i,j] = [0.5-sc/2,0.5-sc/2,1,1]
                    img[j,i] = [0.5-sc/2,0.5-sc/2,1,1]

    return img



def get_ref_img(ref_contact_map):

    N = ref_contact_map.shape[0]
    img = np.ones((N,N,4))

    for i in xrange(N):
        for j in xrange(N):
            if j < (i+5):
                continue
            if len(ref_contact_map) > 0:
                # FN
                if ref_contact_map[i,j] < 8:
                    #img[i,j] = [0.5,0.5,1,1]
                    img[i,j] = [0.5,0.5,0.5,1]
                    img[j,i] = [0.5,0.5,0.5,1]
    return img


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
    M = 0.
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
        M += 1.
    alifile.close()
    #coverage_lst = [1-(coverage[i]/float(N)) for i in range(L)]
    coverage_lst = [M - coverage[i] for i in range(L)]
    return coverage_lst, M


def get_meff_coverage(meff_file):
    meff_lst = []
    with open(meff_file) as f:
        for l in f:
            if 'Meff = ' in l:
                meff = float(l.strip().split()[-1])
            if 'MeffPerPos' in l:
                meff_lst = l.split()[-1].strip('[]').split(',')
                meff_lst = map(float, meff_lst)
    return meff_lst, meff


def is_gap(symbol):
    return symbol == '-'


def embedd_alignment(target_0, target_1, source_0, source_1):
    """ Embedds pairwise alignment (source) in another pairwise alignment (target), 
        where unaligned target_1 is identical to unaligned source_0.
    """
    #print target_1
    #print source_0
    assert target_1.replace('-','') == source_0.replace('-','')

    result_0 = ''
    result_1 = ''

    i = 0
    j = 0
    #sys.stdout.write('\n')
    while i < len(target_0):
        t_0i = target_0[i] 
        t_1i = target_1[i]
        
        if not is_gap(t_0i) and is_gap(t_1i):
            #sys.stdout.write('%s' % t_0i)
            #sys.stdout.write('%s' % '-')
            result_0 += t_0i
            result_1 += '-'
            i += 1
        elif not is_gap(t_0i) and not is_gap(t_1i):
            s_0j = source_0[j]
            s_1j = source_1[j]
            #sys.stdout.write('%s' % s_0j)
            #sys.stdout.write('%s' % s_1j)
            result_0 += s_0j
            result_1 += s_1j
            if not is_gap(s_0j):
                i += 1
            j += 1
        elif is_gap(t_0i):
            s_0j = source_0[j]
            s_1j = source_1[j]
            #sys.stdout.write('%s' % '-')
            #sys.stdout.write('%s' % s_1j)
            result_0 += '-'
            result_1 += s_1j 
            if not is_gap(s_0j):
                i += 1
            j += 1
    
    while j < len(source_0):
        result_0 += '-'
        result_1 += source_1[j] 
        j += 1

    #print result_0
    #print result_1
    #print result_1.replace('-','')
    #print source_1.replace('-','')
    assert result_0.replace('-','') == target_0.replace('-','')
    assert result_1.replace('-','') == source_1.replace('-','')

    return result_0, result_1





def plot_map(fasta_filename, c_filename, factor=1.0, th=-1, f_obs=-1, c2_filename='', psipred_horiz_fname='', psipred_vert_fname='', pdb_filename='', is_heavy=False, chain='', sep=',', outfilename='', ali_filename='',  meff_filename='', name='', start=0, end=-1, pdb_start=0, pdb_end=-1, noalign=False, pdb_alignment='', pdb_id='', binary=False):
  
    #acc = c_filename.split('.')[0]
    #acc = fasta_filename.split('.')[0][:4]
    if name == '': 
        acc = '.'.join(os.path.basename(fasta_filename).split('.')[:-1])
    else:
        acc = name

    ### get sequence
    seq = parse_fasta.read_fasta(open(fasta_filename, 'r')).values()[0][0]
    ref_len = len(seq)

    ### trim sequence according to given positions
    ### default: take full sequence
    if end == -1:
        end = ref_len
    seq = seq[start:end]
    ref_len = len(seq)
    unit = (ref_len/50.0)

    ### get top "factor" * "ref_len" predicted contacts
    contacts = parse_contacts.parse(open(c_filename, 'r'), sep)
    contacts_np = parse_contacts.get_numpy_cmap(contacts, seq_len=ref_len)
    contacts_np = contacts_np[start:end,start:end]

    contacts_x = []
    contacts_y = []
    scores = []
    contact_dict = {}

    count = 0
    for i in range(len(contacts)):
        score = contacts[i][0]
        c_x = contacts[i][1] - 1
        c_y = contacts[i][2] - 1
        
        # only look at contacts within given range
        # default: take full sequence range into account
        if c_x < start or c_x >= end:
            continue
        if c_y < start or c_y >= end:
            continue
        
        pos_diff = abs(c_x - c_y)
        too_close = pos_diff < 5

        if not too_close:
            contacts_x.append(c_x - start)
            contacts_y.append(c_y - start)
            scores.append(score)
            count += 1
           
        if count >= ref_len * factor and th == -1 and f_obs == -1:
            th = score
            break
        
        if score < th and not th == -1 and f_obs == -1:
            factor = count/float(ref_len)
            break

        # if cutoff by fraction of observed contacts:
        # take all contacts and cut list after reading pdb
        #if f_obs != -1:
        #    factor = ref_len
        #    th = -1

 

    ### start plotting
    fig = plt.figure(figsize=(8, 8), dpi=96, facecolor='w')
    ax = fig.add_subplot(111)#, aspect='auto')
    ax.set_adjustable('box-forced')
    ax.tick_params(direction='out', right='off', top='off')
    ax.set_xlim([-unit,ref_len])
    ax.set_ylim([-unit,ref_len])
    
    ### plot alignment coverage if alignemnt given
    if ali_filename or meff_filename:
        # adjust overall canvas  
        ax = plt.subplot2grid((8,8), (1, 1), colspan=7, rowspan=7)#, aspect='auto')
        #ax.set_adjustable('box-forced')
        #ax.set_autoscale_on(False) 
        ax.autoscale(False)
        ax.tick_params(direction='out',labelleft='off', right='off', top='off')
        ax.set_xlim([-unit,ref_len])
        ax.set_ylim([-unit,ref_len])

        if ali_filename:
            coverage_lst, M = get_ali_coverage(ali_filename)
            max_cover = M
        elif meff_filename:
            coverage_lst, Meff = get_meff_coverage(meff_filename)
            max_cover = Meff
        #max_cover = max(coverage_lst)

        coverage_lst = coverage_lst[start:end]
        
        #lt = pow(10, max(1,floor(log10(max_cover)) - 1))
        #upper = int(ceil(max_cover/float(lt)) * lt)
        ax2 = plt.subplot2grid((8,8), (1,0), rowspan=7, sharey=ax)
        #ax2.set_adjustable('box-forced')
        #ax2.set_autoscale_on(False) 
        ax2.autoscale(False)
        #print len([0]+coverage_lst+[0])
        #print len([0]+range(ref_len)+[ref_len-1])

        ax2.plot([0]+coverage_lst+[0], [0]+range(ref_len)+[ref_len-1], 'k', lw=0)
        ax2.axvline(x=max_cover*0.25, lw=0.5, c='black', ls=':')
        ax2.axvline(x=max_cover*0.5, lw=0.5, c='black', ls=':')
        ax2.axvline(x=max_cover*0.75, lw=0.5, c='black', ls=':')
        ax2.fill([0]+coverage_lst+[0], [0]+range(ref_len)+[ref_len-1], facecolor='gray', lw=0, alpha=0.5)
        ax2.set_xticks([0, max_cover])
        ax2.tick_params(axis='x', top='off', direction='out')
        ax2.invert_xaxis()
        #ax2.spines['top'].set_visible(False)
        #ax2.spines['left'].set_visible(False)
        #ax.get_xaxis().tick_bottom()
        #ax.get_yaxis().tick_right()
        ax2.grid()
        ax2.set_ylim([-unit,ref_len])

        ax3 = plt.subplot2grid((8,8), (0,1), colspan=7, sharex=ax)
        #ax3.set_adjustable('box-forced')
        #ax3.set_autoscale_on(False) 
        ax3.autoscale(False)
        ax3.plot([0]+range(ref_len)+[ref_len-1], [0]+coverage_lst+[0], 'k', lw=0)
        ax3.axhline(y=max_cover*0.25, lw=0.5, c='black', ls=':')
        ax3.axhline(y=max_cover*0.5, lw=0.5, c='black', ls=':')
        ax3.axhline(y=max_cover*0.75, lw=0.5, c='black', ls=':')
        ax3.fill([0]+range(ref_len)+[ref_len-1], [0]+coverage_lst+[0], facecolor='gray', lw=0, alpha=0.5)
        #ax3.xaxis.tick_top()
        ax3.set_yticks([0, max_cover])
        ax3.tick_params(labelbottom='off')
        ax2.tick_params(axis='y', right='off', direction='out', left='on')
        #ax3.spines['top'].set_visible(False)
        #ax3.spines['right'].set_visible(False)
        #ax.get_xaxis().tick_top()
        #ax.get_yaxis().tick_left()
        ax3.grid()
        ax3.set_xlim([-unit,ref_len])


    ### plot secondary structure along axis if given
    if psipred_horiz_fname or psipred_vert_fname:
        if psipred_horiz_fname:
            ss = parse_psipred.horizontal(open(psipred_horiz_fname, 'r'))
        else:
            ss = parse_psipred.vertical(open(psipred_vert_fname, 'r'))

        ss = ss[start:end]
        assert len(ss) == ref_len
 
        ax.axhline(y=0, lw=1, c='black')
        ax.axvline(x=0, lw=1, c='black')
        for i in range(len(ss)):
            if ss[i] == 'H':
                #ax.plot(-unit/2, i, 's', c='#8B0043', mec="#8B0043")#, markersize=2)
                #ax.plot(i, -unit/2, 's', c='#8B0043', mec="#8B0043")#, markersize=2)
                #ax.plot(i, -unit/2, 's', c='#8B0043', mec="#8B0043")#, markersize=2)
                ax.add_patch(plt.Rectangle((-unit, i-0.5), unit, 1, edgecolor='#8B0043', facecolor="#8B0043"))
                ax.add_patch(plt.Rectangle((i-0.5, -unit), 1, unit, edgecolor='#8B0043', facecolor="#8B0043"))
            if ss[i] == 'E':
                ax.add_patch(plt.Rectangle((-unit, i-0.5), unit, 1, edgecolor='#0080AD', facecolor="#0080AD"))
                ax.add_patch(plt.Rectangle((i-0.5, -unit), 1, unit, edgecolor='#0080AD', facecolor="#0080AD"))
                #ax.plot(-unit/2, i, 's', c='#0080AD', mec="#0080AD")#, markersize=2)
                #ax.plot(i, -unit/2, 's', c='#0080AD', mec="#0080AD")#, markersize=2)
            if ss[i] == 'C':
                continue


    ### plot reference contacts in the background if given
    if pdb_filename:
        res_lst = parse_pdb.get_coordinates(open(pdb_filename, 'r'), chain)
        cb_lst = parse_pdb.get_cb_coordinates(open(pdb_filename, 'r'), chain)
        atom_seq = parse_pdb.get_atom_seq(open(pdb_filename, 'r'), chain)

        # if not true there is some serious problem with the provided pdb file
        assert len(res_lst) == len(cb_lst) == len(atom_seq)

        ### trim PDB sequence according to given positions
        ### default: take full sequence
        if pdb_end == -1:
            pdb_end = len(res_lst)
        res_lst = res_lst[pdb_start:pdb_end]
        cb_lst = cb_lst[pdb_start:pdb_end]
        atom_seq = atom_seq[pdb_start:pdb_end]

        #print atom_seq
        #print seq

        if noalign:
            dist_mat = get_cb_contacts(cb_lst)
            cb_cutoff = 8
            ref_contact_map = dist_mat < cb_cutoff
            PPV, TP, FP = get_ppvs(contacts_x, contacts_y, ref_contact_map, ref_len, factor)
            tp_colors = get_tp_colors(contacts_x, contacts_y, ref_contact_map)

        else:
            if pdb_alignment and pdb_id:
                #align = parse_hhblits_hhr.parse_alignments(pdb_alignment)[pdb_id]
                #atom_seq_ali = align[0][0]
                #seq_ali = align[0][1]
                seq_ali, atom_seq_ali = parse_a3m.get_pairwise(pdb_alignment, pdb_id)
                seqres_seq = atom_seq_ali.replace('-', '')
                #print seqres_seq
                
                #print atom_seq_ali
                #print seq_ali
                #print ""
                align_seqres = pairwise2.align.globalms(atom_seq, seqres_seq, 2, -1, -0.5, -0.1)
                atom_seq_ali0 = align_seqres[-1][0]
                seqres_seq_ali0 = align_seqres[-1][1]
                #print ""
                #print atom_seq_ali0
                #print seqres_seq_ali0
                #print ""
                atom_seq_ali, seq_ali = embedd_alignment(atom_seq_ali0, seqres_seq_ali0, atom_seq_ali, seq_ali)
                #print atom_seq_ali1
                #print seq_ali1
                #print ""
            else:
                matrix = matlist.blosum62
                #align = pairwise2.align.globalms(atom_seq, seq, 2, -1, -0.5, -0.1)
                #align = pairwise2.align.localds(atom_seq, seq, matrix, -11, -1)
                align = pairwise2.align.globalds(atom_seq, seq, matrix, -25, -1)
                atom_seq_ali = align[-1][0]
                seq_ali = align[-1][1]

            #print atom_seq_ali
            #print seq_ali
            #print len(atom_seq), len(seq), len(res_lst), len(cb_lst)
            #print len(atom_seq_ali), len(seq_ali)

            j = 0
            gapped_res_lst = []
            gapped_cb_lst = []

            for i in xrange(len(atom_seq_ali)):
                if atom_seq_ali[i] == '-':
                    if seq_ali[i] == '-':
                        continue
                    gapped_res_lst.append('-')
                    gapped_cb_lst.append('-')
                elif seq_ali[i] == '-':
                    j += 1
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
            
            # if f_obs given, take top f_obs * num_obs contacts:
            if f_obs != -1:
                num_obs = sum(ref_contacts_x-ref_contacts_y >= 5)
                num_top = int(ceil(f_obs * num_obs))
                contacts_x = contacts_x[:num_top]
                contacts_y = contacts_y[:num_top]
                scores = scores[:num_top]
                th_obs = scores[-1]
            else:
                th_obs = th

            PPVs, TPs, FPs = get_ppvs(contacts_x, contacts_y, ref_contact_map, ref_len, factor, atom_seq_ali=atom_seq_ali)
            tp_colors = get_tp_colors(contacts_x, contacts_y, ref_contact_map, atom_seq_ali=atom_seq_ali)

        if not c2_filename:
            img = get_colors(contacts_np, ref_contact_map=dist_mat, th=th_obs, binary=binary)
            sc = ax.imshow(img, interpolation='none')
        else:
            # plot native contacts in background
            img = get_ref_img(dist_mat)
            sc = ax.imshow(img, interpolation='none')
            
   
        print '%s %s %s %s' % (acc, PPVs[-1], TPs[-1], FPs[-1])
      
        cmap = cm.get_cmap("binary")
        cmap.set_bad([1,1,1,0])
        dist_mat_masked = np.ma.array(dist_mat, mask=np.tri(dist_mat.shape[0], k=-1))
        #sc = ax.imshow(s_score_vec(dist_mat_masked, 5), cmap=cmap, interpolation='none')
        
        ref_contacts_diag_x = []
        ref_contacts_diag_y = []
        for i in range(len(ref_contacts_x)):
            x_i = ref_contacts_x[i]
            y_i = ref_contacts_y[i]
            if not dist_mat_masked.mask[x_i, y_i] and abs(x_i - y_i) >= 5:
                ref_contacts_diag_x.append(x_i)
                ref_contacts_diag_y.append(y_i)
       
        #ax.scatter(ref_contacts_diag_x, ref_contacts_diag_y, marker='+', c='#000000')


    ### plot predicted contacts from second contact map if given
    if c2_filename:
        contacts2 = parse_contacts.parse(open(c2_filename, 'r'))
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
                contacts2_x.append(c_x - start)
                contacts2_y.append(c_y - start)
                scores2.append(score)
                count += 1
               
            if count >= ref_len * factor and th == -1 and f_obs == -1:
                th = score
                break
            
            if score < th and not th == -1 and  f_obs == -1:
                factor = count/float(ref_len)
                break

            # if cutoff by fraction of observed contacts:
            # stop at num_top = f_obs * num_obs 
            if pdb_filename and f_obs != -1:
                if count >= num_top:
                    factor = count/float(ref_len)
                    th = score
                    break

        ### use TP/FP color coding if reference contacts given
        if pdb_filename:
            #PPVs2, TPs2, FPs2 = get_ppvs(contacts2_x, contacts2_y, ref_contact_map, atom_seq_ali, ref_len, factor)
            tp2_colors = get_tp_colors(contacts2_x, contacts2_y, ref_contact_map, atom_seq_ali)
            #print '%s %s %s %s' % (acc, PPVs2[-1], TPs2[-1], FPs2[-1])
            #fig.suptitle('%s\nPPV (upper left) = %.2f | PPV (lower right) = %.2f' % (acc, PPVs[-1], PPVs2[-1]))
            sc = ax.scatter(contacts2_y[::-1], contacts2_x[::-1], marker='s', c=tp2_colors[::-1], s=4, alpha=1, lw=0, edgecolor=tp2_colors[::-1])
            sc = ax.scatter(contacts_x[::-1], contacts_y[::-1], marker='s', c=tp_colors[::-1], s=4, alpha=1, lw=0, edgecolor=tp_colors[::-1])
        else:
            sc = ax.scatter(contacts2_y[::-1], contacts2_x[::-1], marker='0', c='#D70909', edgecolor='#D70909', s=6, linewidths=0.5)
            sc = ax.scatter(contacts_x[::-1], contacts_y[::-1], marker='o', c='#004F9D', edgecolor='#004F9D', s=6, linewidths=0.5)


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
            #cmap = cm.get_cmap("binary")
            #cmap.set_bad([1,1,1,0])
            #contacts_np_masked = np.ma.array(contacts_np, mask=np.tri(contacts_np.shape[0], k=-1))
            #sc = ax.imshow(contacts_np_masked.T, cmap=cmap)
            #sc = ax.imshow(contacts_np, cmap=cmap)
            #sc = ax.imshow(contacts_np + contacts_np.T, cmap=cm.binary, vmin=0.2, vmax=1.0, interpolation='none')
            #sc = ax.scatter(contacts_x[::-1], contacts_y[::-1], marker='o', c=tp_colors[::-1], s=6, alpha=0.75, linewidths=0.0)
            #sc = ax.scatter(contacts_y[::-1], contacts_x[::-1], marker='o', c=tp_colors[::-1], s=6, alpha=0.75, linewidths=0.0)
        else:
            #if c_filename.startswith('data'):
            #    acc = c_filename.split('/')[1]
            #else:
            #    acc = c_filename.split('/')[-1]
            fig.suptitle('%s' % acc)
            #sc = ax.imshow(contacts_np + contacts_np.T, cmap=cm.hot_r)
            #sc = ax.imshow(contacts_np + contacts_np.T,
            #        cmap=cm.binary, vmin=th, vmax=1.0, interpolation='none')
            img = get_colors(contacts_np, th=th)
            sc = ax.imshow(img, interpolation='none')
            #divider1 = make_axes_locatable(ax)
            #cax1 = divider1.append_axes("right", size="2%", pad=0.05)
            #plt.colorbar(sc, cax=cax1)
            #plt.colorbar(sc, ax=ax)
            #sc = ax.scatter(contacts_x[::-1], contacts_y[::-1],
            #        marker='o', c="black", s=6, alpha=0.75,
            #        linewidths=0.1, edgecolors='none')
            #sc = ax.scatter(contacts_y[::-1], contacts_x[::-1], marker='o', c=scores[::-1], s=4, alpha=0.75, cmap=cm.hot_r, linewidths=0.1, edgecolors='none')

    #plt.gca().set_xlim([0,ref_len])
    #plt.gca().set_ylim([0,ref_len])

    ax.grid()
    ax.set_xlim([-unit,ref_len])
    ax.set_ylim([-unit,ref_len])
    #print ax.axis()
    ax.axis([-unit,ref_len, -unit,ref_len])
    #ax.invert_yaxis()
    ax.set_autoscale_on(False) 

    if outfilename:
        if outfilename.endswith('.pdf'):
            pp = PdfPages(outfilename)
            pp.savefig(fig)
            pp.close()
        elif outfilename.endswith('.eps'):
            plt.savefig(outfilename, format='eps', dpi=300)
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
    p.add_argument('-f', '--factor', default=1.0, type=float)
    p.add_argument('-t', '--threshold', default=-1, type=float)
    p.add_argument('--f_obs', default=-1, type=float)
    p.add_argument('--c2', default='')
    p.add_argument('--psipred_horiz', default='')
    p.add_argument('--psipred_vert', default='')
    p.add_argument('--pdb', default='')
    p.add_argument('--heavy', action='store_true')
    p.add_argument('--chain', default='')
    p.add_argument('--alignment', default='')
    p.add_argument('--meff', default='')
    p.add_argument('--name', default='')
    p.add_argument('--start', default=0, type=int)
    p.add_argument('--end', default=-1, type=int)
    p.add_argument('--pdb_start', default=0, type=int)
    p.add_argument('--pdb_end', default=-1, type=int)
    p.add_argument('--noalign', action='store_true')
    p.add_argument('--binary', action='store_true')
    p.add_argument('--pdb_alignment', default='')
    p.add_argument('--pdb_id', default='')

    args = vars(p.parse_args(sys.argv[1:]))

    fasta_filename = args['fasta_file']
    c_filename = args['contact_file']
    psipred_filename = args['psipred_horiz']

    # guessing separator of constraint file
    for line in open(c_filename,'r'):
        # exclude comments
        if line.startswith('#'):
            continue
        # exclude PhyCmap header/tail lines
        if line.strip()[0].isalpha():
            continue
        if len(line.split(',')) != 1:
            sep = ','
        elif len(line.split('\t')) != 1:
            sep = '\t'
        else:
            sep = ' '
    
    plot_map(args['fasta_file'], args['contact_file'], factor=args['factor'], th=args['threshold'], f_obs=args['f_obs'], c2_filename=args['c2'], psipred_horiz_fname=args['psipred_horiz'], psipred_vert_fname=args['psipred_vert'], pdb_filename=args['pdb'], is_heavy=args['heavy'], chain=args['chain'], sep=sep, outfilename=args['outfile'], ali_filename=args['alignment'], meff_filename=args['meff'], name=args['name'], start=args['start'], end=args['end'], pdb_start=args['pdb_start'], pdb_end=args['pdb_end'], noalign=args['noalign'], pdb_alignment=args['pdb_alignment'], pdb_id=args['pdb_id'], binary=args['binary'])


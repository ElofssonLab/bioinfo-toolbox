#!/usr/bin/env python
from __future__ import division
import sys, os, re, string
import argparse
from math import *

# on UPPMAX only
sys.path.append('/sw/apps/bioinfo/biopython/1.59/tintin/lib/python')

from Bio import pairwise2

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
import parse_iupred
import parse_pdb


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
        if str(cb1) == str('-'):
            continue
        for j, cb2 in enumerate(gapped_cb_lst):
            if str(cb2) == str('-'):
                continue
            diff_vec = cb1 - cb2
            #dist_mat[i+offset,j+offset] = np.sqrt(np.sum(diff_vec * diff_vec))
            dist_mat[i,j] = np.sqrt(np.sum(diff_vec * diff_vec))
    return dist_mat



def print_contacts(fasta_filename,score,contacts_x, contacts_y, ref_contact_map, atom_seq_ali, ref_len, factor,disorder):
    #    for num_c in range(min(len(contacts_x), int(ceil(ref_len * factor))) + 1)[1:]:
    TP = 0.0
    FP = 0.0
    disoTP = 0.0
    disoFP = 0.0
    mixTP = 0.0
    mixFP = 0.0
    for i in range(len(contacts_x) ):
        c_x = contacts_x[i]
        c_y = contacts_y[i]
        if atom_seq_ali[c_x] == '-':
            continue
        if atom_seq_ali[c_y] == '-':
            continue
        if len(disorder)>i:
            print "DATA: ",fasta_filename,i,scores[i],c_x,c_y,disorder[c_x],disorder[c_y]
        else:
            print "DATA: ",fasta_filename,i,scores[i],c_x,c_y



def get_ppvs(contacts_x, contacts_y, ref_contact_map, atom_seq_ali, ref_len, factor,disorder):

    PPVs = []
    TPs = []
    FPs = []

    TP=0
    FP=0

    orderPPVs = []
    orderTPs = []
    orderFPs = []
    disoPPVs = []
    disoTPs = []
    disoFPs = []
    mixPPVs = []
    mixTPs = []
    mixFPs = []
    longcount = 1.e-20
    count = 1.e-20
    disocount = 1.e-20
    mixcount = 1.e-20
    ordercount = 1.e-20
    #    for num_c in range(min(len(contacts_x), int(ceil(ref_len * factor))) + 1)[1:]:
    TP = 0.0
    FP = 0.0
    disoTP = 0.0
    disoFP = 0.0
    orderTP = 0.0
    orderFP = 0.0
    mixTP = 0.0
    mixFP = 0.0

    
#    for num_c in range(min(len(contacts_x), int(ceil(ref_len * factor))) + 1)[1:]:
#        TP = 0.0
#        FP = 0.0
#        for i in range(num_c):
#            c_x = contacts_x[i]
#            c_y = contacts_y[i]
#            if atom_seq_ali[c_x] == '-':
#                continue
#            if atom_seq_ali[c_y] == '-':
#                continue
#            if ref_contact_map[c_x, c_y] > 0:
#                TP += 1.0 / (ref_len*factor)
#            else:
#                FP += 1.0 / (ref_len*factor)
#        
#
#        if TP > 0.0:
#            PPVs.append(TP / (TP + FP))
#            TPs.append(TP)
#            FPs.append(FP)

#    for num_c in range(min(len(contacts_x), int(ceil(ref_len * factor))) ):
#    for num_c in xrange(len(atom_seq_ali)):
    for num_c in xrange(len(contacts_x)):
        c_x = contacts_x[num_c]
        c_y = contacts_y[num_c]
        if atom_seq_ali[c_x] == '-':
            continue
        if atom_seq_ali[c_y] == '-':
            continue
        if (disorder[c_x] > 0.5 and disorder[c_y] > 0.5):
            if (disocount < ref_len * factor):
                disocount+=1
            if ref_contact_map[c_x, c_y] > 0:
                disoTP += 1.0 
            else:
                disoFP += 1.0 
            disoPPVs.append(disoTP / (disoTP + disoFP))
            disoTPs.append(disoTP/disocount)
            disoFPs.append(disoFP/disocount)

        elif (disorder[c_x] > 0.5 or disorder[c_y] > 0.5):
            if (mixcount < ref_len * factor):
                mixcount+=1
            if ref_contact_map[c_x, c_y] > 0:
                mixTP += 1.0 
            else:
                mixFP += 1.0 
            mixPPVs.append(mixTP / (mixTP + mixFP))
            mixTPs.append(mixTP/mixcount)
            mixFPs.append(mixFP/mixcount)
        else:
            if (ordercount < ref_len * factor):
                ordercount+=1
            if ref_contact_map[c_x, c_y] > 0:
                orderTP += 1.0 
            else:
                orderFP += 1.0 
            orderPPVs.append(orderTP / (orderTP + orderFP))
            orderTPs.append(orderTP/ordercount)
            orderFPs.append(orderFP/ordercount)
        if (num_c < ref_len * factor):
            if ref_contact_map[c_x, c_y] > 0:
                TP += 1.0 / (ref_len*factor)
            else:
                FP += 1.0 / (ref_len*factor)
            PPVs.append(TP / (TP + FP))
            TPs.append(TP)
            FPs.append(FP)

            

    if len(PPVs) == 0:
        PPVs.append(0.0)
    if len(TPs) == 0:
        TPs.append(0.0)
    if len(FPs) == 0:
        FPs.append(0.0)
    if len(mixPPVs) == 0:
        mixPPVs.append(0.0)
    if len(mixTPs) == 0:
        mixTPs.append(0.0)
    if len(mixFPs) == 0:
        mixFPs.append(0.0)
    if len(orderPPVs) == 0:
        orderPPVs.append(0.0)
    if len(orderTPs) == 0:
        orderTPs.append(0.0)
    if len(orderFPs) == 0:
        orderFPs.append(0.0)
    if len(disoPPVs) == 0:
        disoPPVs.append(0.0)
    if len(disoTPs) == 0:
        disoTPs.append(0.0)
    if len(disoFPs) == 0:
        disoFPs.append(0.0)

    return PPVs, TPs, FPs,orderPPVs, orderTPs, orderFPs,mixPPVs, mixTPs, mixFPs,disoPPVs, disoTPs, disoFPs


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
 

def get_colors(contacts_np, ref_contact_map=[], atom_seq_ali=[], th=1.0,factor=1.0):

    N = contacts_np.shape[0]
    img = np.ones((N,N,4))
    vector = np.zeros(N*N)
    numth = th
    num=0
    for i in xrange(N):
        for j in xrange(N):
            if j < (i+5):
                continue
            else:
                vector[num] = contacts_np[i,j]
                num=num+1
    vector[::-1].sort()
    num=0
    for i in vector:
        if num<=N*factor:
            numth=i
        else:
            break
        num=num+1

    for i in xrange(N):
        for j in xrange(N):
            if j < (i+5):
                continue
            num=num+1
            sc = contacts_np[i,j]
            if len(ref_contact_map) > 0:
                assert N == ref_contact_map.shape[0]
                # FN
                if ref_contact_map[i,j] < 8:
                    if sc <= th:
                    #img[i,j] = [0.5,0.5,1,1]
                        img[i,j] = [0.5,0.5,0.5,1]
                    else:
                        img[i,j] = [0,1,0,1]
                    if sc <= numth:
                        img[j,i] = [0.5,0.5,0.5,1]
                    else:
                        img[j,i] = [0,1,0,1]
                # TP
                #elif contacts_np[i,j] > th and ref_contact_map[i,j] >= 8:
                elif ref_contact_map[i,j] >= 12:
                    if sc > th:
                        img[i,j] = [1,0,0,1]
                    if sc > numth:
                        img[j,i] = [1,0,0,1]
                    #img[j,i] = [1-sc,1-sc,1-sc,1]
                # grey zone between 8 and 12 Angstroem
                elif (ref_contact_map[i,j] < 12 or ref_contact_map[i,j] >= 8):
                    if sc > th:
                        val = (ref_contact_map[i,j] - 8)/(12 - 8)
                        img[i,j] = [0.5+val/2,1-val/2,0,1]
                    if sc > numth:
                        val = (ref_contact_map[i,j] - 8)/(12 - 8)
                        img[j,i] = [0.5+val/2,1-val/2,0,1]
                    #img[j,i] = [1-sc,1-sc,1-sc,1]
            else:
                if sc > th:
                    img[i,j] = [0.5-sc/2,0.5-sc/2,1,1]
                if sc > numth:
                    img[j,i] = [0.5-sc/2,0.5-sc/2,1,1]

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
    #coverage_lst = [1-(coverage[i]/float(N)) for i in range(L)]
    coverage_lst = [N - coverage[i] for i in range(L)]
    return coverage_lst


def get_meff_coverage(meff_file):
    meff_lst = []
    with open(meff_file) as f:
        for l in f:
            if 'MeffPerPos' in l:
                meff_lst = l.split()[-1].strip('[]').split(',')
                meff_lst = map(float, meff_lst)
    return meff_lst


def plot_map(fasta_filename, c_filename, factor=1.0, th=-1, c2_filename='', psipred_horiz_fname='', psipred_vert_fname='',iupred_fname='', pdb_filename='', is_heavy=False, chain='', sep=',', outfilename='', ali_filename='',  meff_filename='', name='', start=0, end=-1):  
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

    if ali_filename:
        coverage_lst = get_ali_coverage(ali_filename)
        max_cover = max(coverage_lst)
        # Current Meff format does not work
        #    elif meff_filename:
        #        coverage_lst = get_meff_coverage(meff_filename)
        #        max_cover = max(coverage_lst)
    else:
        max_cover = 0
    average_disorder=0.
    average_order=0.
    fraction_disorder=0.
    cover_order=0.
    cover_disorder=0.
    if iupred_fname:
        disorder = parse_iupred.pred(open(iupred_fname, 'r'))
    else:
        disorder=np.zeros(ref_len)
    average_disorder = np.sum(disorder)/ref_len
    fraction_disorder = 0.0
    num_disorder=0
    num_order=0
    j=0
    for i in disorder:
        if (i>0.5):
            fraction_disorder +=1/ref_len
            num_disorder+=1
            cover_disorder+=coverage_lst[j]
        else:
            num_order+=1
            cover_order+=coverage_lst[j]
        j+=1

    if (num_disorder > 0):
        cover_disorder=cover_disorder/num_disorder
    if (num_order > 0):
        cover_order=cover_order/num_order
            
    ### get top "factor" * "ref_len" predicted contacts
    contacts = parse_contacts.parse(open(c_filename, 'r'), sep,1)
    contacts_np = parse_contacts.get_numpy_cmap(contacts)
    contacts_np = contacts_np[start:end,start:end]

    contacts_x = []
    contacts_y = []
    scores = []
    mixscores = []
    disoscores = []
    orderscores = []
    tooclose = []    
    contact_dict = {}
    if iupred_fname:
        disorder = parse_iupred.pred(open(iupred_fname, 'r'))
    else:
        disorder = np.zeros(ref_len)

    count = 1.e-20
    mixcount = 1.e-20
    disocount = 1.e-20
    ordercount =1.e-20
    longcount = 1.e-20
    highscore = 0
    numbins=20
    sum=0.0
    longsum=0.0
    disosum=0.0
    ordersum=0.0
    mixsum=0.0
    average=0.0
    longaverage=0.0
    mixaverage=0.0
    disoaverage=0.0
    orderaverage=0.0
    histo=np.zeros(numbins)
    disotop=0
    ordertop=0
    doubletop=0
    mixcount=0
    mixtop=0
    separation = 0.0




    # We actually divide the analysis into three groups (ordered,disordered and mixed)
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
        long_dist = pos_diff > 24
        if not too_close:
            if score > th:
                contacts_x.append(c_x - start)
                contacts_y.append(c_y - start)
                if (disorder[c_x] > 0.5 and disorder[c_y] > 0.5):
                    disocount += 1
                    disoscores.append(score)
                    if (disocount <= ref_len * factor):
                        disosum += score
                        disoaverage=disosum/disocount
                elif (disorder[c_x] > 0.5 or disorder[c_y] > 0.5):
                    mixcount += 1
                    mixscores.append(score)
                    if (mixcount <= ref_len * factor):
                        mixsum += score
                        mixaverage=mixsum/mixcount
                else:
                    ordercount += 1
                    orderscores.append(score)
                    if (ordercount <= ref_len * factor):
                        ordersum += score
                        orderaverage=ordersum/ordercount
                count += 1
                scores.append(score)
                if (count <= ref_len * factor):
                    sum += score
                    average=sum/count
                separation += pos_diff
                if long_dist:
                    longcount += 1
                    if (longcount <= ref_len * factor):
                        longsum += score
                        longaverage=longsum/longcount
        else:
            tooclose.append(score)
                

#    statline="Highs: %.1f (%.1f%%) (%.1f%%)  average:  %.2f (%.2f) (%.2f)  Meff: %.0f  Diso: %.1f%%  " % (count/ref_len,100*mixcount/count,100*disocount/count,average,mixaverage,disoaverage,max_cover,100*fraction_disorder)
#    statline="Highs: %.1f %.3f %.3f  average:  %.2f %.2f %.2f  Meff: %.0f  Diso: %.3f  " % (count/ref_len,mixcount/count,disocount/count,average,mixaverage,disoaverage,max_cover,fraction_disorder)
#    statline="Length: %d NumAli: %d Counts: %d %d %d %.3f %.3f %.3f %.3f\n"  % ( ref_len,max_cover,(count-mixcount-disocount),mixcount,disocount,sum,mixsum,disosum,fraction_disorder)
    statline="NumAli: %d %d %d Length: %d %d %d Counts: %d %d %d %d  RelContacts: %.3f %.3f %.3f Disorder: %.3f Long: %.3f %.3f %.3f \n"  % ( max_cover,cover_order,cover_disorder,ref_len,num_order,num_disorder,count,ordercount,mixcount,disocount,count/(ref_len+1.e-20),ordercount/(1.e-20+num_order),disocount/(1.e-20+num_disorder),fraction_disorder,longcount/count,longcount/(ref_len+1.e-20),separation/count)
    statfig = plt.figure(figsize=(8, 8), dpi=96, facecolor='w')
    plt.hist((tooclose,scores), numbins,range=(0,1), histtype='bar',
             normed=(numbins,numbins), alpha=0.75,
             label=['Too_Close','Contacts'])
    plt.xlabel('Score')
    plt.ylabel('Normalized count')
    statfig.suptitle('%s\n%s\n' %  (c_filename,line))

 

    ### start plotting
    fig = plt.figure(figsize=(8, 8), dpi=96, facecolor='w')
    ax = fig.add_subplot(111)#, aspect='auto')
    ax.set_adjustable('box-forced')
    ax.tick_params(direction='out', right='off', top='off')
    ax.set_xlim([-unit,ref_len])
    ax.set_ylim([-unit,ref_len])
    max_cover=0
    ### plot alignment coverage if alignemnt given (only on Y-axis)
    if ali_filename: # or meff_filename:
        # adjust overall canvas  
        ax = plt.subplot2grid((8,8), (1, 1), colspan=7, rowspan=7)#, aspect='auto')
        #ax.set_adjustable('box-forced')
        #ax.set_autoscale_on(False) 
        ax.autoscale(False)
        ax.tick_params(direction='out',labelleft='off', right='off', top='off')
        ax.set_xlim([-unit,ref_len])
        ax.set_ylim([-unit,ref_len])

        if ali_filename:
            coverage_lst = get_ali_coverage(ali_filename)
        #elif meff_filename:
        #    coverage_lst = get_meff_coverage(meff_filename)
        max_cover = max(coverage_lst)

        
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

        #ax3 = plt.subplot2grid((8,8), (0,1), colspan=7, sharex=ax)
        #ax3.set_adjustable('box-forced')
        #ax3.set_autoscale_on(False) 
        #ax3.autoscale(False)
        #ax3.plot([0]+range(ref_len)+[ref_len-1], [0]+coverage_lst+[0], 'k', lw=0)
        #ax3.axhline(y=max_cover*0.25, lw=0.5, c='black', ls=':')
        #ax3.axhline(y=max_cover*0.5, lw=0.5, c='black', ls=':')
        #ax3.axhline(y=max_cover*0.75, lw=0.5, c='black', ls=':')
        #ax3.fill([0]+range(ref_len)+[ref_len-1], [0]+coverage_lst+[0], facecolor='gray', lw=0, alpha=0.5)
        #ax3.xaxis.tick_top()
        #ax3.set_yticks([0, max_cover])
        #ax3.tick_params(labelbottom='off')
        ax2.tick_params(axis='y', right='off', direction='out', left='on')
        #ax3.spines['top'].set_visible(False)
        #ax3.spines['right'].set_visible(False)
        #ax.get_xaxis().tick_top()
        #ax.get_yaxis().tick_left()
        #ax3.grid()
        #ax3.set_xlim([-unit,ref_len])


    ### plot secondary structure along axis if given
    average_disorder=0.
    fraction_disorder=0.
    #statline = "Highs: %.1f Aver: %.2f  Meff: %.0f" % (count/ref_len,average,max_cover)
    ax.get_xaxis().tick_top()
    ax.get_yaxis().tick_left()
    if iupred_fname:
        ax = plt.subplot2grid((8,8), (1, 1), colspan=7, rowspan=7)#, aspect='auto')
        #ax.set_adjustable('box-forced')
        #ax.set_autoscale_on(False) 
        ax.autoscale(False)
        ax.tick_params(direction='out',labelleft='off', right='off', top='off')
        ax.set_xlim([-unit,ref_len])
        ax.set_ylim([-unit,ref_len])

        average_disorder = np.sum(disorder)/ref_len
        fraction_disorder = 0.0
        for d in disorder:
            if (d>0.5):
                fraction_disorder += 1/ref_len
                
        ax3 = plt.subplot2grid((8,8), (0,1), colspan=7, sharex=ax)
        ax3.set_adjustable('box-forced')
        ax3.set_autoscale_on(False) 
        ax3.autoscale(False)
        ax3.plot([0]+range(ref_len)+[ref_len-1], [0]+disorder+[0], 'b', lw=2)
        ax3.axhline(y=0.5, lw=0.5, c='black', ls=':')
        #ax3.fill([0]+range(ref_len)+[ref_len-1], [0]+disorder+[0], facecolor='gray', lw=0, alpha=0.5)
        ax3.xaxis.tick_top()
        ax3.set_yticks([0, 1])
        ax3.tick_params(labelbottom='off')
        ax3.spines['top'].set_visible(False)
        ax3.spines['right'].set_visible(False)
        ax3.grid()
        ax3.set_xlim([-unit,ref_len])
        #statline = "Highs: %.3f  %.3f   %.3f    Aver: %.2f  Meff: %.0f  Diso: %.3f  " % (count/ref_len,disocount/count,doublecount/count,average,max_cover,fraction_disorder)


    print "STATs: %s %s" % (c_filename,statline)
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
        chain='*'
        # We try to get all chains...
        res_lst = parse_pdb.get_coordinates(open(pdb_filename, 'r'), chain)
#        cb_lst = parse_pdb.get_ca_coordinates(open(pdb_filename, 'r'), chain)
        cb_lst = parse_pdb.get_cb_coordinates(open(pdb_filename, 'r'), chain)
        atom_seq = parse_pdb.get_atom_seq(open(pdb_filename, 'r'), chain)
                
#        align = pairwise2.align.globalms(atom_seq, seq, 2, -1, -0.5, -0.1)
        align = pairwise2.align.globalms(atom_seq, seq, 2, -1, -10.5, -10.1)

        atom_seq_ali = align[-1][0]
        seq_ali = align[-1][1]

        print atom_seq_ali
        print seq_ali

        j = 0
        gapped_res_lst = []
        gapped_cb_lst = []

        for i in xrange(len(atom_seq_ali)):
#            print i,j
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
            
        PPVs, TPs, FPs,orderPPVs, orderTPs, orderFPs,mixPPVs, mixTPs, mixFPs,disoPPVs, disoTPs, disoFPs = get_ppvs(contacts_x, contacts_y, ref_contact_map, atom_seq_ali, ref_len, factor,disorder)
        tp_colors = get_tp_colors(contacts_x, contacts_y, ref_contact_map, atom_seq_ali)
        img = get_colors(contacts_np, ref_contact_map=dist_mat, atom_seq_ali=atom_seq_ali, th=th, factor=factor)
        sc = ax.imshow(img, interpolation='none')
   
        #       
        print 'PPV: %s %s %s %s %s' % (acc, PPVs[-1], orderPPVs[-1], mixPPVs[-1], disoPPVs[-1])
        
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
            PPVs, TPs, FPs,orderPPVs, orderTPs, orderFPs,mixPPVs, mixTPs, mixFPs,disoPPVs, disoTPs, disoFPs = get_ppvs(contacts_x, contacts_y, ref_contact_map, atom_seq_ali, ref_len, factor,disorder)
            tp2_colors = get_tp_colors(contacts2_x, contacts2_y, ref_contact_map, atom_seq_ali)
            print '%s %s %s %s' % (acc, PPVs2[-1], TPs2[-1], FPs2[-1])
            fig.suptitle('%s\nPPV (upper left) = %.2f | PPV (lower right) = %.2f' % (acc, PPVs[-1], PPVs2[-1]))
            sc = ax.scatter(contacts2_y[::-1], contacts2_x[::-1], marker='o', c=tp2_colors[::-1], s=6, alpha=0.75, lw=0)
            sc = ax.scatter(contacts_x[::-1], contacts_y[::-1], marker='o', c=tp_colors[::-1], s=6, alpha=0.75, lw=0)
        else:
            sc = ax.scatter(contacts2_y[::-1], contacts2_x[::-1], marker='0', c='#D70909', edgecolor='#D70909', s=6, linewidths=0.5)
            sc = ax.scatter(contacts_x[::-1], contacts_y[::-1], marker='o', c='#004F9D', edgecolor='#004F9D', s=6, linewidths=0.5)


           
#        ppv='PPV: %.2f %.2f %.2f' % (float(PPVs[-1]), float(TPs[-1]), float(FPs[-1]))
        ppv='PPV: %.2f (%d) %.2f (%d) %.2f (%d) ' % (float(PPVs[-1]),len(PPVs), float(mixPPVs[-1]),len(mixPPVs), float(disoPPVs[-1]),len(disoPPVs))
    else:
        if pdb_filename:
            pdb_acc = parse_pdb.get_acc(open(pdb_filename))
            if pdb_acc:
                if chain:
                    fig.suptitle('%s (PDB: %s, chain %s)\nPPV = %.2f\n%s' % (c_filename , pdb_acc, chain, PPVs[-1],statline),fontsize = 8)
                else:
                    fig.suptitle('%s (PDB: %s)\nPPV = %.2f \n%s' % (c_filename, pdb_acc, PPVs[-1],statline),fontsize = 8)
            else:
                fig.suptitle('%s\nPPV = %.2f\n%s' % (c_filename, PPVs[-1],statline),fontsize = 8)
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
            fig.suptitle('%s\n%s' % (c_filename,statline),fontsize = 8)
            #sc = ax.imshow(contacts_np + contacts_np.T, cmap=cm.hot_r)
            #sc = ax.imshow(contacts_np + contacts_np.T,
            #        cmap=cm.binary, vmin=th, vmax=1.0, interpolation='none')
            img = get_colors(contacts_np, th=th,factor=factor)
            sc = ax.imshow(img, interpolation='none')
            #divider1 = make_axes_locatable(ax)
            #cax1 = divider1.append_axes("right", size="2%", pad=0.05)
            #plt.colorbar(sc, cax=cax1)
            #plt.colorbar(sc, ax=ax)
            #sc = ax.scatter(contacts_x[::-1], contacts_y[::-1],
            #        marker='o', c="black", s=6, alpha=0.75,
            #        linewidths=0.1, edgecolors='none')
            #sc = ax.scatter(contacts_y[::-1], contacts_x[::-1], marker='o', c=scores[::-1], s=10, alpha=0.75, cmap=cm.hot_r, linewidths=0.1, edgecolors='none')

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
            ppstat = PdfPages(outfilename+"_statistics.pdf")
            ppstat.savefig(statfig)
            ppstat.close()
        elif outfilename.endswith('.png'):
            plt.savefig(outfilename)
        else:
            pp = PdfPages('%s.pdf' % outfilename)
            pp.savefig(fig)
            pp.close()
            ppstat = PdfPages(outfilename+"_statistics.pdf")
            ppstat.savefig(statfig)
            ppstat.close()
    else:
        pp = PdfPages('%s_ContactMap.pdf' % c_filename)
        pp.savefig(fig)
        pp.close()
        ppstat = PdfPages('%s_statistics.pdf' % c_filename)
        ppstat.savefig(statfig)
        ppstat.close()


    
if __name__ == "__main__":

    p = argparse.ArgumentParser(description='Plot protein residue contact maps.')
    p.add_argument('fasta_file')#, required=True)
    p.add_argument('contact_file')#, required=True)
    p.add_argument('-o', '--outfile', default='')
    p.add_argument('-f', '--factor', default=1.0, type=float)
    p.add_argument('-t', '--threshold', default=0.4, type=float)
    p.add_argument('--c2', default='')
    p.add_argument('--psipred_horiz', default='')
    p.add_argument('--psipred_vert', default='')
    p.add_argument('--iupred', default='')
    p.add_argument('--pdb', default='')
    p.add_argument('--heavy', action='store_true')
    p.add_argument('--chain', default='')
    p.add_argument('--alignment', default='')
    p.add_argument('--meff', default='')
    p.add_argument('--name', default='')
    p.add_argument('--start', default=0, type=int)
    p.add_argument('--end', default=-1, type=int)

    args = vars(p.parse_args(sys.argv[1:]))

    
    fasta_filename = args['fasta_file']
    c_filename = args['contact_file']
    psipred_filename = args['psipred_horiz']
    iupred_filename = args['iupred']

    # guessing separator of constraint file
    line = open(c_filename,'r').readline()
    if len(line.split(',')) != 1:
        sep = ','
    elif len(line.split(' ')) != 1:
        sep = ' '
    else:
        sep = '\t'

    factor=args['factor']
    cutoff=args['threshold']


    plot_map(args['fasta_file'], args['contact_file'], factor,  th=args['threshold'], c2_filename=args['c2'], psipred_horiz_fname=args['psipred_horiz'], psipred_vert_fname=args['psipred_vert'], iupred_fname=args['iupred'], pdb_filename=args['pdb'], is_heavy=args['heavy'], chain=args['chain'], sep=sep, outfilename=args['outfile'], ali_filename=args['alignment'], meff_filename=args['meff'], name=args['name'], start=args['start'], end=args['end'])

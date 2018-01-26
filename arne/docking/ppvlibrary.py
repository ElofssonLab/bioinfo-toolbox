from __future__ import division
import sys
import string
import re
from math import *
import Bio.PDB
import numpy as np
import operator
# on UPPMAX only
sys.path.append('/sw/apps/bioinfo/biopython/1.59/tintin/lib/python')

from Bio import pairwise2

import numpy as np

import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from os.path import expanduser
home = expanduser("~")
sys.path.append(home + '/git/bioinfo-toolbox/')
from parsing import parse_contacts
from parsing import parse_psipred
from parsing import parse_fasta
from parsing import parse_iupred
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
    for i, cb1 in enumerate(gapped_cb_lst):
        if cb1[0] ==  '-':
            continue
        for j, cb2 in enumerate(gapped_cb_lst):
            if cb2[0] == '-':
                continue
            diff_vec = cb1 - cb2
            dist_mat[i,j] = np.sqrt(np.sum(diff_vec * diff_vec))
    return dist_mat


def print_distances(contacts_x, contacts_y, scores, dist_mat, area, dist, lenA,lenB,seq,ali_lst=[], atom_seq=[], outfile=""):
    num_c = len(contacts_x)
    outstr = ""
    domain = ""
    
    for i in range(num_c):
        c_x = contacts_x[i]
        c_y = contacts_y[i]
        if ali_lst:
            if ali_lst[c_x] < 0:
                continue
            if ali_lst[c_y] < 0:
                continue
        if ali_lst[c_x] < lenA and ali_lst[c_y] < lenA:
            domain="A"
        elif ali_lst[c_x] >= lenA and ali_lst[c_y] >= lenA:
            domain="B"
        else:
            domain="I"
        if ali_lst[c_x] > 0:
            areaX=area[ali_lst[c_x]][1]
            distX=dist[ali_lst[c_x]]
        else:
            areaX=-1
            distX=-1
        if ali_lst[c_y] > 0:
            areaY=area[ali_lst[c_y]][1]
            distY=dist[ali_lst[c_y]]
        else:
            areaY=-1
            distY=-1
        if ali_lst:
            resnoX=ali_lst[c_x]+1
            if resnoX>lenA:
                resnoX-=lenA
            resnoY=ali_lst[c_y]+1
            if resnoX>lenA:
                resnoX-=lenA
            outstr += "%s %s %s %s %s %.2f %.2f %.2f %.2f %s %s %s %s \n" % (domain,c_x,atom_seq[ali_lst[c_x]], c_y,atom_seq[ali_lst[c_y]],
                                                                  areaX,areaY,distX,distY, scores[i],
                                                                             dist_mat[c_x, c_y],resnoX,resnoY)
        else:
            resnoX=c_x+1
            if resnoX>lenA:
                resnoX-=lenA
            resnoY=c_y+1
            if resnoX>lenA:
                resnoX-=lenA
            outstr += "%s %s %s %s %s %.2f %.2f %.2f %.2f %s %s %s %s \n" % (domain,c_x,seq[c_x], c_y,seq[c_y],
                                                                  areaX,areaY,distX,distY, scores[i],
                                                                             dist_mat[c_x, c_y],c_x,c_y)
    if outfile:
        with open(outfile, 'w') as outf:
            outf.write(outstr)
    else:
        print outstr

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


def get_interface_contacts(contacts_x, contacts_y, scores, dist_mat, ref_len, factor, cutoff, atom_seq_ali=[]):
    temp={}
    tempX={}
    tempY={}
    tempcontacts_x=[]
    tempcontacts_y=[]
    tempscores=[]
    for i in range(len(contacts_x)):
        c_x = contacts_x[i]
        c_y = contacts_y[i]
        if atom_seq_ali:
            if atom_seq_ali[c_x] == '-':
                continue
            if atom_seq_ali[c_y] == '-':
                continue
        if dist_mat[c_x, c_y] > cutoff and c_x < ref_len:
            tempcontacts_x.append(c_x)
            tempcontacts_y.append(c_y+ref_len)
            tempscores.append(scores[i])
        else:
            tempcontacts_x.append(c_x)
            tempcontacts_y.append(c_y)
            tempscores.append(scores[i])
            tempcontacts_x.append(c_x+ref_len)
            tempcontacts_y.append(c_y+ref_len)
            tempscores.append(scores[i])
            # Now we need to sort all scores...
    for i in range(len(contacts_x)):
        temp[i]=tempscores[i]
        tempX[i]=tempcontacts_x[i]
        tempY[i]=tempcontacts_y[i]
    temp_sorted = sorted(temp.items(), key=operator.itemgetter(1))
    contacts_x=[]
    contacts_y=[]
    scores=[]
    for i in reversed(temp_sorted):
        contacts_x.append(tempX[i[0]])
        contacts_y.append(tempY[i[0]])
        scores.append(temp[i[0]])
    return (contacts_x,contacts_y,scores)


def get_ppv_helper(contacts_x, contacts_y, ref_contact_map, ref_len, factor, atom_seq_ali=[]):
    num_c = int(ref_len*factor)
    TP = 0.0
    FP = 0.0
    PPV = 0.0
    Zscore = 0.0
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

def get_Zscore(contacts_x, contacts_y, ref_contact_map,  score, atom_seq_ali=[]):
    Zscore = 0.0

    #contact=np.array()
    #noncontact=np.array()
    numcontacts=0.
    numnoncontacts=0.
    sumcontacts=0.
    sumnoncontacts=0.
    sum2contacts=0.
    sum2noncontacts=0.
    for i in range(len(contacts_x)):
        c_x = contacts_x[i]
        c_y = contacts_y[i]
        if atom_seq_ali:
            if atom_seq_ali[c_x] == '-':
                continue
            if atom_seq_ali[c_y] == '-':
                continue
        if ref_contact_map[c_x, c_y] > 0:
            numcontacts+=1
            sumcontacts+=score[i]
            sum2contacts+=score[i]*score[i]
        else:
            numnoncontacts+=1
            sumnoncontacts+=score[i]
            sum2noncontacts+=score[i]*score[i]
                                                                                                                    
    avecontact=sumcontacts/numcontacts
    avenoncontact=sumnoncontacts/numnoncontacts
    sdnoncontact=sqrt((numnoncontacts*sum2noncontacts-sumnoncontacts*sumnoncontacts)/(numnoncontacts*(numnoncontacts-1)))
    Zscore=(avecontact-avenoncontact)/sdnoncontact
    return (Zscore)



def get_Zscore_interface(contacts_x, contacts_y, ref_contact_map, chainlenA,chainlenB, score, atom_seq_ali=[]):
    Zscore = 0.0

    #contact=np.array()
    #noncontact=np.array()
    numcontacts=0.
    numnoncontacts=0.
    sumcontacts=0.
    sumnoncontacts=0.
    sum2contacts=0.
    sum2noncontacts=0.
    for i in range(len(contacts_x)):
        c_x = contacts_x[i]
        c_y = contacts_y[i]
        if ( c_x < chainlenA and c_y > chainlenA) or ( c_x > chainlenA and c_y < chainlenA):
            if atom_seq_ali:
                if atom_seq_ali[c_x] == '-':
                    continue
                if atom_seq_ali[c_y] == '-':
                    continue
            if ref_contact_map[c_x, c_y] > 0:
                numcontacts+=1
                sumcontacts+=score[i]
                sum2contacts+=score[i]*score[i]
            else:
                numnoncontacts+=1
                sumnoncontacts+=score[i]
                sum2noncontacts+=score[i]*score[i]
                                                                                                                    
    avecontact=sumcontacts/numcontacts
    avenoncontact=sumnoncontacts/numnoncontacts
    sdnoncontact=sqrt((numnoncontacts*sum2noncontacts-sumnoncontacts*sumnoncontacts)/(numnoncontacts*(numnoncontacts-1)))
    Zscore=(avecontact-avenoncontact)/sdnoncontact
    return (Zscore)


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

def get_tp_colors(contacts_x, contacts_y, ref_contact_map, atom_seq_ali):

    tp_colors = []

    for i in range(len(contacts_x)):
        c_x = contacts_x[i]
        c_y = contacts_y[i]
        if atom_seq_ali[c_x] == '-':
            #tp_colors.append('green')
            tp_colors.append('grey')
            continue
        if atom_seq_ali[c_y] == '-':
            #tp_colors.append('green')
            tp_colors.append('grey')
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


def plot_map(contacts_x, contacts_y, scores, contact_map, lenA,lenB, factor, numinter, ali_lst=[],atom_seq=[],outfile=""):

    
    tp_colors = get_tp_colors(contacts_x, contacts_y, ref_contact_map, atom_seq_ali)
    img = get_colors(contacts_np, ref_contact_map=dist_mat, th=th_obs, binary=binary)
    sc = ax.imshow(img, interpolation='none')
    
    fig = plt.figure(figsize=(8, 8), dpi=96, facecolor='w')
    ax = fig.add_subplot(111)#, aspect='auto')
    ax.set_adjustable('box-forced')
    ax.tick_params(direction='out', right='off', top='off')
    ax.set_xlim([-unit,ref_len])
    ax.set_ylim([-unit,ref_len])
    max_cover=0
    img = get_colors(contacts_np, ref_contact_map=dist_mat, atom_seq_ali=atom_seq_ali, th=th, factor=factor)
    sc = ax.imshow(img, interpolation='none')
    
    #       
    
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
    

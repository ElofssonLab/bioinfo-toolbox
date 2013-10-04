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
#import parse_psipred
import parse_fasta
#import parse_jones
#import parse_hmm
#import dotter


# scale value x from [min_x, max_x] to [0,1]
def scale(x, min_x, max_x):
    range_x = (max_x - min_x)
    #print range_x
    return (x - min_x) / range_x


# calculates euclidean distance between CB atoms of two reisdues
# in case of glycine (no CB) the CA atom is used instead
def calc_residue_dist(res1, res2):
    na_flag = False
    """
    if 'CB' in res1:
        res1_coord = res1['CB'].coord
    elif 'CA' in res1:
        res1_coord = res1['CA'].coord
    else:
        #print "BLA RES1"
        na_flag = True
    if 'CB' in res2:
        res2_coord = res2['CB'].coord
    elif 'CA' in res2:
        res2_coord = res2['CA'].coord
    else:
        #print "BLA RES2"
        na_flag = True
    """
    if len(res1.get_list()) > 4:
        res1_coord = res1.get_list()[4].coord
    elif len(res1.get_list()) > 1:
        res1_coord = res1.get_list()[1].coord
    else:
        res1_coord = res1.get_list()[0].coord
    
    if len(res2.get_list()) > 4:
        res2_coord = res2.get_list()[4].coord
    elif len(res1.get_list()) > 1:
        res2_coord = res2.get_list()[1].coord
    else:
        res2_coord = res2.get_list()[0].coord

    if not na_flag:
        diff_vector = res1_coord - res2_coord
        return np.sqrt(np.sum(diff_vector * diff_vector))
    else:
        return float('inf')


# calculate the distance between two residues w.r.t. their heavy atoms
def calc_residue_dist_heavy(res1, res2):
    min_dist = float('inf')
    for atm1 in res1.get_list():
        for atm2 in res2.get_list():
            diff_vector = atm1.coord - atm2.coord
            dist = np.sqrt(np.sum(diff_vector * diff_vector))
            if dist < min_dist:
                min_dist = dist
    return min_dist


# creates pairwise eucledian distance matrix of two protein chains
def calc_dist_matrix(chain1, chain2, ref_len):
    #dist_mat = np.zeros((len(chain1), len(chain2)), np.float)
    dist_mat = np.zeros((ref_len, ref_len), np.float)
    dist_mat.fill(float('inf'))
    print len(chain1)
    for row, res1 in enumerate(chain1):
        if row < ref_len:
            for col, res2 in enumerate(chain2):
                if col < ref_len:
                    if res1.id[0] == ' ' and res2.id[0] == ' ':
                        #print res1.get_full_id()
                        #print res2.get_full_id()
                        dist_mat[row, col] = calc_residue_dist(res1, res2)

    return dist_mat


# creates pairwise eucledian distance matrix of two protein chains
def calc_dist_matrix_heavy(chain1, chain2):
    dist_mat = np.zeros((len(chain1), len(chain2)), np.float)
    print len(chain1)
    for row, res1 in enumerate(chain1):
        for col, res2 in enumerate(chain2):
            dist_mat[row, col] = calc_residue_dist_heavy(res1, res2)
    return dist_mat


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



def plot_map(contact_filename, psipred_filename, pdb_filename, rep_len, sep=','):  
    
    #ss = parse_psipred.horizontal(open(psipred_filename, 'r'))
    #seq = parse_psipred.horizontal_seq(open(psipred_filename, 'r'))
    seq = parse_fasta.read_fasta(open(psipred_filename, 'r')).values()[0][0]
    ss = seq
    ref_len = len(seq)

    pdb_code = pdb_filename.split('/')[-1].split('.')[0]
    ref_struct = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filename)
    ref_model = ref_struct[0]
    ref_chain = ref_model.get_list()[0]
    #ref_len = len(ref_chain)
    #print pdb_code
    #print ref_chain

    pdb_offset = 0
    seq_offset = 0
    pdb_peptides = Bio.PDB.PPBuilder().build_peptides(ref_struct)
    print seq
    if len(pdb_peptides) > 0:
        pdb_seq = pdb_peptides[0].get_sequence()
        print str(pdb_seq)
        align = pairwise2.align.globalxs(str(pdb_seq), seq, -10, -0.5)
        print pdb_seq
        print len(pdb_seq)
        print abs(len(pdb_seq) - ref_len)
        print align
        for res in align[-1][0]:
            if res == '-':
                pdb_offset += 1
            else:
                break
        for res in align[-1][1]:
            if res == '-':
                seq_offset += 1
            else:
                break
    print pdb_offset
    print seq_offset
    print ref_len

    #dist_mat = calc_dist_matrix_heavy(ref_chain, ref_chain)
    dist_mat = calc_dist_matrix(ref_chain, ref_chain, ref_len)
    cutoff = 8
    #ref_contact_map = (dist_mat < 8) & (dist_mat > 4)
    ref_contact_map = dist_mat < cutoff
    ref_contact_map += (pdb_offset - seq_offset)
    #ref_contacts = numpy.where((dist_mat < 8) & (dist_mat > 4))
    ref_contacts = np.where(dist_mat < cutoff)

    ref_contacts_x = ref_contacts[0] + (pdb_offset - seq_offset)
    ref_contacts_y = ref_contacts[1] + (pdb_offset - seq_offset)

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

    #print dist_mat
    #print ref_contact_map
    #print 'ref_contacts=' + str(ref_contacts[1])

    contacts = parse_contacts.parse(open(contact_filename, 'r'), sep)
    #contacts_cut = contacts[0:ref_len]

    contacts_nf_x = []
    contacts_nf_y = []
    scores_nf = []
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
            contacts_nf_x.append(c_x)
            contacts_nf_y.append(c_y)
            scores_nf.append(score)

        if not in_same_helix(c_x, c_y, ss) or in_same_helix(c_x, c_y, ss):

            if (not too_close): #and (not too_far):
                #contacts_x.append(c_x)
                #contacts_y.append(c_y)
                #scores.append(score)

                c_key = '%d-%d' % (c_x, c_y)
                if contact_dict.has_key(c_key):
                    contact_dict[c_key].append(score)
                else:
                    contact_dict[c_key] = [score]
                count += 1
             
            elif (not too_close) and too_far:
                factor = round(pos_diff / float(rep_len)) - 1
                #contacts_x.append(c_x)# - (factor * rep_len))
                #contacts_y.append(int(c_y - (factor * rep_len)))
                #scores.append(score)

                c_vkey = '%d-%d' % (c_x, c_y - (factor * rep_len))
                if contact_dict.has_key(c_vkey):
                    contact_dict[c_vkey].append(score)
                else:
                    contact_dict[c_vkey] = [score]
                
                #contacts_x.append(int(c_x + (factor * rep_len)))
                #contacts_y.append(c_y)# - (factor * rep_len))
                #scores.append(score)

                c_hkey = '%d-%d' % (c_x + (factor * rep_len), c_y)
                if contact_dict.has_key(c_hkey):
                    contact_dict[c_hkey].append(score)
                else:
                    contact_dict[c_hkey] = [score]
                count += 1
            
        if count > ref_len * 1.0:
            break
        #if score == 0.0:
        #    break

     
    contacts_x = []
    contacts_y = []
    scores = []
    for key, scs in contact_dict.iteritems():
        c_x = int(key.split('-')[0])
        c_y = int(key.split('-')[1])

        if len(scs) > 1:
            score = np.mean(scs)
            contacts_x.append(c_x)
            contacts_y.append(c_y)
            scores.append(score)
        else:
            score = scs[0]
            contacts_x.append(c_x)
            contacts_y.append(c_y)
            scores.append(score)

    #print numpy.where(contact_map == 1)

    PPVs = []
    TPs = []
    FPs = []

    print len(contacts_x)
    for num_c in range(min(len(contacts_x), ref_len * 1))[1:]:
        TP = 0.0
        FP = 0.0
        for i in range(num_c):
            #c_x = contacts[i][1] - 1
            #c_y = contacts[i][2] - 1
            c_x = contacts_x[i]
            c_y = contacts_y[i]
            if ref_contact_map[c_x, c_y]:
                TP += 1.0
            else:
                FP += 1.0
        PPVs.append(TP / (TP + FP))
        TPs.append(TP / ref_len)
        FPs.append(FP / ref_len)

    #print PPVs[-1]
    #print TPs[-1]
    #print FPs[-1]

    ### get pairwise distances from the alignment
    #alinum = parse_jones.get_numeric('%s.jones' % '.'.join(contact_filename.split('.')[:-1]))
    #alidist = squareform(pdist(alinum.T, 'euclidean'))
    #alidistlog = np.log(alidist)

    """
    dot_matrix = dotter.calc_dot_matrix(seq)
    profile = parse_hmm.read_hmm(open('%s.hmm' % '.'.join(contact_filename.split('.')[:-1])))
    dot_matrix = dotter.calc_dot_matrix_profile(seq, profile)
    tmp_dot_matrix = dot_matrix
    similarities_x = []
    similarities_y = []
    similarities_sc = []
    for (i,j), score in np.ndenumerate(dot_matrix):
        if i == j:
            continue
        if i > j:
            tmp_dot_matrix[i, j] = 0.0
            if score > 0.1:
                similarities_x.append(j)
                similarities_y.append(i)
                similarities_sc.append(score)
        else:
            if score > 0.1:
                similarities_x.append(i)
                similarities_y.append(j)
                similarities_sc.append(score)
    dot_matrix = tmp_dot_matrix
    """


    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    """
    for i in range(len(ss)):
        if ss[i] == 'H':
            plt.plot(i, i, 'o', c='#999999', mec="#444444")
        if ss[i] == 'E':
            plt.plot(i, i, 'D', c='#999999', mec="#444444", markersize=4)
        if ss[i] == 'C':
            plt.plot(i, i, 'D', c='#999999', mec='#999999', markersize=3)
    """
    #ax.imshow(dot_matrix, origin='lower', cmap=cm.binary)
    ax.scatter(ref_contacts_x, ref_contacts_y, marker='o', c='#CCCCCC', lw=0, edgecolor='#CCCCCC')
    ax.scatter(ref_contacts_y, ref_contacts_x, marker='o', c='#CCCCCC', lw=0, edgecolor='#CCCCCC')
    #sc = ax.scatter(contacts_x[::-1], contacts_y[::-1], marker='o', c=scores[::-1], s=8, alpha=0.75, cmap=cm.jet, linewidths=0.5)
    ax.scatter(contacts_x[::-1], contacts_y[::-1], marker='o', c='#', s=6, linewidths=0)
    ax.scatter(contacts2_x[::-1], contacts2_y[::-1], marker='o', c='#', s=6, linewidths=0)
    #sc = ax.scatter(contacts_x[::-1], contacts_y[::-1], marker='o', c='#0080AD', edgecolor='#0080AD', s=7, linewidths=0.5)
    #sc = ax.scatter(contacts_x[::-1], contacts_y[::-1], marker='o', c='#8B0043', edgecolor='#8B0043', s=7, linewidths=0.5)
    #sc = ax.scatter(contacts_y[::-1], contacts_x[::-1], marker='o', c='#0080AD', edgecolor='#0080AD', s=7, linewidths=0.5)
    #sc = ax.scatter(similarities_y[::-1], similarities_x[::-1], marker='x', edgecolor=similarities_sc[::-1], s=7, cmap=cm.binary, linewidths=0.5)
    #sc = ax.scatter(similarities_y[::-1], similarities_x[::-1], marker='x', edgecolor='#8B0043', s=7, linewidths=0.5)
    fig.suptitle(contact_filename)
    print ref_len

    plt.gca().set_xlim([0,ref_len])
    plt.gca().set_ylim([0,ref_len])
    #plt.colorbar(sc)

    pp = PdfPages('%s_ContactMap.pdf' % contact_filename)
    pp.savefig(fig)
    pp.close()

    """
    outfile = open('%s.contacts' % '.'.join(contact_filename.split('.')[0:-1]),'w')
    for i in range(len(scores)):
        outfile.write('%s,%s,%s\n' % (int(contacts_x[i] + 1), int(contacts_y[i] + 1), scores[i]))
    """
    """
    #pylab.matshow(numpy.transpose(contact_map), cmap=pylab.cm.binary, extent=(0, 143, 143, 0))
    #pylab.colorbar()in_same_helix(1,3,ss)

    pylab.plot(ref_contacts[0], ref_contacts[1], 'o', mec='0.8', mfc='0.8', mew=0.0, markersize=4)
    pylab.plot(contacts_x, contacts_y, 'o', mec='red', mfc='None', ls='None', markersize=3)

    pylab.xlabel('Residue index')
    pylab.ylabel('Residue index')
    pylab.title('%s - %s' % (contact_filename.split('/')[-1], pdb_code))
    pylab.savefig('%s_ContactMap.pdf' % contact_filename)
    pylab.close()
    """
    fig.clf()
    ax2 = fig.add_subplot(111)
    ax2.plot(PPVs)
    #pp = PdfPages('%s_PPVs.pdf' % contact_filename)
    #pp.savefig(fig)
    #pp.close() 


def plot_map_nopdb(contact_filename, psipred_filename, rep_len, sep=','):  

    ss = parse_psipred.horizontal(open(psipred_filename, 'r'))
    #print ss
    ref_len = len(ss)

    contacts = parse_contacts.parse(open(contact_filename, 'r'), sep)
    #contacts_cut = contacts[0:ref_len]

    contacts_x = []
    contacts_y = []
    scores = []

    count = 0
    #for i in range(ref_len * 1):
    for i in range(len(contacts)):
        score = contacts[i][0]
        c_x = contacts[i][1] - 1
        c_y = contacts[i][2] - 1
        pos_diff = abs(c_x - c_y)
        too_close = pos_diff < 5
        too_far = pos_diff > rep_len * 1.5
        if not in_same_helix(c_x, c_y, ss):
            if (not too_close) and (not too_far):
                contacts_x.append(c_x)
                contacts_y.append(c_y)
                scores.append(score)
                count += 1
             
            elif (not too_close) and too_far:
                factor = round(pos_diff / float(rep_len)) - 1
                contacts_x.append(c_x)# - (factor * rep_len))
                contacts_y.append(int(c_y - (factor * rep_len)))
                scores.append(score)
                contacts_x.append(int(c_x + (factor * rep_len)))
                contacts_y.append(c_y)# - (factor * rep_len))
                scores.append(score)
                count += 2
            
        if count > ref_len * 1.0:
            break



    scores_scaled = []
    for score in scores:
        scores_scaled.append(scale(score, min(scores), max(scores)))

    #print scores_scaled
    #print scores
    #print numpy.where(contact_map == 1)

    #pylab.matshow(numpy.transpose(contact_map), cmap=pylab.cm.binary, extent=(0, 143, 143, 0))
    #pylab.colorbar()
    #pylab.plot(ref_contacts[0], ref_contacts[1], 'o', mec='0.8', mfc='0.8', mew=0.0, markersize=4)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    for i in range(len(ss)):
        if ss[i] == 'H':
            plt.plot(i + 1, i + 1, 'o', c='#999999', mec="#444444")
        if ss[i] == 'E':
            plt.plot(i + 1, i + 1, 'D', c='#999999', mec="#444444", markersize=4)
        if ss[i] == 'C':
            plt.plot(i + 1, i + 1, '.', c='#999999', mec='#999999')#, markersize=3)

    sc = ax.scatter(contacts_x[::-1], contacts_y[::-1], marker='o', c=scores[::-1], s=15, alpha=0.75, cmap=cm.jet, linewidths=0.5)

    plt.gca().set_xlim([0,ref_len])
    plt.gca().set_ylim([0,ref_len])
    plt.colorbar(sc)
    #cbar = plt.colorbar(ax, ticks=[min(scores), max(scores)])
    #cbar.ax.set_yticklabels(['Low', 'High'])
    #cbar.set_label(r'Contact Score')

    pp = PdfPages('%s_ContactMap.pdf' % contact_filename)
    pp.savefig(fig)
    pp.close()

    outfile = open('%s.contacts' % '.'.join(contact_filename.split('.')[0:-1]),'w')
    for i in range(len(scores)):
        outfile.write('%s,%s,%s\n' % (int(contacts_x[i] + 1), int(contacts_y[i] + 1), scores[i]))

    """
    pylab.plot(contacts_x, contacts_y, 'o', mec=scale(contacts[contacts_x][contacts_y]), mfc='None', ls='None', markersize=3)

    # plot secondary structure elements
    for i in range(len(ss)):
        if ss[i] == 'H':
            pylab.plot(i, i, 'o', mec='blue', mfc='blue', ls='None', markersize=4)
        if ss[i] == 'E':
            pylab.plot(i, i, 'o', mec='green', mfc='blue', ls='None', markersize=2)
        if ss[i] == 'C':
            pylab.plot(i, i, '.', mec='gray', mfc='gray', ls='None', markersize=2)

    pylab.xlabel('Residue index')
    pylab.ylabel('Residue index')
    pylab.title('%s' % contact_filename.split('/')[-1])
    pylab.savefig('%s_ContactMap.pdf' % contact_filename)
    pylab.close()
    """

if __name__ == "__main__":

    if len(sys.argv) < 3:
        sys.stderr.write('Usage: python plot_contact_map.py <contact_filename> <psipred_filename> <pdb_filename> <rep_len>\n')

    contact_filename = sys.argv[1]
    psipred_filename = sys.argv[2]

    # guessing separator of constraint file
    line = open(contact_filename,'r').readline()
    if len(line.split(',')) != 1:
        sep = ','
    elif len(line.split(' ')) != 1:
        sep = ' '
    else:
        sep = '\t'

    if len(sys.argv) == 5:
        pdb_filename = sys.argv[3]
        rep_len = int(sys.argv[4])
        plot_map(contact_filename, psipred_filename, pdb_filename, rep_len, sep)
    else:
        rep_len = int(sys.argv[3])
        print "sep=" + sep
        plot_map_nopdb(contact_filename, psipred_filename, rep_len, sep)


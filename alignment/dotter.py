import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform

from Bio.SubsMat import MatrixInfo

from os.path import expanduser
home = expanduser("~")
sys.path.append(home + '/bioinfo-toolbox')

from parsing import parse_hmm


"""
TODO: MSP length approximation by Karlin and Altschul
def calc_MSP(seq):

    b62 = MatrixInfo.blosum62
    N = len(seq)
    MSP_global = np.log(N * N)
    MSP_res =
    
    return MSP_global / MSP_res
"""


def calc_dot_matrix(seq):
    """
     calculates similarity score matrix for given sequence to itself
     according to the DOTTER dot-plot program
     Reference: Erik L.L. Sonnhammer, Richard Durbin - "A dot-matrix program 
            with dynamic threshold control suited for genomic DNA and protein
            sequence analysis" - Gene 167 (1996) 1-10
    """

    aa_list = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M',
           'F', 'P', 'S', 'T', 'W', 'Y', 'V'] #, '-', 'B', 'Z', 'J']
    b62 = MatrixInfo.blosum62

    N = len(seq)
    alpha = 20
    W = 10
    
    score_vec = {}
    for aa_1 in aa_list:
        score_list = []
        for aa_2  in seq:
            if (aa_1, aa_2) in b62:
                score_list.append(b62[(aa_1, aa_2)])
            else:
                score_list.append(b62[(aa_2, aa_1)])
        score_vec[aa_1] = np.array(score_list)

    newsum = np.zeros(N)
    oldsum = np.zeros(N)
    zero_vec = np.zeros(N)

    dot_matrix = np.zeros((N, N))

    for i in range(N):
        tmp = oldsum
        oldsum = newsum
        newsum = tmp

        add_vec = score_vec[seq[i]]
        if i > W:
            del_vec = score_vec[seq[i - W]]
        else:
            del_vec = zero_vec

        newsum[0] = add_vec[0]
        for j in range(1,W):
            newsum[j] = oldsum[j - 1] + add_vec[j]
        for j in range(W, N):
            newsum[j] = oldsum[j - 1] + add_vec[j] - del_vec[j - W]
            if newsum[j] > 0 and i > W:
                score = newsum[j] / float(W)
                dot_matrix[i - W/2, j - W/2] = score
     
    plt.imshow(dot_matrix, origin='lower', cmap=cm.binary)
    plt.show()
    return dot_matrix
    #print np.max(dot_matrix)
    #print np.min(dot_matrix)



def calc_dot_matrix_profile(seq, profile):
    """ calculates similarity score matrix for given profile hmm to itself
        according to the DOTTER dot-plot program
        Reference: Erik L.L. Sonnhammer, Richard Durbin - "A dot-matrix program 
            with dynamic threshold control suited for genomic DNA and protein
            sequence analysis" - Gene 167 (1996) 1-10
    """

    aa_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P',
           'Q', 'R', 'S', 'T', 'V', 'W', 'Y'] #, '-', 'B', 'Z', 'J']
    b62 = MatrixInfo.blosum62

    N = len(seq)
    alpha = 20
    W = 15
    
    score_vec = {}
    for aa_1 in aa_list:
        score_list = []
        for (aa_2, freq_list) in profile:
            tmp_score = 0.0
            for i in range(len(freq_list[:-1])):
                aa_3 = aa_list[i]
                hmm_prob = freq_list[i]
                if hmm_prob == '*':
                    hmm_prob = 10000
                # hmm probabilities are stored in a log scale:
                # -1000 * log2(frequency)
                # see hhsuite usequide for more information
                # here we nee the reverse function to obtain probabilities:
                freq = 2 ** (float(hmm_prob) / -1000)
                if (aa_1, aa_3) in b62:
                    blosum_score = b62[(aa_1, aa_3)]
                    tmp_score += blosum_score * freq
                else:
                    blosum_score = b62[(aa_3, aa_1)]
                    tmp_score += blosum_score * freq
            score_list.append(tmp_score)
        score_vec[aa_1] = np.array(score_list)

    newsum = np.zeros(N)
    oldsum = np.zeros(N)
    zero_vec = np.zeros(N)

    dot_matrix = np.zeros((N, N))

    for i in range(N):
        tmp = oldsum
        oldsum = newsum
        newsum = tmp

        add_vec = score_vec[seq[i]]
        if i > W:
            del_vec = score_vec[seq[i - W]]
        else:
            del_vec = zero_vec

        newsum[0] = add_vec[0]
        for j in range(1,W):
            newsum[j] = oldsum[j - 1] + add_vec[j]
        for j in range(W, N):
            newsum[j] = oldsum[j - 1] + add_vec[j] - del_vec[j - W]
            if newsum[j] > 0 and i > W:
                score = newsum[j] / float(W)
                dot_matrix[i - W/2, j - W/2] = score
     
    plt.imshow(dot_matrix, origin='lower', cmap=cm.binary)
    plt.show()
    #print np.max(dot_matrix)
    #print np.min(dot_matrix)
    return dot_matrix



if __name__ == '__main__':

    #print 'Doing nothing.'
    seqfile = sys.argv[1]

    seq = ''
    with open(seqfile) as f:
        for l in f:
            if l.startswith('>'):
                continue
            else:
                seq += l.strip()

    if len(sys.argv) == 3:
        hmmfile = sys.argv[2]
        with open(hmmfile) as f:
            profile = parse_hmm.read_hmm(f)
        calc_dot_matrix_profile(seq, profile)
    else:
        calc_dot_matrix(seq)













import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform

sys.path.append('/bubo/sw/apps/bioinfo/biopython/1.59/tintin/lib/python')
from Bio.SubsMat import MatrixInfo


def get_char(filename):
    
    alifile = open(filename, 'r')
    aliarr = np.array([list(line.strip()) for line in alifile])

    return aliarr
    


def get_numeric(filename):

    aa_num_dict = {'A': 1, 'R': 2, 'N': 3, 'D': 4, 'C': 5, 'E': 6, 'Q': 7,
                   'G': 8, 'H': 9, 'I': 10, 'L': 11, 'K': 12, 'M': 13,
                   'F': 14, 'P': 15, 'S': 16, 'T': 17, 'W': 18, 'Y': 19,
                   'V': 20, '-': 21, 'B': 7, 'Z': 3, 'J': 11}

    alifile = open(filename, 'r')
    aliarr = np.array([list(line.strip()) for line in alifile])
    alinum = np.empty(aliarr.shape)

    for (i, j), letter in np.ndenumerate(aliarr):
        alinum[i, j] = aa_num_dict[letter]

    alifile.close()

    """
    print aliarr
    print alinum
    alicov = np.cov(alinum, rowvar=False)
    print alicov
    print alicov.shape
    alicor = np.corrcoef(alinum, rowvar=False)
    
    alidist = squareform(pdist(alinum.T, 'euclidean'))
    print alidist.shape

    alidistlog = np.log(alidist)
    
    plt.imshow(alicov, origin='lower')
    plt.savefig('%s_cov.png' % filename, bbox_inches='tight', dpi=300)
    plt.close()

    plt.imshow(alicor, origin='lower')
    plt.savefig('%s_cor.png' % filename, bbox_inches='tight', dpi=300)
    plt.close()

    plt.imshow(alidist, origin='lower')
    plt.savefig('%s_dist.png' % filename, bbox_inches='tight', dpi=300)
    plt.close()

    plt.imshow(alidistlog, origin='lower')
    plt.savefig('%s_distlog.png' % filename, bbox_inches='tight', dpi=300)
    plt.close()
    
    plt.imshow(alidistlog, origin='lower', cmap=cm.Spectral)
    plt.show()
    """

    return alinum


def calc_blosum(aliarr):

    b62 = MatrixInfo.blosum62
    
    (N, L) = aliarr.shape
    
    print (N, L)


    seqblos = np.zeros((L, L))

    local_set = set

    for i in range(L):
        for j in range(i):
            f_ij = {}
            ali_ij = zip(aliarr[:,i], aliarr[:,j])
            #ali_ij.sort(key=lambda x: tuple(x))
            tmp = map(local_set, ali_ij)

            """
            for n in range(N):
                aa_i = aliarr[n,i]
                aa_j = aliarr[n,j]
                if aa_i == '-' or aa_j == '-':
                    continue
                if (aa_i, aa_j) in b62:
                    seqblos[i,j] += b62[(aa_i, aa_j)]
                else:
                    seqblos[i,j] += b62[(aa_j, aa_i)]
            """

    plt.imshow(seqblos, origin='lower', cmap=cm.binary)
    plt.show()
   

def calc_blosum_seq(aliarr):

    b62 = MatrixInfo.blosum62
    seq = aliarr[0,:]
    L = len(seq)
    seqblos = np.zeros((L, L))
    seqblos2 = np.zeros((L, L))

    for i in range(L):
        aa_i = seq[i]
        for j in range(i+1):
            aa_j = seq[j]
            if aa_i == '-' or aa_j == '-':
                continue
            if (aa_i, aa_j) in b62:
                seqblos[i,j] += b62[(aa_i, aa_j)]
                seqblos2[i,j] += seqblos[i,j] * seqblos[i,j]
            else:
                seqblos[i,j] += b62[(aa_j, aa_i)]
                seqblos2[i,j] += seqblos[i,j] * seqblos[i,j]
       
    seqblos_var = np.zeros((L, L))
    for i in range(L):
        aa_i = seq[i]
        for j in range(i+1):
            aa_j = seq[j]
    

    plt.imshow(seqblos, origin='lower', cmap=cm.binary)
    plt.show()


""" 
TODO: MSP length approximation by Karlin and Altschul
def calc_MSP(seq):

    b62 = MatrixInfo.blosum62
    N = len(seq)
    MSP_global = np.log(N * N)
    MSP_res =
    
    return MSP_global / MSP_res
"""


# calculates similarity score matrix for given sequence to itself
# according to the DOTTER dot-plot program
# Reference: Erik L.L. Sonnhammer, Richard Durbin - "A dot-matrix program 
#            with dynamic threshold control suited for genomic DNA and protein
#            sequence analysis" - Gene 167 (1996) 1-10
def calc_dotter(seq):

    aa_list = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M',
           'F', 'P', 'S', 'T', 'W', 'Y', 'V'] #, '-', 'B', 'Z', 'J']
    b62 = MatrixInfo.blosum62

    N = len(seq)
    alpha = 20
    W = 25
    
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
        for j in range(1,N):
            newsum[j] = oldsum[j - 1] + add_vec[j]
        for j in range(W, N):
            newsum[j] = oldsum[j - 1] + add_vec[j] - del_vec[j - W]
            if newsum[j] > 0 and i > W:
                score = newsum[j] / float(W)
                if score < 1:
                    score = 0
                if score > 2:
                    score = 2
                dot_matrix[i - W/2, j - W/2] = score
     
    return dot_matrix
    #plt.imshow(dot_matrix, origin='lower', cmap=cm.binary)
    #plt.show()
    print np.max(dot_matrix)
    print np.min(dot_matrix)



if __name__ == '__main__':

    aliarr = get_char(sys.argv[1])
    #calc_blosum_seq(aliarr)
    calc_dotter(aliarr[0])














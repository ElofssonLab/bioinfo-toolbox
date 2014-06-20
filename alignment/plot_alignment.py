import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def plot(filename):
    alifile = open(filename, 'r')
    #tmp = open('.'.join(filename.split('.')[:-1]) + '.fas', 'w')
    N = 0
    L = 0
    for line in alifile:
        N += 1
        #tmp.write('>seq:%s\n' % N)
        #tmp.write(line)
        if L == 0:
            for chr in line:
                L += 1
    alifile.close()
    #tmp.close()
    aliarr = np.zeros((N, L))
    alifile = open(filename, 'r')
    i = 0
    for line in alifile:
        j = 0
        for chr in line:
            if chr != '-':
                aliarr[i,j] = 1.0
            else:
                aliarr[i,j] = 0.0
            j += 1
        i += 1
    alifile.close()
    print aliarr.shape
    print aliarr
    plt.imshow(aliarr, cmap=cm.binary, aspect='auto')
    #plt.show()
    plt.savefig('%s.png' % filename, bbox_inches='tight', dpi=300)
    #plt.savefig('%s.pdf' % filename)


def get_frac_gaps(filename):
    alifile = open(filename, 'r')
    N = 0.
    for line in alifile:
        if not line.startswith('>'):
            N += 1.
    alifile.close()
    print N
    alifile = open(filename, 'r')
    frac_gaps = 0.
    for line in alifile:
        if not line.startswith('>'):
            ngaps = line.count('-')
            frac_gaps += ngaps/N
    alifile.close()
    return frac_gaps


if __name__ == "__main__":
    
    if sys.argv[1].endswith('.jones'):
        plot(sys.argv[1])
    elif sys.argv[1].endswith('.a3m'):
        print sys.argv[1] + ' ' + str(get_frac_gaps(sys.argv[1]))

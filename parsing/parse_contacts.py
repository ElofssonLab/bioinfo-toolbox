#!/usr/bin/env python
import numpy as np


def parse(afile, sep=' ', min_dist=4):
    
    """Parse contact file.
    @param  afile   contact file
    @param  sep     separator of contact file (default=' ')
    Ensures: Output is sorted by confidence score.
    @return [(score, residue a, residue b)]
    """

    contacts = []
    for aline in afile:
        if aline.strip() != '':
            line_arr = aline.strip().split(sep)
            if line_arr[0].startswith('E'):
                continue
            i = int(line_arr[0])
            j = int(line_arr[1])
            score = float(line_arr[-1])
            if abs(i - j) > min_dist:
                contacts.append((score, i, j))
    afile.close()

    contacts.sort(key=lambda x: x[0], reverse=True)
    return contacts


def get_numpy_cmap(contacts, seq_len=-1, sep=' '):

    """Convert contacts into numpy matrix.
    @param  contacts    contact list as obtained from "parse"
    @param  seq_len     sequence length
    @param  sep         separator of contact file (default=' ')
    @return np.array((seq_len, seq_len), score)
    """

    max_i = max(contacts, key=lambda(x):x[1])[1]
    max_j = max(contacts, key=lambda(x):x[2])[2]
    n = int(max(seq_len, max_i, max_j))
    cmap = np.zeros((n,n))

    for c in contacts:
        i = c[1] - 1
        j = c[2] - 1
        cmap[i,j] = c[0]
    
    return cmap



def write(contacts, outfile, sep=' '):

    """Write contact file.
    @param  contacts    contact list
    @param  outfile     output contact file
    @param  sep     separator of contact file (default=' ')
    """

    for c in contacts:
        outfile.write('%d%s%d%s%f\n' % (c[1], sep, c[2], sep, c[0]))


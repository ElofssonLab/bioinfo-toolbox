#!/usr/bin/env python

def pred(hfile):

    """Reads psipred output .ss2 file.
    @param  hfile   psipred .ss2 file
    @return disorder list
    """

    result = []
    for l in hfile:
        if l.startswith('#'):
            continue
        if not l.strip():
            continue
        l_arr = l.strip().split()
        result.append(float(l_arr[2]))
    return result



def seq(file):

    """Reads iupred output .horiz file.
    @param  file   iupred (long or short) output file
    @return sequencex
    """

    result = ''
    for line in file:
        if line.startswith('#'):
            continue
        if not line.strip():
            continue
        line_arr = line.strip().split(' ')
        result += line_arr[1]
    return result

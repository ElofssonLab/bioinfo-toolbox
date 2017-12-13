#!/usr/bin/env python

from __future__ import print_function

import sys

def extract_names(fname):
    with open(fname) as f:
        for l in f.readlines():
            l_arr = l.strip().split(',')
            uni1 = l_arr[2]
            uni2 = l_arr[5]
            pdb = l_arr[7]
            chain = l_arr[8]
            chain1 =pdb+"_"+chain[0:1]
            chain2 =pdb+"_"+chain[1:2]
            print(uni1,uni2,chain1,chain2)

if __name__ == '__main__':
    filename=sys.argv[1]
#    pdbname=sys.argv[2]
#    uniname=sys.argv[3]
    extract_names(filename)

#!/usr/bin/env python

import sys, os, argparse
import numpy as np
import re


parser = argparse.ArgumentParser(description="Converting CASP-RR file t gramm constrants")
#parser.add_argument('file', metavar='file', type=argparse.FileType('r'), nargs=1, help='filename')
parser.add_argument('--threshold',"-t", type=float, help='distance skip',default=2.)
#parser.add_argument("--sepseq","-sep","-S",required=False, help='Separation sequence between protein in MSA' ,default="GGGGGGGGGGGGGGGGGGGG")
parser.add_argument('-seq','--sequence','-s', required= False, help='sequence file of first sequene to identify domain borders')
parser.add_argument('file', metavar='file', type=str, nargs=1, help='filename')
parser.add_argument('-max','--maxhits','-m', required= False, help='maximum number of sequences to include')
args = parser.parse_args()


#print args


for infilef in args.file:
#    print infilef
    infile = open(infilef)

#if args.maxhits:
#    maxhits=args.maxhits
#else:
#    maxhits=1.e20
dompos=0
borders=[]
if args.sequence:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    #import A3MIO
    # We can read it as fasta as we only care about the first sequence withouth gaps
    import re
    with open(args.sequence, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq=record
    borders+=[len(seq)]
    
#print ("test",borders[0])
for l in infile:
    #if (counter > maxhits): break
    #if '>' in l and not counter == 0:
    #print (l)
    try:
        A,B,mini,maxi,prob,stdev=l.split()
        A=int(A)
        B=int(B)
        maxi=float(maxi)
        #print ("TEST",A,B,borders[0],maxi,args.threshold,)
        if (B>borders[0] and A<=borders[0]):
            print (A-1,B-borders[0]-1,maxi+args.threshold)
    except:
        #print (l)
        continue
    

#!/usr/bin/env python

import sys, os, argparse
import numpy as np
import re

def skipgaps(string,gaps):
    newstring=''
    for j in range(len(string)):
        if (gaps[j]==1):
            newstring+=string[j]
    return newstring


parser = argparse.ArgumentParser(description="Trimming extra characters in aligned sequence from an a3m file")
parser.add_argument('-o','--orgname', help='Keep original filenames', action="store_true")
#parser.add_argument('file', metavar='file', type=argparse.FileType('r'), nargs=1, help='filename')
parser.add_argument('-name', type=str, help='name')
parser.add_argument('--cutoff',"-c", type=float, help='cutoff',default=0.75)
parser.add_argument("--sepseq","-sep","-S",required=False, help='Separation sequence between protein in MSA' ,default="GGGGGGGGGGGGGGGGGGGG")
parser.add_argument('-seq','--sequence','-s', required= False, help='sequence file to identify domain baorders')
parser.add_argument('file', metavar='file', type=str, nargs=1, help='filename')
parser.add_argument('-max','--maxhits','-m', required= False, help='maximum number of sequences to include')
args = parser.parse_args()


cutoff=args.cutoff
#print args

#if args.orgname:
#        print "verbosity turned on"

for infilef in args.file:
#    print infilef
    infile = open(infilef)

#if args.maxhits:
#    maxhits=args.maxhits
#else:
#    maxhits=1.e20
dompos=0
if args.sequence:
    seqfile = open (args.sequence)
    for l in seqfile:
        if not ">" in l:
            l = l.strip()
            dompos+=len(l)
    
# Added functionality to remove gaps in first sequence..


maxlen=100000
counter = 0
gaps=np.zeros(maxlen)

seqname=""
for l in infile:
    #if (counter > maxhits): break
    if '>' in l and not counter == 0:
        if args.orgname:
            #sys.stdout.write('\n'+l)
            seqname=('\n'+l)
        else:
            #sys.stdout.write('\n>sequence{0:07d}/1-100\n'.format(counter))
            seqname='\n>sequence{0:07d}/1-100\n'.format(counter)
        #sys.stdout.write('>sequence{0:07d}/1-100\n'.format(counter))
        counter += 1
    elif '>' not in l:
        l = l.strip()
        upperseq = ''.join([c for c in l if not c.islower()])
        upperseq = upperseq.replace('X', '-')
        if counter == 1:
            for i in range(len(upperseq)):
                #print i,upperseq[i]
                if (upperseq[i]!="-"):
                    gaps[i]=1
                else:
                    gaps[i]=0
            counter += 1
        #print upperseq,gaps
        new=skipgaps(upperseq,gaps)
        count=new.count('-')
        frac=float(count)/float(len(new))
        if (frac<(1-cutoff)):
            #        if re.search("[A-Z]",new):
            #sys.stdout.write(upperseq)
            try:
                sys.stdout.write(seqname)
            except (BrokenPipeError, IOError):
                sys.exit()
            if dompos>0:
                try:
                    sys.stdout.write(new[:dompos]+args.sepseq+new[dompos:])
                except (BrokenPipeError, IOError):
                    sys.exit()
            else:
                try:
                    sys.stdout.write(new)
                except (BrokenPipeError, IOError):
                    sys.exit()
    elif '>' in l and counter == 0:
        if args.name:
            #sys.stdout.write(">"+args.name+" "+l+"\n")
            seqname=">"+args.name+" "+l+"\n"
        elif args.orgname:
            #sys.stdout.write(l)
            seqname=l
        else:
            #sys.stdout.write('>target/1-100\n')
            seqname='>target/1-100\n'

        counter += 1

try:
    sys.stdout.write('\n')
except (BrokenPipeError, IOError):
    sys.exit()

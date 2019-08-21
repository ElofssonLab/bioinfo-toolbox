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

cutoff=0.75
parser = argparse.ArgumentParser(description="Trimming extra characters in aligned sequence from an a3m file")
parser.add_argument('-o','--orgname', help='Keep original filenames', action="store_true")
#parser.add_argument('file', metavar='file', type=argparse.FileType('r'), nargs=1, help='filename')
parser.add_argument('-name', type=str, help='name')
parser.add_argument('file', metavar='file', type=str, nargs=1, help='filename')
args = parser.parse_args()

#print args

#if args.orgname:
#        print "verbosity turned on"

for infilef in args.file:
#    print infilef
    infile = open(infilef)


# Added functionality to remove gaps in first sequence..


maxlen=100000
counter = 0
gaps=np.zeros(maxlen)

seqname=""
for l in infile:
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
            sys.stdout.write(seqname)
            sys.stdout.write(new)
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
sys.stdout.write('\n')

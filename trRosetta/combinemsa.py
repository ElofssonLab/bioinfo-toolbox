#!/usr/bin/python3

import sys, getopt,re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#from Bio.SeqFeature import SeqFeature, FeatureLocation


file1=sys.argv[1]
file2=sys.argv[2]

type="fasta"
handle1 = open(file1, 'rU')
handle2 = open(file2, 'rU')


for record in SeqIO.parse(handle1, type ) :
    print (record.name)
    break
seq1=record


for record in SeqIO.parse(handle2, type ) :
    print (record.name)
    break
seq2=record


seq=seq1+seq2
print (seq)

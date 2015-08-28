#!/usr/bin/python

import sys, getopt,re
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#from Bio.SeqFeature import SeqFeature, FeatureLocation


file=sys.argv[1]

# Avoid having too many sequences in a single dir..
#subdir=

# Open GenBank file
handle = open(file, 'rU')

j=0
cutoff=100000
#print "opening "+ file +"\n"
# For each record (mitochrodrial genome, in this case)...
num=0
OutFile = open(file + ".split-" + str(j), 'w')
for record in SeqIO.parse(handle, 'fasta') :
   # Grab the entire sequence
   num=num+1
   if (num>cutoff):
      OutFile.close()
      OutFile = open(file + ".split-" + str(j), 'w')
      j=j+1
      num=0
   SeqIO.write(record, OutFile, "fasta")
OutFile.close()

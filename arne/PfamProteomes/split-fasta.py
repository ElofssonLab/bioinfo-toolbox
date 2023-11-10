#!/usr/bin/env python

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
handle = open(file, 'r')

j=0
cutoff=1
num=0
OutFile = open(file + ".split-" + str(cutoff) +"-"  + str(j) +".fasta", 'w')
InfoFile = open(file + ".split-" + str(cutoff) +"-"  + str(j) +".names", 'w')
for record in SeqIO.parse(handle, 'fasta') :
   # Grab the entire sequence
   num=num+1
   if (num>=cutoff):
      OutFile.close()
      InfoFile.close()
      OutFile = open(file + ".split-" + str(cutoff) +"-"  + str(j) +".fasta", 'w')
      InfoFile = open(file + ".split-" + str(cutoff) +"-"  + str(j) +".names", 'w')
      j=j+1
      num=0
   SeqIO.write(record, OutFile, "fasta")
   name= record.id
   subdirA=re.sub(r'.*\_','',name)[:2]
   subdirB=re.sub(r'.*\_','',name)[2:4]
   subdirC=re.sub(r'.*\_','',name)[4:6]
   subdirD=re.sub(r'.*\_','',name)[6:8]
   subdirE=re.sub(r'.*\_','',name)[8:10]
   subdirF=re.sub(r'.*\_','',name)[10:12]
   subdirG=re.sub(r'.*\_','',name)[12:14]
   subdir = subdirA+"/"+subdirB+"/"+subdirC+"/"+subdirD+"/"+subdirE+"/"+subdirF+"/"+subdirG
   InfoFile.write(subdir+"/"+name+".fa\n")
   
OutFile.close()
InfoFile.close()

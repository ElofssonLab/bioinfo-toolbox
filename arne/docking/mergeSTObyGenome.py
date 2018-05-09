#!/usr/bin/env python

import sys, getopt,re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#from Bio.SeqFeature import SeqFeature, FeatureLocation


fileA=sys.argv[1]
fileB=sys.argv[2]


handleA = open(fileA, 'rU')

dataA={}
dataB={}

#print "opening "+ fileA +"\n"
# For each record 
first=True
for record in SeqIO.parse(handleA, 'stockholm') :
   if first:
      seqA=record
      first=False
   else:
      organism= re.sub(r'.*TaxID=','',record.description)
      organism= re.sub(r'\s.*','',organism)
      if (not organism in dataA.keys()):
         #print (record.name,record.description,organism)
         dataA[organism]=record

handleB = open(fileB, 'rU')
#print "opening "+ fileB +"\n"        
first=True
for record in SeqIO.parse(handleB, 'stockholm') :
   if first:
      seqB=record
      first=False
   else:
      organism= re.sub(r'.*TaxID=','',record.description)
      organism= re.sub(r'\s.*','',organism)
      if (not organism in dataB.keys()):
         #        print record.name,organism
         dataB[organism]=record

# First we shoudl always use sequecne 1 in both files...

print "> " + seqA.name + " " + seqB.name
print seqA.seq+seqB.seq

for key in dataA.keys():
   if (key in dataB.keys()):
      print "> " + key 
      print dataA[key].seq+dataB[key].seq


#!/usr/bin/env python

from __future__ import print_function

import sys, getopt,re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

#from Bio.SeqFeature import SeqFeature, FeatureLocation

#sepseq="GGGGGGGGGGGGGGGGGGGG" # fir old reasons I keep polyA
sepseq="AAAAAAAAAAAAAAAAAAAA"


fileA=sys.argv[1]
fileB=sys.argv[2]
virustsv=sys.argv[3]

df=pd.read_csv(virustsv,sep="\t")

#print (df.loc[df["refseq id"]=="EU371560"]["host tax id"])
#sys.exit()

handleA = open(fileA, 'r')

dataA={}
dataB={}

# We modify to just use virushostdb from japan


#print ("opening "+ fileA +"\n")
# For each record 
first=True
for record in SeqIO.parse(handleA, 'stockholm') :
   if first:
      seqA=record
      first=False
   else:
      #foo=re.split(r'\|',record.description)
      #print (foo)
      geneid,name,organism,virus,host,refid,region,subregion=re.split(r'\|',record.description)

      try:
         organism=str(int((df.loc[df["refseq id"]==refid]["host tax id"].mean())))
      except:
         continue
      if (not organism in dataA.keys()):
         #print ("A",record.description,organism)
         dataA[organism]=record

handleB = open(fileB, 'r')
#print ("opening "+ fileB +"\n"        )
first=True
for record in SeqIO.parse(handleB, 'stockholm') :
   if first:
      seqB=record
      first=False
   else:
      organism= re.sub(r'.*OX=','',record.description)
      organism= re.sub(r'\s.*','',organism)

      if (not organism in dataB.keys()):
         #print ("B",record.description,organism)
         dataB[organism]=record

# First we shoudl always use sequecne 1 in both files...

print ("> " + seqA.name + " " + seqB.name)
print (seqA.seq+sepseq+seqB.seq)

#print (dataA.keys())
#print (dataB.keys())
for key in dataA.keys():
   if (key in dataB.keys()):
      print ("> REFID=" + key )
      print (dataA[key].seq+sepseq+dataB[key].seq)


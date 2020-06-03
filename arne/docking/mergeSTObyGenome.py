#!/usr/bin/env python
from __future__ import print_function

import sys, getopt,re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#from Bio.SeqFeature import SeqFeature, FeatureLocation

sepseq="AAAAAAAAAAAAAAAAAAAA"

fileA=sys.argv[1]
fileB=sys.argv[2]


handleA = open(fileA, 'r')

dataA={}
dataB={}

#print ("opening "+ fileA +"\n")
# For each record 
first=True
for record in SeqIO.parse(handleA, 'stockholm') :
   if first:
      seqA=record
      first=False
   else:
      if re.match(r'.*TaxID=',record.description):
         organism= re.sub(r'.*TaxID=','',record.description)
         organism= re.sub(r'\s.*','',organism)
      elif re.match(r'.*Tax=',record.description):
         organism= re.sub(r'.*Tax=','',record.description)
         organism= re.sub(r'\s.*','',organism)
      elif re.match(r'.*OS=',record.description):
         organism= re.sub(r'.*OS=','',record.description)
         organism= re.sub(r'OX=.*','',organism)
      elif re.match(r'.*RepID=',record.description):
         organism= re.sub(r'.*RepID=','',record.description)
         organism= re.sub(r'.*\_','',organism)
         organism= re.sub(r'\s.*','',organism)
      else: # This is for NCBI annotations with
         temp=re.split(r'\|',record.description)
         try:
            geneid,name,host,organism,region=temp
         except:
            continue
            #print (re.split(r'\|',record.description))
            #
         #print (organism)
         #try: 
         #   geneid,name,host,organism=re.split(r'\|',record.description)
         #except:
         #   print ("skipping: ",record.name)
         #   continue
         
      if (not organism in dataA.keys()):
         #print ((record.name,record.description,organism))
         dataA[organism]=record


handleB = open(fileB, 'r')
#print ("opening "+ fileB +"\n"        )
first=True
for record in SeqIO.parse(handleB, 'stockholm') :
   if first:
      seqB=record
      first=False
   else:
      if re.match(r'.*TaxID=',record.description):
         organism= re.sub(r'.*TaxID=','',record.description)
         organism= re.sub(r'\s.*','',organism)
      elif re.match(r'.*Tax=',record.description):
         organism= re.sub(r'.*Tax=','',record.description)
         organism= re.sub(r'\s.*','',organism)
      elif re.match(r'.*OS=',record.description):
         organism= re.sub(r'.*OS=','',record.description)
         organism= re.sub(r'OX=.*','',organism)
      elif re.match(r'.*RepID=',record.description):
         organism= re.sub(r'.*RepID=','',record.description)
         organism= re.sub(r'.*\_','',organism)
         organism= re.sub(r'\s.*','',organism)
      else: # This is for NCBI annotations with 
         try: 
            geneid,name,host,organism,region=re.split(r'\|',record.description)
         except:
            #print ("skipping: ",record.name,re.split(r'\|',record.description))
            continue
      if (not organism in dataB.keys()):
         #        print (record.name,organism)
         dataB[organism]=record

# First we shoudl always use sequecne 1 in both files...

print ("> " + seqA.name + " AND " + seqB.name)
print (seqA.seq+sepseq+seqB.seq)

for key in dataA.keys():
   if (key in dataB.keys()):
      print ("> " + key )
      print (dataA[key].seq+sepseq+dataB[key].seq)


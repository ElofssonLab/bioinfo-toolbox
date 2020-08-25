#!/usr/bin/env python3
from __future__ import print_function

import sys, getopt,re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
from argparse import RawTextHelpFormatter
        
##args = docopt.docopt(__doc__)
#out_dir = args['--output_folder']
 

p = argparse.ArgumentParser(description = '- plotting trRosetta maps-',
                            formatter_class=RawTextHelpFormatter)
p.add_argument('-fileA','--inputA','-i', required= True, help='Input File A')
p.add_argument('-fileB','--inputB','-j', required= True, help='Input file B')
#p.add_argument('-dom','--domain','-d', required= False, help='positions of domain borders', nargs='+')
p.add_argument('-seq','--sepseq','-s', required= False, help='sequence file to identify domain baorders (20 Gly)',default="GGGGGGGGGGGGGGGGGGGG")
p.add_argument('--genus','-g', required= False, help='Use only genus names', action='store_true')
p.add_argument('--host','-l', required= False, help='use host to match fileB', action='store_true')
#p.add_argument('-out','--output','-o', required= False, help='output image')
#parser.add_argument('--nargs', nargs='+')
ns = p.parse_args()

#from Bio.SeqFeature import SeqFeature, FeatureLocation

sepseq=ns.sepseq

fileA=ns.inputA
fileB=ns.inputB

if (ns.genus):
   use_genus=True  # THis is a flag that we set to only match on Genus (i.e. first word in genus)
else:
   use_genus=False

if (ns.host):
   use_host=True # If this is set we use the hostname for virus proteins (from NCBI virues
else:
   use_host=False

#print ("TEST:",use_host,use_genus)
   
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
         if use_host:
            organism=host
            #print (organism)
         #try: 
         #   geneid,name,host,organism=re.split(r'\|',record.description)
         #except:
         #   print ("skipping: ",record.name)
         #   continue
      if use_genus:
         organism=re.sub(r'\s+','',organism)
      if (not organism in dataA.keys()):
         #print ("A:",organism,(record.name,record.description,organism))
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
         if use_host:
            organism=host
      if use_genus:
         organism=re.sub(r'\s+','',organism)
      if (not organism in dataB.keys()):
         #print ("B",record.name,organism)
         dataB[organism]=record

# First we shoudl always use sequecne 1 in both files...

print ("> " + seqA.name + " AND " + seqB.name)
print (seqA.seq+sepseq+seqB.seq)

for key in dataA.keys():
   if (key in dataB.keys()):
      print ("> " + key )
      print (dataA[key].seq+sepseq+dataB[key].seq)


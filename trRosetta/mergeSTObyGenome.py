#!/usr/bin/env python3
from __future__ import print_function

import sys, getopt,re

from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
from argparse import RawTextHelpFormatter
from Bio.Align import MultipleSeqAlignment        

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
p.add_argument('--skipped','-S', required= False, help='file to print skipped files (endings will be added)')
p.add_argument('--paralogs','-p', required= False, help='file to print paralogous files (endinfs will be added)')
p.add_argument('--output','-out','-o', required= False, help='output file (if missing stdout)')
#p.add_argument('-out','--output','-o', required= False, help='output image')
#parser.add_argument('--nargs', nargs='+')
p.add_argument('--max','-m', required= False, help='Maximum number of sequences to include',default=1000000)
#p.add_argument('-id', required= False, help='Maximum number of sequences to include',default=1.0,type=float) # This is way too slow..

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


#skippingA=MultipleSeqAlignment([])
#skippingB=MultipleSeqAlignment([])
#paralogA=MultipleSeqAlignment([])
#paralogB=MultipleSeqAlignment([])
skippingA=[]
skippingB=[]
paralogA=[]
paralogB=[]

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
         #if use_host:
         #   organism=host
         #   #print (organism)
         #try: 
         #   geneid,name,host,organism=re.split(r'\|',record.description)
         #except:
         ##   print ("skipping: ",record.name)
         #   skippingA+=[record]
         #   continue
      #gaps=len(re.findall(r'[A-Z]',str(record.seq)))/len(record.seq)
      print (record.seq,gaps)
      sys.exit()
      if use_genus:
         organism=re.sub(r'\s+','',organism)
      if (not organism in dataA.keys()):
         #print ("A:",organism,(record.name,record.description,organism))
         dataA[organism]=record
      else:
         paralogA+=[record]

         
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
            skippingB+=[record]
            continue
         if use_host:
            organism=host
      if use_genus:
         organism=re.sub(r'\s+','',organism)

      if (not organism in dataB.keys()):
         #print ("B",record.name,organism)
         dataB[organism]=record
      else:
         paralogB+=[record]

# First we should always use sequecne 1 in both files...

sequences=[]

sequences+=[SeqRecord(seq=seqA.seq+sepseq+seqB.seq,id=seqA.name+" "+seqB.name,name=seqA.name+" "+seqB.name,description=seqA.name+" "+seqB.name,dbxrefs=[])]

count=0
if (not ns.output):
   print ("> " + seqA.name + " AND " + seqB.name)
   print (seqA.seq+sepseq+seqB.seq)
for key in dataA.keys():
   if ((key in dataB.keys()) and  (dataA[key].id!=dataB[key].id)):
      sequences+=[SeqRecord(seq=dataA[key].seq+sepseq+dataB[key].seq,id=key,name=key,description=key,dbxrefs=[])]

      if count < ns.max:
         if (not ns.output):
            print ("> " + key )
            print (dataA[key].seq+sepseq+dataB[key].seq)
            count+=1
      else:
         skippingA+=[SeqRecord(seq=dataA[key].seq,id=key,name=key,description=key,dbxrefs=[])]
         skippingB+=[SeqRecord(seq=dataB[key].seq,id=key,name=key,description=key,dbxrefs=[])]
   else:
      skippingA+=[SeqRecord(seq=dataA[key].seq,id=key,name=key,description=key,dbxrefs=[])]
for key in dataB.keys():
   if not (key in dataA.keys()):
      skippingB+=[SeqRecord(seq=dataB[key].seq,id=key,name=key,description=key,dbxrefs=[])]

      
selected=[]      
if (ns.output):
   msa=MultipleSeqAlignment(sequences)
   for record in msa:
      # Check protein id
      #id=0
      #selmsa=MultipleSeqAlignment(selected)
      #for sel in selmsa:
      #   j=0
      #   k=0
      #   #print ("test",sel)
      #   for i in range(len(sel.seq)):
      #      #print ("FOO",i,j,k,sel.seq[i],record.seq[i])
      #      if (sel.seq[i]!="-"):
      #         k+=1
      #         if (sel.seq[i]==record.seq[i]):
      #            j+=1
      #   ###print (sel.id,record.id,len(sel.seq),j,k,j/k)
      #   ##a = sel.seq
      #   ##b = record.seq
      #   ##k = len(sel.seq)
      #   ##print (sel.seq)
      #   ##k = sel.seq.findall('[A-Z]')
      #   ##print (k)
      #   #k = len(re.findall('[A-Z]',str(sel.seq)))
      #   #j = k-len(re.findall('[A-Z]',''.join(x for x, y in zip(sel.seq, record.seq) if x != y)))
      #   ##print (sel.seq)
      #   ##print (record.seq)
      #   ##print (j/k)
      #   id=max(id,j/k)
      if (count<ns.max):   
         selected+=[record]
         count+=1
      ##print ("max",record.name,id)
      ##print ("selected",selected)
   selmsa=MultipleSeqAlignment(selected)
   AlignIO.write(selmsa, ns.output+".fasta", "fasta")
   #AlignIO.write(msa, ns.output+".fasta", "fasta")


   
if (ns.skipped):
   msaA=MultipleSeqAlignment(skippingA)
   msaB=MultipleSeqAlignment(skippingB)
   AlignIO.write(msaA, ns.skipped+"-A.fasta", "fasta")
   AlignIO.write(msaB, ns.skipped+"-B.fasta", "fasta")
if (ns.paralogs):
   msaA=MultipleSeqAlignment(paralogA)
   msaB=MultipleSeqAlignment(paralogB)
   AlignIO.write(msaA, ns.paralogs+"-A.fasta", "fasta")
   AlignIO.write(msaB, ns.paralogs+"-B.fasta", "fasta")

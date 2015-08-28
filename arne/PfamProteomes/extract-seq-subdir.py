#!/usr/bin/python

import sys, getopt,re
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#from Bio.SeqFeature import SeqFeature, FeatureLocation


file=sys.argv[1]
dir=sys.argv[2]

# Avoid having too many sequences in a single dir..
#subdir=

# Open GenBank file
handle = open(file, 'rU')
count=0
#print "opening "+ file +"\n"
# For each record (mitochrodrial genome, in this case)...
for record in SeqIO.parse(handle, 'fasta') :
#   print record.name
   # Grab the entire sequence
#   seq = str(record.seq)
   name= re.sub(r'\|','-',str(record.name))

#   record.id=str(record.name)
   # Look at all features for this record
   #   for feature in record.features:
   subdirA=re.sub(r'.*\_','',name)[:2]
   subdirB=re.sub(r'.*\_','',name)[2:5]
   subdirC=re.sub(r'.*\_','',name)[5:8]

#   print subdirA,subdirB,subdirC

   subdir = subdirA+"/"+subdirB+"/"+subdirC
   if (not os.path.isdir(dir + "/" +subdir)): 
      os.makedirs(dir + "/" +subdir)
   if (not os.path.isfile(dir + "/" + subdir + "/" + name +  '.fa')): 
      try:
         OutFile = open(dir + "/" + subdir + "/" + name +  '.fa', 'w')
#      print "FILE: " + name + "\n"
         SeqIO.write(record, OutFile, "fasta")
         count=count+1
      except:
         print "error" + dir + "/" + subdir + "/" + name +  '.fa'
#   else:
#      print "skipping: " + name + "\n"
print "Wote number of files: ",file,count

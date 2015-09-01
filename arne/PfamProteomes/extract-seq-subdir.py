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

infofile = open(file + ".names", 'w')

count=0
#print "opening "+ file +"\n"
# For each record (mitochrodrial genome, in this case)...
for record in SeqIO.parse(handle, 'fasta') :
##   print record.name
   # Grab the entire sequence
#   seq = str(record.seq)
   name= re.sub(r'\|','-',str(record.name))

#   record.id=str(record.name)
   # Look at all features for this record
   #   for feature in record.features:
   subdirA=re.sub(r'.*\_','',name)[:2]
   subdirB=re.sub(r'.*\_','',name)[2:4]
   subdirC=re.sub(r'.*\_','',name)[4:6]
   subdirD=re.sub(r'.*\_','',name)[6:8]
   subdirE=re.sub(r'.*\_','',name)[8:10]
   subdirF=re.sub(r'.*\_','',name)[10:12]
   subdirG=re.sub(r'.*\_','',name)[12:14]

#   print subdirA,subdirB,subdirC

   subdir = subdirA+"/"+subdirB+"/"+subdirC+"/"+subdirD+"/"+subdirE+"/"+subdirF+"/"+subdirG
   if (not os.path.isdir(dir + "/" +subdir)): 
      os.makedirs(dir + "/" +subdir)
   infofile.write(dir + "/" + subdir + "/" + name +  '.fa' + "\n")
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


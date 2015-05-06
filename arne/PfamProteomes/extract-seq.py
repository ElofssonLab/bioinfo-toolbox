#!/usr/bin/python

import sys, getopt

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#from Bio.SeqFeature import SeqFeature, FeatureLocation


file=sys.argv[1]
dir=sys.argv[2]

# Open GenBank file
handle = open(file, 'rU')

print "opening "+ file +"\n"
# For each record (mitochrodrial genome, in this case)...
for record in SeqIO.parse(handle, 'fasta') :
   print record.name
   # Grab the entire sequence
#   seq = str(record.seq)
   name= str(record.name)
#   record.id=str(record.name)
   # Look at all features for this record
   #   for feature in record.features:
   OutFile = open(dir + "/" + name +  '.fa', 'w')
#   print "FILE: " + name + "\n"
   SeqIO.write(record, OutFile, "fasta")

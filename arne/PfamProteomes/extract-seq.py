#!/usr/bin/python

import sys, getopt,re

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

type='fasta'
if re.search(".gb$",file):
   type='gb'
elif re.search(".fa",file):   
   type='fasta'
   

print ("opening "+ file +"\n")
# For each record (mitochrodrial genome, in this case)...
for record in SeqIO.parse(handle, type ) :
   print (record.name)
   # Grab the entire sequence
#   seq = str(record.seq)
   name= re.sub(r'[\<\>\/\\\|]','-',str(record.name))
   print (record.name,name)
#   record.id=str(record.name)
   # Look at all features for this record
   #   for feature in record.features:
   OutFile = open(dir + "/" + name +  '.fa', 'w')
#   print ("FILE: " + name + "\n")
   SeqIO.write(record, OutFile, "fasta")
   if record.features:
      for feature in record.features:
         if feature.type == "CDS":
            #print ("Feature:",feature.qualifiers)
            #print (feature.qualifiers["translation"])
            #print (feature.gene)
            #print (feature.location)
            #print (feature.qualifiers["protein_id"])
            print (feature.qualifiers["label"])
            #print (feature.location.extract(record).seq)
            # label seems to be the only existig case
            try:
               name=feature.qualifiers["label"][0]
            except:
               name="Unknown Name"
            try:
               id=feature.qualifiers["protein_id"][0]
            except:
               id=name
            try:
               product=feature.qualifiers["product"][0]
            except:
               product=""
               
            newseq = SeqRecord(Seq(feature.qualifiers["translation"][0]),
                               name=name,
                               id=id,
                               description=product)
            name=feature.qualifiers["label"] [0]
            #print (newseq)
            OutFile = open(dir + "/" + name +  '.fa', 'w')
            SeqIO.write(newseq, OutFile, "fasta")                


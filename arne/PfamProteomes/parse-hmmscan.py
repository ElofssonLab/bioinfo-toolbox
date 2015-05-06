#!/usr/bin/python

import sys, getopt

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SearchIO
#from Bio.SearchIO import HmmerIO
#from Bio.SeqFeature import SeqFeature, FeatureLocation
from pprint import pprint


cutoff=0.001

PfamFile=sys.argv[1]
# directory with fasta files.
SeqFile=sys.argv[2]


# Open GenBank file
handle = open(SeqFile, 'rU')

print "opening "+ SeqFile +"\n"
# For each record (mitochrodrial genome, in this case)...
SeqRecord=SeqIO.parse(handle, 'fasta')


#print SeqRecord

#print "opening "+ PfamFile +"\n"
# For each record (mitochrodrial genome, in this case)...
for record in SearchIO.parse(PfamFile, 'hmmscan3-domtab') :
#for record in SearchIO.parse(file, 'hmmer3-text') :
   for target in record:
      print "\n\n **** target ****\n\n"
      print target
      print "\n\n **** info ****\n\n"

      print target.id
      print target.accession
      print target.evalue
      print target.description


      pprint (vars(target))

      for hit in target:
         if (hit.evalue < cutoff):
            print "\n\n **** hit ****\n\n"
            print hit
            print "\n\n **** info ****\n\n"
            pprint (vars(hit))


      #sys.exit()
RECORD={}
for entry in SeqRecord:
#   print entry.seq
#  print entry.id
   RECORD[entry.id]=entry

#for entry in RECORD:
#   print entry.seq
#   print entry.id


print RECORD["sp|P0ACS2|SOXR_ECOLI"]
#print RECORD



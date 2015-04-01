#!/usr/bin/python

import sys, getopt

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SearchIO
#from Bio.SearchIO import HmmerIO
#from Bio.SeqFeature import SeqFeature, FeatureLocation




file=sys.argv[1]
# directory with fasta files.
#dir=sys.argv[2]


print "opening "+ file +"\n"
# For each record (mitochrodrial genome, in this case)...
for record in SearchIO.parse(file, 'hmmscan3-domtab') :
#for record in SearchIO.parse(file, 'hmmer3-text') :
   for target in record:
      print "\n\n **** target ****\n\n"
      print target
      print target.id
      print target.accession
      print target.evalue
      
#      for qresult in record:
#         print qresult
#      for frag in record:
#         print frag
#   sys.exit()

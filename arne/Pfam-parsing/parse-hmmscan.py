#!/usr/bin/python

import sys, getopt

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SearchIO
#from Bio.SearchIO import HmmerIO
#from Bio.SeqFeature import SeqFeature, FeatureLocation


file=sys.argv[1]
# directory with fata files.
dir=sys.argv[2]

# Open GenBank file
handle = open(file, 'rU')

print "opening "+ file +"\n"
# For each record (mitochrodrial genome, in this case)...
for record in SearchIO(handle, 'hmmer3-domtab') :
   print record


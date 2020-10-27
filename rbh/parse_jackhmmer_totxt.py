import pandas as pd
from Bio import SearchIO
import sys

'''
parse jackhmmer tblout file to three-column txt file: queryID,targetID,evalue
'''

infile=sys.argv[1]

rows=[]
for qresult in SearchIO.parse(infile,'hmmer3-tab'):
    for i in range(len(qresult.hits)): 
        rows.append(qresult.id+' '+qresult.hits[i].id+' '+str(qresult.hits[i].evalue))

##overwrite input file with string
with open(infile,'w+') as f:
     for i in rows: 
         f.write(i+'\n')



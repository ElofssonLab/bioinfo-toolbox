#!/usr/bin/env python
import os, sys
import pandas as pd

data_file="data/ecod.latest.domains.tsv"
#data_file="foo.tsv"
data = pd.read_csv(data_file, sep='\t',na_values=['na'])
#n = data[data['pdb'] == '4ag1']

# Filter for resolution etc etc ..


for p in set(data['pdb']):
    pdb=data[data['pdb'] == p]
    numchain=len(set(pdb['chain']))
    for c in set(pdb['chain']):
        chain=pdb[pdb['chain'] == c]
        numdom=len(set(chain['ecod_domain_id']))
        if (numchain == 1 and numdom > 1):
            print "Domain: ",p,c,numchain,numdom
#        for d in set(chain['ecod_domain_id']):
#            domain=chain[chain['ecod_domain_id'] == d ]
#                #print domain
                # print "TEST: ",p,c,d,numchain,numdom


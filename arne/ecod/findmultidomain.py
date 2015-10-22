#!/usr/bin/env python
import os, sys, re
import pandas as pd

data_file="data/ecod.latest.domains.tsv"
#data_file="foo.tsv"
data = pd.read_csv(data_file, sep='\t',na_values=['na'])
#n = data[data['pdb'] == '4ag1']

# Filter for resolution etc etc ..
fold={}

for p in set(data['pdb']):
    pdb=data[data['pdb'] == p]
    numchain=len(set(pdb['chain']))
    for c in set(pdb['chain']):
        chain=pdb[pdb['chain'] == c]
        numdom=len(set(chain['ecod_domain_id']))
# we only want entries with a single chain and multiple domains
        if (numchain == 1 and numdom > 1):
            new=True
            newfolds=[]
#    check if any fold has been used before
            for d in set(chain['ecod_domain_id']):
                domain=chain[chain['ecod_domain_id'] == d ]
#                print domain
                for f in set(domain['f_id']):
                    ecod=re.split('\.',f)
                    FOLD=ecod[0]+"."+ecod[1]+"."+ecod[2]
#                    (A,X,H,T)=re.split('\.',f)
#                    FOLD=A+"."+X+"."+T
#                    print FOLD
                    if FOLD in fold.keys():
                        new=False
                    newfolds.append(FOLD)
            for f in newfolds:
                fold[f]=f
            if (new):
                print "Domain: ",p,c,numchain,numdom,newfolds
                

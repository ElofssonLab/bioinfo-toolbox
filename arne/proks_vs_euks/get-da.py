#!/usr/bin/env python

import re
import sys
name=''
string=''
with open('../data/pfam/Pfam-A.regions.uniprot.selected_fields.sorted.tsv') as fp:
    for line in fp:
        #print (line)
        p = re.compile(r'\W+')
        s=p.split(line)
        #print (s[0],name)
        if s[0] != name:
            print (string)
            string=s[0]+" "+s[1]
            name=s[0]
        else:
            string+=" "+s[1]
print string

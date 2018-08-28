#!/usr/bin/env python
import os
import re
from multiprocessing import Pool
from Bio import SeqIO
from Bio import SwissProt
from Bio.SwissProt import KeyWList



handle = open("data/proteomes/UP000009229_760568.txt")
handle = open("data/proteomes/UP000001554_7739.txt")
#record = SwissProt.parse(handle)
#descriptions = [record.description for record in SwissProt.parse(handle)]
#len(descriptions)
for record in SwissProt.parse(handle):
    print (record.__dict__)
    print (record.accessions)
    #print (record.cross_references)
    for  db in record.cross_references:
        if ( db[0] == "Pfam"):
            print db


keywordfile="/pfs/nobackup/home/w/wbasile/annotate_uniprot_proteomes/bin/keywlist.txt"
handle = open(keywordfile)
for keyword in  KeyWList.parse(handle):
    print (keyword)

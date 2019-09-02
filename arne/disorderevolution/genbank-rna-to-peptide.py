#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
import re
from multiprocessing import Pool
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqUtils.ProtParamData import kd
from Bio.SeqUtils import ProtParam
from Bio.SeqUtils import GC
import subprocess
from Bio import SwissProt
import sys

def index_genbank_features(gb_record, feature_type, qualifier) :
    answer = dict()
    for (index, feature) in enumerate(gb_record.features) :
        if feature.type==feature_type :
            if qualifier in feature.qualifiers :
                #There should only be one locus_tag per feature, but there
                #are usually several db_xref entries
                for value in feature.qualifiers[qualifier] :
                    if value in answer :
                        print ("WARNING - Duplicate key %s for %s features %i and %i" \
                           % (value, feature_type, answer[value], index))
                    else :
                        answer[value] = index
    #return answer
    return value,index


filename=sys.argv[1]
record = next(SeqIO.parse(filename, "genbank"))
from Bio import SeqIO
for seq_record in SeqIO.parse(filename, "genbank"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    #print(len(seq_record))
    (aaseq,feature)=index_genbank_features(record,"CDS","translation")
    print (aaseq)
    start=record.features[feature]
    for (index, feature) in enumerate(start) :
        print (index,feature)
    print (start)
    #print (record.features[2])
    #xsys.exit()
    for (index, feature) in enumerate(record.features) :
        print (index,feature)
    #print(seq_record.features)
    sys.exit()

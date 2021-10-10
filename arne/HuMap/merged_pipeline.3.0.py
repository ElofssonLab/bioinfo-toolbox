#!/usr/bin/env python
# coding: utf-8

import sys
## add bioinfo-toolbox
#sys.path.append('/scratch2/wzhu/bioinfo-toolbox/')
#from parsing import parse_pdb
#from parsing import fasta
import parse_PDB_B

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import pandas as pd
import statistics
import pickle
import csv
from prody import *
from matplotlib.pylab import *
from Bio.SubsMat import MatrixInfo as matlist

peptide = sys.argv[1]
model = sys.argv[2]
uniprot = sys.argv[3]

##### Parse fasta file ########
fasta_sequences = SeqIO.parse(open(peptide),'fasta')
for fasta in fasta_sequences:
    petide_name, peptide_seq = fasta.id, str(fasta.seq)

##### Parse PDB file ########
with open(model, 'r') as pdb_file:
    for record in SeqIO.parse(pdb_file, 'pdb-atom'):
        models_seq =record.seq
        
res_dict = parse_PDB_B.get_res_dict(open(model, 'r'), "*")

#####get_ca_coordinates#######

ca_lst = []

# need to sort to get the sequence correct
sorted_keys = sorted(res_dict.keys())
    
for i in sorted_keys:
    ca_lst.append(res_dict[i][0])


##### Align ########
## ------alignment method: Julie use local alignments instead of global aln.
#alignments = pairwise2.align.globalms(peptide_seq, models_seq, 2, -1, -0.5, -0.1)
 
matrix = matlist.blosum62
alignments = pairwise2.align.localds(peptide_seq,models_seq,matrix,-10,-0.5)
atom_seq_ali = alignments[-1][0]
seq_ali = alignments[-1][1]

######## Prepare dssp #########
p = PDBParser()
structure = p.get_structure(str(model), model)
model1 = structure[0]
dssp = DSSP(model1, model, dssp='mkdssp')
a_key = list(dssp.keys())[2]

list_dssp = []

for a in dssp:
    list_dssp.append(a[2])

list_sasa = []

for a in dssp:
    list_sasa.append(a[3])
    
list_rsa = []    
list_ss = []

#### Prepar Pfam ########
try:
    Dict_pfam = prody.searchPfam(uniprot)
except (KeyError, UnboundLocalError,ValueError):
    pass

####### get aligned B-factor #####
j = 0
B_lst = []

for i in range(len(atom_seq_ali)):
    if atom_seq_ali[i] == '-':
        continue
        j += 1
    elif seq_ali[i] == '-':
        j += 1
        continue
    else:
        B_lst.append(ca_lst[i][3])
        list_ss.append(list_dssp[i])
        list_rsa.append(list_sasa[i])
        j += 1
        try:
            for pfamID in Dict_pfam:
                startPFAM = (Dict_pfam[pfamID]['locations'][0]['start'])
                endPFAM = (Dict_pfam[pfamID]['locations'][0]['end'])
                list_PFAM = []
                Pfam_range =range(int(startPFAM),int(endPFAM)+1)
                if int(startPFAM) in Pfam_range:
                    list_PFAM.append(Dict_pfam[pfamID]['accession'])
        except NameError:
            pass            

   
PFAM = "none"   
try:          
    PFAM = max(set(list_PFAM), key=list_ss.count)            
except NameError:
    pass

SS = max(set(list_ss), key=list_ss.count)
RSA_mean = statistics.mean(list_rsa)
mean = statistics.mean(B_lst)
#print( "%s , %s, %s , %s" % (peptide, mean, SS, RSA_mean))

data_dict = {"PeptideID":peptide,"AFscore":mean,"SS":SS, "RSA":RSA_mean, "PFAM":PFAM}

with open('Data_pipeline_TPP.2.csv', 'a') as file:
        w = csv.DictWriter(file, data_dict.keys())

        if file.tell() == 0:
            w.writeheader()

        w.writerow(data_dict)


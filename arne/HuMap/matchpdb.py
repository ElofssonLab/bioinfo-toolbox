#!/usr/bin/env python3
# coding: utf-8

import sys
#import parse_PDB_B

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqIO
from Bio.PDB import PDBParser
from Bio.PDB import Selection
from Bio.PDB import PDBIO
#from Bio.PDB import StructureBuilder
from Bio.PDB.Polypeptide import PPBuilder
import re
#from Bio.PDB.DSSP import DSSP
#import pandas as pd
#import statistics
#import pickle
#import csv
#from prody import *
#from matplotlib.pylab import *
from Bio.SubsMat import MatrixInfo as matlist

parser = PDBParser(PERMISSIVE=1)
name1=sys.argv[1]
name2=sys.argv[2]
outname=re.sub(".pdb","_reorder.pdb",name1)
outname1=re.sub(".pdb","_matching.pdb",name1)
outname2=re.sub(".pdb","_matching.pdb",name2)
structure1 = parser.get_structure("protein1", name1)
structure2 = parser.get_structure("protein2", name2)
seq1 = Selection.unfold_entities(structure1, "R")
seq2 = Selection.unfold_entities(structure2, "R")

ppb=PPBuilder()
seq1=[]
seq2=[]
for pp in ppb.build_peptides(structure1):
    seq1+=[pp.get_sequence()]
for pp in ppb.build_peptides(structure2):
    seq2+=[pp.get_sequence()]

    
    
##### Parse fasta file ########
#fasta_sequences = SeqIO.parse(open(peptide),'fasta')
#for fasta in fasta_sequences:
#    petide_name, peptide_seq = fasta.id, str(fasta.seq)

print (len(seq1[0]),len(seq1[1]),seq1)
print (len(seq2[0]),len(seq2[1]),seq2)
                
#print (record1.residue[0], record2)
#res_dict1 = parse_PDB_B.get_res_dict(open(protein1, 'r'), "*")
#res_dict2 = parse_PDB_B.get_res_dict(open(protein2, 'r'), "*")

#####get_ca_coordinates#######


#print (chain1,chain2)
    
matrix = matlist.blosum62


residue_to_remove1 = [[],[]]
residue_to_remove2 = [[],[]]

newpdb1=structure1
newpdb2=structure2
chain=0

alignments1 = pairwise2.align.localds(seq1[0],seq2[0],matrix,-10,-0.5,score_only=True)
alignments2 = pairwise2.align.localds(seq1[1],seq2[1],matrix,-10,-0.5,score_only=True)
alignments3 = pairwise2.align.localds(seq1[0],seq2[1],matrix,-10,-0.5,score_only=True)
alignments4 = pairwise2.align.localds(seq1[1],seq2[0],matrix,-10,-0.5,score_only=True)

# Check if the files are in reverse order
if (alignments3>alignments1 and alignments4>alignments2):
    print ("Renaming chain ids",alignments3,alignments4,alignments1,alignments2)
    chain1=[]
    tempid=["TMP1","TMP2","TMP3","TMP4","TMP5","TMP6"]
    for chain in structure1.get_chains():
        chain1+=[chain.id]
    i=0
    for chain in structure1.get_chains():
        chain.id=tempid[i]
        i+=1
    i=1
    for chain in structure1.get_chains():
        chain.id=chain1[i]
        i-=1
    # Get a list of the chains in a structure
    chains = list(structure1.get_chains())
    print (chains)
    newstruct=structure1.copy()
    chains = list(newstruct.get_chains())
    print (chains)
    for c in chains:
        c.detach_parent()
        newstruct[0].detach_child(c.id)
    #chains = list(newstruct.get_chains())
    #print (chains)
    newstruct.id=0
    newstruct.serial_num=0
    newstruct[0].add(chains[1])
    newstruct[0].add(chains[0])
    chains = list(newstruct.get_chains())
    print (chains)
    # Add it onto structure1
    #newstruct.add(chains[0])
    
    io = PDBIO()
    io.set_structure(newstruct)
    io.save(outname)
        
    structure1 = parser.get_structure("protein1", outname)
    seq1 = Selection.unfold_entities(structure1, "R")
    seq1=[]
    for pp in ppb.build_peptides(structure1):
        seq1+=[pp.get_sequence()]

    
chain1=[]
chain2=[]
for chain in structure1.get_chains():
    chain1+=[chain]
for chain in structure2.get_chains():
    chain2+=[chain]
print (chain1,chain2)
chain=0
for c in structure1.get_chains():
    k=0
    j=0
    print ("CHAIN:",chain)
    alignments = pairwise2.align.localds(seq1[chain],seq2[chain],matrix,-10,-0.5)
    ali1 = alignments[-1][0]
    ali2 = alignments[-1][1]
    
    print (ali1)
    print (ali2)
    for i in range(len(ali1)):
        if ali1[i] == '-':
            residue_to_remove2[chain]+=[k]
            k += 1
            continue
        elif ali2[i] == '-':
            residue_to_remove1[chain].append(j)
            #print ("TestRem",j,residue_to_remove1)
            j += 1
            continue
        else:
            j += 1
            k += 1
            #print (i,j,k,ali1[i],ali2[i])
        #print (structure1[0]["A"][i])
    chain+=1

c=0
print (residue_to_remove1)
print (residue_to_remove2)
remove1=[[],[]]
remove2=[[],[]]
for chain in structure1[0]:
    i=0
    j=0
    if (len(residue_to_remove1[c])>0):
        for residue in chain:
            #res=residue.get_full_id()
            #resid=int(res[3][1])
            print ("test",c,i,j)
            if j>=len(residue_to_remove1[c]): continue
            if (i==residue_to_remove1[c][j]):
                print ("Removing1",c,i,j,chain,residue)
                #newpdb1[0][chain.id].detach_child(residue.id)
                #chain.detach_child(residue.id)
                remove1[c]+=[residue]
                j+=1
            i+=1
    c+=1
c=0
for chain in structure2[0]:
    i=0
    j=0
    if (len(residue_to_remove2[c])>0):
        for residue in chain:
            #res=residue.get_full_id()
            #resid=int(res[3][1])
            if j>=len(residue_to_remove2[c]): continue
            if (i==residue_to_remove2[c][j]):
                print ("Removing2",c,i,j,chain,residue)
                #newpdb2[0][chain.id].detach_child(residue.id)
                #chain.detach_child(residue.id)
                remove2[c]+=[residue]
                j+=1
            i+=1
    c+=1 
print (remove1)
print (remove2)
c=0
for chain in structure1[0]:
    for res in remove1[c]:
        chain.detach_child(res.id)
    c+=1
c=0
for chain in structure2[0]:
    for res in remove2[c]:
        chain.detach_child(res.id)
    c+=1
io = PDBIO()

#structure1.renumber_residues()
#structure2.renumber_residues()
        
io.set_structure(structure1)
io.save(outname1)
io.set_structure(structure2)
io.save(outname2)


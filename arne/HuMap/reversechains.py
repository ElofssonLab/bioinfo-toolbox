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
#name2=sys.argv[2]
outname=re.sub(".pdb","_reorder.pdb",name1)
#outname1=re.sub(".pdb","_matching.pdb",name1)
#outname2=re.sub(".pdb","_matching.pdb",name2)
structure1 = parser.get_structure("protein1", name1)
#structure2 = parser.get_structure("protein2", name2)
seq1 = Selection.unfold_entities(structure1, "R")
#seq2 = Selection.unfold_entities(structure2, "R")

ppb=PPBuilder()
seq1=[]
#seq2=[]
for pp in ppb.build_peptides(structure1):
    seq1+=[pp.get_sequence()]
#for pp in ppb.build_peptides(structure2):
#    seq2+=[pp.get_sequence()]

    

#
#chains = list(structure1.get_chains())
#print (chains)
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
id1=chains[0].id
id2=chains[1].id
chains[1].id=id1
chains[0].id=id2
newstruct[0].add(chains[1])
newstruct[0].add(chains[0])

chains = list(newstruct.get_chains())
print (chains)
# Add it onto structure1
#newstruct.add(chains[0])
    
io = PDBIO()
io.set_structure(newstruct)
io.save(outname)


#!/usr/bin/env python3

import sys
#import parse_PDB_B

from Bio import SeqIO
from Bio.PDB import MMCIFParser
from Bio.PDB import Selection
from Bio.PDB import MMCIFIO
from Bio.PDB import PDBIO
from Bio.PDB import Select
from Bio.PDB.Polypeptide import PPBuilder
import re

parser = MMCIFParser()
name=sys.argv[1]
outname=re.sub(".cif","",name)
structure = parser.get_structure("protein", name)

class ChainSelect(Select):
    def __init__(self,chain):
        self.chain=chain
    def accept_chain(self,chain):
        if (re.search(chain.get_id(),self.chain)):
            return 1
        else:
            return 0

io=PDBIO()
io.set_structure(structure)
for chain1 in structure.get_chains():
    for chain2 in structure.get_chains():
        if chain2>chain1:
            chains=str(chain1.id)+str(chain2.id)
            io.save(outname+"_"+chains+".pdb",ChainSelect(chains))

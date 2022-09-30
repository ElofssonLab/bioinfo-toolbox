#!/usr/bin/env python3
import argparse
import re
import io
import sys

def print_to_string(*args, **kwargs):
    output = io.StringIO()
    print(*args, file=output, **kwargs)
    contents = output.getvalue()
    output.close()
    return contents

#from Bio.PDB.vectors import rotaxis, calc_angle, calc_dihedral
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import Selection
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB import Selection
from Bio.PDB.PDBIO import PDBIO
d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
 'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


chains=["A","B","C","D","E","F","G","H","I"]
chainsrev=["B","A"]
#chainsrev=chains.reverse()

if __name__ == "__main__":
    arg_parser = argparse.\
        ArgumentParser(
                description="Convert a pdb/mcif to trRosetta distances/angles")

    in_group = arg_parser.add_mutually_exclusive_group(required=True)
    in_group.add_argument("-p", "--pdb_file", type=argparse.FileType('r'))
    in_group.add_argument("-m", "--mmCIF_file", type=argparse.FileType('r'))
    arg_parser.add_argument("-r","--reverse",required=False,action="store_true",help="Reverse the chain order")
    arg_parser.add_argument("-g","--gap",required=False,default="GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG",type=str,help="sequence separating the two chains")
    args=arg_parser.parse_args()

    if args.pdb_file:
        from Bio.PDB.PDBParser import PDBParser
        bio_parser = PDBParser(PERMISSIVE=1)
        structure_file = args.pdb_file
        structure_id = args.pdb_file.name[:-4]
    else:
        from Bio.PDB.MMCIFParser import MMCIFParser
        bio_parser = MMCIFParser()
        structure_file = args.mmCIF_file
        structure_id = args.mmCIF_file.name[:-4]

    # Load structure
    structure = bio_parser.get_structure(structure_id, structure_file)

    # Get residues and length of protein
#    residues = []
#    resnum = []
#    for chain in structure[0]:
#        for residue1 in structure[0][chain.id]:
#            if not is_aa(residue1):
#                continue
#            residues.append(residue1.get_resname())
#            resnum.append(residue1.get_resname())
#    plen = len(residues)

polyinsert=args.gap
seq = ''
for model in structure:
    for chain in model:
        for residue in chain:
            seq+=d3to1[residue.resname]
#print('>some_header\n',''.join(seq))
#print(seq,polyinsert)
#print(re.search(polyinsert, seq, re.IGNORECASE))
#seq=[]
sequence=seq.split(polyinsert)
numchains=len(sequence)
#print(seq,numchains)
if (numchains == 1):
    sys.exit()

#print (seq[0],seq2)
#sys.exit()

if args.reverse: # Only works for dimers
    CHAIN=chainsrev[0]
    firstchain=False
else:
    CHAIN=chains[0]
    firstchain=True
    #skiplen=200
resid=0


i=0
c=1

leninsert=len(polyinsert)
skip=0

last=0
lastres=[-1*leninsert]
#print (numchains,sequence)
for i in range(1,numchains):
    lastres.append(len(sequence[i])+last)
    last=lastres[i]+leninsert
    #print (i,lastres,last)
lastres.append(len(seq))
lastres.append(len(seq))
#print (lastres)
for model in structure:
    for chain in model:
        #print (chain)
        for residue in chain:
            resid+=1
            #print (residue,residue.get_id()[1],lastres[c],leninsert,chain.id,CHAIN,skip)
            if (residue.get_id()[1]==(lastres[c]+1)):
                #print ("TEST1")
                skip=residue.get_id()[1] # +leninsert # -lastres
                i=0
                #skiplen+=20000
                if args.reverse:
                    CHAIN=chainsrev[c]
                    firstchain=True
                else:
                    print ("TER")
                    CHAIN=chains[c]
                resid=1
                c=c+1
            elif (residue.get_id()[1]>=(lastres[c-1]+leninsert) and residue.get_id()[1]<lastres[c+1]):
                #print ("TEST3")
                if (firstchain):
                    for atom in residue:
                        i+=1
                        print("{:6s}{:5d}  {:4s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}".format("ATOM",i,atom.id,residue.get_resname(),CHAIN,residue.get_id()[1]-skip,"",atom.get_coord()[0],atom.get_coord()[1],atom.get_coord()[2],1.,atom.get_bfactor()))
            else:
                #print ("TEST2",c,residue.get_id()[1],lastres[c],leninsert,lastres[c-1])
                skip

print ("TER")
if (args.reverse):
    i=0
    resid=0
    CHAIN="B"
    skip=0
    for model in structure:
        for chain in model:
            #print (chain)
            for residue in chain:
                resid+=1
                #print (residue,residue.get_id()[1],skiplen,lastres,chain.id,CHAIN)
                if (residue.get_id()[1]==lastres[c]+1):
                    skip=residue.get_id()[1]-1+leninsert # -lastres
                    i=0
                    print ("TER")
                    #skiplen+=20000
                    break
                else:
                    if (firstchain):
                        for atom in residue:
                            i+=1
                            print("{:6s}{:5d}  {:4s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}".format("ATOM",i,atom.id,residue.get_resname(),CHAIN,residue.get_id()[1]-skip,"",atom.get_coord()[0],atom.get_coord()[1],atom.get_coord()[2],1.,atom.get_bfactor()))

print ("END")

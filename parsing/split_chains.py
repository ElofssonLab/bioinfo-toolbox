#!/usr/bin/env python3
import argparse
import re
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
polygly="GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"

if __name__ == "__main__":
    arg_parser = argparse.\
        ArgumentParser(
                description="Convert a pdb/mcif to trRosetta distances/angles")

    in_group = arg_parser.add_mutually_exclusive_group(required=True)
    in_group.add_argument("-p", "--pdb_file", type=argparse.FileType('r'))
    in_group.add_argument("-m", "--mmCIF_file", type=argparse.FileType('r'))
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


seq = ''
for model in structure:
    for chain in model:
        for residue in chain:
            seq+=d3to1[residue.resname]
#print('>some_header\n',''.join(seq))
#print(seq,polygly)
#print(re.search(polygly, seq, re.IGNORECASE))
seq1,seq2=seq.split(polygly)
#print (seq1,seq2)
#sys.exit()

i=0
lastres=len(seq1)
lengly=len(polygly)
skip=0
CHAIN="A"
skiplen=200
resid=0
for model in structure:
    for chain in model:
        #print (chain)
        for residue in chain:
            resid+=1
            #print (residue,residue.get_id()[1],skiplen,lastres,chain.id,CHAIN)
            if (residue.get_id()[1]==lastres+1):
                skip=residue.get_id()[1]-1+lengly # -lastres
                i=0
                print ("TER")
                skiplen+=20000
                CHAIN="B"
                resid=1
            elif (residue.get_id()[1]<=lastres+lengly and residue.get_id()[1]>lastres):
                skip
            else:
                for atom in residue:
                    i+=1
                    print("{:6s}{:5d}  {:4s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}".format("ATOM",i,atom.id,residue.get_resname(),CHAIN,residue.get_id()[1]-skip,"",atom.get_coord()[0],atom.get_coord()[1],atom.get_coord()[2],1.,atom.get_bfactor()))

print ("TER")
print ("END")

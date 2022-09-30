#!/usr/bin/env python3
import argparse
#from Bio.PDB.vectors import rotaxis, calc_angle, calc_dihedral
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB import Selection
from Bio.PDB.PDBIO import PDBIO

chains=["A","B","C","D","E","F","G","H","I"]
chainsrev=["B","A"]
#chainsrev=chains.reverse()


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
                description="Convert a pdb/mcif to trRosetta distances/angles")

    in_group = arg_parser.add_mutually_exclusive_group(required=True)
    in_group.add_argument("-p", "--pdb_file", type=argparse.FileType('r'))
    in_group.add_argument("-m", "--mmCIF_file", type=argparse.FileType('r'))
    arg_parser.add_argument("-s","--skiplen",required=False,default=200,type=int,help="Skip length separating chains (default=200)")
    arg_parser.add_argument("-r","--reverse",required=False,action="store_true",help="Reverse the chain order")
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
    #print (structure_id,structure_file)
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

c=0
i=0
lastres=0
skip=0
if args.reverse:
    CHAIN="B"
    firstchain=False
    FIRSTCHAIN="B"
else:
    CHAIN=chains[c]
    FIRSTCHAIN="A"
    firstchain=True
#skiplen=200
skiplen=args.skiplen
resid=0
for model in structure:
    for chain in model:
        #print (chain)
        for residue in chain:
            #print (residue,residue.get_id()[1],skiplen,lastres,chain.id,CHAIN)
            if (residue.get_id()[1]-skiplen>lastres): # or chain.id!=CHAIN):
                if (residue.get_id()[1]-skiplen>lastres):
                    print ("TER")
                    #skiplen+=20000
                    firstchain=True
                    #resid=1
                    skip=residue.get_id()[1]-1 # -lastres
                i=0
                c+=1
                if args.reverse:
                    CHAIN="A"
                else:
                    CHAIN=chains[c]
                
            lastres=residue.get_id()[1]
            if (firstchain):
                resid+=1
                for atom in residue:
                    i+=1
                    print("{:6s}{:5d}  {:4s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}".format("ATOM",i,atom.id,residue.get_resname(),CHAIN,residue.get_id()[1]-skip,"",atom.get_coord()[0],atom.get_coord()[1],atom.get_coord()[2],1.,atom.get_bfactor()))
            #if (args.reverse):
            #    lastres=0
                

print ("TER")
                        
if (args.reverse):
    i=0
    resid=0
    CHAIN="B"
    FIRSTCHAIN="A"
    skip=0
    skiplen=args.skiplen
    for model in structure:
        for chain in model:
            #print (chain)
            for residue in chain:
                resid+=1
                #print (residue,residue.get_id()[1],skiplen,lastres,chain.id,CHAIN)
                if (chain.id!=FIRSTCHAIN or residue.get_id()[1]-skiplen>lastres):
                    break
                for atom in residue:
                    i+=1
                    print("{:6s}{:5d}  {:4s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}".format("ATOM",i,atom.id,residue.get_resname(),CHAIN,residue.get_id()[1]-skip,"",atom.get_coord()[0],atom.get_coord()[1],atom.get_coord()[2],1.,atom.get_bfactor()))
                lastres=residue.get_id()[1]

    print ("TER")
print ("END")

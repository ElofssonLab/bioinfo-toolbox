#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import sys
import re
import argparse
from argparse import RawTextHelpFormatter
from Bio.PDB import *
from Bio import SeqIO
from Bio.Seq import Seq

##args = docopt.docopt(__doc__)
#out_dir = args['--output_folder']
 

p = argparse.ArgumentParser(description = 'Extracting pdb chains and sequences of those from a PDB file',
                            formatter_class=RawTextHelpFormatter)
p.add_argument('-pdb','--pdb','-p', required= True, help='Pdb file for analysis od distances')
#p.add_argument('-seqdir','--seqdir','-s', required= False, help='Output directory for sequences (default ./)')
#p.add_argument('-out','--output','-o', required= False, help='output dir for pdb fdiles (default ./)')
#parser.add_argument('--nargs', nargs='+')
ns = p.parse_args()

p = PDBParser()
structure = p.get_structure('', ns.pdb)

code=ns.pdb.replace(".pdb","_")

three2one = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D',
             'CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
             'HIS':'H','ILE':'I','LEU':'L','LYS':'K',
             'MET':'M','PHE':'F','PRO':'P','SER':'S',
             'THR':'T','TRP':'W','TYR':'Y','VAL':'V',
             'MSE':'M'}


class NotDisordered(Select):
    def accept_atom(self, atom):
        return not atom.is_disordered() or atom.get_altloc() == "A"

class NotAlt(Select):
    def accept_atom(self, atom):
        return not atom.get_altloc() == "A"

class NotHet(Select):
    def accept_res(self, res):
        return not res.get_resname()[0:1] == "H_"
    
for chain in structure[0]: 
    print (chain,len(chain))
    
    io=PDBIO() 
    io.set_structure(chain) 
    io.save(code+chain.id+".pdb", select=NotAlt())
    seq=Seq('')
    for s in chain.get_residues():
        #print (s,s.get_full_id()[3][0])
        if s.get_full_id()[3][0]==' ':
            seq+=three2one[s.get_resname()]
    print (seq)
    sequence=Seq(seq)
    sequence.id=code+chain.id
    # Using CA-CA 
    ppb = CaPPBuilder() 
    #with open(code+chain.id+".fasta", "w") as output_handle:
    #    SeqIO.write(seq, output_handle, "fasta")

    #for pp in ppb.build_peptides(chain): 
    #    seq=Seq(pp.get_sequence())
    #    seq.id=code+chain.id
    #    print (seq,seq.id)
    #    #with open(code+chain.id+".fasta", "w") as output_handle:
    #    #    SeqIO.write(seq, output_handle, "fasta")
    output_handle=open(code+chain.id+".fasta", "w")
    output_handle.write(">"+code+chain.id+"\n")
    output_handle.write(str(seq)+"\n")
    output_handle.close()

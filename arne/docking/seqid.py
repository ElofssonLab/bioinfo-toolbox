#!/usr/bin/env python2

import sys
import argparse
import Bio.PDB
from Bio import pairwise2
from os.path import expanduser
home = expanduser("~")
sys.path.append(home + '/git/bioinfo-toolbox')

from parsing import parse_contacts
from parsing import parse_fasta
from parsing import parse_pdb


if __name__ == "__main__":
    p = argparse.ArgumentParser(description='Get sequence identity from two fasta files.')
    p.add_argument('fasta_fileA')
    p.add_argument('fasta_fileB')
    args = vars(p.parse_args(sys.argv[1:]))
    fasta_filenameA = args['fasta_fileA']
    fasta_filenameB = args['fasta_fileB']
    seqA = parse_fasta.read_fasta(open(fasta_filenameA, 'r')).values()[0][0]
    seqB = parse_fasta.read_fasta(open(fasta_filenameB, 'r')).values()[0][0]
    align = pairwise2.align.localms(seqA,seqB , 1, -1, -0.5, -0.1)
    minlen=len(seqA)
    if len(seqB)<minlen:
        minlen=len(seqB)
#    print "Identity: ",fasta_filenameA,fasta_filenameB,float(align[0][2])/float(align[0][4]-align[0][3])
    print "Identity: ",fasta_filenameA,fasta_filenameB,float(align[0][2])/float(minlen)


#!/usr/bin/env python3
from Bio import SeqIO
import argparse
from argparse import RawTextHelpFormatter

p = argparse.ArgumentParser(description = '- Sto to fasta format-',
                            formatter_class=RawTextHelpFormatter)
p.add_argument('-file','--input','-i', required= True, help='Input File A')
#p.add_argument('-out','--output','-o', required= False, help='output image')
#parser.add_argument('--nargs', nargs='+')
ns = p.parse_args()

file=ns.input
handle = open(file, 'r')
for record in SeqIO.parse(handle,"stockholm"): 
    print (">",record.id,"\n",record.seq)


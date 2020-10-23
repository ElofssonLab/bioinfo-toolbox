#!/usr/bin/env python3
from Bio import SeqIO
import argparse
from argparse import RawTextHelpFormatter

p = argparse.ArgumentParser(description = '- Convert Stockholm to trimmed a3m format -',
                            formatter_class=RawTextHelpFormatter)
p.add_argument('-i', required=True, help='input Stockholm file')
p.add_argument('-o', required=True, help='output trimmed a3m file')
ns = p.parse_args()

msa = []
for record in SeqIO.parse(open(ns.i),"stockholm"):
    msa.append(['>'+record.id.strip()+'\n', record.seq.strip()])

qmask = [pos for pos, ch in enumerate(msa[0][1]) if ch != '-']

trimmed = []
for code, seq in msa:
    tseq = ''
    for pos in qmask: tseq += seq[pos]
    tseq += '\n'

    trimmed.append([code, tseq])

outfile = open(ns.o, 'w')
for code, seq in trimmed: outfile.write(code+seq)
outfile.close()
    

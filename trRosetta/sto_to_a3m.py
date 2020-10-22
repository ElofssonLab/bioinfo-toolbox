#!/usr/bin/env python3
from Bio import SeqIO
import argparse
from argparse import RawTextHelpFormatter

p = argparse.ArgumentParser(description = '- Convert Stockholm to trimmed a3m format -',
                            formatter_class=RawTextHelpFormatter)
p.add_argument('-i', required=True, help='input Stockholm file')
p.add_argument('-top', required=False, default=None, help='maximum number of hits to include in a3m MSA')
p.add_argument('-cov', required=False, default=None, help='minimum coverage (0-100) to include hits')
p.add_argument('-o', required=True, help='output trimmed a3m file')
ns = p.parse_args()

msa = []
for record in SeqIO.parse(open(ns.i),"stockholm"):
    msa.append(['>'+record.id.strip()+'\n', record.seq.strip()])

qmask = [pos for pos, ch in enumerate(msa[0][1]) if ch != '-']

if ns.cov != None: cov = int(ns.cov)
else: cov = None
if ns.top != None: top = int(ns.top)
else: top = None

trimmed = []
for code, seq in msa:
    tseq = ''
    for pos in qmask: tseq += seq[pos]
    tseq += '\n'

    trimmed.append([code, tseq])

outfile = open(ns.o, 'w')
for code, seq in trimmed: outfile.write(code+seq)
outfile.close()
    

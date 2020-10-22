#!/usr/bin/env python3
from Bio import SeqIO
import argparse
from argparse import RawTextHelpFormatter

def parse_msa(msa_file, top, cov):
    query = ''
    hits = {}
    proteomes = {}
    for record in SeqIO.parse(open(msa_file), "fasta"):
        if query == '': query = [record.id.strip(), record.seq.strip()]
        else:
            frag = record.id.split('|')[1]+'_'+record.id.split('/')[-1].rstrip()
            prot = record.id.split('|')[-1].split('/')[0]
            seq = record.seq.strip()

            if cov != None:
                if (len(seq.strip('-'))/len(query[1]))*100 < cov: continue
            if len(hits.keys) == top: break

            if prot not in proteomes: proteomes[prot] = []
            proteomes[prot].append(frag)
            hits[frag] = seq

    return query, hits, proteomes

def merge_msa(q1, h1, p1, q2, h2, p2, sep='A'*20):
    merged_msa = []
    qmcode = q1[0]+'-'+q2[0]+'\n'
    qmseq = q1[1]+sep+q2[1]+'\n'
    merged_msa.append([qmcode, qmseq])

    proteomes1 = set(p1.keys())
    proteomes2 = set(p2.keys())
    common = proteomes1.intersection(proteomes2)
    for proteome in common:
        frag1 = p1[proteome][0]
        frag2 = p2[proteome][0]
        mcode = '>'+frag1+'-'+frag2+'\n'
        mseq = h1[frag1]+sep+h2[frag2]+'\n'
        merged_msa.append([mcode, mseq])

    return merged_msa


p = argparse.ArgumentParser(description = '- Merge hits belonging to the same proteome from a couple of a3m MSA -',
                            formatter_class=RawTextHelpFormatter)
p.add_argument('-i1', required=True, help='input a3m file A')
p.add_argument('-i2', required=True, help='input a3m file B')
p.add_argument('-top', required=False, default=None, help='number of top-hits to include from input MSAs')
p.add_argument('-cov', required=False, default=None, help='minimum coverage (0-100) to include hits')
ns = p.parse_args()


if ns.cov != None: cov = int(ns.cov)
else: cov = None
if ns.top != None: top = int(ns.top)
else: top = None

q1, h1, p1 = parse_msa(ns.i1, top, cov)
q2, h2, p2 = parse_msa(ns.i2, top, cov)

merged12 = merge_msa(q1, h1, p1, q2, h2, p2)
merged21 = merge_msa(q2, h2, p2, q1, h1, p1)

out = ns.i1.rstrip(ns.i1.split('.')[-1])
out += '-'+ns.i2.rstrip(ns.i2.split('.')[-1])
out += '.trimmed'
with open(out, 'w') as f: 
    for code, seq in merged12: f.write(code+seq)

out = ns.i2.rstrip(ns.i2.split('.')[-1])
out += '-'+ns.i1.rstrip(ns.i1.split('.')[-1])
out += '.trimmed'
with open(out, 'w') as f:
    for code, seq in merged21: f.write(code+seq)

#!/usr/bin/env python
import os, sys, myfunc
from math import ceil

file_pairalnfile="/data3/wk/MPTopo/pfamAna_refpro/cellular_filter_all/pairwise/withinClan/Pfam-A-full.perTM75_nseq20.nr100.filtered.withinclan.max30000.kalignP.pairaln"

(idList, annoList, seqList) = myfunc.ReadFasta(file_pairalnfile)

numseq = len(idList)

outpath = "splitted"

os.system("mkdir -p %s"%outpath)

nsplit = 10

numPair = numseq / 2
pairPerSplit = int(ceil(float(numPair) / nsplit))

bp = 0
for i in xrange(nsplit):
    outfile=outpath + os.sep + "split_%d" %i + ".fa"
    fpout = open(outfile, "w")
    for p in range(bp, bp + pairPerSplit):
        if p < numPair:
            anno1 = annoList[2*p]
            anno2 = annoList[2*p+1]
            seq1 = seqList[2*p]
            seq2 = seqList[2*p+1]
            fpout.write(">%s\n"%anno1)
            fpout.write("%s\n"%seq1)
            fpout.write(">%s\n"%anno2)
            fpout.write("%s\n"%seq2)
    bp += pairPerSplit
    fpout.close()
    print "%s output" % outfile


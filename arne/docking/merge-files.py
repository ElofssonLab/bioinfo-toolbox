#!/usr/bin/env python

from __future__ import print_function
import sys
import os


def mergefile(fileA,fileB,outfile):
    MERGE="~/git/bioinfo-toolbox/arne/docking/mergeMSAbyGenome.py"
    os.system(MERGE+"  "+fileA+"  "+fileB+ "   > " + outfile)  
if __name__ == '__main__':
    dirname=sys.argv[1]
    filename=sys.argv[2]
    fasta=sys.argv[3]
    extension=".trimmed"
    cutoffs=["HH1","HH0.0001","HH1.e-10"]
    with open(filename) as f:
        for l in f.readlines():
            l_arr = l.strip().split()
            A=l_arr[0]
            B=l_arr[1]
            for i in cutoffs:
                fileA=dirname+"/"+A+"."+fasta+"."+i+extension
                fileB=dirname+"/"+B+"."+fasta+"."+i+extension
                outfile=dirname+"/"+A+"-"+B+"."+i+extension
                if (not os.path.isfile(outfile)):
                    print(fileA,fileB,outfile)
                    mergefile(fileA,fileB,outfile)

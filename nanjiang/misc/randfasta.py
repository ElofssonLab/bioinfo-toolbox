#!/usr/bin/env python
# Filename: randfasta.py
# Description:
#   output a number of randomly selected sequences from a fasta file
# Author:
#   Nanjiang Shu  nanjiang.shu@scilifelab.se
import os
import sys
import random
import myfunc

BLOCK_SIZE = 100000

usage="""
Usage:  randfasta.py [-i] fastafile  [-n INT ] 

Description: 
    output a number of randomly selected sequences from a fasta file 

OPTIONS:
  -i    FILE  Set the input file
  -n     INT  Output N sequences instead of all
  -o    FILE  Output the result to file
  -seed  INT  Set random seed, (default: set by time)
  -bs, --block-size INT
              Size for blocks when reading file, (default: 100000)
  -h, --help  Print this help message and exit

Created 2011-04-08, updated 2011-10-31, Nanjiang Shu

Examples:
    randfasta.py test.fa -n 20 -o test.sel20.fa # output 20 random sequence from test.fa
"""

def PrintHelp():#{{{
    print usage
#}}}
def RandFasta(inFile,N,rand_seed, fpout):#{{{
    (idList,annotationList, seqList) = myfunc.ReadFasta(inFile, BLOCK_SIZE)
    if idList == None:
        print >> sys.stderr, "Failed to read fastafile %s. Exit."%inFile
        return -1
    random.seed(rand_seed)
    Nseq=len(idList)
    if N > Nseq:
        N=Nseq
    idxArray=range(Nseq)
    idxSample=random.sample(idxArray,N)
    for i in xrange(N):
        idx = idxSample[i]
        fpout.write(">%s\n"% annotationList[idx])
        fpout.write("%s\n"% seqList[idx])
    return 0
#}}}

def main():#{{{
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outFile=""
    inFile=""
    N=999999999
    rand_seed=None

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            isNonOptionArg=False
            i = i + 1
        elif sys.argv[i] == "--":
            isNonOptionArg=True
            i = i + 1
        elif sys.argv[i][0] == "-":
            if sys.argv[i] in [ "-h", "--help"]:
                PrintHelp()
                return 1
            elif sys.argv[i] in ["-i",  "--infile"]:
                inFile, i = myfunc.my_getopt_str(sys.argv, i)
            elif sys.argv[i] in [ "-n" ,  "--n"]:
                N,i = myfunc.my_getopt_int(sys.argv, i)
            elif sys.argv[i] in [ "-seed" , "--seed"]:
                rand_seed, i = myfunc.my_getopt_int(sys.argv, i)
            elif sys.argv[i] in [ "-o" , "--outfile"]:
                outFile,i  = myfunc.my_getopt_str(sys.argv, i)
            elif sys.argv[i] in [ "-bs" ,  "--block-size" ,  "-block-size"]:
                BLOCK_SIZE, i = myfunc.my_getopt_int(sys.argv, i)
                if BLOCK_SIZE < 0:
                    print >> sys.stderr,"Error! BLOCK_SIZE should >0"
                    return 1
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i])
                return 1
        else:
            inFile = sys.argv[i]
            i+=1

    if myfunc.checkfile(inFile, "Input file") != 0:
        return 1
    fpout = myfunc.myopen(outFile, sys.stdout, "w", False)
    RandFasta(inFile, N, rand_seed,  fpout)
    myfunc.myclose(fpout)
#}}}
if __name__ == '__main__' :
    main()

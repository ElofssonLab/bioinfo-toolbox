#!/usr/bin/env python
# Filename: catfasta.py
# Description:
#   Extract sequences from a fasta file with specified starting index and
#   ending index
# Author:
#   Nanjiang Shu  nanjiang.shu@scilifelab.se
import sys
import re
import os
import myfunc

BLOCK_SIZE = 100000

usage="""
Usage: catfasta.py  [-i] fastafile

Description:
    Extract sequences from a fasta file with specified starting index and
    ending index

OPTIONS:
  -i    FILE      Set the input file
  -b    INT       Set the index of the sequence starts to output, (default: 0)
  -e    INT       Set the index of the sequence stops to output, (default: END)
  -o    FILE      Output the result to file, (default: stdout)
  -h, --help      print this help message and exit

Created 2010-08-23, updated 2012-06-10, Nanjiang Shu

Examples:
    catfasta.py -b 0 -e 20 test.fa  #extract the first 20 sequences from test.fa
"""

def PrintHelp():
    print usage

def CatFasta(inFile,begin,end,fpout):#{{{
    cntSeq=0
    fpin = open(inFile, "r")
    line = fpin.readline()
    while line:
        line = line.rstrip('\n').strip()
        if line and line[0] == ">":
            if cntSeq > end + 1:
                break
            cntSeq+=1
            if cntSeq > begin and cntSeq <= end:
                print >> fpout,"%s"%line
            line = fpin.readline()
            while line and line[0] != ">":
                line = line.rstrip('\n').strip()
                if cntSeq > begin and cntSeq <= end:
                    print >> fpout,"%s"%line
                line = fpin.readline()
        else:
            line = fpin.readline()

    fpin.close()
    return 0
#}}}
def CatFasta2(inFile,beginSeqIndex,endSeqIndex,fpout):#{{{
    """
    Read fasta file by blocks
    """
    cntSeq=0
    fpin = open(inFile, "r")
    buff = fpin.read(BLOCK_SIZE)
    brokenseq=""; ##for the seq broken by BLOCK
    while buff:
        if cntSeq > endSeqIndex:
            break
        beg=0
        end=0
        while 1:
            if brokenseq:
                end=buff.find("\n>")
                if end >= 0:
                    seq=brokenseq+buff[0:end]
                    brokenseq=""
                    beg=end
                    if cntSeq > beginSeqIndex and cntSeq <= endSeqIndex:
                        fpout.write("%s\n"%seq)
                else:
                    brokenseq += buff
                    break

            beg=buff.find(">",beg)
            end=buff.find("\n>",beg+1)
            if beg >= 0:
                cntSeq+=1
                if end >=0:
                    seq=buff[beg:end]
                    beg=end
                    if cntSeq > beginSeqIndex and cntSeq <= endSeqIndex:
                        fpout.write("%s\n"%seq)
                else:
                    brokenseq=buff[beg:]
                    break
            else:
                brokenseq+=buff
                break
        buff = fpin.read(BLOCK_SIZE)
    if brokenseq:
        if cntSeq > beginSeqIndex and cntSeq <= endSeqIndex:
            fpout.write("%s\n"%brokenseq)

    fpin.close()
    return 0
#}}}

def main(g_params):
    # Check argv
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outFile = ""
    inFile = ""
    begin = 0
    end = 999999999

    method = 2

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
            if sys.argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif sys.argv[i] in ["-i", "--i", "--infile"]:
                inFile, i = myfunc.my_getopt_str(sys.argv, i)
            elif sys.argv[i] in ["-b", "--b", "--begin"]:
                begin, i = myfunc.my_getopt_int(sys.argv, i)
            elif sys.argv[i] in ["-e", "--e", "--end"]:
                end, i = myfunc.my_getopt_int(sys.argv, i)
            elif sys.argv[i] in ["-o" , "--o","--outfile"]:
                outFile, i = myfunc.my_getopt_str(sys.argv,i )
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i])
                return 1
        else:
            inFile=sys.argv[i]
            i+=1

    if myfunc.checkfile(inFile, "Input file") != 0:
        return 1

    fpout = myfunc.myopen(outFile, sys.stdout, "w", False)

    if method == 1:
        CatFasta(inFile,begin, end, fpout)
    else:
        CatFasta2(inFile,begin, end, fpout)

    myfunc.myclose(fpout)
    return 0

if __name__ == '__main__' :
    g_params = {}
    sys.exit(main(g_params));

#!/usr/bin/env python
import sys
import re
import os
import myfunc

progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

usage="""
usage:  %s msafile_in_fasta_format [-o OUTFILE]

Description: input the multiple sequence alignment, output sequence (with gaps
             removed) in fasta format

  -o  OUTFILE   Output the result to outfile
  -f  STR       Set format of the input file, default: auto
                if auto, the format will be detected automatically
  -h,--help     Print this help message and exit

Created 2013-01-10, updated 2013-01-10, Nanjiang Shu 
"""%(progname)

GAP='-'
BLOCK_SIZE = 100000

def PrintHelp():
    print usage
def AutoDetectFormat(infile):#{{{
    return "fasta"
#}}}
def MSA2Seq_fasta(infile, outfile):#{{{
    try: 
        fpout = myfunc.myopen(outfile, sys.stdout, "w", False)

        fpin = open (infile, "rb")
        unprocessedBuffer=""
        isEOFreached = False
        while 1:
            buff = fpin.read(BLOCK_SIZE)
            if len(buff) < BLOCK_SIZE:
                isEOFreached=True
            buff = unprocessedBuffer + buff
            recordList = []
            unprocessedBuffer = myfunc.ReadFastaFromBuffer(buff,recordList, isEOFreached)
            for rd in recordList:
                anno = rd[1]
                seq = rd[2].replace("-", "").replace(".", "").replace(" ", "")
                fpout.write(">%s\n"%anno)
                fpout.write("%s\n"%seq)
            if isEOFreached == True:
                break
        fpin.close()
        myfunc.myclose(fpout)
    except IOError:
        print >> sys.stderr, "Failed to read file", infile
        return 1
#}}}
def MSA2Seq(infile, input_format, outfile):#{{{
    if input_format == "fasta":
        return MSA2Seq_fasta(infile, outfile)
    else:
        return 1
#}}}
def main():#{{{
    numArgv = len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile = ""
    infile = ""
    input_format = "auto"

    i = 1
    isNonOptionArg = False
    while i < numArgv:
        if isNonOptionArg == True:
            infile = sys.argv[i]
            isNonOptionArg = False
            i = i + 1
        elif sys.argv[i] == "--":
            isNonOptionArg=True
            i = i + 1
        elif sys.argv[i][0] == "-":
            if sys.argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif sys.argv[i] in  ["-i", "--i", "--infile", "--infile"]:
                infile = sys.argv[i+1]
                i = i + 2
            elif sys.argv[i] in ["-f", "--f", "-format", "--format"]:
                input_format = sys.argv[i+1]
                i = i + 2
            elif sys.argv[i] in ["-o", "--o", "-outfile", "--outfile"]:
                outfile = sys.argv[i+1]
                i = i + 2
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i])
                return 1
        else:
            infile = sys.argv[i]
            i+=1

    if infile == "":
        print >> sys.stderr,"Error! msafile not set."
        return 1

    #sizeAASeqFile = os.path.getsize(fastaFile)
    if input_format == "auto":
        input_format = AutoDetectFormat(infile)


    MSA2Seq(infile, input_format, outfile)

#}}}
if __name__ == '__main__' :
    sys.exit(main())

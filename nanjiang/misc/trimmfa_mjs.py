#!/usr/bin/env python
# Given multiple sequence alignment, remove gaps in the target sequence
# and corresponding positions of the rest sequences in the alignment.
# Sequence descriptions will also be changed according to Marcin
import sys ,re,os
import myfunc

progname =  os.path.basename(sys.argv[0])

usage="""
Usage: %s MSA-FASTA-FILE [-o OUTFILE]

Description: Given multiple sequence alignment, remove gaps of the target
             sequence and the corresponding positions in the MSA. Replace 'X'
             to GAP as well.

  -o  OUTFILE   Output the result to OUTFILE
  -h, --help    Print this help message and exit

Created 2013-03-07, updated 2013-03-07, Nanjiang Shu
"""%(progname)

GAP='-'
BLOCK_SIZE=100000

def PrintHelp():
    print usage


def GetTrimmedPosition(alignedseq):#{{{
    seqlen = len(alignedseq)
    trimmedSeqPosList = []
    gapSegList = myfunc.GetSegPos(alignedseq, GAP)
    num = len(gapSegList)
    if num < 1:
        trimmedSeqPosList.append((0, seqlen))
    else:
        if gapSegList[0][0] > 0:
            trimmedSeqPosList.append((0, gapSegList[0][0]))
        for i in xrange(num-1):
            trimmedSeqPosList.append((gapSegList[i][1], gapSegList[i+1][0]))
        if gapSegList[num-1][1] < seqlen:
            trimmedSeqPosList.append((gapSegList[num-1][1], seqlen))
    return trimmedSeqPosList
#}}}

def GetPositionTermGapless(alignedseq):#{{{
    seqlen = len(alignedseq)
    i = 0
    while i < seqlen:
        if alignedseq[i] == GAP:
            i+=1
        else:
            break
    beg = i
    i = seqlen - 1
    while i >= 0:
        if alignedseq[i] == GAP:
            i -= 1
        else:
            break
    end = i+1
    return (beg, end)
#}}}

def main():#{{{
    numArgv = len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile=""
    infile=""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            infile = sys.argv[i]
            isNonOptionArg=False
            i = i + 1
        elif sys.argv[i] == "--":
            isNonOptionArg=True
            i = i + 1
        elif sys.argv[i][0] == "-":
            if sys.argv[i] ==  "-h" or  sys.argv[i] == "--help":
                PrintHelp()
                return 1
            elif sys.argv[i] in [ "-o", "--o"]:
                outfile=sys.argv[i+1]
                i = i + 2
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i])
                return 1
        else:
            infile=sys.argv[i]
            i+=1

    if infile == "":
        print >> sys.stderr,"Error! MSA file not set."
        return 1
    elif not os.path.exists(infile):
        print >> sys.stderr,"Error! MSA file %s does not exist."%(infile)
        return 1

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    isFirstSeqSet = False
    targetSeqID = ""
    targetSeq = ""
    trimmedSeqPosList = [] #segment list after removal of gaps in 
                           #the target sequence
    trimmedseqLength = 0
    fpin = open (infile, "rb")
    if not fpin:
        print >> sys.stderr, "Failed to open input file %s"%(infile)
        return -1
    unprocessedBuffer=""
    isEOFreached = False
    processedTopoIDSet = set([])
    cntseq = 0
    while 1:
        buff = fpin.read(BLOCK_SIZE)
        if len(buff) < BLOCK_SIZE:
            isEOFreached=True
        buff = unprocessedBuffer + buff
        recordList = []
        unprocessedBuffer = myfunc.ReadFastaFromBuffer(buff,recordList, isEOFreached)
        if len(recordList) > 0: 
            if not isFirstSeqSet:
                targetSeqID = recordList[0][0]
                targetSeq = recordList[0][2]
                if not targetSeqID == "target":
                    print >> sys.stderr, "Error, the first sequence is not target sequence in file %s"%(infile)
                    return 1
                trimmedSeqPosList = GetTrimmedPosition(targetSeq)
                lengthSegList = [(e-b) for (b,e) in trimmedSeqPosList]
                trimmedseqLength = sum(lengthSegList)
                trimmedseqLength = 100
                isFirstSeqSet = True
            for rd in recordList:
                desp = ""
                if rd[0] == "target":
                    desp = "target/1-%d"%(trimmedseqLength)
                else:
                    cntseq += 1
                    desp = "sequence%07d/1-%d"%(cntseq, trimmedseqLength)
                slist = [rd[2][p[0]:p[1]] for p in trimmedSeqPosList]
                trimmedseq = "".join(slist)
                trimmedseq = trimmedseq.replace("X", GAP)
                print >> fpout, ">%s"%(desp)
                print >> fpout, "%s"%(trimmedseq)
        if isEOFreached == True:
            break
    fpin.close()
    myfunc.myclose(fpout)


#}}}
if __name__ == '__main__' :
    sys.exit(main())

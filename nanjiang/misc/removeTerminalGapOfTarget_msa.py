#!/usr/bin/env python
# Given multiple sequence alignment, remove the terminal gaps of the target
# sequence and the corresponding positions in the MSA
import sys,re,os;
import myfunc ;
import cProfile;


usage="""
Usage: removeTerminalGapOfTarget_msa.py MSA-FASTA-FILE [-o OUTFILE]

Description: Given multiple sequence alignment, remove the terminal 
             gaps of the target sequence and the corresponding positions 
             in the MSA

  -o  OUTFILE   Output the result to OUTFILE
  -h, --help       Print this help message and exit

Created 2013-03-06, updated 2013-03-06, Nanjiang Shu 
"""

GAP='-';
BLOCK_SIZE=100000;

def PrintHelp():
    print usage;

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
    numArgv = len(sys.argv);
    if numArgv < 2:
        PrintHelp();
        return 1;

    outfile="";
    infile="";

    i = 1;
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            infile = sys.argv[i];
            isNonOptionArg=False;
            i = i + 1;
        elif sys.argv[i] == "--":
            isNonOptionArg=True;
            i = i + 1;
        elif sys.argv[i][0] == "-":
            if sys.argv[i] ==  "-h" or  sys.argv[i] == "--help":
                PrintHelp();
                return 1;
            elif sys.argv[i] in [ "-o", "--o"]:
                outfile=sys.argv[i+1];
                i = i + 2;
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i]);
                return 1;
        else:
            infile=sys.argv[i];
            i+=1;

    if infile == "":
        print >> sys.stderr,"Error! MSA file not set.";
        return 1;
    elif not os.path.exists(infile):
        print >> sys.stderr,"Error! MSA file %s does not exist."%(infile);
        return 1;

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False);
    isFirstSeqSet = False
    targetSeqID = ""
    targetSeqAnno = ""
    targetSeq = ""
    targetSeqBeg = 0; #starting position in alignment without terminal gaps
    targetSeqEnd = 0;

    fpin = open (infile, "rb");
    if not fpin:
        print >> sys.stderr, "Failed to open input file %s"%(infile);
        return -1;
    unprocessedBuffer="";
    isEOFreached = False;
    processedTopoIDSet = set([]);
    while 1:
        buff = fpin.read(BLOCK_SIZE);
        if len(buff) < BLOCK_SIZE:
            isEOFreached=True;
        buff = unprocessedBuffer + buff;
        recordList = [];
        unprocessedBuffer = myfunc.ReadFastaFromBuffer(buff,recordList, isEOFreached);
        if len(recordList) > 0: 
            if not isFirstSeqSet:
                targetSeqID = recordList[0][0]
                targetSeqAnno = recordList[0][1]
                targetSeq = recordList[0][2]
                if not targetSeqID == "target":
                    print >> sys.stderr, "Error, the first sequence is not target sequence in file %s"%(infile)
                    return 1
                (targetSeqBeg, targetSeqEnd) = GetPositionTermGapless(
                        targetSeq);
                isFirstSeqSet = True;
            for rd in recordList:
                print >> fpout, ">%s"%(rd[1]);
                print >> fpout, "%s"%(rd[2][targetSeqBeg:targetSeqEnd]);
        if isEOFreached == True:
            break;
    fpin.close();
    myfunc.myclose(fpout);


#}}}
if __name__ == '__main__' :
    # Check argv
    isPrintSeqID=True;
    sys.exit(main());
#    cProfile.run("main()");

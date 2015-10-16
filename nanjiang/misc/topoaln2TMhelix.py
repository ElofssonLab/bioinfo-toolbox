#!/usr/bin/env python
# Description:
import os
import sys
import re
import myfunc
usage = """
usage: topoaln2TMhelix.py  topofile  [-o OUTFILE] 

Description: given topology file in fasta format
             output a list of TM helices
             (gaps are not removed)
Options:
  -o FILE    Output the result to FILE, (default: stdout)
  -h, --help      Print this help message and exit

Created 2012-08-30, updated 2012-08-30, Nanjiang Shu 
"""

BLOCK_SIZE=100000
def PrintHelp():
    print usage

def GetSeqDict(fastafile):#{{{
    seqDict = {}
    (idList, seqList) = myfunc.ReadFasta_without_annotation(fastafile)
    for i in xrange(len(idList)):
        seqDict[idList[i]] = seqList[i]
    return seqDict
#}}}

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile = ""
    seqdbfile = ""
    infile = ""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            infile = argv[i]
            isNonOptionArg = False
            i += 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-outfile", "--outfile"]:
                outfile = argv[i+1]
                i += 2
            elif argv[i] in ["-seqdb", "--seqdb"]:
                seqdbfile = argv[i+1]
                i += 2
            elif argv[i] in ["-q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            infile = argv[i]
            i += 1
    if infile == "":
        print >> sys.stderr, "topofile not set"
        return 1
    elif not os.path.exists(infile):
        print >> sys.stderr, "topofile %s does not exist"%(infile)
        return 1
#     if seqdbfile == "":
#         print >> sys.stderr, "seqdbfile file not set"
#         return 1
#     elif not os.path.exists(seqdbfile):
#         print >> sys.stderr, "seqdbfile file %s does not exist"%(seqdbfile)
#         return 1
#     seqDict = GetSeqDict(seqdbfile)
#     if seqDict == {}:
#         print >> sys.stderr, "Failed to read seqdbfile %s"%(seqdbfile)
#         return 1

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    fpin = open (infile, "rb");
    if not fpin:
        print >> sys.stderr, "Failed to open input file %s"%(infile)
        return 1
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
            idListTopo = [r[0] for r in recordList];
            topoList = [r[2] for r in recordList];
            for i in xrange(len(idListTopo)):
                seqid = idListTopo[i]
                topo = topoList[i]

                posTM = myfunc.GetTMPosition(topo)
                if len(posTM) > 0:
                    cnt = 0
                    for (b,e) in posTM:
                        seg = topo[b:e]
                        fpout.write("%s\t%4d\t%s\n"%(seqid, cnt+1, seg))
                        cnt += 1

        if isEOFreached == True:
            break;
    fpin.close();

    myfunc.myclose(fpout)

#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

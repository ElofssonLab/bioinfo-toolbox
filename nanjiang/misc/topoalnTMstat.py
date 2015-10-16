#!/usr/bin/env python
# Description:
import os
import sys
import re
import myfunc
usage = """
usage: topoalnTMstat.py  topofile  [-o OUTFILE] 

Description: given pairwise topology alignment file in fasta format
             output statistics
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

def GetTopoAlignStat(topo1, topo2):
    stat = []
    posTM1 = myfunc.GetTMPosition(topo1)
    if len(posTM1) > 0:
        for (b, e) in posTM1:
            segList1 = []
            segList2 = []
            cntTM = 0
            cntGap = 0
            cntSeq = 0
            for j in xrange(b,e):
                if topo1[j] == 'M':
                    segList2.append(topo2[j])
                    if topo2[j] == 'M':
                        cntTM += 1
                    elif topo2[j] == '-':
                        cntGap += 1
                    else:
                        cntSeq += 1
            rd = {}
            sizeSeg = len(segList2)
            freqTM = myfunc.FloatDivision(cntTM, sizeSeg)
            freqGap = myfunc.FloatDivision(cntGap, sizeSeg)
            freqSeq = myfunc.FloatDivision(cntSeq, sizeSeg)

            rd['seg2'] = ''.join(segList2)
            rd['freqTM'] = freqTM
            rd['freqGap'] = freqGap
            rd['freqSeq'] = freqSeq
            stat.append(rd)
    return stat

def WriteStat(id1, id2, stat, fpout):
    for rd in stat:
        fpout.write("%s %s %s %6.2f %6.2f %6.2f\n" % (
            id1, id2, rd['seg2'], rd['freqTM'], rd['freqGap'], rd['freqSeq']))

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
    unprocessedRecordList = []
    isEOFreached = False;
    processedTopoIDSet = set([]);
    while 1:
        buff = fpin.read(BLOCK_SIZE);
        if len(buff) < BLOCK_SIZE:
            isEOFreached=True;
        buff = unprocessedBuffer + buff;
        recordList = unprocessedRecordList
        unprocessedBuffer = myfunc.ReadFastaFromBuffer(buff,recordList, isEOFreached);
        numRecord = len(recordList)
        if numRecord > 0: 
            numPair = numRecord/2
            for i in xrange(numPair):
                rd1 = recordList[2*i]
                rd2 = recordList[2*i+1]
                id1 = rd1[0]
                topo1 = rd1[2]
                id2 = rd2[0]
                topo2 = rd2[2]
                stat1 = GetTopoAlignStat(topo1, topo2)
                stat2 = GetTopoAlignStat(topo2, topo1)
                WriteStat(id1, id2, stat1, fpout)
                WriteStat(id2, id1, stat1, fpout)

        if numRecord%2 == 1:
            unprocessedRecordList = [recordList[-1]]
        else:
            unprocessedRecordList = []

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

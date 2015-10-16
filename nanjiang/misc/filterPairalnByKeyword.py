#!/usr/bin/env python
# Description:
# Filter the allinfo file by topcons_single
# filter uniref id by uniprot id and repalce the annotation line
import os
import sys
import myfunc
usage = """
usage: filterPairalnByKeyword.py pairalnfile -idwithanno FILE -o OUTFILE
Description:
# Filter pairaln file by keyword matching in idwithanno

Options:
  -q              Quiet mode
  -case           Case sensitive, (default: no)
  -h, --help      Print this help message and exit

Created 2012-06-07, updated 2012-06-07, Nanjiang Shu  
"""
BLOCK_SIZE=100000;

def PrintHelp():
    print usage
def ReadSeqIDWithAnnoFile(infile):#{{{
    try:
        seqid2AnnoDict = {}
        fpin = open(infile, "r")
        line = fpin.readline()
        while line:
            line = line.strip()
            if line and line[0] != "#":
                strs = line.split("\t")
                if len(strs) == 2:
                    seqid2AnnoDict[strs[0]] = strs[1]
            line = fpin.readline()
        fpin.close()
        return seqid2AnnoDict
    except IOError:
        print >> sys.stderr, "Failed to read file", infile
        return {}
#}}}
def FilterRecord(recordList, seqid2AnnoDict, fpout):#{{{
    numRecord = len(recordList)
    numPair = numRecord/2
    keywordList = g_params['keywordList']
    isCaseSensitive = g_params['isCaseSensitive']
    for i in xrange(numPair):
        rd1 = recordList[2*i]
        rd2 = recordList[2*i+1]
        isFilter = False
        seqid1 = rd1[0]
        seqid2 = rd2[0]
        anno1 = rd1[1]
        anno2 = rd2[1]
        for seqid in [seqid1, seqid2]:
            if seqid in seqid2AnnoDict:
                seqdef = seqid2AnnoDict[seqid]
                if not isCaseSensitive:
                    seqdef = seqdef.lower()
                for keyword in keywordList:
                    if not isCaseSensitive:
                        keyword = keyword.lower()
                        if seqdef.find(keyword) != -1:
                            isFilter = True
                            break
            if isFilter:
                break
        if not isFilter:
            fpout.write(">%s" % (anno1))
            if seqid1 in seqid2AnnoDict:
                fpout.write(" %s"%(seqid2AnnoDict[seqid1]))
            fpout.write("\n")
            fpout.write("%s\n"%(rd1[2]))
            fpout.write(">%s" % (anno2))
            if seqid2 in seqid2AnnoDict:
                fpout.write(" %s"%(seqid2AnnoDict[seqid2]))
            fpout.write("\n")
            fpout.write("%s\n"%(rd2[2]))
#}}}
def FilterPairalnByKeyword(infile, seqid2AnnoDict, fpout):#{{{
    try:
        fpin = open(infile, "r")
        unprocessedBuffer="";
        isEOFreached = False;
        processedTopoIDSet = set([]);
        remainedRd = None
        while 1:
            buff = fpin.read(BLOCK_SIZE);
            if len(buff) < BLOCK_SIZE:
                isEOFreached=True;
            buff = unprocessedBuffer + buff;
            recordList = [];
            if remainedRd != None:
                recordList.append(remainedRd)
            unprocessedBuffer = myfunc.ReadFastaFromBuffer(
                    buff, recordList, isEOFreached);
            numRecord = len(recordList)
            if numRecord > 0: 
                FilterRecord(recordList, seqid2AnnoDict, fpout)
                if numRecord%2 == 1:
                    remainedRd = recordList[numRecord-1]
                else:
                    remainedRd = None

            if isEOFreached == True:
                break;
        fpin.close();
    except IOError:
        print >> sys.stderr, "Failed to open file %s for read"%(infile);
        return 1
#}}}
def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile = ""
    infile = ""
    idwithannofile = ""

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
            elif argv[i] in ["-o", "--o"]:
                outfile = argv[i+1]
                i += 2
            elif argv[i] in ["-case", "--case"]:
                g_params['isCaseSensitive'] = True
                i += 1
            elif argv[i] in ["-idwithanno", "--idwithanno"]:
                idwithannofile = argv[i+1]
                i += 2
            elif argv[i] in ["-key", "--key"]:
                g_params['keywordList'].append(argv[i+1])
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
        print >> sys.stderr, "infile not set.exit"
        return 1
    if idwithannofile == "":
        print >> sys.stderr, "idwithannofile not set.exit"
        return 1

    seqid2AnnoDict = ReadSeqIDWithAnnoFile(idwithannofile)
    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    FilterPairalnByKeyword(infile, seqid2AnnoDict, fpout)
    myfunc.myclose(fpout)
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['keywordList'] = []
    g_params['isCaseSensitive'] = False
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

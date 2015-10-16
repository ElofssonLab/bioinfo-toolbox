#!/usr/bin/env python
# Description:
import os
import sys
import myfunc
import libtopologycmp as lcmp
usage = """
usage:  stat_dup_paircmp.py  -dup DUPFILE -paircmp paircmpFile
Description:

Options:
  -outpath DIR    Set ouput path
  -q              Quiet mode
  -h, --help      Print this help message and exit

Created 2012-11-30, updated 2012-11-30, Nanjiang Shu  
"""
BLOCK_SIZE = 100000

def PrintHelp():
    print usage
def ReadDupPairList(infile):
    try:
        fpin = open(infile, "r")
        lines = fpin.readlines()
        fpin.close()
        pairlist = []
        for line in lines:
            if line:
                try:
                    pair = line.split()[0].split("-")
                    pairlist.append((pair[0],pair[1]))
                except (IndexError):
                    pass
        return pairlist
    except IOError:
        return 1
def StatDupPaircmp(recordList, dupPairSet, fpout):
    for record in recordList:
        if record == {}:
            continue
        cmpclass = record['cmpclass']
        seqid1 = record['id1']
        seqid2 = record['id2']
        seqLength1 = record['seqLength1']
        seqLength2 = record['seqLength2']
        numTM1 = record['numTM1']
        numTM2 = record['numTM2']
        seqidt = record['seqidt']
        if (seqid1, seqid2) in dupPairSet or (seqid2, seqid1) in dupPairSet:
            fpout.write("%s %s %20s %2d %2d %4d %4d %6.1f\n"%(
                seqid1, seqid2, cmpclass, numTM1, numTM2, seqLength1,
                seqLength2, seqidt))

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = "./"
    dupfile = ""
    paircmpfile = ""
    outfile = ""
    
    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            isNonOptionArg = False
            i += 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-outpath", "--outpath"]:
                outpath = argv[i+1]
                i += 2
            elif argv[i] in ["-o", "--o"]:
                outfile = argv[i+1]
                i += 2
            elif argv[i] in ["-l", "--l"] :
                idListFile = argv[i+1]
                i += 2
            elif argv[i] in ["-dup", "--dup"] :
                dupfile = argv[i+1]
                i += 2
            elif argv[i] in ["-paircmp", "--paircmp"] :
                paircmpfile = argv[i+1]
                i += 2
            elif argv[i] in ["-q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            print >> sys.stderr, "Error! Wrong argument:", argv[i]
            return 1
    
    if paircmpfile == "":
        return 1
    if dupfile == "":
        return 1

    dupPairList = ReadDupPairList(dupfile)
    dupPairSet = set(dupPairList)
    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)

    fpin = open(paircmpfile, "r")
    unprocessedBuffer=""
    cntTotalReadInRecord = 0
    cntTotalOutputRecord = 0
    isEOFreached = False
    while 1:
        buff = fpin.read(BLOCK_SIZE)
        if buff == "":
            isEOFreached = True
        buff = unprocessedBuffer + buff
        pairCmpRecordList=[]
        unprocessedBuffer = lcmp.ReadPairCmpResultFromBuffer(buff,pairCmpRecordList)
            
        if len(pairCmpRecordList) > 0: 
            StatDupPaircmp(pairCmpRecordList, dupPairSet, fpout)
            cntTotalReadInRecord += len(pairCmpRecordList)
        if isEOFreached == True:
            break
    fpin.close()
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

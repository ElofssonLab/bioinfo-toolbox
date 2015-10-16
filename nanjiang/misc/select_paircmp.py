#!/usr/bin/env python
# /data3/wk/MPTopo/src/select_paircmp.py
import os
import sys
import libtopologycmp as lcmp
import myfunc
import copy
import subprocess

DEBUG_UNMAPPED_TM_POSITION = 0

BLOCK_SIZE = 100000
progname =  os.path.basename(sys.argv[0])
usage="""
Usage: %s paircmp-file [-o OUTFILE]

Description:
    Select paircmp file based on filtering schemes
OPTIONS:
  -o OUTFILE    Output the result to file
  -q            Quiet mode
  -h, --help    Print this help message and exit

  -tableinfo FILE    Set pairwise alignment table info, get more pairwise
                     statistics
  -seqidttype INT    Set sequence identity type, (default: 1)
                     0: seqIDT = numIDTRes /alnLength
                     1: seqIDT = numIDTRes / min(len1, len2)
                     2: seqIDT = numIDTRes / (alnLength - NumGAP)
                     Note: if seqidttype = 1 or 2, tableinfo file must be set
  -alignrange STR    Select alignment with different alignment ranges
                     all, full, part, (default: all)
  -cmpclass   STR    cmpclass 
  -min-seqidt FLOAT  (default: 0)
  -max-seqidt FLOAT  (default: 100)
  -signalp    FILE   signalp file
  -rmsp              Remove pairs with signalp
  -restrictidlist FILE  Set restriction seq idlist

Created 2013-06-26, updated 2013-08-13, Nanjiang Shu 

Example:
    %s t1.paircmp -tableinfo t1.tableinfo -seqidttype 1 -alignrange full -cmpclass TM2SEQ -cmpclass TM2GAP_AND_TM2SEQ
"""%(progname, progname)

def PrintHelp():#{{{
    print usage
#}}}

def AddTableInfo(recordList, pairalnStat):#{{{
    if pairalnStat != {}:
        for record in recordList:
            pairid = record['id1'] + '-' + record['id2']
            if pairid in pairalnStat:
                record['seqidt1'] = pairalnStat[pairid]['seqidt1']
                record['seqidt2'] = pairalnStat[pairid]['seqidt2']
#                 print pairid, record['seqidt1']
#}}}
def AddSignalPInfo(recordList, signalpDict): #{{{
    if signalpDict != {}:
        for record in recordList:
            id1 =  record['id1']
            id2 =  record['id2']
            try:
                record['sp1'] = signalpDict[id1]
            except KeyError:
                record['sp1'] = -1
            try:
                record['sp2'] = signalpDict[id2]
            except KeyError:
                record['sp2'] = -1

#}}}
def AddDupInfo(recordList, dupPairSet):#{{{
    for record in recordList:
        id1 =  record['id1']
        id2 =  record['id2']
        if (id1,id2) in dupPairSet or (id2,id1) in dupPairSet:
            record['isDup'] = True
        else:
            record['isDup'] = False
#}}}

def FilterPairCmpResult(recordList, cmpclassList, rltyDict, #{{{
        restrictIDSet):
    """
    Filter paircmp result by g_params
    return newList
    """
    newList = []
    pairListSet = set([])
    numInputRecord = len(recordList)
    seqidttype = g_params['seqidttype']
    isRemoveSignalP = g_params['isRemoveSignalP']
    isRemoveDup = g_params['isRemoveDup']
    for record in recordList:
#         print record['id1'], record['seqidt1'] 
#         print record['id2'], record['seqidt1']
        if record == {}:
            continue


        id1 = record['id1']
        id2 = record['id2']

        if ((g_params['isRestrictIDListSet'] == True) and 
                ((not id1 in restrictIDSet) 
                or (not id2 in restrictIDSet))):
            continue



        if isRemoveSignalP:
            if 'sp1' in record and record['sp1'] != -1:
                continue
            if 'sp2' in record and record['sp2'] != -1:
                continue

        if isRemoveDup:
            if 'isDup' in record and record['isDup']:
                continue
        if record['isLocalAlignment'] and g_params['alignrange'] != 'all':
            if record['alignrange'] != g_params['alignrange']:
                continue

        if rltyDict != {}:
            if id1 in rltyDict:
                rlty = rltyDict[id1]
#                 print "rlty[%s]=%.1f"%(id1, rlty)
                if rlty < g_params['minRLTY'] or rlty > g_params['maxRLTY']:
                    continue
            if id2 in rltyDict:
#                 print "rlty[%s]=%.1f"%(id1, rlty)
                if rlty < g_params['minRLTY'] or rlty > g_params['maxRLTY']:
                    continue
        seqidt = lcmp.GetSeqIDT(record, seqidttype)
        if (seqidt < g_params['minSeqIDT'] or seqidt >=
                g_params['maxSeqIDT']):
            continue
        if len(cmpclassList) == 0 or record['cmpclass'] in  cmpclassList:
            newList.append(record)
    if g_params['isDEBUG']:
        if numOutputRecord < numInputRecord:
            print "%d pairs dropped" % (numInputRecord-numOutputRecord)
    return newList
#}}}
def GetSeqIDTGroupIndex(seqidt, seqIDTGroupList):#{{{
    numGroup = len(seqIDTGroupList)/2
    for i in xrange(numGroup):
        if seqidt >= seqIDTGroupList[i*2] and seqidt < seqIDTGroupList[i*2+1]:
            return i
    return numGroup
#}}}


def main(g_params):#{{{
    argv = sys.argv
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    infile = ""
    outpath = "./"
    isQuiet = False
    tableinfoFile = ""
    cmpclassList = []
    restrictIDListFile = ""

    signalpFile = ""
    dupFile = ""
    outfile = ""
    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            infile = sys.argv[i]
            isNonOptionArg = False
            i += 1
        elif sys.argv[i] == "--":
            isNonOptionArg=True
            i += 1
        elif sys.argv[i][0] == "-":
            if sys.argv[i] in [ "-h", "--help"]:
                PrintHelp()
                sys.exit()
            elif argv[i] in [ "-o", "--o"]:
                (outfile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in [ "-cmpclass", "--cmpclass"]:
                (tmpstr, i) = myfunc.my_getopt_str(argv, i)
                cmpclassList.append(tmpstr)
            elif argv[i] in [ "-signalp", "--signalp"]:
                (signalpFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in [ "-restrictidlist", "--restrictidlist"]:
                (restrictIDListFile, i) = myfunc.my_getopt_str(argv, i)
                g_params['isRestrictIDListSet'] = True
            elif argv[i] in [ "-dup", "--dup", "-dupfile", "--dupfile"]:
                (dupFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in [ "-rmsp", "--rmsp"]:
                g_params['isRemoveSignalP']  = True; i+=1
            elif argv[i] in [ "-rmdup", "--rmdup"]:
                g_params['isRemoveDup']  = True; i+=1
            elif argv[i] in ["-seq2fammap", "--seq2fammap"]:
                (seq2famMapfile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-seqidttype", "--seqidttype"]:
                g_params['seqidttype'], i = myfunc.my_getopt_int(argv,i)
            elif argv[i] in ["-tableinfo", "--tableinfo"]:
                tableinfoFile, i = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-min-seqidt", "--min-seqidt"]:
                g_params['minSeqIDT'], i  = myfunc.my_getopt_float(argv, i)
            elif argv[i] in ["-max-seqidt", "--max-seqidt"]:
                g_params['maxSeqIDT'], i  = myfunc.my_getopt_float(argv, i)
            elif argv[i] in ["-evodist", "--evodist"]:
                g_params['isEvodist'] = True
                i += 1
            elif argv[i] in ["-alignrange", "--alignrange"]:
                g_params['alignrange'],i  = myfunc.my_getopt_str(argv,i)
                if not g_params['alignrange'] in ['all', 'full', 'part']:
                    print >> sys.stderr, "alignrange must be one of [all, full, part]"
                    return 1
                else:
                    if g_params['alignrange'] == 'full':
                        g_params['alignrange'] = 'FULL_ALIGNED'
                    elif g_params['alignrange'] == 'part':
                        g_params['alignrange'] = 'PART_ALIGNED'
            elif argv[i] in ["-debug", "--debug"]:
                if argv[i+1][0].lower() == 'y':
                    g_params['isDEBUG'] = True
                else:
                    g_params['isDEBUG'] = False
                i += 2
            elif argv[i] in ["-debug-unmapped-position", "--debug-unmapped-position"]:
                DEBUG_UNMAPPED_TM_POSITION = 1
                i += 2
            elif sys.argv[i] == "-q":
                isQuiet=True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", sys.argv[i]
                return -1
        else:
            infile = sys.argv[i]
            i += 1
    if infile == "":
        print >> sys.stderr, "infile not set. Exit."
        return -1
    elif not os.path.exists(infile):
        print >> sys.stderr, "infile %s does not exists. Exit."%infile

    try:
        fpin = open(infile, "rb")
    except IOError:
        print >> sys.stderr, "Failed to open input file %s"%(infile)
        return -1


    pairalnStat = {}
    if g_params['seqidttype'] != 0:
        if tableinfoFile == "" or not os.path.exists(tableinfoFile):
            print >> sys.stderr, "tableinfoFile must be set when seqidttype is set to 1 or 2"
            print >> sys.stderr, "but seqidttype = %d is set. Exit."%g_params['seqidttype']
            return -1
        pairalnStat = lcmp.ReadPairAlnTableInfo(tableinfoFile)


    rootname = os.path.basename(os.path.splitext(infile)[0])

    binpath = os.path.dirname(sys.argv[0])

    signalpDict = {}
    if signalpFile != "":
        signalpDict = lcmp.ReadSignalPDict(signalpFile)
    if signalpDict != {}:
        g_params['isSignalPSet'] = True

    dupPairList = []
    if dupFile != "":
        dupPairList = lcmp.ReadDupPairList(dupFile)
    if len(dupPairList) > 0:
        g_params['isDupSet'] = True
    dupPairSet = set(dupPairList)

    restrictIDSet = set([])
    if restrictIDListFile != "":
        restrictIDSet = set(myfunc.ReadIDList(restrictIDListFile))

    rltyDict = {}
    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
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

        AddTableInfo(pairCmpRecordList, pairalnStat)
        AddSignalPInfo(pairCmpRecordList, signalpDict)
        AddDupInfo(pairCmpRecordList, dupPairSet)

        cntTotalReadInRecord += len(pairCmpRecordList)
        pairCmpRecordList = FilterPairCmpResult(pairCmpRecordList, cmpclassList, rltyDict, restrictIDSet)

        if len(pairCmpRecordList) > 0: 
            lcmp.WritePairCmpRecord(pairCmpRecordList, cntTotalOutputRecord, fpout)
            cntTotalOutputRecord += len(pairCmpRecordList)
        if isEOFreached == True:
            break
    fpin.close()

    print "cntTotalReadInRecord =", cntTotalReadInRecord
    print "cntTotalOutputRecord =", cntTotalOutputRecord
    myfunc.myclose(fpout)
    return 0
#}}}
def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isDEBUG'] = False
    g_params['selecttype'] = 'all'
    g_params['outpath'] = ""

    g_params['minGapFraction'] = 0.0
    g_params['maxGapFraction'] = 1.0

    g_params['minRLTY'] = 0.0 # minimal reliability score
    g_params['maxRLTY'] = 100.0 # maximal reliability score
    g_params['seqidttype'] = 0

    g_params['minDGvalue'] = -999999.0
    g_params['maxDGvalue'] = 100
    g_params['minSeqIDT'] = 0.0
    g_params['maxSeqIDT'] = 100.0
    g_params['isShowProgress'] = True

    g_params['isPrintDIFFPair'] = False
    g_params['DIFFPairList'] = []

    g_params['isPrintCountPairInFam'] = False
    g_params['countPairInFam'] = []
    g_params['isPrintFileRltyCmpclass'] = False
    g_params['isPrintFileRltyHelixCmpclass'] = False
    g_params['isPrintNumTMHeatMap'] = True

    g_params['isEvodist'] = False
    g_params['pairwise_comparison_method']  = 0
    g_params['thHigh_seqidt'] = 30.0
    g_params['thHigh_evodist'] = 1.0

    g_params['isFilterPredictedSeq'] = False
    g_params['isRLTYSupplied'] = False
    g_params['isSignalPSet'] = False
    g_params['isRemoveSignalP'] = False
    g_params['isRemoveDup'] = False
    g_params['isDupSet'] = False
    g_params['numTMHeatMapMode'] = "full"
    g_params['isRestrictIDListSet'] = False
    g_params['alignrange'] = 'all'
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

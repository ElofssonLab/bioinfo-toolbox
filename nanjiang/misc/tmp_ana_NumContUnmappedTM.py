#!/usr/bin/env python

import os
import sys
import re
import myfunc
import libtopologycmp as lcmp

DEBUG_UNMAPPED_TM_POSITION = 0
DEBUG_13 = 1

BLOCK_SIZE = 100000
NUMTMSET = set(range(100))
NUMTMSET = set([2, 4, 6, 8, 10, 12])

def FilterSegPos(posList, string, neighbour_char):#{{{
### return only list of "0110"
    newList = []
    N = len(string)
    for (b,e) in posList:
        if b == 0 or e == N:
            continue
        if b>0 and string[b-1] != neighbour_char:
            continue
        if e < N-1 and string[e] != neighbour_char:
            continue
        newList.append((b,e))
    return newList#}}}
def GetSegPos(string, keyC):#{{{
    """
    Get segment of a continue keyC state
    e.g. given a string "0001100022000111100"
    and keyC = '1'
    return [(3,5), (13,17)]
    """
    posList = []
    ex = "(%s+)"%(keyC)
    m = re.finditer(ex,string)
    for i in m:
        posList.append((i.start(0), i.end(0)))
    return posList
#}}}
def CountContinuousUnmappedTM_0(recordList, freqDict, isFilterNeighbour=True):#{{{
### one of the topology should be all mapped and the other one should only have TM2GAP
    if isFilterNeighbour:
        str_filter = "True"
    else:
        str_filter = "False"
    for record in recordList:
        if record == {}:
            continue
        if not "mapArray" in  record:
            continue
        id1 = record['id1']
        id2 = record['id2']
        key = "%s-%s"%(id1, id2)


        mapArrayList = [record['mapArray'][i] for i in range(2)]

        mapArrayList[0] = filter(lambda a:  a!= -1, mapArrayList[0])

        mapArrayList[1] = filter(lambda a:  a!= -1, mapArrayList[1])

        seqidList = [record['id%d'%(i+1)] for i in range(2)]
        numTMList  = [len(mapArray) for mapArray in mapArrayList]
        if min(numTMList) <= 1: #ignore single spanning TM protein in location analysis
            continue
        max_st_list = [max(mapArray) for mapArray in mapArrayList]
        if set(max_st_list) == set([0,1]):
            if set(mapArrayList[0]) == set([0, 1]):
                pMapArray = mapArrayList[0]
            else:
                pMapArray = mapArrayList[1]

            str_maparray_list = ["%d"%x for x in pMapArray]
            str_maparray = "".join(str_maparray_list)
            neighbour_char = "0"
            st = 1

            posContList0 = GetSegPos(str_maparray, "%d"%st)
            if isFilterNeighbour:
                posContList = FilterSegPos(posContList0, str_maparray,
                        neighbour_char)
            else:
                posContList = posContList0
            if len(posContList) > 0:
                if DEBUG_UNMAPPED_TM_POSITION:
                    print >> sys.stderr, "Method 0 isFilterNeighbour = %s"%(str_filter)
                    print >> sys.stderr, "Unmapped TM position for %s - %s" %(id1, id2)
                    print >> sys.stderr, "mapArray1:", record['mapArray'][0]
                    print >> sys.stderr, "mapArray2:", record['mapArray'][1]
                    print >> sys.stderr, "posContList (nonFilt ):", posContList0
                    print >> sys.stderr, "posContList (Filtered):", posContList
                    print >> sys.stderr
                for (b,e) in posContList:
                    numTMSeg = e-b
                    if DEBUG_13 and numTMSeg == 13:
                        print >> sys.stderr
                        print >> sys.stderr, "Unmapped TM position for %s - %s" %(id1, id2)
                        print >> sys.stderr, "mapArray1:", record['mapArray'][0]
                        print >> sys.stderr, "mapArray2:", record['mapArray'][1]

                    if not numTMSeg in freqDict:
                        freqDict[numTMSeg] = 0
                    freqDict[numTMSeg] += 1

#}}}
def CountContinuousUnmappedTM_1(recordList, freqDict, isFilterNeighbour=True):#{{{
### one of the topology should have specified numTM
    if isFilterNeighbour:
        str_filter = "True"
    else:
        str_filter = "False"
    for record in recordList:
        if record == {}:
            continue
        if not "mapArray" in  record:
            continue
        id1 = record['id1']
        id2 = record['id2']
        key = "%s-%s"%(id1, id2)

        mapArrayList = [record['mapArray'][i] for i in range(2)]
        seqidList = [record['id%d'%(i+1)] for i in range(2)]
        numTMList  = [len(mapArray) for mapArray in mapArrayList]
        if min(numTMList) <= 1: #ignore single spanning TM protein in location analysis
            continue
        if ((not numTMList[0] in NUMTMSET) and (not numTMList[1] in
            NUMTMSET)):
            continue

        for i in range(2):
            pMapArray = mapArrayList[i]

            str_maparray_list = ["%d"%x for x in pMapArray]
            str_maparray = "".join(str_maparray_list)
            neighbour_char = "0"
            st = 1

            posContList0 = GetSegPos(str_maparray, "%d"%st)
            if isFilterNeighbour:
                posContList = FilterSegPos(posContList0, str_maparray,
                        neighbour_char)
            else:
                posContList = posContList0
            if len(posContList) > 0:
                if DEBUG_UNMAPPED_TM_POSITION:
                    print >> sys.stderr, "Method 1 isFilterNeighbour=%s"%(str_filter)
                    print >> sys.stderr, "Unmapped TM position for %s - %s" %(id1, id2)
                    print >> sys.stderr, "mapArray1:", record['mapArray'][0]
                    print >> sys.stderr, "mapArray2:", record['mapArray'][1]
                    print >> sys.stderr, "posContList (nonFilt ):", posContList0
                    print >> sys.stderr, "posContList (Filtered):", posContList
                    print >> sys.stderr
                for (b,e) in posContList:
                    numTMSeg = e-b
                    if DEBUG_13 and numTMSeg == 13:
                        print >> sys.stderr
                        print >> sys.stderr, "Unmapped TM position for %s - %s" %(id1, id2)
                        print >> sys.stderr, "mapArray1:", record['mapArray'][0]
                        print >> sys.stderr, "mapArray2:", record['mapArray'][1]
                    if not numTMSeg in freqDict:
                        freqDict[numTMSeg] = 0
                    freqDict[numTMSeg] += 1
#}}}
def CountContinuousUnmappedTM_2(recordList, freqDict, isFilterNeighbour=True):#{{{
### TM2SEQ are filtered
    if isFilterNeighbour:
        str_filter = "True"
    else:
        str_filter = "False"
    for record in recordList:
        if record == {}:
            continue
        if not "mapArray" in  record:
            continue
        id1 = record['id1']
        id2 = record['id2']
        key = "%s-%s"%(id1, id2)

        mapArrayList = [record['mapArray'][i] for i in range(2)]
        seqidList = [record['id%d'%(i+1)] for i in range(2)]
        numTMList  = [len(mapArray) for mapArray in mapArrayList]
        if min(numTMList) <= 1: #ignore single spanning TM protein in location analysis
            continue
        max_st_list = [max(mapArray) for mapArray in mapArrayList]

        if max(mapArrayList[0]+mapArrayList[1]) == 2:
            continue

        for i in range(2):
            pMapArray = mapArrayList[i]

            str_maparray_list = ["%d"%x for x in pMapArray]
            str_maparray = "".join(str_maparray_list)
            neighbour_char = "0"
            st = 1

            posContList0 = GetSegPos(str_maparray, "%d"%st)
            if isFilterNeighbour:
                posContList = FilterSegPos(posContList0, str_maparray,
                        neighbour_char)
            else:
                posContList = posContList0
            if len(posContList) > 0:
                if DEBUG_UNMAPPED_TM_POSITION:
                    print >> sys.stderr, "Method 2 isFilterNeighbour=%s"%(str_filter)
                    print >> sys.stderr, "Unmapped TM position for %s - %s" %(id1, id2)
                    print >> sys.stderr, "mapArray1:", record['mapArray'][0]
                    print >> sys.stderr, "mapArray2:", record['mapArray'][1]
                    print >> sys.stderr, "posContList (nonFilt ):", posContList0
                    print >> sys.stderr, "posContList (Filtered):", posContList
                    print >> sys.stderr
                for (b,e) in posContList:
                    numTMSeg = e-b
                    if DEBUG_13 and numTMSeg == 13:
                        print >> sys.stderr
                        print >> sys.stderr, "Unmapped TM position for %s - %s" %(id1, id2)
                        print >> sys.stderr, "mapArray1:", record['mapArray'][0]
                        print >> sys.stderr, "mapArray2:", record['mapArray'][1]
                    if not numTMSeg in freqDict:
                        freqDict[numTMSeg] = 0
                    freqDict[numTMSeg] += 1
#}}}
def CountContinuousUnmappedTM_3(recordList, freqDict, isFilterNeighbour=True):#{{{
### one of the topology should have specifed numTM and TM2SEQ are filtered
    if isFilterNeighbour:
        str_filter = "True"
    else:
        str_filter = "False"
    for record in recordList:
        if record == {}:
            continue
        if not "mapArray" in  record:
            continue
        id1 = record['id1']
        id2 = record['id2']
        key = "%s-%s"%(id1, id2)

        mapArrayList = [record['mapArray'][i] for i in range(2)]
        seqidList = [record['id%d'%(i+1)] for i in range(2)]
        numTMList  = [len(mapArray) for mapArray in mapArrayList]
        if min(numTMList) <= 1: #ignore single spanning TM protein in location analysis
            continue
        max_st_list = [max(mapArray) for mapArray in mapArrayList]

        if ((not numTMList[0] in NUMTMSET) and (not numTMList[1] in
            NUMTMSET)):
            continue

        if max(mapArrayList[0]+mapArrayList[1]) == 2:
            continue

        for i in range(2):
            pMapArray = mapArrayList[i]

            str_maparray_list = ["%d"%x for x in pMapArray]
            str_maparray = "".join(str_maparray_list)
            neighbour_char = "0"
            st = 1

            posContList0 = GetSegPos(str_maparray, "%d"%st)
            if isFilterNeighbour:
                posContList = FilterSegPos(posContList0, str_maparray,
                        neighbour_char)
            else:
                posContList = posContList0
            if len(posContList) > 0:
                if DEBUG_UNMAPPED_TM_POSITION:
                    print >> sys.stderr, "Method 3 isFilterNeighbour=%s"%(str_filter)
                    print >> sys.stderr, "Unmapped TM position for %s - %s" %(id1, id2)
                    print >> sys.stderr, "mapArray1:", record['mapArray'][0]
                    print >> sys.stderr, "mapArray2:", record['mapArray'][1]
                    print >> sys.stderr, "posContList (nonFilt ):", posContList0
                    print >> sys.stderr, "posContList (Filtered):", posContList
                    print >> sys.stderr
                for (b,e) in posContList:
                    numTMSeg = e-b
                    if DEBUG_13 and numTMSeg == 13:
                        print >> sys.stderr
                        print >> sys.stderr, "Unmapped TM position for %s - %s" %(id1, id2)
                        print >> sys.stderr, "mapArray1:", record['mapArray'][0]
                        print >> sys.stderr, "mapArray2:", record['mapArray'][1]
                    if not numTMSeg in freqDict:
                        freqDict[numTMSeg] = 0
                    freqDict[numTMSeg] += 1
#}}}

def Ana_NumContUnmappedTM(infile):
    methodList = [0,1,2,3]
    outpath = os.path.dirname(infile)
    if outpath == "":
        outpath = "."
    try:
        freqDict = {}
        for method in methodList:
            freqDict[2*method] = {}
            freqDict[2*method+1] = {}
        unprocessedBuffer=""
        cntTotalReadInRecord = 0
        cntTotalOutputRecord = 0
        isEOFreached = False
        fpin = open(infile)
        while 1:
            buff = fpin.read(BLOCK_SIZE)
            if buff == "":
                isEOFreached = True
            buff = unprocessedBuffer + buff
            pairCmpRecordList=[]
            unprocessedBuffer = lcmp.ReadPairCmpResultFromBuffer(buff,pairCmpRecordList)
            if len(pairCmpRecordList) > 0: 
                CountContinuousUnmappedTM_0(pairCmpRecordList, freqDict[0],
                        isFilterNeighbour=True)
                CountContinuousUnmappedTM_0(pairCmpRecordList, freqDict[1],
                        isFilterNeighbour=False)
                CountContinuousUnmappedTM_1(pairCmpRecordList, freqDict[2],
                        isFilterNeighbour=True)
                CountContinuousUnmappedTM_1(pairCmpRecordList, freqDict[3],
                        isFilterNeighbour=False)
                CountContinuousUnmappedTM_2(pairCmpRecordList, freqDict[4],
                        isFilterNeighbour=True)
                CountContinuousUnmappedTM_2(pairCmpRecordList, freqDict[5],
                        isFilterNeighbour=False)
                CountContinuousUnmappedTM_3(pairCmpRecordList, freqDict[6],
                        isFilterNeighbour=True)
                CountContinuousUnmappedTM_3(pairCmpRecordList, freqDict[7],
                        isFilterNeighbour=False)
                cntTotalReadInRecord += len(pairCmpRecordList)
            if isEOFreached == True:
                break
        fpin.close()
        for method in methodList:
            for idx in [2*method, 2*method+1]:
                if idx == 2*method:
                    str_filter = "True"
                    outfile = (outpath + os.sep +
                            "tmp_ana_numContTM_method%d_filternb.txt" %
                            (method))
                else:
                    str_filter = "False"
                    outfile = (outpath + os.sep +
                            "tmp_ana_numContTM_method%d_nonfilternb.txt" %
                            (method))
                fpout = open(outfile, "w")
                print
                print >> fpout, "#numTM count Method_%d isFilterNeighbour=%s"%(method,
                        str_filter)
                for i in range(1, 21):
                    msg = "%-5d %5d"
                    try:
                        print >> fpout, msg%(i, freqDict[idx][i])
                    except KeyError:
                        print >> fpout, msg%(i, 0)
                fpout.close()
                cmd = "/data3/wk/MPTopo/src/tmp_plot_histogram_logscale.sh %s"
                os.system(cmd%(outfile))
    except IOError:
        return 1

try:
    infile = sys.argv[1]
    print "#Analyzing %s"%(infile)
    Ana_NumContUnmappedTM(infile)
except IndexError:
    progname =  os.path.basename(sys.argv[0])
    print >> sys.stderr, "Usage: %s paircmpFile"%(progname)
    sys.exit()

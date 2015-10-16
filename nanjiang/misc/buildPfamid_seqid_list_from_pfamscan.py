#!/usr/bin/env python
# Description:
import os
import sys
import tempfile
import myfunc
usage = """
usage: buildPfamid_seqid_list_from_pfamscan.py pfamscanfile -outpath DIR
    Description:
        Create two lists
    list 1:
        Pfamid1 seqid1 seqid2 seqid3 ...
        pfamid2 seqid1 seqid2 seqid3 ...
    list 2:
        seqid1 pfamid1 pfamid2 ...
        seqid2 pfamid1 pfamid2 ...     
Options:
  -q              Quiet mode
  -evalue FLOAT   Set e-value threshold, (default: 1e-3)
  -m INT          Set method, (default: 1)
                  method 0: old method, read the whole input file into memory
                  method 1: memory friendly
  -h, --help      Print this help message and exit

Created 2012-05-24, updated 2012-05-24, Nanjiang Shu 

Examples:
"""

rundir = os.path.dirname(sys.argv[0])
binpath = rundir

# ChangeLog 2014-11-24 
#   pfam-B hits has lots of over laps, not very good, for the domain
#   arrangement, show only for Pfam-A
#   output also domainlistperseq

# Note, this version use lots of memory when the pfamscan file is large.
# need to be modified 

def PrintHelp():
    print usage

def IsOverlappingDomain(thisRangeList, totalRangeList, seqid, hit):#{{{
    if len(totalRangeList) < 1 or len(thisRangeList) < 1:
        return False
    else:
        sumCover = 0
        thisSumLen = 0
        totalSumLen = 0
        for i in xrange(len(thisRangeList)):
            (b1, e1) = thisRangeList[i]
            thisSumLen += (e1-b1)
        for j in xrange(len(totalRangeList)):
            (b2, e2) = totalRangeList[j]
            totalSumLen += (e2-b2)
        for i in xrange(len(thisRangeList)):
            (b1, e1) = thisRangeList[i]
            for j in xrange(len(totalRangeList)):
                (b2, e2) = totalRangeList[j]
                sumCover += max(myfunc.coverage(b1,e1,b2,e2),0)
        percentCoverage = (sumCover / (float(thisSumLen+totalSumLen)/2.0))*100
        if percentCoverage < 30.0:
            return False
        else:
            print >> sys.stderr, "Overlapping found:", "%s %s %4d %4d" % (
                seqid, hit['pfamid'], hit['alnBeg'], hit['alnEnd']
                ), "%4d %4d %4d %.1f" % ( sumCover, thisSumLen, totalSumLen,
                        percentCoverage)
            return True
#}}}
def ReadPfamScan(infile):#{{{
    try:
        evalue_threshold = g_params['evalue_threshold']
        seqIDPfamScanDict = {}
        fpin = open(infile, "r")
        line = fpin.readline()
        while line:
            if line != "" and line[0] != "#":
                strs = line.split()
                if len(strs) >= 15:
                    seqid = strs[0]
                    tmpdict = {}
                    tmpdict['alnBeg'] = int (strs[1])
                    tmpdict['alnEnd'] = int (strs[2])
                    tmpdict['pfamid'] = strs[5].split('.')[0]
#                     tmpdict['pfamname'] = strs[6]
                    tmpdict['evalue'] = float(strs[12])
                    evalue = tmpdict['evalue']
                    tmpdict['clanid'] = strs[14]
                    if evalue <= evalue_threshold:
                        if seqid in seqIDPfamScanDict:
                            seqIDPfamScanDict[seqid].append(tmpdict)
                        else:
                            seqIDPfamScanDict[seqid] = []
                            seqIDPfamScanDict[seqid].append(tmpdict)
            line = fpin.readline()
        fpin.close()
        return seqIDPfamScanDict
    except IOError:
        print >> sys.stderr, "Failed to read file %s" %infile
        return {}
#}}}
def ReadPfamScan2(infile):#{{{
# a quick solution, to same a little memory
    evalue_threshold = g_params['evalue_threshold']
    seqIDPfamScanDict = {}
    hdl = myfunc.ReadLineByBlock(infile)
    if hdl.failure:
        return {}
    lines = hdl.readlines()

    while lines != None:
        for line in lines:
            if line != "" and line[0] != "#":
                strs = line.split()
                if len(strs) >= 15:
                    seqid = strs[0]
                    alnBeg = int (strs[1])
                    alnEnd = int (strs[2])
                    pfamid = strs[5].split('.')[0]
#                     tmpdict['pfamname'] = strs[6]
                    evalue = float(strs[12])
                    clanid = strs[14]
                    tup_info = (alnBeg, alnEnd, pfamid, clanid)
                    if evalue <= evalue_threshold:
                        if seqid in seqIDPfamScanDict:
                            seqIDPfamScanDict[seqid].append(tup_info)
                        else:
                            seqIDPfamScanDict[seqid] = []
                            seqIDPfamScanDict[seqid].append(tup_info)
        lines = hdl.readlines()
    if hdl:
        hdl.close()
    return seqIDPfamScanDict
#}}}
def CleanPfamScan(inDict):#{{{
# return cleaned dict
    newDict = {}
    for seqid in inDict:
        thisRecord = inDict[seqid]
        if len(thisRecord) == 1:
            newDict[seqid] = thisRecord
        else:
            #grouped by PFAMID
            tmppfamiddict = {}
            numHit = len(thisRecord)
            for j in range(numHit):
                pfamid = thisRecord[j]['pfamid']
                if pfamid in tmppfamiddict:
                    tmppfamiddict[pfamid].append(j)
                else:
                    tmppfamiddict[pfamid] = []
                    tmppfamiddict[pfamid].append(j)
            numPfamOfHit = len(tmppfamiddict)
            if numPfamOfHit <= 0:
                print >> sys.stderr, "Error pfamscan hit", thisRecord
            elif numPfamOfHit == 1:
                newDict[seqid] = []
                newDict[seqid].append(thisRecord[0])
            else:
                seqCoverRangeList = []
                for pfamid in tmppfamiddict:
                    if pfamid[0:2] == "PB": # ignore pfam-B for hit overlapping check
                        continue
                    hitIndexList = tmppfamiddict[pfamid]
                    thisSeqCoverRangeList = []
                    for j in hitIndexList:
                        thisSeqCoverRangeList.append((thisRecord[j]['alnBeg'],
                            thisRecord[j]['alnEnd']))
                    if not IsOverlappingDomain(thisSeqCoverRangeList,
                            seqCoverRangeList, seqid, thisRecord[j]):
                        newDict[seqid] = []
                        newDict[seqid].append(thisRecord[hitIndexList[0]])
                        seqCoverRangeList += thisSeqCoverRangeList
    return newDict
#}}}
def WriteIDMapFile(outfile, idmapdict):#{{{
    try:
        fpout = open(outfile, "w")
        for idd in idmapdict:
            fpout.write("%s %d"%(idd, len(idmapdict[idd])))
            for idd2 in idmapdict[idd]:
                fpout.write(" %s"%(idd2))
            fpout.write("\n")
        fpout.close()
        return 0
    except IOError:
        print >> sys.stderr, "Failed to write to file %s" %outfile
        return 1

#}}}
def BuildPfamSeqDB_method_0(pfamscanfile, outpath):#{{{
    """
    Quick dirty method, consume lots of memory when the pfamscanfile is large (e.g. 10GB)

    """
    rootname = os.path.basename(os.path.splitext(pfamscanfile)[0]);
#     seqIDPfamScanDict = ReadPfamScan(pfamscanfile)
    seqIDPfamScanDict = ReadPfamScan2(pfamscanfile)
#     return 0
    seqid2pfamidDict = {}
    seqid2clanidDict = {}
    pfamid2seqidDict = {}
    clanid2seqidDict = {}
    outFile_domainpattern = "%s%s%s.domainlistperseq"%(outpath, os.sep, rootname)
    fpout_domain = myfunc.myopen(outFile_domainpattern, None, "w", True) #only for pfam-A

    for seqid in seqIDPfamScanDict:
        seqid2pfamidDict[seqid] = []
        seqid2clanidDict[seqid] = []
        addedClanIDSet = set([])
        domainList = []
        domainRangeList = []
        for hit in seqIDPfamScanDict[seqid]:
            (alnBeg, alnEnd, pfamid, clanid) = hit
            if clanid == "No_clan":
                clanid = pfamid

            if pfamid[0:2] != "PB":
                domainList.append(pfamid)
                domainRangeList.append((alnBeg, alnEnd))

            seqid2pfamidDict[seqid].append(pfamid)

            if not clanid in addedClanIDSet:
                seqid2clanidDict[seqid].append(clanid)
                addedClanIDSet.add(clanid)

            if pfamid in pfamid2seqidDict:
                pfamid2seqidDict[pfamid].append(seqid)
            else:
                pfamid2seqidDict[pfamid] = []
                pfamid2seqidDict[pfamid].append(seqid)

            if clanid in clanid2seqidDict:
                clanid2seqidDict[clanid].append(seqid)
            else:
                clanid2seqidDict[clanid] = []
                clanid2seqidDict[clanid].append(seqid)

        if len(domainList) > 0:
            fpout_domain.write("%s %d"%(seqid, len(domainList)))
            for j in xrange(len(domainList)):
                fpout_domain.write(" %s,%d,%d"%(domainList[j], domainRangeList[j][0], domainRangeList[j][1]))
            fpout_domain.write("\n")


    myfunc.myclose(fpout_domain)

    outListFile1 = outpath + os.sep + rootname + ".seqid2pfamid"
    outListFile2 = outpath + os.sep + rootname + ".pfamid2seqid"
    outListFile3 = outpath + os.sep + rootname + ".seqid2clanid"
    outListFile4 = outpath + os.sep + rootname + ".clanid2seqid"

    WriteIDMapFile(outListFile1, seqid2pfamidDict)
    WriteIDMapFile(outListFile2, pfamid2seqidDict)
    WriteIDMapFile(outListFile3, seqid2clanidDict)
    WriteIDMapFile(outListFile4, clanid2seqidDict)
#     clean_seqIDPfamScanDict = CleanPfamScan(seqIDPfamScanDict)
#     #CleanPfamScan(seqIDPfamScanDict) # just check, do not change the content of the seqIDPfamScanDict
# 
#     outfile5 = outpath + os.sep + rootname + ".cleanedpfamscan"
#     fpout = open (outfile5, "w")
#     for seqid in clean_seqIDPfamScanDict: # one hit per pfamid
#         for hit in clean_seqIDPfamScanDict[seqid]:
#             print >> fpout, "%s %s %20s %8s %4d %4d %8g" % (seqid,
#                     hit['pfamid'], hit['pfamname'], hit['clanid'],
#                     hit['alnBeg'], hit['alnEnd'], hit['evalue'])
#     fpout.close()
    print "result output to"
    print "\t", outListFile1
    print "\t", outListFile2
    print "\t", outListFile3
    print "\t", outListFile4
#     print "\t", outfile5
    print "\t", outFile_domainpattern
#}}}
def BuildPfamSeqDB_method_1(pfamscanfile, outpath):#{{{
    """
    Method do not read the whole file into memory
    """
    rootname = os.path.basename(os.path.splitext(pfamscanfile)[0]);
#     seqIDPfamScanDict = ReadPfamScan(pfamscanfile)
    seqIDPfamScanDict = ReadPfamScan2(pfamscanfile)
#     return 0
    seqid2pfamidDict = {}
    seqid2clanidDict = {}
    pfamid2seqidDict = {}
    clanid2seqidDict = {}
    outFile_domainpattern = "%s%s%s.domainlistperseq"%(outpath, os.sep, rootname)
    fpout_domain = myfunc.myopen(outFile_domainpattern, None, "w", True) #only for pfam-A

    for seqid in seqIDPfamScanDict:
        seqid2pfamidDict[seqid] = []
        seqid2clanidDict[seqid] = []
        addedClanIDSet = set([])
        domainList = []
        domainRangeList = []
        for hit in seqIDPfamScanDict[seqid]:
            (alnBeg, alnEnd, pfamid, clanid) = hit
            if clanid == "No_clan":
                clanid = pfamid

            if pfamid[0:2] != "PB":
                domainList.append(pfamid)
                domainRangeList.append((alnBeg, alnEnd))

            seqid2pfamidDict[seqid].append(pfamid)

            if not clanid in addedClanIDSet:
                seqid2clanidDict[seqid].append(clanid)
                addedClanIDSet.add(clanid)

            if pfamid in pfamid2seqidDict:
                pfamid2seqidDict[pfamid].append(seqid)
            else:
                pfamid2seqidDict[pfamid] = []
                pfamid2seqidDict[pfamid].append(seqid)

            if clanid in clanid2seqidDict:
                clanid2seqidDict[clanid].append(seqid)
            else:
                clanid2seqidDict[clanid] = []
                clanid2seqidDict[clanid].append(seqid)

        if len(domainList) > 0:
            fpout_domain.write("%s %d"%(seqid, len(domainList)))
            for j in xrange(len(domainList)):
                fpout_domain.write(" %s,%d,%d"%(domainList[j], domainRangeList[j][0], domainRangeList[j][1]))
            fpout_domain.write("\n")


    myfunc.myclose(fpout_domain)

    outListFile1 = outpath + os.sep + rootname + ".seqid2pfamid"
    outListFile2 = outpath + os.sep + rootname + ".pfamid2seqid"
    outListFile3 = outpath + os.sep + rootname + ".seqid2clanid"
    outListFile4 = outpath + os.sep + rootname + ".clanid2seqid"

    WriteIDMapFile(outListFile1, seqid2pfamidDict)
    WriteIDMapFile(outListFile2, pfamid2seqidDict)
    WriteIDMapFile(outListFile3, seqid2clanidDict)
    WriteIDMapFile(outListFile4, clanid2seqidDict)
#     clean_seqIDPfamScanDict = CleanPfamScan(seqIDPfamScanDict)
#     #CleanPfamScan(seqIDPfamScanDict) # just check, do not change the content of the seqIDPfamScanDict
# 
#     outfile5 = outpath + os.sep + rootname + ".cleanedpfamscan"
#     fpout = open (outfile5, "w")
#     for seqid in clean_seqIDPfamScanDict: # one hit per pfamid
#         for hit in clean_seqIDPfamScanDict[seqid]:
#             print >> fpout, "%s %s %20s %8s %4d %4d %8g" % (seqid,
#                     hit['pfamid'], hit['pfamname'], hit['clanid'],
#                     hit['alnBeg'], hit['alnEnd'], hit['evalue'])
#     fpout.close()
    print "result output to"
    print "\t", outListFile1
    print "\t", outListFile2
    print "\t", outListFile3
    print "\t", outListFile4
#     print "\t", outfile5
    print "\t", outFile_domainpattern
#}}}

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = "./"
    pfamscanfile = ""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            isNonOptionArg = False
            pfamscanfile = argv[i]
            i += 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-evalue", "--evalue"] :
                g_params['evalue_threshold'] = float(argv[i+1])
                i += 2
            elif argv[i] in ["-m", "--m", "-method", "--method"] :
                g_params['method'] = int(argv[i+1])
                i += 2
            elif argv[i] in ["-outpath", "--outpath"] :
                outpath = argv[i+1]
                i += 2
            elif argv[i] in ["-q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            pfamscanfile = argv[i]
            i += 1

    if pfamscanfile == "":
        print >> sys.stderr, "Error! pfamscanfile not set. Exit."
        return 1


    if not os.path.exists(outpath):
        os.system("mkdir -p %s"%outpath)

    if g_params['method'] == 0:
        BuildPfamSeqDB_method_0(pfamscanfile, outpath)
    elif g_params['method'] == 1:
        BuildPfamSeqDB_method_1(pfamscanfile, outpath)
    else:
        print >> sys.stderr, "Wrong method %d"%(g_params['method'])
        return 1

    return 0

#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['method'] = 1
    g_params['evalue_threshold'] = 1e-3
    return g_params
#}}}
if __name__ == '__main__' :

    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

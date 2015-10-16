#!/usr/bin/env python
# Description:
# run hhsearch given the pairwise tableinfo file
import os
import sys
import myfunc
import itertools
progname =  os.path.basename(sys.argv[0])
usage = """
Usage: %s [-l] tableinfoFile -hhprofile DIR

Description: Search the shorter sequence in longer sequence by hhsearch,
             predict whether one is a duplicated form of the other

OPTIONS:
  -outpath    DIR   Set ouput path, default = ./
  -dupfile   FILE   Output duplication file
  -hhprofile  DIR   Multiple path can be set by evoking multiple times, ordered
  -hhsearch   DIR   When this is supplied, hhr files will first be searched in
                    the path before running
  -topofile  FILE   Topology file
  -forcewrite       force write the already exist file (default: not)
  -h, --help        Print this help message and exit

Created 2012-07-05, updated 2013-03-14, Nanjiang Shu

Examples:

# running hhsearch for pairs listed in tableinfo file
    %s -l test.tableinfo -hhprofile hhprofile -outpath hhsearch 

# predict duplications using already prepared hhsearch result
    %s -l test.tableinfo -hhprofile hhprofile -hhsearch hhsearch -outpath hhsearch -topofile data.topo -dupfile out.dup
"""%(progname, progname, progname)

EVALUE_THRESHOLD = 1e-3

def PrintHelp():
    print usage
def ReadSeqPathMapDict(infile):#{{{
    hdl = myfunc.ReadLineByBlock(infile)
    if hdl.failure:
        return {}
    dt = {}
    lines = hdl.readlines()
    while lines != None:
        for line in lines:
            if not line or  line[0] == "#":
                continue
            strs = line.split()
            if len(strs) == 2:
                dt[strs[0]] = strs[1]
        lines = hdl.readlines()
    return dt
#}}}
def ExtractHit(line):#{{{
    hit = {}
    if len(line) < 37:
        return hit
    strs = line[35:].split()
    if len(strs) >= 8:
        evalue = float(strs[1])

        strs11 = strs[6].split("-")
        posQueryList = []
        for ss in strs11:
            if ss.find("(") == -1:
                posQueryList.append(ss)
            else:
                posQueryList.append(ss.split("(")[0])

        posQuery = [int(x) for x in posQueryList]
        posQuery[0] -= 1 # let the sequence position start from zero

        strs11 = strs[7].split("-")
        posTemplateList = []
        for ss in strs11:
            if ss.find("(") == -1:
                posTemplateList.append(ss)
            else:
                posTemplateList.append(ss.split("(")[0])

        posTemplate = [int(x) for x in posTemplateList]
        posTemplate[0] -= 1

        try:
            lengthTemplate = int(line[95:].split(")")[0])
        except (IndexError, ValueError, TypeError):
            print >> sys.stderr, "hhrfile error in line \"%s\""%(line)
            return {}

        hit['evalue'] = evalue
        hit['posQuery'] = posQuery
        hit['posTemplate'] = posTemplate
        hit['lengthTemplate'] = lengthTemplate
    return hit
#}}}
def IsDuplicatedByHHSearch_obsolete(hhrfile):#{{{
    try:
        fpin = open(hhrfile,"r")
        lines = fpin.readlines()
        fpin.close()
        hitList = []
        numLine = len(lines)
        i = 0
        while i < numLine:
            line = lines[i]
            if line.find(" No Hit") == 0:
                j = 1
                while i+j < numLine and lines[i+j] != "":
                    hit = ExtractHit(lines[i+j])
                    if hit != {}:
                        hitList.append(hit)
                    else:
                        break
                    j += 1
                break
            else:
                i += 1
        if len(hitList) < 2:
            return False
        else:
            sortedHitList = sorted(hitList, key=lambda x:x['evalue'], reverse=False)
            hit1 = hitList[0]
            hit2 = hitList[1]
            if hit2['evalue'] > 1e-3:
                return False
            else:
                (b1, e1) = hit1['posTemplate']
                (b2, e2) = hit2['posTemplate']
                overlap = max(0, myfunc.coverage(b1, e1, b2, e2))
                if overlap / float(e1-b1) < 0.2 and overlap / float(e2-b2) < 0.2:
                    return True
                else:
                    return False
    except IOError:
        print >> sys.stderr, "Failed to read hhrfile %s"%hhrfile
        return False
#}}}
def IsDuplicatedByHHSearch(hhrfile, seqid1="", seqid2="", cnt=0):#{{{
    try:
        # Read in hhsearch hits
        fpin = open(hhrfile,"r")
        lines = fpin.readlines()
        fpin.close()
    except IOError:
        print >> sys.stderr, "Failed to read hhrfile %s"%hhrfile
        return False

    lengthQuery = 0
    lengthTemplate = 0
    hitList = []
    numLine = len(lines)
    i = 0
    while i < numLine:
        line = lines[i]
        if line.find("Match_columns") == 0:
            try:
                lengthQuery = int(line.split()[1])
            except (IndexError, ValueError):
                print >> sys.stderr, "Error in hhrfile %s. Ignore"%(hhrfile)
                return False
            i += 1
        elif line.find(" No Hit") == 0:
            j = 1
            while i+j < numLine and lines[i+j] != "":
                hit = ExtractHit(lines[i+j])
                if hit != {}:
                    hitList.append(hit)
                else:
                    break
                j += 1
            break
        else:
            i += 1

    isDup = False
# checking whether the template is a duplicated form of the query
    numHit = len(hitList)
    numGoodHit = 0
    if numHit >= 2: # there should be at least two hits
        sortedHitList = sorted(hitList, key=lambda x:x['evalue'], reverse=False)
        if hitList[1]['evalue'] <= 1e-3: # there should be at leave two hits with evalue < th
            lengthTemplate = hitList[0]['lengthTemplate']
            countGoodHit = 0 # there should be at least two good hits.
            posListCoverageInTemplate = [] # covered segment list in template
            idxGoodHitList = []

            for i in xrange(len(sortedHitList)):
                hit = hitList[i]
                (b1, e1) = hit['posQuery']
                if (e1-b1)/float(lengthQuery) > 0.5:
                    idxGoodHitList.append(i)
            numGoodHit = len(idxGoodHitList)
            if numGoodHit >= 2: # there should be >= 2 good Hits
# if any of two hits are not overlapping, consider it as a duplication
                for pair in itertools.combinations(idxGoodHitList, 2):
                    hit1 = hitList[pair[0]]
                    hit2 = hitList[pair[1]]
                    (b1, e1) = hit1['posTemplate']
                    (b2, e2) = hit2['posTemplate']
                    overlap = max(0, myfunc.coverage(b1, e1, b2, e2))
                    if (overlap / float(e1-b1) < 0.2) and (overlap / float(e2-b2) < 0.2):
                        isDup = True
                # if non pair are not overlapping, return false
    if isDup:
        ss_isdup = 'y'
    else:
        ss_isdup = 'n'
    sys.stdout.write("%d: %s-%s %s numHit=%d numGoodHit=%d\n" %(cnt, seqid1,
        seqid2, ss_isdup, numHit, numGoodHit))
    return isDup
#}}}
def IsDuplicated(hitList, seqlen1, seqlen2):#{{{
    """
    Check the template is a duplicated form of the query
    The input is a list of hits retrieved by HHsearch
    each hit in the hit list contains information of
    posQuery: position of the hit in the query sequence
    posTemplate: position of the hit in the template sequence
    numTM1: numTM of the hit in the query sequence
    numTM2: numTM of the hit in the template sequence
    """
# checking whether the template is a duplicated form of the query
    numHit = len(hitList)
# if any of two hits are not overlapping, consider it as a duplication
    for pair in itertools.combinations(range(numHit), 2):
        hit1 = hitList[pair[0]]
        hit2 = hitList[pair[1]]
        (b1_query, e1_query) = hit1['posQuery']
        (b2_query, e2_query) = hit2['posQuery']
        (b1_temp, e1_temp) = hit1['posTemplate']
        (b2_temp, e2_temp) = hit2['posTemplate']
        overlap_query = max(0, myfunc.coverage(b1_query, e1_query, b2_query, e2_query))
        overlap_temp = max(0, myfunc.coverage(b1_temp, e1_temp, b2_temp, e2_temp))
# if one query segmnet find two hits in two template segments and both with
# similar number of TM helices
        if (    (overlap_query / float(e1_query-b1_query) >= 0.75
                    or overlap_query / float(e2_query-b2_query) >= 0.75)
                and (overlap_temp / float(e1_temp-b1_temp) < 0.25) 
                and (overlap_temp / float(e2_temp-b2_temp) < 0.25) 
                and hitList[pair[0]]['numTM1'] > 0
                and hitList[pair[0]]['numTM2'] > 0
                and hitList[pair[1]]['numTM1'] > 0
                and hitList[pair[1]]['numTM2'] > 0
                and abs(hitList[pair[0]]['numTM2']-hitList[pair[1]]['numTM2']) <=2
                ):
            return True
    return False
#}}}
def CheckDuplication(hhrfile, seqid1, seqid2, topoDict, cnt):#{{{
    hitinfo = {}
    try:
        # Read in hhsearch hits
        fpin = open(hhrfile,"r")
        lines = fpin.readlines()
        fpin.close()
    except IOError:
        print >> sys.stderr, "Failed to read hhrfile %s"%hhrfile
        return {}

    try:
        topo1 = topoDict[seqid1]
    except KeyError:
        topo1 = ""
    try:
        topo2 = topoDict[seqid2]
    except KeyError:
        topo2 = ""


    lengthQuery = 0
    lengthTemplate = 0
    hitList = []
    numLine = len(lines)
    i = 0
    while i < numLine:
        line = lines[i]
        if line.find("Match_columns") == 0:
            try:
                lengthQuery = int(line.split()[1])
                hitinfo['seqLen1'] = lengthQuery
            except (IndexError, ValueError):
                print >> sys.stderr, "Error in hhrfile %s. Ignore"%(hhrfile)
                return {}
            i += 1
        elif line.find(" No Hit") == 0:
            j = 1
            while i+j < numLine and lines[i+j] != "":
                hit = ExtractHit(lines[i+j])
                if hit != {} and hit['evalue'] <= EVALUE_THRESHOLD:
                    posQuery = hit['posQuery']
                    posTemplate = hit['posTemplate']
                    if topo1 != "":
                        hit['numTM1'] = len(myfunc.GetTMPosition(topo1[posQuery[0]:posQuery[1]]))
                    else:
                        hit['numTM1'] = 0
                    if topo2 != "":
                        hit['numTM2'] = len(myfunc.GetTMPosition(topo2[posTemplate[0]:posTemplate[1]]))
                    else:
                        hit['numTM2'] = 0
                    hitList.append(hit)
                else:
                    break
                j += 1
            break
        else:
            i += 1

    hitList = sorted(hitList, key=lambda x:x['evalue'], reverse=False)
    hitinfo['hit'] = hitList
    numHit = len(hitList)
    if numHit < 2: # there should be at least two hits
        isDup = False

    else:
        isDup = IsDuplicated(hitList, len(topo1), len(topo2))

    if isDup:
        ss_isdup = 'y'
        hitinfo['isDup'] = 'y'
    else:
        ss_isdup = 'n'
        hitinfo['isDup'] = 'n'
    sys.stdout.write("%d: %s-%s %s numHit=%d\n" %(cnt, seqid1,
        seqid2, ss_isdup, numHit))
    return hitinfo
#}}}
def GetProfileFileName(hhprofilepathList, hhprofilepathMapDictList, seqid,#{{{
        ext):
    numList = len(hhprofilepathList)
    for i in xrange(numList):
        hhprofilepath = hhprofilepathList[i]
        hhprofilepathMapDict = hhprofilepathMapDictList[i]
        try:
            datafile = (hhprofilepath + os.sep + hhprofilepathMapDict[seqid]+
                    os.sep + seqid + ext)
            return  datafile
        except (KeyError): 
            pass
    return ""
#}}}
def RunHHSearchPairwise(tableinfoFile,  #{{{
        hhprofilepathList, hhprofilepathMapDictList,
        hhsearchpathList, hhsearchpathMapDictList,
        topoDict, outpath, dupfile):
    fpoutDup = None
    if dupfile != "":
        fpoutDup = myfunc.myopen(dupfile, sys.stdout, "w", False)

    hdl = myfunc.ReadLineByBlock(tableinfoFile)
    if hdl.failure:
        return 1
    cnt = 0
    lines = hdl.readlines()
    while lines != None:
        for line in lines:
            if not line or line[0] == "#":
                continue
            strs = line.split()
            try:
                seqid1 = strs[0]
                seqid2 = strs[1]
            except (IndexError, ValueError):
                print >> sys.stderr, "Bad record line \"%s\""%(line)
                continue

            try:
                topo1 = topoDict[seqid1]
            except KeyError:
                topo1 = ""
            try:
                topo2 = topoDict[seqid2]
            except KeyError:
                topo2 = ""

            seqlen1 = len(topo1)
            seqlen2 = len(topo2)

            pairlist = [(seqid1, seqlen1), (seqid2, seqlen2)]
            pairlist = sorted(pairlist, key=lambda x:x[1], reverse=False) # short - long

            hhrfile = "%s%s%s_%s.hhr"%(outpath, os.sep, seqid1, seqid2)
            if g_params['isUsePreBuildHHSearchResult']:
                keystr = "%s_%s"%(seqid1, seqid2)
                tmp_hhrfile = GetProfileFileName(hhsearchpathList,
                        hhsearchpathMapDictList, keystr, ".hhr")
                if os.path.exists(tmp_hhrfile):
                    hhrfile = tmp_hhrfile
                else:
                    print >> sys.stderr, "hhrfile %s does not exist in"\
                            " the prebuilt path"%(hhrfile)


            # update seqid1 and seqid2 (shorter - longer)
            seqid1 = pairlist[0][0] # shorter sequence
            seqid2 = pairlist[1][0] # longer sequence

            try:
                topo1 = topoDict[seqid1]
            except KeyError:
                topo1 = ""
            try:
                topo2 = topoDict[seqid2]
            except KeyError:
                topo2 = ""

            seqlen1 = len(topo1)
            seqlen2 = len(topo2)
            numTM1 = len(myfunc.GetTMPosition(topo1))
            numTM2 = len(myfunc.GetTMPosition(topo2))


            if not os.path.exists(hhrfile) or g_params['isForceOverWrite']:
                a3mfile = GetProfileFileName(hhprofilepathList, #query
                        hhprofilepathMapDictList, pairlist[0][0], ".a3m")
                hhmfile = GetProfileFileName(hhprofilepathList,  #template
                        hhprofilepathMapDictList, pairlist[1][0], ".hhm")
                if a3mfile == "" or not os.path.exists(a3mfile):
                    print >> sys.stderr, "a3mfile not found for %s. Ignore." %(pairlist[0][0])
                elif hhmfile == "" or not os.path.exists(hhmfile):
                    print >> sys.stderr, "hhmfile not found for %s. Ignore." %(pairlist[1][0])
                else:
                    tmp_hhrfile = "%s.tmp"%(hhrfile)
                    cmd = "hhsearch -i %s -d %s -o %s -v 0 -nocons -nopred -nodssp" % (
                            a3mfile, hhmfile, tmp_hhrfile)
                    os.system(cmd)
                    if os.path.exists(tmp_hhrfile):
                        os.system("/bin/mv -f %s %s"%(tmp_hhrfile, hhrfile))
                        print hhrfile, "output"
            if fpoutDup and os.path.exists(hhrfile):
                ss_isdup = ""
                hitinfo = {}
#                 if IsDuplicatedByHHSearch(hhrfile, seqid1, seqid2, cnt):
#                     ss_isdup = 'y'
#                 else:
#                     ss_isdup = 'n'
                hitinfo = CheckDuplication(hhrfile, seqid1, seqid2, topoDict, cnt)
                if hitinfo != {}:
                    fpoutDup.write("%s-%s %s %4d %4d %4d %4d" %(
                        seqid1, seqid2, hitinfo['isDup'],
                        len(topo1), len(topo2), numTM1, numTM2))
                    if 'hit' in hitinfo:
                        for j in xrange(len(hitinfo['hit'])):
                            hit = hitinfo['hit'][j]
                            ss_hit = "%d-%d(nTM=%d) %d-%d(nTM=%d)"%(
                                    hit['posQuery'][0], hit['posQuery'][1], hit['numTM1'],
                                    hit['posTemplate'][0], hit['posTemplate'][1], hit['numTM2'])
                            fpoutDup.write(" | %35s"%(ss_hit))
                    fpoutDup.write("\n")
            cnt += 1

        lines = hdl.readlines()
    hdl.close()
    myfunc.myclose(fpoutDup)
    return 0
#}}}
def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = "./"
    tableinfoFile = ""
    hhprofilepathList = []
    hhsearchpathList = []
    dupfile = ""
    topofile = ""
# /data3/wk/MPTopo/pfamAna_refpro/pred_topcons_single_method4/refpro20120604-celluar.selmaxlength-m1.topcons-single_topcons_single.m1.agree-44.topo

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            tableinfoFile = argv[i]
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
                (outpath, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-hhprofile", "--hhprofile", "-hhprofile1", "--hhprofile1"] :
                (ss, i) = myfunc.my_getopt_str(argv, i)
                hhprofilepathList.append(ss)
            elif argv[i] in ["-hhsearch", "--hhsearch"] :
                (ss, i) = myfunc.my_getopt_str(argv, i)
                hhsearchpathList.append(ss)
            elif argv[i] in ["-dupfile", "--dupfile"] :
                (dupfile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-topofile", "--topofile"] :
                (topofile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-l", "--l"] :
                (tableinfoFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True; i += 1
            elif argv[i] in ["-overwrite", "-forcewrite", "--forcewrite"]:
                g_params['isForceOverWrite'] = True; i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            tableinfoFile = argv[i]
            i += 1

    if tableinfoFile == "":
        print >> sys.stderr, "tableinfoFile not set. exit"
        return 1
    if len(hhprofilepathList) < 1:
        print >> sys.stderr, "hhprofilepath not set. exit"
        return 1
    if not os.path.exists(outpath):
        os.makedirs(outpath)
        if not os.path.exists(outpath):
            print >> sys.stderr, "failed to created outpath %s"%(outpath)
            return 1
    topoDict = {}
    if topofile != "":
        (idList, topoList) = myfunc.ReadFasta_without_annotation(topofile)
        for i in xrange(len(idList)):
            topoDict[idList[i]] = topoList[i]

    if dupfile != "" and topoDict == {}:
        print >> sys.stderr, "Error! topoDict is empty when dupfile"\
                "is set. Exit"
        return 1

# read in hhprofile dict
    hhprofilepathMapDictList = []
    for hhprofilepath in hhprofilepathList:
        hhprofilemapfile = hhprofilepath + os.sep + "id2pathmap.txt"
        if not os.path.exists(hhprofilemapfile):
            print >> sys.stderr, "hhprofilemapfile not exist. exit"
            hhprofilepathMapDictList.append({})
        else:
            hhprofilepathMapDictList.append(ReadSeqPathMapDict(hhprofilemapfile))
            #print
            #print hhprofilemapfile
            #print ReadSeqPathMapDict(hhprofilemapfile)

# read in index dictionary for hhsearch result file
    hhsearchpathMapDictList = []
    if len(hhsearchpathList) > 0:
        g_params['isUsePreBuildHHSearchResult'] = True
        for hhsearchpath in hhsearchpathList:
            hhsearchmapfile = hhsearchpath + os.sep + "id2pathmap.txt"
            if not os.path.exists(hhsearchmapfile):
                print >> sys.stderr, "hhsearchmapfile not exist. exit"
                hhsearchpathMapDictList.append({})
            else:
                hhsearchpathMapDictList.append(ReadSeqPathMapDict(hhsearchmapfile))

    RunHHSearchPairwise(tableinfoFile, 
            hhprofilepathList, hhprofilepathMapDictList, 
            hhsearchpathList, hhsearchpathMapDictList, topoDict,
            outpath, dupfile)

    return 0
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['hhprofilepath'] = ""
    g_params['hhprofilepathMapDict'] = {}
    g_params['isForceOverWrite'] = False
    g_params['isUsePreBuildHHSearchResult'] = False
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

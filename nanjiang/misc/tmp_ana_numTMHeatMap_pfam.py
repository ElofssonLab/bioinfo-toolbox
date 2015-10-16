#!/usr/bin/env python
# Description: analyzing the pfam families distribution of major duplicated
# numTM pairs, e.g. 6-12, 8-4, 
# also analyze the topology variation for each protein family

# 2013-05-03  works only for norm_diag

# ChangeLog 2013-12-03 
#   add option "-hhsearch" when this is supplied, hhr file will first be
#   searched in hhsearch path by using the index file id2pathmap.txt
# ChangeLog 2015-10-13 
#   add the option -mindiffpair
import os
import sys
import re
import myfunc
import libtopologycmp as lcmp
import subprocess
import inspect

DEBUG_UNMAPPED_TM_POSITION = 0
DEBUG_13 = 1

DATADIR3 = os.environ['DATADIR3']

cmpclassDefinitionList = [ "IDT", "INV", "TM2GAP", "TM2SEQ", "TM2GAP_AND_TM2SEQ" , "AllTM2GAP", "AllTM2SEQ"]
cmpClassList_method1 = ["IDT","INV","TM2GAP", "TM2SEQ", "TM2GAP_AND_TM2SEQ"]
cmpClassList_method3 = ["IDT","INV","TM2GAP", "TM2SEQ", "TM2SP", "Mixed"]
cmpClassList_mp3_cmpdup = ["IDT","INV","DUP", "TM2GAP", "TM2SEQ", "TM2SP", "Mixed"]


BLOCK_SIZE = 100000

binpath = os.path.dirname(sys.argv[0])
if binpath  == "":
    binpath = "."

progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))
usage_short="""
Usage: %s paircmpfile [-o OUTFILE]
"""%(progname)

usage_ext="""
Description:
    Analyze numTMHeatMap and frequency of pfam families of special pairs

OPTIONS:
  -o OUTFILE    Output the result to OUTFILE
  -clanidlist   FILE 
  -pfamidlist   FILE
  -seqid2clanid FILE
  -seqid2pfamid FILE
  -clanid2seqid FILE
  -pfamid2seqid FILE
  -predTMdbname STR   database name for the predicted topology
  -pfamdef      FILE
  -signalp      FILE
  -pfamtype     STR   can be all, family, domain, repeat, motif, (default: "")
  -outname STR        Set Outname, defalt: tmp_ana_numTMHeatMap_
  -pairlist FILE      Analyze only for those pairs
  -pairlistwithpfamid  FILE
       supply file for selected pairs
  -mp INT             pairwise comparion method, 0, 1, or 3 (default: 0)
                      3: helix level, TM2TM, TM2GAP, TM2SEQ, TM2SP
                      protein level, IDT, INV, TM2GAP, TM2SEQ, TM2SP, Mixed
  -winsize  INT       window size for the running average, (default: 50)
  -prokar             Do analysis for only prokaryotic proteins
  -eukar              Do analysis for only eukaryotic proteins
  -mindiffpair  INT   Set minimal number of different topology pairs to be considered as different, (default: 1)
  -q                  Quiet mode
  -h, --help          Print this help message and exit

Created 2013-05-02, updated 2015-10-13, Nanjiang Shu
"""
usage_exp="""
Examples:
"""

def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}
def ReadPairListWithFamID(infile):#{{{
    try:
        pairlistDict = {}
        pairlistSet = set([])
        fpin = open(infile,"r")
        lines = fpin.readlines()
        fpin.close()
        for line in lines:
            if line and line[0] != "#":
                strs = line.split()
                if len(strs) >= 3:
                    id1 = strs[0]
                    id2 = strs[1]
                    if (id1,id2) in pairlistSet or (id2,id1) in pairlistSet:
                        continue
                    famid = strs[2]
                    pairlistSet.add((id1,id2))
                    if not famid in pairlistDict:
                        pairlistDict[famid] = []
                    pairlistDict[famid].append([id1,id2])  
        return pairlistDict
    except IOError:
        print >> sys.stderr, "Failed to read pairlistwithclanid file ", infile
        return {}
#}}}
def ReadPfamDefFile(infile):#{{{
    try:
        fpin = open(infile, "r")
        lines = fpin.readlines()
        fpin.close()
        dtClan = {}
        dtPfam = {}
        for line in lines:
            strs = line.split("\t")
            try:
                pfamid = strs[0]
                pfamDefShort = strs[3]
                dtPfam[pfamid] = pfamDefShort

                clanid = strs[1]
                clanDefShort = strs[2]
                if clanid != "\N":
                    dtClan[clanid] = clanDefShort
                else:
                    dtClan[pfamid] = pfamDefShort
            except IndexError:
                pass
        return (dtPfam, dtClan)
    except IOError:
        print >> sys.stderr, "Failed to read file %s"%infile
        return ({}, {})
#}}}
def InitTableNumTMHeatMap(dataTable, classList, MAX_NUMTM, SPE_PAIR_LIST):#{{{
    for cls in classList:
        dataTable[cls] = {}
        dt = dataTable[cls]
        dt['data'] = []
        dt['maxNumTM'] = 0
        dt['numPair'] = 0
        dt['pairInfoLists'] = []
        for i in xrange(len(SPE_PAIR_LIST)):
            dt['pairInfoLists'].append([])   #record special pair lists
        for i in xrange(MAX_NUMTM):
            dt['data'].append([0]*MAX_NUMTM)
#}}}
def WriteNumTMHeatMap(data, maxNumTM, count, mode, outfile):#{{{
    try:
        fpout = open(outfile, "w")
        #maxNumTM = dataTable['maxNumTM']
        #count = dataTable['numPair']
        #data = dataTable['data']
        scale_norm_col_list = []
# normalized so that the sum of diagonal = 100
        if mode == "norm_diag":
            diag = [data[i][i] for i in xrange(maxNumTM+1)]
            scale_norm_diag = myfunc.FloatDivision(count,sum(diag))
        elif mode == "norm_col":
            for j in xrange(0, maxNumTM+1):
                li = [data[i][j] for i in xrange(0, maxNumTM+1)]
                scale_norm_col_list.append(myfunc.FloatDivision(count,sum(li)))
        elif mode == "no_norm":
            scale_norm_diag = 1.0

        for i in xrange(0, maxNumTM+1):
            if mode == "norm_diag":
                scale = scale_norm_diag
            elif mode == "norm_row":
                scale = myfunc.FloatDivision(count,sum(data[i]))
            for j in xrange(0, maxNumTM+1):
                if mode in ["norm_col", "norm_diag"]:
                    if mode == "norm_col":
                        scale = scale_norm_col_list[j]
                    fpout.write(" %6.3g"%(myfunc.FloatDivision(data[i][j],count)*scale*100))
                else:
                    fpout.write(" %6d"%(data[i][j]))
            fpout.write("\n")
        fpout.close()
        return 0
    except IOError:
        print >> sys.stderr, "Failed to write to file %s"%outfile
        return 1
#}}}
def IncrementSumTableWithDiffTopo(idx, countAllList, sumWindowWithDiffTopoList):#{{{
    """
    increment the sum table if there is difference in topology
    """
    mindiffpair = g_params['mindiffpair']
    cntDiffPair = 0
    isHaveDiffTopo = False
    cntList = countAllList[idx][3]
    for j in xrange(1,len(cntList)):
        cntDiffPair += cntList[j]
        if cntList[j] >= mindiffpair:
            sumWindowWithDiffTopoList[j] += 1
    if cntDiffPair >= mindiffpair:
        isHaveDiffTopo = True
    if isHaveDiffTopo:
        sumWindowWithDiffTopoList[0] += 1
#}}}
def DecrementSumTableWithDiffTopo(idx, countAllList, sumWindowWithDiffTopoList):#{{{
    """
    increment the sum table if there is difference in topology
    """
    mindiffpair = g_params['mindiffpair']
    cntDiffPair = 0
    isHaveDiffTopo = False
    cntList = countAllList[idx][3]
    for j in xrange(1,len(cntList)):
        cntDiffPair += cntList[j]
        if cntList[j] >= mindiffpair:
            sumWindowWithDiffTopoList[j] -= 1
    if cntDiffPair >= mindiffpair:
        isHaveDiffTopo = True
    if isHaveDiffTopo:
        sumWindowWithDiffTopoList[0] -= 1
#}}}

def WriteFamPairCount(freqList, pairInfoList, famDefDict, #{{{
        cmpclassList, pairwise_comparison_method, isCmpDup, outfile):
    """
    Write the number of pairs for each protein family as well as the frequency 
    of topology variations in different classes for each family
    Input:
        freqList:   [(pfamid, [numpair,numseq, numseq_TMpro, (id1,id2,cmpclass), ()...]), ...]
        pairInfoList: a list of tuples, [(id1,id2,cmpclass)]
    """
    numpair_total = len(pairInfoList)

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
# write header line
    fpout.write("#%-7s %18s %7s %7s %9s %6s %8s"%("PfamID", "PfamDef",
        "NumPair", "NumSeq", "NumSeq_TM", "%#", "NumTotal"))
    for ss in cmpclassList:
        fpout.write(" %6s"%(ss))
    fpout.write("\n")

    CountListAll = [] # list of [(numpair, numseq, numseq_TMpro, [#IDT, #INV, #...]), ([])]

    for tup in freqList:
        famid = tup[0]
        pairInfoList_fam = tup[1] #(pfamid, [numpair, numseq, numseq_TMpro, (id1,id2,cmpclass),(id1,id2,cmpclass)])
        numpair = tup[1][0]
        numseq = tup[1][1]
        numseq_TMpro = tup[1][2]
        try:
            famdef = famDefDict[famid]
        except KeyError:
            famdef = ""

        cmpclasslist_fam = [] # a list of [cmpclass, cmpclass, ...]
        for tt in pairInfoList_fam[3:]:
            cmpclass = tt[2]
            if isCmpDup:
                if cmpclass.find("TM2GAP|DUP") == 0:
                    cmpclass = "DUP"
                else:
                    cmpclass = cmpclass.split('|')[0]
            else:
                cmpclass = cmpclass.split('|')[0]


            cmpclasslist_fam.append(cmpclass)
        cntList = []
        for cls in cmpclassList:
            cntList.append(cmpclasslist_fam.count(cls))
        fpout.write("%-8s %18s %7d %7d %9d %6.2f %8d"%(famid, famdef, numpair,numseq,
            numseq_TMpro, float(numpair)/numpair_total*100.0, numpair_total))

        for cnt in cntList:
            fpout.write(" %6d"%cnt)
        fpout.write("\n")

        CountListAll.append((numpair, numseq, numseq_TMpro, cntList))
    myfunc.myclose(fpout)
    print "file %s output"%(outfile)

# output the fraction of topN largest families that have topology variations
# sorted in descending order by "numpair", "numseq", numseq_TMpro #{{{
    mindiffpair = g_params['mindiffpair']
    itemList = ["numpair", "numseq", "numseq_TMpro"]
    outfileList = []
    for pivIdx in xrange(len(itemList)): #[0,1,2]

        outList = sorted(CountListAll, key=lambda x:x[pivIdx], reverse=True)

        freqTopNList = [] # [ [topN, min, frac_DIFF, frac_INV, frac_TM2GAP]]
        sumTopNWithDiffTopoList = [0]*(len(cmpclassList))
        for i in xrange(len(outList)):
            isHaveDiffTopo = False
            cntDiffPair = 0
            cntList = outList[i][3]
            for j in xrange(1,len(cntList)):
                cntDiffPair += cntList[j]
                if cntList[j] >= mindiffpair :
                    sumTopNWithDiffTopoList[j]+=1

            if cntDiffPair >= mindiffpair:
                isHaveDiffTopo = True

            if isHaveDiffTopo:
                sumTopNWithDiffTopoList[0] += 1
            fracList = []
            for j in xrange(len(sumTopNWithDiffTopoList)):
                fracList.append(myfunc.FloatDivision(sumTopNWithDiffTopoList[j], i+1))

            freqTopNList.append( [i+1, outList[i][pivIdx]]+ fracList)
        outfile1 = outfile + ".topNdifffam.sortby_%s.mindiffpair_%d.txt"%(itemList[pivIdx], mindiffpair)
        outfileList.append(outfile1)
        fpout = myfunc.myopen(outfile1, sys.stdout, "w", False)

        ss_sort_item = "min_%s"%(itemList[pivIdx])
        fpout.write("#%-7s %*s %7s"%("topN", len(ss_sort_item), ss_sort_item, "DIFF"))
        for ss in cmpclassList[1:]:
            fpout.write(" %7s"%(ss))
        fpout.write("\n")
        for i in xrange(len(freqTopNList)):
            d = freqTopNList[i]  #  [topN, min, frac_DIFF, frac_INV, frac_TM2GAP]
            fpout.write("%-8d %*d"%(d[0], len(ss_sort_item), d[1]))
            fracList = d[2:]
            for tt in fracList:
                fpout.write(" %7.2f"%(tt*100))
            fpout.write("\n")
        myfunc.myclose(fpout)
        print "file %s output"%(outfile1)
# make plot
    cmd = ["%s/plotMaxFracFamilyWithTopoVariation.sh"%(binpath)] + outfileList
    try:
        subprocess.check_output(cmd)
    except subprocess.CalledProcessError, e:
        print e
#}}}

# using running average
# sorted in descending order by "numpair", "numseq", numseq_TMpro #{{{
    itemList = ["numpair", "numseq", "numseq_TMpro"]
    winsize = g_params['winsize']
    outfileList = []
    for pivIdx in xrange(len(itemList)): #[0,1,2]

        outList = sorted(CountListAll, key=lambda x:x[pivIdx], reverse=True)

        freqTopNList = [] # [ [idxWindow, min, frac_DIFF, frac_INV, frac_TM2GAP]]
        sumWindowWithDiffTopoList = [0]*(len(cmpclassList)) # [#DIFF, #INV, ...]
        isFirstWindow = True # whether it is the first window
        for i in xrange(len(outList)-winsize+1): #iterating over windows
            if isFirstWindow:
                for iw in xrange(winsize):
                    IncrementSumTableWithDiffTopo(i+iw, outList, sumWindowWithDiffTopoList)
                isFirstWindow = False
            else: # if not first window, minus the previous one, and plus the next one
                #iw_previous = i-1
                #iw_next = i+winsize-1
                DecrementSumTableWithDiffTopo(i-1, outList, sumWindowWithDiffTopoList)
                IncrementSumTableWithDiffTopo(i+winsize-1, outList, sumWindowWithDiffTopoList)

            fracList = []
            for j in xrange(len(sumWindowWithDiffTopoList)):
                fracList.append(myfunc.FloatDivision(sumWindowWithDiffTopoList[j], winsize))
            freqTopNList.append( [i+1, outList[i+winsize/2][pivIdx]]+ fracList)
        outfile1 = outfile + ".difffam_win%d.sortby_%s.mindiffpair_%d.txt"%(winsize, itemList[pivIdx], mindiffpair)
        outfileList.append(outfile1)
        fpout = myfunc.myopen(outfile1, sys.stdout, "w", False)

        ss_sort_item = "min_%s"%(itemList[pivIdx])
        fpout.write("#%-7s %*s %7s"%("idxWin", len(ss_sort_item), ss_sort_item, "DIFF"))
        for ss in cmpclassList[1:]:
            fpout.write(" %7s"%(ss))
        fpout.write("\n")
        for i in xrange(len(freqTopNList)):
            d = freqTopNList[i]  #  [topN, min, frac_DIFF, frac_INV, frac_TM2GAP]
            fpout.write("%-8d %*d"%(d[0], len(ss_sort_item), d[1]))
            fracList = d[2:]
            for tt in fracList:
                fpout.write(" %7.2f"%(tt*100))
            fpout.write("\n")
        myfunc.myclose(fpout)
        print "file %s output"%(outfile1)
# make plot
    cmd = ["%s/plotMaxFracFamilyWithTopoVariation.sh"%(binpath)] + outfileList
    try:
        subprocess.check_output(cmd)
    except subprocess.CalledProcessError, e:
        print e
#}}}

#}}}
def WriteSpecialPair(dataTable, all_pairInfoList,  #{{{
        seqid2pfamidDict, seqid2clanidDict,
        tm_pfamidSet, tm_clanidSet, 
        pfamidDefDict, clanidDefDict, 
        SPE_PAIR_LIST, outfile):
    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    numpair_total = len(all_pairInfoList)
    pairInfoLists = dataTable['pairInfoLists']
    (freqListPfam, freqListClan) = AnaFamFrequency(pairInfoLists, seqid2pfamidDict,
            seqid2clanidDict, tm_pfamidSet, tm_clanidSet)
    for i in xrange(len(SPE_PAIR_LIST)):
        cnt_this_pair_pfam = 0
        cnt_this_pair_clan = 0
        for tup in freqListPfam[i]:
            cnt_this_pair_pfam += tup[1][0]
        for tup in freqListClan[i]:
            cnt_this_pair_clan += tup[1][0]

        print >> fpout
        pair = SPE_PAIR_LIST[i]
        print >> fpout, pair, "PfamID", "%5d %5.1f %8d %6.2f"%(cnt_this_pair_pfam,
                cnt_this_pair_pfam*g_params['scale_count'], numpair_total,
                myfunc.FloatDivision(cnt_this_pair_pfam, numpair_total)*100)
        print >> fpout
        for tup in freqListPfam[i]:
            try:
                pfamdef = pfamidDefDict[tup[0]]
            except KeyError:
                pfamdef = ""
            fpout.write("%-8s %20s %5d %5.1f %6.2f  "%(tup[0], pfamdef,
                tup[1][0], tup[1][0]*g_params['scale_count'],
                float(tup[1][0])/cnt_this_pair_pfam*100))
            for pp in tup[1][1]:
                fpout.write("(%s %s) "%(pp[0],pp[1]))
            fpout.write("\n")

        print  >> fpout
        pair = SPE_PAIR_LIST[i]
        print >> fpout, pair, "ClanID", "%5d %5.1f %8d %6.2f"%(cnt_this_pair_clan,
                cnt_this_pair_clan*g_params['scale_count'], numpair_total,
                myfunc.FloatDivision(cnt_this_pair_clan, numpair_total)*100)
        print >> fpout
        for tup in freqListClan[i]:
            try:
                clandef = clanidDefDict[tup[0]]
            except KeyError:
                clandef = ""
            fpout.write("%-8s %20s %5d %5.1f %6.2f  "%(tup[0], clandef, tup[1][0],
                tup[1][0]*g_params['scale_count'],
                float(tup[1][0])/cnt_this_pair_clan*100))
            for pp in tup[1][1]:
                fpout.write("(%s %s) "%(pp[0],pp[1]))
            fpout.write("\n")
        print >> fpout, "#====================================================="
    myfunc.myclose(fpout)
    return 0
#}}}
def FillSymmetricDataTableNumTMHeatMap(dataTable, classList):#{{{
    for cls in classList:
        dt = dataTable[cls]
        mtx = dt['data']
        maxNumTM = dt['maxNumTM']
        for i in xrange(maxNumTM):
            for j in xrange(i+1, maxNumTM):
                mtx[j][i] = mtx[i][j]
#}}}
def AnaFamFrequency_onelist_sub(id1, id2, cmpclass, seqid2pfamidDict,#{{{
        pfamid2seqidDict, tm_pfamidSet, idSet_TMpro, usedPfamIDSet,
        freqDictPfam):
    """
    sub-procedure for the function AnaFamFrequency_onelist
    Input:
        id1
        id2
        seqid2pfamidDict
        pfamid2seqidDict
        tm_pfamidSet        Set of ids for TM protein families
        idSet_TMpro         Set of seqids for TM proteins
        usedPfamIDSet       pfamid used in pair selection
    Output:
        freqDictPfam
    """
    try:
        pfamidlist1 = seqid2pfamidDict[id1]
        pfamidlist2 = seqid2pfamidDict[id2]
        common_pfamidlist = list(set(pfamidlist1) & set(pfamidlist2))

        for pfamid in common_pfamidlist:
            if pfamid.find("CL") == -1 and pfamid not in usedPfamIDSet:
# ignore pfamid that are not used in pair selection
                continue
            try:
                idlist = pfamid2seqidDict[pfamid]
                numseq = len(idlist)
                numseq_TMpro = len(idSet_TMpro & set(idlist))
            except KeyError:
                print >> sys.stderr, "%s not in pfamid2seqidDict"%(pfamid)
                numseq = -1
                numseq_TMpro = -1
            if pfamid in tm_pfamidSet:
                if not pfamid in freqDictPfam:
                    freqDictPfam[pfamid] = [0, numseq, numseq_TMpro]
                freqDictPfam[pfamid].append((id1,id2,cmpclass))
                freqDictPfam[pfamid][0] += 1
    except KeyError:
        msg = "%s - %s not found in seqid2pfamidDict"
        print >> sys.stderr, msg%(id1, id2)
        pass
#}}}
def AnaFamFrequency_onelist(pairInfoList,  #{{{
        seqid2pfamidDict, seqid2clanidDict, 
        pfamid2seqidDict, clanid2seqidDict, 
        tm_pfamidSet, tm_clanidSet, idSet_TMpro, usedPfamIDSet):
    """
    Get the frequency of topology variation in different classes for each
    protein family

    Input:
        pairInfoList:       a list of tuples of (id1, id2, cmpclass)
        seqid2pfamidDict:   dictionary, seqid2pfamidDict[id] = [pfamid]
        seqid2clanidDict:   dictionary, seqid2clanidDict[id] = [clanid]
        tm_pfamidSet:       set of pfamids for TM proteins
        tm_clanidSet:       set of clanids for TM proteins
        idSet_TMpro:        set of seqids for all predicted TM proteins
    Output:
        (list1, list2)
        list1:  list of frequency for families each item is [numPair, (id1, id2, cmpclass)]
        list1:  list of frequency for clans, each item is [numPair, (id1, id2, cmpclass)]
    """
    freqDictPfam = {}
    freqDictClan = {}
    # data structure of freqDictPfam
    #  {pfamid: [numPair, numSeq, numTMseq, (id1_1, id2_1, cmpclass_1), (id1_2,
    #  id2_2, cmpclass_2), ...]}
# Note that numTMseq is not for all number of TM proteins, but defined in the
    # topofile, for example, those predicted with TOPCONS-single 4/4. This
    # number can be much smaller than the actual number of TM proteins within
    # the family
# format of freqDictPfam
# dt['pfamid'] = [5, pairinfolist]
    for tup in pairInfoList:
        id1 = tup[0]
        id2 = tup[1]
        cmpclass  = tup[2]

        AnaFamFrequency_onelist_sub(id1, id2, cmpclass, seqid2pfamidDict,
                pfamid2seqidDict, tm_pfamidSet, idSet_TMpro, usedPfamIDSet,
                freqDictPfam)
        AnaFamFrequency_onelist_sub(id1, id2, cmpclass, seqid2clanidDict,
                clanid2seqidDict, tm_clanidSet, idSet_TMpro, usedPfamIDSet,
                freqDictClan)
    # sort the dictionary in descending order by the number of pairs in the
    # family
    list1 = sorted(freqDictPfam.items(), key=lambda x:x[1][2],
        reverse=True)
    list2 = sorted(freqDictClan.items(), key=lambda x:x[1][2],
        reverse=True)
    # after sorting, the data structure of the list is
    # [(pfamid, [numPair, numseq, numseq_TMpro, (id1,id2,cmpclass), (id1,id2,cmpclass)...]), (pfamid...)]
    return (list1, list2)
#}}}
def AnaFamFrequency(pairInfoLists, seqid2pfamidDict, seqid2clanidDict, #{{{
        tm_pfamidSet, tm_clanidSet):
    N = len(pairInfoLists)
    freqListPfam = []
    freqListClan = []
    for i in xrange(N):
        freqDictPfam = {}
        freqDictClan = {}
        for tup in pairInfoLists[i]:
            id1 = tup[0]
            id2 = tup[1]
            try:
                pfamidlist1 = seqid2pfamidDict[id1]
                pfamidlist2 = seqid2pfamidDict[id2]
                common_pfamidlist = list(set(pfamidlist1) & set(pfamidlist2))

                for pfamid in common_pfamidlist:
                    if pfamid in tm_pfamidSet:
                        if not pfamid in freqDictPfam:
                            freqDictPfam[pfamid] = [0, []]
                        freqDictPfam[pfamid][0] += 1
                        freqDictPfam[pfamid][1].append((id1, id2))
            except KeyError:
                print >> sys.stderr, "%s - %s not found in seqid2pfamidDict"%(id1, id2)
                pass

            try:
                clanidlist1 = seqid2clanidDict[id1]
                clanidlist2 = seqid2clanidDict[id2]
                common_clanidlist = list(set(clanidlist1) & set(clanidlist2))
                for clanid in common_clanidlist:
                    if clanid in tm_clanidSet:
                        if not clanid in freqDictClan:
                            freqDictClan[clanid] = [0, []]
                        freqDictClan[clanid][0] += 1
                        freqDictClan[clanid][1].append((id1, id2))
            except KeyError:
                print >> sys.stderr, "%s - %s not found in seqid2clanidDict"%(id1, id2)
                pass
        freqListPfam.append(sorted(freqDictPfam.items(), key=lambda x:x[1][0],
            reverse=True))
        freqListClan.append(sorted(freqDictClan.items(), key=lambda x:x[1][0],
            reverse=True))

    return (freqListPfam, freqListClan)
#}}}
def AnaPairCmpResultNumTMHeatMap(recordList, dataTable, pairInfoListDict, #{{{
        classList, signalpDict, SPE_PAIR_LIST, alignrange):
    """
    Input:
        recordList:     list of paircmp
        alignrange:     FULL_ALIGNED or PART_ALIGNED
        classList:      ALL or RMSP
        signalpDict:    dictionary of signal peptide. signalpDict[id] = INT
        SPE_PAIR_LIST:  list of tuples of (numTM, numTM) 
    Output:
        dataTable:        data for storing statistis of counts for (numTM, numTM)
        pairInfoListDict: pairInfoListDict["ALL"] is a list of tuples of (id1,
                          id2, cmpclass)
    """
    isOnlyAnaProkar = g_params['isOnlyAnaProkar']
    isOnlyAnaEukar = g_params['isOnlyAnaEukar']
    prokarSeqIDSet = g_params['prokarSeqIDSet']
    eukarSeqIDSet = g_params['eukarSeqIDSet']

    for record in recordList:
        if record == {}:
            continue

        cmpclass = record['cmpclass']

        if (cmpclass.find('UNALIGNED') == 0 or cmpclass.find('AMBIGUOUS') == 0):
            continue
        if record['isLocalAlignment'] and alignrange != 'all':
            if record['alignrange'] != alignrange:
                continue
        try:
            cmpclass = record['cmpclass']
            numTM1 = record['numTM1']
            numTM2 = record['numTM2']
            id1 = record['id1']
            id2 = record['id2']
        except KeyError:
            print >> sys.stderr, "KeyError for the record %s - %s" %(
                    record['id1'], record['id2'])
            continue

        if isOnlyAnaProkar:
            if not (id1 in prokarSeqIDSet and id2 in prokarSeqIDSet):
                continue
        if isOnlyAnaEukar:
            if not (id1 in eukarSeqIDSet and id2 in eukarSeqIDSet):
                continue

        for cls in classList:
            dt = dataTable[cls]
            data = dt['data']
            if cls == "ALL":
                isAdd = True
            elif cls == "RMSP":
                if id1 in signalpDict:
                    isAdd = False
                elif id2 in signalpDict:
                    isAdd = False
                else:
                    isAdd = True
            elif cls == "RMDUP":
                isAdd = False
            else:
                isAdd = False



            if isAdd:
                pairInfoListDict[cls].append((id1, id2, cmpclass))
                data[min(numTM1,numTM2)][max(numTM1, numTM2)] += 1
                if max(numTM1, numTM2) > dt['maxNumTM']:
                    dt['maxNumTM'] = max(numTM1, numTM2)

                pair = (min(numTM1, numTM2), max(numTM1,numTM2))
                try:
                    idx = SPE_PAIR_LIST.index(pair)
                    dt['pairInfoLists'][idx].append((id1,id2))
                except (ValueError, IndexError):
                    pass
                dt['numPair'] += 1
    return 0
#}}}
def CountSpecialPair(recordList, pairInfoLists, SPE_PAIR_LIST):#{{{
### one of the topology should be all mapped and the other one should only have TM2GAP
    for record in recordList:
        if record == {}:
            continue
        id1 = record['id1']
        id2 = record['id2']
        numTM1 = record['numTM1']
        numTM2 = record['numTM2']
        pair = (min(numTM1, numTM2), max(numTM1,numTM2))
        try:
            idx = SPE_PAIR_LIST.index(pair)
            pairInfoLists[idx].append((id1,id2))
        except (ValueError, IndexError):
            pass
#}}}
def Ana_NumTMHeatMap(infile, seqid2pfamidDict, seqid2clanidDict,  #{{{
        tm_pfamidSet, tm_clanidSet, pfamidDefDict, clanidDefDict,
        signalpDict,classList_TableNumTMHeatMap, SPE_PAIR_LIST,
        pfamid2seqidDict, clanid2seqidDict, idSet_TMpro, usedPfamIDSet,
        alignrange):

    dataTableNumTMHeatMap = {}
    InitTableNumTMHeatMap(dataTableNumTMHeatMap, classList_TableNumTMHeatMap,
            100, SPE_PAIR_LIST)
    pairInfoListDict = {}
    for cls in classList_TableNumTMHeatMap:
        pairInfoListDict[cls] = []

    if g_params['outpath'] != "":
        outpath = g_params['outpath']
    else:
        outpath = os.path.dirname(infile)
        if outpath == "":
            outpath = "."

    try:
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
                #CountSpecialPair(pairCmpRecordList, pairInfoLists, SPE_PAIR_LIST)
                AnaPairCmpResultNumTMHeatMap(pairCmpRecordList,
                        dataTableNumTMHeatMap, pairInfoListDict,
                        classList_TableNumTMHeatMap, signalpDict,
                        SPE_PAIR_LIST, alignrange)
                cntTotalReadInRecord += len(pairCmpRecordList)
                print "cntTotalReadInRecord = ", cntTotalReadInRecord
            if isEOFreached == True:
                break
        fpin.close()

        #print "count 6,12 = ", dataTableNumTMHeatMap['RMSP']['data'][6][12]
        #print "count = ", dataTableNumTMHeatMap['RMSP']['data']
        #print "count special", dataTableNumTMHeatMap['RMSP']['pairInfoLists']
        try:
            for i in xrange(len(dataTableNumTMHeatMap['RMSP']['pairInfoLists'])):
                print SPE_PAIR_LIST[i], len(dataTableNumTMHeatMap['RMSP']['pairInfoLists'][i])
        except KeyError:
            pass


#         if g_params['numTMHeatMapMode'] == "full":
#             FillSymmetricDataTableNumTMHeatMap(dataTableNumTMHeatMap,
#                     classList_TableNumTMHeatMap)

        for cls in classList_TableNumTMHeatMap: # ["ALL", "RMSP"]
            (freqListPfam, freqListClan) = AnaFamFrequency_onelist(
                    pairInfoListDict[cls], seqid2pfamidDict, seqid2clanidDict,
                    pfamid2seqidDict, clanid2seqidDict, tm_pfamidSet,
                    tm_clanidSet, idSet_TMpro, usedPfamIDSet)
            #print "tm_clanidSet", tm_clanidSet
            if g_params['pairwise_comparison_method'] == 1:
                cmpclassList = cmpClassList_method1
            elif g_params['pairwise_comparison_method'] == 3:
                cmpclassList = cmpClassList_method3
            isCmpDup = False
            outFileFamPairCount = "%s%s%s.%s.%s.pfam.paircount.txt"%(
                    outpath, os.sep, g_params['outname'], alignrange, cls)
            WriteFamPairCount(freqListPfam, pairInfoListDict[cls], pfamidDefDict,
                    cmpclassList, g_params['pairwise_comparison_method'], 
                    isCmpDup, outFileFamPairCount)
            outFileFamPairCount = "%s%s%s.%s.%s.clan.paircount.txt"%(
                    outpath, os.sep, g_params['outname'], alignrange, cls)
            WriteFamPairCount(freqListClan, pairInfoListDict[cls], clanidDefDict,
                    cmpclassList, g_params['pairwise_comparison_method'],
                    isCmpDup, outFileFamPairCount)

            if g_params['pairwise_comparison_method'] == 3: 
                # if mp=3, write another statistics with cmpdup
                isCmpDup = True
                cmpclassList = cmpClassList_mp3_cmpdup
                outFileFamPairCount = "%s%s%s.%s.%s.cmpdup.pfam.paircount.txt"%(
                        outpath, os.sep, g_params['outname'], alignrange, cls)
                WriteFamPairCount(freqListPfam, pairInfoListDict[cls], pfamidDefDict,
                        cmpclassList, g_params['pairwise_comparison_method'],
                        isCmpDup, outFileFamPairCount)
                outFileFamPairCount = "%s%s%s.%s.%s.cmpdup.clan.paircount.txt"%(
                        outpath, os.sep, g_params['outname'], alignrange, cls)
                WriteFamPairCount(freqListClan, pairInfoListDict[cls], clanidDefDict,
                        cmpclassList, g_params['pairwise_comparison_method'],
                        isCmpDup, outFileFamPairCount)

            for mode_norm in ["norm_diag", "no_norm"]:
                if mode_norm in ["norm_diag", "no_norm"]:
                    heatmapmode = 'half'
                else:
                    heatmapmode = 'full'
                outFileNumTMHeatMap = "%s%s%s.%s.%s.%s.%s.txt"%(outpath, os.sep,
                        g_params['outname'], alignrange, heatmapmode, cls, mode_norm)
                if heatmapmode == 'full':
                    mtx = myfunc.FillSymmetricMatrix(
                            dataTableNumTMHeatMap[cls]['data'],
                            dataTableNumTMHeatMap[cls]['maxNumTM'])
                else:
                    mtx = dataTableNumTMHeatMap[cls]['data']

                if mode_norm == "no_norm":
                    for i in range(dataTableNumTMHeatMap[cls]['maxNumTM']):
                        mtx[i][i] = 0

                if WriteNumTMHeatMap(mtx,
                        dataTableNumTMHeatMap[cls]['maxNumTM'],
                        dataTableNumTMHeatMap[cls]['numPair'], mode_norm,
                        outFileNumTMHeatMap) == 0:
                    print "heatmap %s output"%(outFileNumTMHeatMap) 
                    cmd = "%s/plotNumTMHeatMap.sh %s" %(binpath,
                            outFileNumTMHeatMap)
                    os.system(cmd)
                outFileSpecialPairAna = "%s%s%s.%s.%s.%s.%s.specialpairana.txt"%(
                        outpath, os.sep, g_params['outname'],alignrange,
                        g_params['numTMHeatMapMode'], cls, mode_norm)
                WriteSpecialPair(dataTableNumTMHeatMap[cls], pairInfoListDict[cls],
                        seqid2pfamidDict, seqid2clanidDict, tm_pfamidSet,
                        tm_clanidSet, pfamidDefDict, clanidDefDict, 
                        SPE_PAIR_LIST, outFileSpecialPairAna)
                print "Anafile %s output"%(outFileSpecialPairAna)



#         for i in xrange(len(SPE_PAIR_LIST)):
#             print
#             pair = SPE_PAIR_LIST[i]
#             print pair
#             print len(pairInfoLists[i])
#             print pairInfoLists[i]

    except IOError:
        return 1
#}}}


def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    SPE_PAIR_LIST = [
            (2,1),
            (2,4),
            (2,6),
            (2,8),
            (3,6),
            (3,7),
            (4,6),
            (4,8),
            (4,10),
            (5,7),
            (5,10),
            (6,8),
            (6,10),
            (6,12),
            (7,14),
            (8,10),
            (8,12),
            (10,12),
            (10,13),
            (11,13),
            (12,14)
            ]

    outfile = ""

    infile = ""
    pfamDefFile = "%s/data/pfam/pfam26.0/Pfam-A.clans.tsv"%(DATADIR3)
    signalpFile = "%s/wk/MPTopo/pfamAna_refpro/pred_signalp/refpro20120604-celluar.selmaxlength-m1.nr100.signalp_list"%(DATADIR3)

    #seqid2clanidMapFile = "%s/wk/MPTopo/pfamAna_refpro/pfammap_from_uniprot/refpro20120604-celluar.selmaxlength-m1.nr100.filter.fragmented.seqid2clanid"%(DATADIR3)
    #seqid2pfamidMapFile = "%s/wk/MPTopo/pfamAna_refpro/pfammap_from_uniprot/refpro20120604-celluar.selmaxlength-m1.nr100.filter.fragmented.seqid2pfamid"%(DATADIR3)
    seqid2clanidMapFile = ""
    seqid2pfamidMapFile = ""
    tm_pfamidListFile = ""
    tm_clanidListFile = ""
    pfamid2seqidMapFile = ""
    clanid2seqidMapFile = ""
    dbname_predTM = ""
    pairlistwithpfamidFile = ""

    pfamtype = ""

    pairListFile = ""

    #classList_TableNumTMHeatMap = ["ALL", "RMSP"] 
    classList_TableNumTMHeatMap = ["ALL"] 

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
            elif argv[i] in ["-o", "--o", "-outfile"]:
                (outfile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-outpath", "--outpath"]:
                (g_params['outpath'], i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-l", "--l"] :
                (fileListFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-pfamdef", "--pfamdef"] :
                (pfamDefFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-signalp", "--signalp"] :
                (signalpFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-mp", "--mp"]:
                g_params['pairwise_comparison_method'], i = myfunc.my_getopt_int(argv, i)
            elif argv[i] in ["-mindiffpair", "--mindiffpair"]:
                g_params['mindiffpair'], i = myfunc.my_getopt_int(argv, i)
            elif argv[i] in ["-pfamtype", "--pfamtype"]:
                pfamtype, i = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-clanidlist", "--clanidlist"] :
                (tm_clanidListFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-pfamidlist", "--pfamidlist"] :
                (tm_pfamidListFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-seqid2clanid", "--seqid2clanid"] :
                (seqid2clanidMapFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-seqid2pfamid", "--seqid2pfamid"] :
                (seqid2pfamidMapFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-pfamid2seqid", "--pfamid2seqid"] :
                (pfamid2seqidMapFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-clanid2seqid", "--clanid2seqid"] :
                (clanid2seqidMapFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-pairlistwithpfamid", "--pairlistwithpfamid"] :
                (pairlistwithpfamidFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-predTMdbname", "--predTMdbname"] :
                (dbname_predTM, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-pairlist", "--pairlist"] :
                (pairListFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-winsize", "--winsize"] :
                (g_params['winsize'], i) = myfunc.my_getopt_int(argv, i)
            elif argv[i] in ["-outname", "--outname"] :
                (g_params['outname'], i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True; i += 1
            elif argv[i] in ["-prokar", "--prokar"]:
                g_params['isOnlyAnaProkar'] = True; i += 1
            elif argv[i] in ["-eukar", "--eukar"]:
                g_params['isOnlyAnaEukar'] = True; i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            infile = argv[i]
            i += 1

    if myfunc.checkfile(infile, "%s (line %d): infile"%( __file__,
        inspect.currentframe().f_lineno)) != 0:
        return 1

    dirpath = myfunc.my_dirname(infile)


    # try to obtain Pfam family tag
    tag = ""
    if pfamtype != "":
        if pfamtype.upper().find("FAM") != -1:
            tag = ".Family"
        elif pfamtype.upper().find("DOM") != -1:
            tag = ".Domain"
        elif pfamtype.upper().find("REP") != -1:
            tag = ".Repeat"
        elif pfamtype.upper().find("MOT") != -1:
            tag = ".Motif"
        else:
            tag = ""
    else:
        if infile.find(".Family.") != -1:
            tag = ".Family"
        elif infile.find(".Domain.") != -1:
            tag = ".Domain"
        elif infile.find(".Repeat.") != -1:
            tag = ".Repeat"
        elif infile.find(".Motif.") != -1:
            tag = ".Motif"
        else:
            tag = ""

    if seqid2clanidMapFile == "":
        seqid2clanidMapFile = "%s/wk/MPTopo/pfamAna_refpro/pfammap_from_uniprot/Pfam-A-full.seqfrompfamfasta.percentTMpro_scampi.perTM75_nseq20.nr100.filter.fragmented.seqid2clanid"%(DATADIR3)
    if myfunc.checkfile(seqid2clanidMapFile, "%s (line %d): seqid2clanidMapFile"%( __file__,
        inspect.currentframe().f_lineno)):
        return 1

    if seqid2pfamidMapFile == "":
        seqid2pfamidMapFile = "%s/wk/MPTopo/pfamAna_refpro/pfammap_from_uniprot/Pfam-A-full.seqfrompfamfasta.percentTMpro_scampi.perTM75_nseq20%s.nr100.filter.fragmented.seqid2pfamid"%(DATADIR3, tag)
    if myfunc.checkfile(seqid2pfamidMapFile, "%s (line %d): seqid2pfamidMapFile"%( __file__,
        inspect.currentframe().f_lineno)):
        return 1

    if pfamid2seqidMapFile == "":
        pfamid2seqidMapFile = "%s/wk/MPTopo/pfamAna_refpro/pfammap_from_uniprot/Pfam-A-full.seqfrompfamfasta.percentTMpro_scampi.perTM75_nseq20.nr100.filter.fragmented.pfamid2seqid"%(DATADIR3)
    if myfunc.checkfile(pfamid2seqidMapFile, "%s (line %d): pfamid2seqidMapFile"%( __file__,
        inspect.currentframe().f_lineno)):
        return 1

    if clanid2seqidMapFile == "":
        clanid2seqidMapFile = "%s/wk/MPTopo/pfamAna_refpro/pfammap_from_uniprot/Pfam-A-full.seqfrompfamfasta.percentTMpro_scampi.perTM75_nseq20%s.nr100.filter.fragmented.clanid2seqid"%(DATADIR3, tag)
    if myfunc.checkfile(clanid2seqidMapFile, "%s (line %d): clanid2seqidMapFile"%( __file__,
        inspect.currentframe().f_lineno)):
        return 1


    if tm_pfamidListFile == "":
        tm_pfamidListFile = "%s/data/pfam/pfam26.0/Pfam-A-full.seqfrompfamfasta.percentTMpro_scampi.perTM75_nseq20%s.pfamidlist"%(DATADIR3, tag)
    if myfunc.checkfile(tm_pfamidListFile, "%s (line %d): tm_pfamidListFile"%( __file__,
        inspect.currentframe().f_lineno)):
        return 1

    if tm_clanidListFile == "":
        tm_clanidListFile = "%s/data/pfam/pfam26.0/Pfam-A-full.seqfrompfamfasta.percentTMpro_scampi.perTM75_nseq20.clanidlist"%(DATADIR3)
    if myfunc.checkfile(tm_clanidListFile, "%s (line %d): tm_clanidListFile"%( __file__,
        inspect.currentframe().f_lineno)):
        return 1

    if dbname_predTM == "":
        dbname_predTM = "%s/wk/MPTopo/pfamAna_refpro/pred_topcons_single_method4/refpro20120604-celluar.selmaxlength-m1.topcons-single_topcons_single.m1.agree-44.RMSP"%(DATADIR3)
    if myfunc.checkfile("%s0.db"%(dbname_predTM), "%s (line %d): dbname_predTM"%( __file__,
        inspect.currentframe().f_lineno)):
        return 1

    if g_params['isOnlyAnaProkar']:
        prokarseqidfile = "%s/data/uniprot/reference_proteome/refpro20120604-celluar.selmaxlength-m1.nr100.filter.fragmented.Prokaryota.seqidlist"%(DATADIR3)
        g_params['prokarSeqIDSet'] = set(myfunc.ReadIDList(prokarseqidfile))
        if len(g_params['prokarSeqIDSet']) < 1:
            return 1
    if g_params['isOnlyAnaEukar']:
        eukarseqidfile = "%s/data/uniprot/reference_proteome/refpro20120604-celluar.selmaxlength-m1.nr100.filter.fragmented.Eukaryota.seqidlist"%(DATADIR3)
        g_params['eukarSeqIDSet'] = set(myfunc.ReadIDList(eukarseqidfile))
        if len(g_params['eukarSeqIDSet']) < 1:
            return 1

    if pairlistwithpfamidFile == "":
        pairlistwithpfamidFile = "%s/../../Pfam-.maxpair100.pairlistwithpfamid"%(dirpath)
    if myfunc.checkfile(pairlistwithpfamidFile, 
            "%s (line %d): pairlistwithpfamidFile"%( __file__,
                inspect.currentframe().f_lineno)):
        return 1

    pfamid_2_seqidpair_Dict = ReadPairListWithFamID(pairlistwithpfamidFile)
    usedPfamIDSet = set(pfamid_2_seqidpair_Dict.keys()) # pfamids used in pair selection

    if pairListFile != "":
        li = myfunc.ReadPairList(pairListFile)
        SPE_PAIR_LIST = []
        for tup in li:
            SPE_PAIR_LIST.append((int(tup[0]), int(tup[1])))


    (pfamidDefDict, clanidDefDict) = ReadPfamDefFile(pfamDefFile)
    signalpDict = lcmp.ReadSignalPDict(signalpFile)

    seqid2clanidDict = myfunc.ReadFam2SeqidMap(seqid2clanidMapFile)
    seqid2pfamidDict = myfunc.ReadFam2SeqidMap(seqid2pfamidMapFile)

    clanid2seqidDict = myfunc.ReadFam2SeqidMap(clanid2seqidMapFile)
    pfamid2seqidDict = myfunc.ReadFam2SeqidMap(pfamid2seqidMapFile)

    tm_pfamidList = myfunc.ReadIDList(tm_pfamidListFile)
    tm_clanidList = myfunc.ReadIDList(tm_clanidListFile)

    tm_pfamidSet = set(tm_pfamidList)
    tm_clanidSet = set(tm_clanidList)

    hdl_predTM = myfunc.MyDB(dbname_predTM)
    if not hdl_predTM.failure:
        idSet_TMpro = set(hdl_predTM.indexedIDList)
    else:
        idSet_TMpro = set([])

    #classList_TableNumTMHeatMap = ["ALL", "RMSP", "RMDUP"] 
    #alignrangeList = ['FULL_ALIGNED', 'all', 'PART_ALIGNED']
    alignrangeList = ['FULL_ALIGNED']

    if g_params['outpath'] != "" and not os.path.exists(g_params['outpath']):
        cmd = ["mkdir", "-p", g_params['outpath']]
        try:
            subprocess.check_call(cmd)
        except subprocess.CalledProcessError, e:
            print e
            return 1

    print "#Analyzing %s"%(infile)
    for alignrange in alignrangeList:
        Ana_NumTMHeatMap(infile, seqid2pfamidDict, seqid2clanidDict,
                tm_pfamidSet, tm_clanidSet, pfamidDefDict, clanidDefDict,
                signalpDict, classList_TableNumTMHeatMap, SPE_PAIR_LIST, 
                pfamid2seqidDict, clanid2seqidDict, idSet_TMpro, usedPfamIDSet,
                alignrange)

    if not hdl_predTM.failure:
        hdl_predTM.close()

#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['numTMHeatMapMode'] = "half"
    g_params['scale_count'] = 1.0
    g_params['outname'] = "tmp_ana_numTMHeatMap"
    g_params['pairwise_comparison_method'] = 1
    g_params['winsize'] = 50
    g_params['mindiffpair'] = 1
    g_params['outpath'] = ""
    g_params['isOnlyAnaProkar'] = False
    g_params['isOnlyAnaEukar'] = False
    g_params['prokarSeqIDSet'] = set([])
    g_params['eukarSeqIDSet'] = set([])
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

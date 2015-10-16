#!/usr/bin/env python

# Given the data produced by compareMSATopo.py -mode 0
# output the following tables and also draw the figure
# 1. Fraction of pairwise comparison classes -vs- sequence identity
# 2. Normalized frequency of occurrence of unmapped TM regions
#    happened at N, C-terminals and internal regions  -vs-
#    sequence identity
# 3. Distribution of number of TM helices of unmapped TM regions at
#    N, C-terminals and internal regions categorized by sequence
#    identity

# 2014-06-23
#   Add a new statistics, same nTM and same orientation as SA
#                         and others as DI
# ChangeLog 2012-08-13#{{{
#   cmpclass "DIFF" further classified as "DIFF1" and "DIFF2"
#   DIFF1: with the same number of TM helices
#   DIFF2: with different number of TM helices
# ChangeLog 2012-10-11
# ChangeLog 2012-11-19:
#   statistics with rlty bins output
# ChangeLog 2012-11-29 
#   added two options: -print-rlty-cmpclass
#                      -print-rlty-helixcmpclass
#
# ChangeLog 2013-03-15
#   -signalp added
# ChangeLog 2013-03-20
#    -rmsp function removed. signal peptide filtering is done in compareMSATopo
# ChangeLog 2013-05-30
#    1. fix seqidt problem, pairs with 100% sequence identity are filtered away
#    2. option -alignrange is added
#    3. for heatmap, when norm_mode == diag, output only half of the heatmap,
#    since it is symmetric
# ChangeLog 2013-06-06 
#    1. use new scheme for determing of N and C terminal
#    anything before aligned TM region is N-terminal and anything after
#    aligned TM region is C-termianl
#    e.g. for the alignment
#           0   0   0   0  1  1 
#    2 1    0   0   0   0
#    N-Term                C-Term
#    then one can also check those families with most variation at both N and C
#    terminal versas those with most at N or at C
# ChangeLog 2013-09-18
#    For the five topology variation category, INV|noSP, TM2GAP|noSP, Mixed,
#    TM2SEQ|noSP, SP2TM. "Other" is changed to "Mixed". Besides, INV|noSP,
#    TM2GAP|noSP and Mixed are considered as real topology variations. And
#    SP2TM is considered as mainly due to mispredictions while TM2SEQ|noSP can
#    be errors or not errors.
# 
#}}}
import os
import sys
import libtopologycmp as lcmp
from operator import itemgetter
import numpy
from math import ceil
import myfunc
import re
import copy
import subprocess

rundir = os.path.dirname(os.path.realpath(__file__))

DEBUG_UNMAPPED_TM_POSITION = 0
#define TM map value, method 1
_UNALIGNED = -1
_TM2TM = 0
_TM2GAP = 1
_TM2SEQ = 2
_TM2SP = 3

GAP = "-"

DATADIR3 = os.environ['DATADIR3']
binpath = os.path.dirname(sys.argv[0])

SEQIDT_GROUP_ALL = [
        0,100]
SEQIDT_GROUP_1 = [
        0,20,
        20,30,
        30,100
        ]
SEQIDT_GROUP_2 = [
        0 ,10,
        10 ,20,
        20 ,30,
        30 ,40,
        40 ,50,
        50 ,60,
        60 ,70,
        70 ,80,
        80 ,90,
        90 ,100 
        ]
SEQIDT_GROUP_3 = [
        0,25,
        25,100
        ]
SEQIDT_GROUP_4 = [
        0,  10,
        10, 20,
        20, 30,
        30, 80,
        80, 100
    ]
SEQIDT_GROUP_5 = [
        0,  20,
        20, 30,
        30, 40,
        40, 100
    ]

BLOCK_SIZE = 100000
progname =  os.path.basename(sys.argv[0])
usage="""
Usage: %s paircmp-file [-outpath outpath]

Description: Analyze the paircmp file, output three tables
             1. Fraction of pairwise comparison classes -vs- sequence identity
             2. Normalized frequency of occurrence of unmapped TM regions
                happened at N, C-terminals and internal regions  -vs-
                sequence identity
             3. Distribution of number of TM helices of unmapped TM regions at
                N, C-terminals and internal regions categorized by sequence
                identity
OPTIONS:
  -outpath DIR  Output the result to DIR, (default: ./)   
  -diffseqidtgroup 0|1|2|3|9
                Select SeqIDT classifications
                9: for evodist
  -mp INT       pairwise comparion method, 0 or 1, (default: 0)
                3: helix level, TM2TM, TM2GAP, TM2SEQ, TM2SP
                    protein level, IDT, INV, TM2GAP, TM2SEQ, TM2SP, Mixed
  -q            Quiet mode
  -h, --help    Print this help message and exit

  -type STR          selecting type, nterm, cterm, internal, all
                     (default: all)
  -tableinfo FILE    Set pairwise alignment table info, get more pairwise
                     statistics
  -seqidttype INT    Set sequence identity type, (default: 0)
                     0: seqIDT = numIDTRes /alnLength
                     1: seqIDT = numIDTRes / min(len1, len2)
                     2: seqIDT = numIDTRes / (alnLength - NumGAP)
                     Note: if seqidttype = 1 or 2, tableinfo file must be set
  -seqdef    FILE    Set Sequence definition file
  -seq2fammap FILE   Set seqid to famid mapping file
  -printdiffseq      Write sequence pair with different topology to outfile
  -printcountpair    Write count all pairs in every family
  -rltyinfo  FILE    Set reliability information
  -filter-predseq    y | n 
                     Filter predicted sequences
  -dgscore FILE      Set dgscore list file
                     when dgscorelist file is supplied, output also the statistics
                     of dgscores of TM2TM, TM2GAP and TM2SEQ
  -topoaln  FILE     Set aligned topology file, this is used when analyzing the
                     location of unmapped TM helices
  -seqaln   FILE     Supply sequence alignment file
  -debug y|n         Whether show debugging information, (default: no)
  -testlevel INT     Just run test level, debug
  -evodist           Whether it is evolutionary distance instead of seqidt
  -print-rlty-cmpclass
                     Print file with content \"rlty cmpclass\"
  -print-rlty-helixcmpclass
                     Print file with content \"rlty helixcmpclass\"
  -heatmap STR       Set heatmap matrix mode, STR = half or full
                     (default: full)
  -signalp FILE      Supply signal peptide definition file, if set, proteins with 
                     signal peptide will be removed
  -dupfile FILE      Supply duplication definition file
  -rmdup             Remove duplications
  -ncinterdef STR    Set Nterm, Cterm, Internal definition method list, (default: 4) 
                     e.g.  -ncinterdef \"0 3 4\", can be a subset of [0, 1, ,2 ,3 ,4]
  -thrltylist STR    Set reliability threshold list, (default: 0.0)
                     e.g.  -thrltylist \"0.0 80.0\"
  -seqidtmode STR    Set sequence identity subgroup mode list, (default: all)e.g.
                     -seqidtmode \"low high all\"
  -alignrange STR    Select alignment with different alignment ranges
                     all, full, part, (default: all)
  -selidlist  FILE   Analyze only for supplied ids
  -seqid2pfamid FILE

Selection control options:
  -gap  FLOAT  Select only the TM with gap fraction >= threshold, 
               (default: 0.5)
  -dg   FLOAT  Select only the TM with DG values <= threshold, 
               (default: 1.0)
  -ps   FLOAT  Set reliability score 0-100. (default: 0)
               select only the TM with topology prediction reliability >=
               threshold
  -min-seqidt  FLOAT, (default: 0)
  -max-seqidt  FLOAT, (default: 100)
               Select only pairs with global sequence identity within [minSeqIDT, maxSeqIDT]

Created 2011-11-08, updated 2015-03-23, Nanjiang Shu 
"""%(progname)

def PrintHelp():#{{{
    print usage
#}}}
def OpenFileRltyCmpclass(g_params, outpath, rootname):#{{{
    for key in ['min_ps', 'avg_ps', 'max_ps']:
        outfile =  (outpath + os.sep + rootname + '_' 
                    + key +'_rlty_cmpclass.txt')
        try: 
            g_params['fpout_rlty_cmpclass_' + key] = open(outfile, "w")
            print "Write to file %s" %(outfile)
        except IOError:
            print >> sys.stderr, "Failed to write to file"%(outfile)
            g_params['fpout_rlty_cmpclass_' + key] = None
#}}}
def OpenFileRltyHelixCmpclass(g_params, outpath, rootname):#{{{
    for key in ['min_ps', 'avg_ps', 'max_ps']:
        outfile =  (outpath + os.sep + rootname + '_' 
                    + key +'_rlty_helixcmpclass.txt')
        try: 
            g_params['fpout_rlty_helixcmpclass_' + key] = open(outfile, "w")
            print "Write to file %s" %(outfile)
        except IOError:
            print >> sys.stderr, "Failed to write to file"%(outfile)
            g_params['fpout_rlty_helixcmpclass_' + key] = None
#}}}
def CloseFileRltyCmpclass():#{{{
    for key in ['min_ps', 'avg_ps', 'max_ps']:
        try: 
            g_params['fpout_rlty_cmpclass_' + key].close()
        except (IOError, KeyError):
            pass
#}}}
def CloseFileRltyHelixCmpclass():#{{{
    for key in ['min_ps', 'avg_ps', 'max_ps']:
        try: 
            g_params['fpout_rlty_helixcmpclass_' + key].close()
        except (IOError, KeyError):
            pass
#}}}

def GOAnalysis_mp3(dataTableCmpClass, cmpClassList, outpath, rootname, addname,
        gomapfile, gotermfile): #{{{
    """
    Topology variation analysis for proteins with different molecular functions
    """
# write pairinfo list
    return 
    gofreqfileList = []
    pairinfofileList = []
    if "pairinfo" in dataTableCmpClass:
        idListAll = []
        for i in range(len(cmpClassList)):
            outfilepairinfo = (outpath + os.sep + rootname + '_' + addname +
                    ".%s.pairinfo.txt"%cmpClassList[i])
            pairinfofileList.append(outfilepairinfo)
            WritePairInfo(dataTableCmpClass['pairinfo'][i], outfilepairinfo)
# write seqidlist
            outfileseqidlist = (outpath + os.sep + rootname + '_' + addname +
                    ".%s.seqidlist"%(cmpClassList[i]))
            idList = []
            for tup in dataTableCmpClass['pairinfo'][i]:
                idList.append(tup[0])
                idList.append(tup[1])
            if myfunc.WriteIDList(idList, outfileseqidlist) == 0:
# Do GO analysis
                outfile_goanalysis = (outpath + os.sep + rootname + '_' +
                        addname + ".%s.gofreq.txt"%cmpClassList[i])
                gofreqfileList.append(outfile_goanalysis)
                cmd = "%s/anaGOterm_uniprotid.py %s $file"\
                        " -gomap %s -goterm %s  -o %s"%(binpath,
                                outfileseqidlist, gomapfile, gotermfile,
                                outfile_goanalysis)
                os.system(cmd)

            idListAll += idList
# Write all seqidlist
        outfileseqidlist = (outpath + os.sep + rootname + '_' + addname +
                ".%s.seqidlist"%("All"))
        if myfunc.WriteIDList(idListAll, outfileseqidlist) == 0:
# Do GO analysis
            outfile_goanalysis = (outpath + os.sep + rootname + '_' +
                    addname + ".%s.gofreq.txt"%("All"))
            gofreqfileList.append(outfile_goanalysis)
            cmd = "%s/anaGOterm_uniprotid.py %s $file"\
                    " -gomap %s -goterm %s  -o %s"%(binpath,
                            outfileseqidlist, gomapfile, gotermfile,
                            outfile_goanalysis)
            os.system(cmd)
    cmd = "Rscript %s/mergeGOfreqdata.R  %s"%(binpath, " ".join(gofreqfileList))
    os.system(cmd)
    cmd = "%s/tmp_cmpclass_for_GO.py -mp %d %s"%(binpath,
            g_params['pairwise_comparison_method'],
            " ".join(pairinfofileList))

# plot GO analysis figures
    mergedWithConfFileList = []
    mergedWithPvalueFileList = []
    t_bsname =  (outpath + os.sep + rootname + '_' + addname )
    for tt in cmpClassList:
        mergedWithPvalueFileList.append("%s.%s.merged.withpvalue.txt"%(t_bsname,
            tt))
        mergedWithConfFileList.append("%s.%s.merged.withconfint.txt"%(t_bsname,
            tt))
    cmd = "%s/plotGOFreqAna.sh -plotwithpvalue %s"%(binpath, 
            " ".join(mergedWithPvalueFileList))
    os.system(cmd)
    cmd = "%s/plotGOFreqAna.sh -plotwithconf %s"%(binpath, 
            " ".join(mergedWithConfFileList))
    os.system(cmd)
#}}}
def ReadIDMap2(infile):#{{{
    try:
        seq2famDict = {}
        fpin = open(infile,"r")
        lines = fpin.readlines()
        for line in lines:
            if line:
                strs = line.split()
                if len(strs) > 2:
                    seq2famDict[strs[0]] = strs[2:]
                else:
                    print >> sys.stderr, "broken item in file %s: line \"%s\"" \
                            % (infile, line)
        fpin.close()
        return seq2famDict
    except IOError:
        print "Failed to read listfile ", infile
        return {}
#}}}
def IsIdenticalOrIsoformProtein(seqdef1, seqdef2): #{{{
    seqdef1.lstrip("PREDICTED: ")
    seqdef2.lstrip("PREDICTED: ")
    if seqdef1 == seqdef2:
        return True
    else:
        if (seqdef1.find('isoform') != -1 and seqdef2.find('isoform') != -1 and
            seqdef1[0:10] == seqdef2[0:10]):
            return True
        else:
            return False
#}}}
def ReadSeqDefInfo_refseq(infile):#{{{
    try:
        seqInfoDict = {}
        fpin = open(infile,"r")
        line = fpin.readline()
        line = fpin.readline()
        while line:
            strs = line.split('|')
            if len(strs) == 4:
                gid = strs[0].strip()
                refseqid = strs[1].strip()
                pfamid = strs[2].strip()
                seqdef = strs[3].strip()
                seqInfoDict[gid] = {}
                seqInfoDict[gid]['pfamid'] = pfamid 
                seqInfoDict[gid]['refseqid'] = refseqid 
                seqInfoDict[gid]['seqdef'] = seqdef 
            line = fpin.readline()
        fpin.close()
        return seqInfoDict
    except IOError:
        print >> sys.stderr, "Error! file seqDefFile (%s) does not exist." %infile
        return {}
#}}}
def ReadSeqDefInfo_idwithanno(infile): #{{{
    try:
        seqInfoDict = {}
        fpin = open(infile,"r")
        line = fpin.readline()
        while line:
            if line and line[0] != "#":
                strs = line.split("\t")
                if len(strs) >= 2:
                    seqid = strs[0].strip()
                    ss = strs[1].split()
                    seqdef = ""
                    if len(ss) > 1:
                        seqdef = " ".join(ss[1:])
                    if not seqid in seqInfoDict:
                        seqInfoDict[seqid] = {}
                        seqInfoDict[seqid]['seqdef'] = seqdef
            line = fpin.readline()
        fpin.close()
        return seqInfoDict
    except IOError:
        print >> sys.stderr, "Error! file seqDefFile (%s) does not exist." %infile
        return {}

#}}}
def ReadSeqDefInfo_uniref(infile):#{{{
    try:
        seqInfoDict = {}
        fpin = open(infile,"r")
        line = fpin.readline()
        while line:
            if line and line[0] != "#":
                strs = line.split("\t")
                if len(strs) >= 2:
                    seqid = strs[0].strip()
                    seqname = strs[1].strip()
                    if not seqid in seqInfoDict:
                        seqInfoDict[seqid] = {}
                        seqInfoDict[seqid]['seqdef'] = seqname 
            line = fpin.readline()
        fpin.close()
        return seqInfoDict
    except IOError:
        print >> sys.stderr, "Error! file seqDefFile (%s) does not exist." %infile
        return {}
#}}}
def ReadDGScore_old(infile):#{{{
    try:
        dgScoreDict = {}
        fpin = open(infile,"r")
        line = fpin.readline()
        while line:
            if line and line[0] != "#":
                strs = line.split()
                numStr = len(strs)
                if numStr >= 2:
                    seqid = strs[0]
                    if numStr == 2:
                        dgscore = float(strs[1])
                    elif numStr == 3:
                        dgscore = float(strs[2])
                    if not seqid in dgScoreDict:
                        dgScoreDict[seqid] = []
                    dgScoreDict[seqid].append(dgscore)
            line = fpin.readline()
        fpin.close()
        return dgScoreDict
    except IOError:
        print >> sys.stderr, "Error! file seqDefFile (%s) does not exist." %infile
        return {}
#}}}
def ReadTopoAln(infile):#{{{
    try:
        topoalnDict = {}
        (idList, annoList, seqList) = myfunc.ReadFasta(infile)
        numSeq = len(idList)
        if numSeq > 1:
            numPair = numSeq / 2
            for i in xrange(numPair):
                id1 = idList[2*i]
                id2 = idList[2*i+1]
                ss = "%s-%s"%(id1,id2) 
                topoalnDict[ss] = [seqList[2*i], seqList[2*i+1]]
        return topoalnDict
    except IOError:
        print >> sys.stderr, "Error! file seqDefFile (%s) does not exist." %infile
        return {}
#}}}
def GetRlty(record, rltyDict, method):#{{{
    try:
        id1 = record['id1']
        id2 = record['id2']
        # ps is using the minimum of the pair
        ps1 = rltyDict[id1]
        ps2 = rltyDict[id2]

        if method == 'min':
            ps = min(ps1,ps2)
        elif method == 'max':
            ps = max(ps1,ps2)
        elif method == 'avg':
            ps = (ps1+ps2)/2.0
        return ps
    except KeyError:
        return -100.0

#}}}

def AnaUnalignedTerminal(recordList, dataTable):#{{{
    for record in recordList:
        if record == {}:
            continue
        if not record['isLocalAlignment']:
            continue
        if record['alignrange'] == "PART_ALIGNED":
            dataTable['numPair'] += 1
            dataTable['numTM_Nterm1'].append(record['numTM_unaligned_Nterm1'])
            dataTable['numTM_Nterm2'].append(record['numTM_unaligned_Nterm2'])
            dataTable['numTM_Cterm1'].append(record['numTM_unaligned_Cterm1'])
            dataTable['numTM_Cterm2'].append(record['numTM_unaligned_Cterm2'])
            dataTable['diff_numTM_Nterm'].append(abs(record['numTM_unaligned_Nterm1']-record['numTM_unaligned_Nterm2']))
            dataTable['diff_numTM_Cterm'].append(abs(record['numTM_unaligned_Cterm1']-record['numTM_unaligned_Cterm2']))
#}}}

def AnaTM2SEQSegment(recordList, dataTable, topoalnDict, seqalnDict): #{{{
# only for helix level TM2SEQ
    data = dataTable['data']
    data_shift = dataTable['data_shift']
    data_noshift = dataTable['data_noshift']
    seqidttype = g_params['seqidttype']
    for record in recordList:
        if record == {}:
            continue
        cmpclass = record['cmpclass']
        if "mapArray" not in record:
            continue
        seqidt = lcmp.GetSeqIDT(record, seqidttype)
        idList = [record['id1'], record['id2']]
        ss = "-".join(idList)
        try:
            topoaln = topoalnDict[ss]
        except:
            print >> sys.stderr, "topoaln not found for %s"%(ss)
            continue
        try:
            seqaln = seqalnDict[ss]
        except KeyError:
            print >> sys.stderr, "seqaln not found for %s"%(ss)
            continue

        alignedMapArray = GetAlignedMapArray_mp1(record['mapArray'])

        for pivIdx in range(2):
            idx1 = pivIdx
            idx2 = (pivIdx+1)%2
            id1 = idList[idx1]
            id2 = idList[idx2]
            mapArray = record['mapArray'][pivIdx]
            posTM = myfunc.GetTMPosition(topoaln[pivIdx])
            numTM1 = len(record['mapArray'][idx1])
            numTM2 = len(record['mapArray'][idx2])
            for i in xrange(len(alignedMapArray[pivIdx])):
                tup = alignedMapArray[pivIdx][i]
                mp = tup[0]
                iTM = tup[1]
                if mp == 2: #isTM2SEQ
                    segTM = seqaln[idx1][posTM[iTM][0]:posTM[iTM][1]]
                    segNonTM = seqaln[idx2][posTM[iTM][0]:posTM[iTM][1]]
                    segTM_topo = topoaln[idx1][posTM[iTM][0]:posTM[iTM][1]]
                    segNonTM_topo = topoaln[idx2][posTM[iTM][0]:posTM[iTM][1]]
                    alignfactor = lcmp.GetAlignmentFactorFromPairAlignment(segTM, segNonTM, True)
                    seqidt_seg = lcmp.GetSeqIDT(alignfactor, seqidttype)
                    rd = (id1, id2, cmpclass, seqidt, posTM[iTM][0], posTM[iTM][1],
                            segTM, segNonTM, seqidt_seg, segTM_topo, segNonTM_topo)
#                     print rd

                    if (((i-1)>=0 and alignedMapArray[(pivIdx+1)%2][i-1][0] in [2]) 
                        or ((i+1)<len(alignedMapArray[pivIdx]) and
                            alignedMapArray[(pivIdx+1)%2][i+1][0] in [2])
                        ):
                        isShiftedTM = True
                    else:
                        isShiftedTM = False
                    data.append(rd)
                    if not isShiftedTM:
                        data_noshift.append(rd)
                    else:
                        data_shift.append(rd)

#}}}
def AnaTM2GAPSegment(recordList, dataTable, topoalnDict, seqalnDict): #{{{
# only for helix level TM2GAP
    data = dataTable['data']
    seqidttype = g_params['seqidttype']
    for record in recordList:
        if record == {}:
            continue
        cmpclass = record['cmpclass']
        if "mapArray" not in record:
            continue
        seqidt = lcmp.GetSeqIDT(record, seqidttype)
        idList = [record['id1'], record['id2']]
        ss = "-".join(idList)
        try:
            topoaln = topoalnDict[ss]
        except:
            print >> sys.stderr, "topoaln not found for %s"%(ss)
            continue
        try:
            seqaln = seqalnDict[ss]
        except KeyError:
            print >> sys.stderr, "seqaln not found for %s"%(ss)
            continue

        alignedMapArray = GetAlignedMapArray_mp1(record['mapArray'])

        for pivIdx in range(2):
            idx1 = pivIdx
            idx2 = (pivIdx+1)%2
            id1 = idList[idx1]
            id2 = idList[idx2]
            mapArray = record['mapArray'][pivIdx]
            posTM = myfunc.GetTMPosition(topoaln[pivIdx])
            numTM1 = len(record['mapArray'][idx1])
            numTM2 = len(record['mapArray'][idx2])
            for i in xrange(len(alignedMapArray[pivIdx])):
                tup = alignedMapArray[pivIdx][i]
                mp = tup[0]
                iTM = tup[1]
                if mp == 1: #isTM2GAP
                    segTM = seqaln[idx1][posTM[iTM][0]:posTM[iTM][1]]
                    segNonTM = seqaln[idx2][posTM[iTM][0]:posTM[iTM][1]]
                    segTM_topo = topoaln[idx1][posTM[iTM][0]:posTM[iTM][1]]
                    segNonTM_topo = topoaln[idx2][posTM[iTM][0]:posTM[iTM][1]]
                    alignfactor = lcmp.GetAlignmentFactorFromPairAlignment(segTM, segNonTM, True)
                    seqidt_seg = lcmp.GetSeqIDT(alignfactor, seqidttype)
                    rd = (id1, id2, cmpclass, seqidt, posTM[iTM][0], posTM[iTM][1],
                            segTM, segNonTM, seqidt_seg, segTM_topo, segNonTM_topo)
                    data.append(rd)

#}}}
def AnaTM2TMSegment(recordList, dataTable, topoalnDict, seqalnDict): #{{{
# only for helix level TM2TM
    data = dataTable['data']
    seqidttype = g_params['seqidttype']
    for record in recordList:
        if record == {}:
            continue
        cmpclass = record['cmpclass']
        if "mapArray" not in record:
            continue
        seqidt = lcmp.GetSeqIDT(record, seqidttype)
        idList = [record['id1'], record['id2']]
        ss = "-".join(idList)
        try:
            topoaln = topoalnDict[ss]
        except:
            print >> sys.stderr, "topoaln not found for %s"%(ss)
            continue
        try:
            seqaln = seqalnDict[ss]
        except KeyError:
            print >> sys.stderr, "seqaln not found for %s"%(ss)
            continue

        alignedMapArray = GetAlignedMapArray_mp1(record['mapArray'])

        for pivIdx in range(2):
            idx1 = pivIdx
            idx2 = (pivIdx+1)%2
            id1 = idList[idx1]
            id2 = idList[idx2]
            mapArray = record['mapArray'][pivIdx]
            posTM = myfunc.GetTMPosition(topoaln[pivIdx])
            numTM1 = len(record['mapArray'][idx1])
            numTM2 = len(record['mapArray'][idx2])
            for i in xrange(len(alignedMapArray[pivIdx])):
                tup = alignedMapArray[pivIdx][i]
                mp = tup[0]
                iTM = tup[1]
                if mp == 0: #isTM2TM
                    segTM = seqaln[idx1][posTM[iTM][0]:posTM[iTM][1]]
                    segNonTM = seqaln[idx2][posTM[iTM][0]:posTM[iTM][1]]
                    segTM_topo = topoaln[idx1][posTM[iTM][0]:posTM[iTM][1]]
                    segNonTM_topo = topoaln[idx2][posTM[iTM][0]:posTM[iTM][1]]
                    alignfactor = lcmp.GetAlignmentFactorFromPairAlignment(segTM, segNonTM, True)
                    seqidt_seg = lcmp.GetSeqIDT(alignfactor, seqidttype)
                    rd = (id1, id2, cmpclass, seqidt, posTM[iTM][0], posTM[iTM][1],
                            segTM, segNonTM, seqidt_seg, segTM_topo, segNonTM_topo)
                    data.append(rd)

#}}}
def AnaTM2SPSegment(recordList, dataTable, topoalnDict, seqalnDict): #{{{
# only for helix level TM2SP
    data = dataTable['data']
    seqidttype = g_params['seqidttype']
    for record in recordList:
        if record == {}:
            continue
        cmpclass = record['cmpclass']
        if "mapArray" not in record:
            continue
        seqidt = lcmp.GetSeqIDT(record, seqidttype)
        idList = [record['id1'], record['id2']]
        ss = "-".join(idList)
        try:
            topoaln = topoalnDict[ss]
        except:
            print >> sys.stderr, "topoaln not found for %s"%(ss)
            continue
        try:
            seqaln = seqalnDict[ss]
        except KeyError:
            print >> sys.stderr, "seqaln not found for %s"%(ss)
            continue

        alignedMapArray = GetAlignedMapArray_mp1(record['mapArray'])

        for pivIdx in range(2):
            idx1 = pivIdx
            idx2 = (pivIdx+1)%2
            id1 = idList[idx1]
            id2 = idList[idx2]
            mapArray = record['mapArray'][pivIdx]
            posTM = myfunc.GetTMPosition(topoaln[pivIdx])
            numTM1 = len(record['mapArray'][idx1])
            numTM2 = len(record['mapArray'][idx2])
            for i in xrange(len(alignedMapArray[pivIdx])):
                tup = alignedMapArray[pivIdx][i]
                mp = tup[0]
                iTM = tup[1]
                if mp == 3: #isTM2SP
                    segTM = seqaln[idx1][posTM[iTM][0]:posTM[iTM][1]]
                    segNonTM = seqaln[idx2][posTM[iTM][0]:posTM[iTM][1]]
                    segTM_topo = topoaln[idx1][posTM[iTM][0]:posTM[iTM][1]]
                    segNonTM_topo = topoaln[idx2][posTM[iTM][0]:posTM[iTM][1]]
                    alignfactor = lcmp.GetAlignmentFactorFromPairAlignment(segTM, segNonTM, True)
                    seqidt_seg = lcmp.GetSeqIDT(alignfactor, seqidttype)
                    rd = (id1, id2, cmpclass, seqidt, posTM[iTM][0], posTM[iTM][1],
                            segTM, segNonTM, seqidt_seg, segTM_topo, segNonTM_topo)
                    data.append(rd)

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
def AddSeqDefInfo(recordList, seqInfoDict):#{{{
    if seqInfoDict != {}:
        for record in recordList:
            id1 =  record['id1']
            id2 =  record['id2']
            if id1 in seqInfoDict:
                record['seqdef1'] = seqInfoDict[id1]['seqdef']
            else:
                record['seqdef1'] = ""
            if id2 in seqInfoDict:
                record['seqdef2'] = seqInfoDict[id2]['seqdef']
            else:
                record['seqdef2'] = ""
#}}}
def AddIDMap2Info(recordList, seq2famDict):#{{{
    if seq2famDict != {}:
        for record in recordList:
            id1 =  record['id1']
            id2 =  record['id2']
            record['pfamid2'] = []
            if id1 in seq2famDict:
                record['pfamid1'] = seq2famDict[id1]
            else:
                record['pfamid1'] = []
            if id2 in seq2famDict:
                record['pfamid2'] = seq2famDict[id2]
            else:
                record['pfamid2'] = []
            record['pfamid-inter'] = list(set(record['pfamid1']) &
                    set(record['pfamid2']))
            record['pfamid-union'] = list(set(record['pfamid1']) |
                    set(record['pfamid2']))
#}}}
def AddDGScore(recordList, dgScoreDict):#{{{
    if dgScoreDict != {}:
        for record in recordList:
            id1 =  record['id1']
            id2 =  record['id2']
            record['dgscore'] = []
            try:
                record['dgscore'].append(dgScoreDict[id1])
            except KeyError:
                print >> sys.stderr, "no dgscore for id %s"%id1
                record['dgscore'].append([])
            try:
                record['dgscore'].append(dgScoreDict[id2])
            except KeyError:
                record['dgscore'].append([])
                print >> sys.stderr, "no dgscore for id %s"%id2
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

def IsAnaHasCtermVariation(ana):#{{{
    if ana == {}:
        return False
    if 'Cterm' in ana and ana['Cterm'] != {} > 0:
        return True
    else:
        return False
        #}}}
def IsHasCtermVariation(record):#{{{
    if record == {}:
        return False
    if record['cmpclass'] != 'DIFF':
        return False
    if (IsAnaHasCtermVariation(record['ana1']) or
        IsAnaHasCtermVariation(record['ana2'])):
        return True
    else:
        return False
#}}}
def IsAnaHasNtermVariation(ana):#{{{
    if ana == {}:
        return False
    if 'Nterm' in ana and ana['Nterm'] != {} > 0:
        return True
    else:
        return False
        #}}}
def IsHasNtermVariation(record):#{{{
    if record == {}:
        return False
    if record['cmpclass'] != 'DIFF':
        return False
    if (IsAnaHasNtermVariation(record['ana1']) or
        IsAnaHasNtermVariation(record['ana2'])):
        return True
    else:
        return False
#}}}

def IsAnaHasInternalVariation(ana):#{{{
    if ana == {}:
        return False
    if 'internal' in ana and len(ana['internal']) > 0:
        return True
    else:
        return False
        #}}}
def IsHasInternalVariation(record):#{{{
    if record == {}:
        return False
    if record['cmpclass'] != 'DIFF':
        return False
    if (IsAnaHasInternalVariation(record['ana1']) or
        IsAnaHasInternalVariation(record['ana2'])):
        return True
    else:
        return False
#}}}

def IsIndelTMDuplicated(posTM, idxTMList, duphit, idx_in_duppair):#{{{
    """
    Check whether the indels of TM helices is a result of duplication
    posTM:      a list of positions (begin,end) of TM helices
    idxTMList:  a list of indeces to TM helices
    duphit   :  hit of duplicated domains
    idx_in_duppair:  index of the sequence in the duplication pair, 0 or 1
    """
    li = []
    for idx in idxTMList:
        li.append(posTM[idx][0])
        li.append(posTM[idx][1])
    posTMregion = (min(li), max(li))
    (b1, e1) = posTMregion
    for i in xrange(len(duphit)):
        tup = duphit[i][idx_in_duppair]
        posDupregion = (tup[0], tup[1])
        (b2, e2) = posDupregion
        overlap = max(0, myfunc.coverage(b1, e1, b2, e2))
        if (myfunc.FloatDivision(overlap, e1-b1) > 0.5 or
                myfunc.FloatDivision(overlap, e2-b2) > 0.5):
            return True
    return False
#}}}

def FilterPairCmpResult(recordList, rltyDict, selectIDListSet):#{{{
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
#        isSatisfied1 = False #debug

        if record == {}:
            continue

#         if record['alignrange'] == g_params['alignrange']:
#             isSatisfied1 = True

        key = "%s-%s"%(record['id1'], record['id2'])

        cmpclass = record['cmpclass']
        if (cmpclass.find('UNALIGNED') == 0 
                or cmpclass.find('AMBIGUOUS') == 0):
            #if isSatisfied1: print key, "cmpclass=", cmpclass, "UNALIGNED or AMBIGUOUS" #debug
            continue

        if not (len(selectIDListSet) <= 0 or key in selectIDListSet):
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
            # controlling selection of pairs with FULL_ALIGNED or all
            if record['alignrange'] != g_params['alignrange']: 
                continue

        id1 = record['id1']
        id2 = record['id2']

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
        if (seqidt < g_params['minSeqIDT'] or seqidt >= g_params['maxSeqIDT']):
            continue

        if record['cmpclass'] != 'DIFF':
            newList.append(record)
        else:
            if g_params['selecttype'] == 'internal':
                if not IsHasInternalVariation(record):
                    continue
            elif g_params['selecttype'] == 'nterm':
                if not IsHasNtermVariation(record):
                    continue
            elif g_params['selecttype'] == 'cterm':
                if not IsHasCtermVariation(record):
                    continue
#             if 'isFilterPredictedSeq' in g_params and g_params['isFilterPredictedSeq']:
#                 if (('seqdef1' in record and record['seqdef1'].find('PREDICTED') >= 0)
#                         or 'seqdef2' in record and
#                         record['seqdef2'].find('PREDICTED') >= 0):
#                     continue

#             if 'seqdef1' in record and 'seqdef2' in record:
#                 if IsIdenticalOrIsoformProtein(record['seqdef1'], record['seqdef2']):
#                     continue

            id1 = record['id1']
            id2 = record['id2']
            if ((id1+'-'+id2 in pairListSet) or (id2+'-'+id1 in pairListSet)):
                continue

            newRecord = {}
            newana1  = lcmp.SelectAnaDIFF(record['ana1'], g_params)
            newana2  = lcmp.SelectAnaDIFF(record['ana2'], g_params)
            if newana1 != {} or newana2 != {}:
                lcmp.CopyGeneralInfo_pairwise(newRecord, record)
                newRecord['ana1'] = newana1
                newRecord['ana2'] = newana2
                newList.append(newRecord)
                pairListSet.add(id1+'-'+id2)
            elif g_params['isDEBUG']:
                print "pair %s - %s dropped by SelectAnaDIFF"%(id1,id2)

    numOutputRecord = len(newList)
    if g_params['isDEBUG']:
        if numOutputRecord < numInputRecord:
            print "%d pairs dropped" % (numInputRecord-numOutputRecord)

    return newList
#}}}
def GetIndexOfBins(x, binList, numBin):#{{{
    for i in xrange(numBin):
        if x >= binList[i*2] and x < binList[i*2+1]:
            return i
    return numBin
#}}}
def GetSeqIDTGroupIndex(seqidt, seqIDTGroupList):#{{{
    numGroup = len(seqIDTGroupList)/2
    for i in xrange(numGroup):
        if seqidt >= seqIDTGroupList[i*2] and seqidt < seqIDTGroupList[i*2+1]:
            return i
    return numGroup
#}}}
def GetClassIndex(cls, classList):#{{{
    try:
        return classList.index(cls)
    except ValueError:
        return len(classList)
#}}}
def GetUnmappedTMPositionIndex(mp, cmpclass):#{{{
    """
    mp      : state of TM helix mapping
    cmpclass: state of protein level topology variation classification
    """
# Groups for pairwise_comparison_method 1:
#   0. TM2GAP in All
#   1. TM2SEQ in All
#   2. TM2GAP in TM2GAP
#   3. TM2GAP in TM2GAP_and_TM2SEQ
#   4. TM2SEQ in TM2SEQ
#   5. TM2SEQ in TM2GAP_and_TM2SEQ
# Groups for pairwise_comparison_method 3:
#   0. TM2GAP in all
#   1. TM2SEQ in all
#   2. TM2GAP in TM2GAP
#   3. TM2GAP in Mixed
#   4. TM2GAP in DUP
#   5. TM2SEQ in TM2SEQ
#   6. TM2SEQ in Mixed
    if g_params['pairwise_comparison_method'] == 1:
        if mp == _TM2GAP:
            if cmpclass == "TM2GAP":
                return 2
            elif cmpclass == "TM2GAP_AND_TM2SEQ":
                return 3
            else:
                return -1
        elif mp == _TM2SEQ:
            if cmpclass == "TM2SEQ":
                return 4
            elif cmpclass == "TM2GAP_AND_TM2SEQ":
                return 5
            else:
                return -1
        else:
            return -1
    elif g_params['pairwise_comparison_method'] == 3:
# TM2GAP in all      0
# TM2SEQ in all      1
# TM2GAP in TM2GAP   2
# TM2GAP in Mixed    3
# TM2GAP in DUP      4
# TM2SEQ in TM2SEQ   5
# TM2SEQ in Mixed    6
        if mp == _TM2GAP:
            if cmpclass == "TM2GAP":
                return 2
            elif cmpclass == "Mixed":
                return 3
            elif cmpclass == "DUP":
                return 4
            else:
                return -1
        elif mp == _TM2SEQ:
            if cmpclass == "TM2SEQ":
                return 5
            elif cmpclass == "Mixed":
                return 6
            else:
                return -1
        else:
            return -1
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

def InitTableHelixCmpClass(dataTable, numGroup, numClass):#{{{
    dataTable['freq'] = []
    dataTable['subsum'] = []
    for i in xrange(numGroup+1):
        dataTable['freq'].append([0]*numClass)
        dataTable['subsum'].append(0)
#}}}
def InitTableHelixDGScore(dataTable, numGroup, numClass):#{{{
    dataTable['data'] = []
    dataTable['data_noshift'] = [] # shifted TM2SEQ is not included
    dataTable['data_shift'] = [] # shifted TM2SEQ is not included
    for i in xrange(numClass):
        dataTable['data'].append([])
        dataTable['data_noshift'].append([])
        dataTable['data_shift'].append([])
        for j in xrange(numGroup):
            dataTable['data'][i].append([])
            dataTable['data_noshift'][i].append([])
            dataTable['data_shift'][i].append([])
#}}}
def InitTableNumTMHeatMap(dataTable, classList, MAX_NUMTM):#{{{
    for cls in classList:
        dataTable[cls] = {}
        dt = dataTable[cls]
        dt['data'] = []
        dt['maxNumTM'] = 0
        dt['numPair'] = 0
        for i in xrange(MAX_NUMTM):
            dt['data'].append([0]*MAX_NUMTM)
#}}}
def InitTableUnmappedTMPosition(dataTable, itemList):#{{{
# Groups for pairwise_comparison_method 1:
#   0. TM2GAP in Only TM2GAP
#   1. TM2GAP in Both TM2GAP and TM2SEQ
#   2. TM2SEQ in Only TM2SEQ
#   3. TM2SEQ in Both TM2GAP and TM2SEQ
#   4. TM2GAP in All
#   5. TM2SEQ in All
# Groups for pairwise_comparison_method 3:
#   0. TM2GAP in all
#   1. TM2SEQ in all
#   2. TM2GAP in TM2GAP
#   3. TM2GAP in Mixed
#   4. TM2GAP in DUP
#   5. TM2SEQ in TM2SEQ
#   6. TM2SEQ in Mixed
# item:
#   different methods for determining N- & C-terminal and internal regions
    if g_params['pairwise_comparison_method'] == 1:
        numGroup = 6
    elif g_params['pairwise_comparison_method'] == 3:
        numGroup = 7

    for item in itemList:
        dataTable[item] = {}
        dt = dataTable[item]
        dt['dataPosition'] = []
        dt['dataNumContinuousTM'] = []
        dt['dataTransition'] = []

        # difference in the number of N-terminal helices, data stored as a list
        # of tuples [(diff_numTM, 1), (diff_numTM, -1)]
# dataDiffNumTMNterm_NtermStatus: diff_numTM_at_Nterminal vs Nterminal status
        dt['dataDiffNumTMNterm_NtermStatus'] = []
        dt['dataDiffNumTMNterm_CtermStatus'] = []
        dt['dataDiffNumTMCterm_NtermStatus'] = []
        dt['dataDiffNumTMCterm_CtermStatus'] = []
        for i in xrange(numGroup):
            dt['dataPosition'].append([])
            dt['dataTransition'].append([])
            dt['dataNumContinuousTM'].append([])
#}}}
def InitTableTM2SEQ(dataTable):#{{{
    dataTable['data'] = []
    dataTable['data_shift'] = []
    dataTable['data_noshift'] = []
#}}}
def InitTableTM2TM(dataTable):#{{{
    dataTable['data'] = []
#}}}
def InitTableTM2GAP(dataTable):#{{{
    dataTable['data'] = []
#}}}
def InitTableTM2SP(dataTable):#{{{
    dataTable['data'] = []
#}}}
def InitTableTM2GAP_add_terminal_numTM(dataTable, numGroup, numClass):#{{{
    dataTable['freq'] = []
    dataTable['subsum'] = []
    for i in xrange(numGroup):
        dataTable['freq'].append([0]*numClass)
        dataTable['subsum'].append(0)
#}}}
def InitTableUnaligned(dataTable):#{{{
    dataTable['numPair'] = 0
    dataTable['numTM_Nterm1'] = []
    dataTable['numTM_Nterm2'] = []
    dataTable['numTM_Cterm1'] = []
    dataTable['numTM_Cterm2'] = []
    dataTable['diff_numTM_Nterm'] = []
    dataTable['diff_numTM_Cterm'] = []
    dataTable['diff_numTM_Nterm_hist'] = {}
    dataTable['diff_numTM_Cterm_hist'] = {}
#}}}
def InitTableCmpClass(dataTable, numGroup, numClass):#{{{
# classList  = cmpClassList ["OK","SHIFT","INV","INV_SHIFT","DUP", "SIGNALP", "DIFF"]
# group is seqIDTgroup
#freq 2d array
#subsum 1d array
#      class1 class2
# g1    cnt   cnt 
# g2    cnt   cnt  
# g3    cnt   cnt  
#subsum sum1 sum2
    dataTable['freq'] = []
    dataTable['subsum'] = []
    dataTable['pairinfo'] = []
    for i in xrange(numGroup):
        dataTable['freq'].append([0]*numClass)
        dataTable['subsum'].append(0)
    for i in xrange(numClass):
        dataTable['pairinfo'].append([])
#}}}
def InitTableCmpClass_ps_bin(dataTable, numGroup, numClass):#{{{
    dataTable['min_ps'] = {}
    dataTable['avg_ps'] = {}
    dataTable['max_ps'] = {}
    numBinsSeqIDT = len(SEQIDT_GROUP_4)/2
    for key in ['min_ps', 'avg_ps', 'max_ps']:
        dt = dataTable[key]
        dt['freq'] = []
        dt['freq_by_seqidt'] = []
        dt['subsum'] = []
        for i in xrange(numGroup):
            dt['freq'].append([0]*numClass)
            dt['freq_by_seqidt'].append([0]*numBinsSeqIDT)
            dt['subsum'].append(0)
#}}}
def InitTableNCTermInter(dataTable, numGroup, numClass):#{{{
# classList  = diffClassList=["Nterm", "Cterm", "Internal"]
# numGroup = 3; groupList = [0-20, 20-30, 30-100]
#freq 2d array
#subsum 1d array
#      class1 class2
# g1    cnt   cnt 
# g2    cnt   cnt  
# g3    cnt   cnt  
#subsum sum1 sum2
#
# freq1 is counting for all occurrences of unmapped TMs at N, C-terminals and
# internal regions.
#
# freq2 is counting for pairs with unmapped TMs at N,C-terminals and internal
# regions. In this case,  for each pair, occurrence of unmapped TMs at internal
# regions is counted once even there are several.
    dataTable['freq1'] = []
    dataTable['freq2'] = []
    dataTable['subsum1'] = []
    dataTable['subsum2'] = []
    for i in xrange(numGroup):
        dataTable['freq1'].append([0]*numClass)
        dataTable['freq2'].append([0]*numClass)
        dataTable['subsum1'].append(0)
        dataTable['subsum2'].append(0)
#}}}
def InitTableNumTMDistribution(dataTable, numGroup, numClass, MAX_NUMTM):#{{{
# classList  = diffClassList=["Nterm", "Cterm", "Internal"]
# groupList = seqIDTGroupList
# dataTable is a 3D array
    dataTable['freq'] = []
    dataTable['subsum'] = []
    for i in xrange(numGroup):
        dataTable['freq'].append([])
        dataTable['subsum'].append([])
        dataTable['subsum'][i] = [0] * numClass
        for j in xrange(numClass):
            dataTable['freq'][i].append([])
            dataTable['freq'][i][j] = [0] * MAX_NUMTM
#}}}

def FilterSegPos(posList, string, neighbour_char):#{{{
### return only list of "0110"
    newList = []
    N = len(string)
    for (b,e) in posList:
        if b>0 and string[b-1] != neighbour_char:
            continue
        if e < N-1 and string[e] != neighbour_char:
            continue
        newList.append((b,e))
    return newList#}}}
def GetAlignedRegion_pairwise(pairaln):#{{{
    """
    Get aligned region after removal of terminal gaps
    """
    try:
        lengthAln = len(pairaln[0])
        b = [0]*2
        e = [lengthAln]*2
        for i in xrange(2):
            while pairaln[i][b[i]] == "-":
                b[i] += 1
            while pairaln[i][e[i]-1] == "-":
                e[i] -= 1
        return (max(b),min(e))
    except IndexError:
        return (-1, -1)
#}}}
def GetInterRegion_pairwise(pairaln, method):#{{{
    """
    Get internal region of the pairwise topology alignment
    method: (method 0 and 2 are for gapless seq)
    1. just exclude terminal gaps
    3. (max(firstTMEnd1, firstTMEnd2), min(lastTMBeg1, lastTMBeg2))
    4. after removel of terminal gaps, 
        (min(firstTMEnd1, firstTMEnd2), max(lastTMBeg1, lastTMBeg2))
    """
    try:
        if method in [1,4]:
            lengthAln = len(pairaln[0])
            bAlnList = [0]*2
            eAlnList = [lengthAln]*2
            for i in xrange(2):
                while pairaln[i][bAlnList[i]] == "-":
                    bAlnList[i] += 1
                while pairaln[i][eAlnList[i]-1] == "-":
                    eAlnList[i] -= 1
            bAln = max(bAlnList)
            eAln = min(eAlnList)
            if method == 1:
                return (bAln, eAln)
            elif method == 4:
                posTMList = []
                posTMList.append(myfunc.GetTMPosition(pairaln[0][bAln:eAln]))
                posTMList.append(myfunc.GetTMPosition(pairaln[1][bAln:eAln]))
                endFirstTM = min([posTMList[i][0][1] for i in range(2)])+bAln
                beginLastTM = max([posTMList[i][len(posTMList[i])-1][0] for i
                    in range(2)])+bAln
                return (endFirstTM, beginLastTM)
        elif method == 3:
            posTMList = []
            for i in xrange(2):
                posTMList.append(myfunc.GetTMPosition(pairaln[i]))
            endFirstTM = max([posTMList[i][0][1] for i in range(2)])
            beginLastTM = min([posTMList[i][len(posTMList[i])-1][0] for i
                in range(2)])
            return (endFirstTM, beginLastTM)
    except IndexError:
        return (-1, -1)
#}}}
def GetFirstLastTM2TM(mapArray):#{{{
    n = len(mapArray)
    try:
        b = mapArray.index(0)
    except ValueError:
        return (-1,-1)
    i = n-1
    while i >= 0:
        if mapArray[i] == 0:
            e = i
            break
        else:
            i -= 1
    return (b,e)
#}}}

def GetAlignedMapArray_mp1(mapArrayList):
# alignedMapArray is a list of tuples
# each tuple contains two values
# (helixcmpclass, idxTM)
# helixcmpclass: -1, 0, 1, 2 and None for gap or non TM region
# idxTM idx of the original mapArray

    mapArray1 = mapArrayList[0]
    mapArray2 = mapArrayList[1]
    numTM1 = len(mapArrayList[0])
    numTM2 = len(mapArrayList[1])
    alignedMapArray1 = []
    alignedMapArray2 = []


    cnt1=0
    cnt2=0
    isHead1 = True
    isHead2 = True
    while cnt1 < numTM1 and cnt2 < numTM2:
        if mapArray1[cnt1] == _UNALIGNED or mapArray2[cnt2] == _UNALIGNED:
            if  mapArray2[cnt2] != _UNALIGNED:
                if isHead2:
                    alignedMapArray1.append((mapArray1[cnt1], cnt1))
                    alignedMapArray2.append((None, None))
                    cnt1 += 1
                else:
                    alignedMapArray2.append((mapArray2[cnt2], cnt2))
                    alignedMapArray1.append((None, None))
                    cnt2 += 1
            elif mapArray1[cnt1] != _UNALIGNED:
                if isHead1:
                    alignedMapArray2.append((mapArray2[cnt2], cnt2))
                    alignedMapArray1.append((None, None))
                    cnt2 += 1
                else:
                    alignedMapArray1.append((mapArray1[cnt1],cnt1))
                    alignedMapArray2.append((None,None))
                    cnt1 += 1
            else:
                alignedMapArray1.append((mapArray1[cnt1], cnt1))
                alignedMapArray2.append((mapArray2[cnt2], cnt2))
                cnt1 += 1
                cnt2 += 1
        else: # none of them are unaligned
            if mapArray1[cnt1] != _TM2TM:
                alignedMapArray1.append((mapArray1[cnt1], cnt1))
                alignedMapArray2.append((None, None))
                cnt1 += 1
            elif mapArray2[cnt2] != _TM2TM:
                alignedMapArray2.append((mapArray2[cnt2], cnt2))
                alignedMapArray1.append((None, None))
                cnt2 += 1
            else:
                alignedMapArray1.append((mapArray1[cnt1], cnt1))
                alignedMapArray2.append((mapArray2[cnt2], cnt2))
                cnt1 += 1
                cnt2 += 1
            isHead1 = False
            isHead2 = False

    #tail
    while cnt1 < numTM1:
        alignedMapArray1.append((mapArray1[cnt1],cnt1))
        alignedMapArray2.append((None, None))
        cnt1 +=1
    while cnt2 < numTM2:
        alignedMapArray1.append((None, None))
        alignedMapArray2.append((mapArray2[cnt2], cnt2))
        cnt2 +=1
    return [alignedMapArray1, alignedMapArray2]
#}}}

def AnaPairCmpResultCmpClass(recordList, dataTable, classList, #{{{
        seqIDTGroupList, topoalnDict, 
        isCmpSP, isCmpDup, isRmUnalignedSP, isAna5,isAnaSADI):
    # isCmpSP: whether the comparison of signal peptide is also used
    # isCmpDup: whether duplications are also classified
    # isRmUnalignedSP: whether do not included those with unaligned SP
    # isAna5: whether using the new 5 category for different topology
    # isAnaSADI: whether for cmpdup3, analysis also SAME (same numTM and same orientation) and DIFF (others)
    freq = dataTable['freq']
    subsum = dataTable['subsum']
    numGroup = len(seqIDTGroupList)/2
    numClass = len(classList)
    seqidttype = g_params['seqidttype']
    for record in recordList:
        if record == {}:
            continue

        seqidt = lcmp.GetSeqIDT(record, seqidttype)
        idxGroup = GetSeqIDTGroupIndex(seqidt, seqIDTGroupList)
        cmpclass = record['cmpclass']


        if isRmUnalignedSP and cmpclass.find("unalignedSP") != -1:
            continue

        if isCmpDup:
            if cmpclass.find("TM2GAP|DUP") == 0:
                cmpclass = "DUP"
            else:
                cmpclass = cmpclass.split('|')[0]




        elif isAna5:
            if cmpclass.find("IDT") == 0:
                cmpclass = "IDT"
            elif cmpclass.find("INV|noSP") == 0:
                cmpclass = "INV"
            elif cmpclass.find("TM2GAP|noSP") == 0:
                cmpclass = "TM2GAP"
            elif cmpclass.find("TM2SEQ|noSP") == 0:
                cmpclass = "TM2SEQ"
            elif cmpclass.find("TM2SEQ|SP2TM") != -1 or (
                    cmpclass.find("TM2GAP|SP2TM") != -1 
                        and abs(record['numTM1'] - record['numTM2']) == 1):
                cmpclass = "SP2TM"
            else:
                cmpclass = "Mixed"
        elif (not isCmpSP) or (not isCmpDup):
            cmpclass = cmpclass.split('|')[0]

        id1 = record['id1']
        id2 = record['id2']

        key = "%s-%s"%(id1, id2)
        try:
            pairaln = topoalnDict[key]
        except KeyError:
            print >> sys.stderr, "topology alignment does not find for %s"%(key)
            continue
        NtermStatus1 = lcmp.GetNtermState(pairaln[0])
        NtermStatus2 = lcmp.GetNtermState(pairaln[1])
        seqlen1 = record['seqLength1']
        seqlen2 = record['seqLength2']
        numTM1 = record['numTM1']
        numTM2 = record['numTM2']


        if isAnaSADI: # added 2014-06-23
            if cmpclass in ["TM2GAP", "TM2SEQ", "TM2SP", "Mixed"]:
                if numTM1 == numTM2 and NtermStatus1 == NtermStatus2:
                    cmpclass = "%s|%s"%(cmpclass, "SAME")
                else:
                    cmpclass = "%s|%s"%(cmpclass, "DIFF")


        if record['cmpclass'] == "DIFF":
            if record['numTM1']  == record['numTM2']:
                cmpclass += "1"
            else:
                cmpclass += "2"
        idxClass = GetClassIndex(cmpclass, classList)
        if idxClass == -1:
            continue

        if idxGroup != numGroup and idxClass != numClass:

            freq[idxGroup][idxClass] += 1
            subsum[idxGroup] += 1
            #add record['pfamid-inter'] to pairinfo 2014-09-23, pfamid-inter
            # is the intersection of pfamidlist_of_seqid1 and
            # pfamidlist_of_seqid2

            pfamidlist_itersection = []
            try:
                pfamidlist_itersection = record['pfamid-inter']
            except:
                pfamidlist_itersection = []

            dataTable['pairinfo'][idxClass].append((id1, id2, NtermStatus1,
                NtermStatus2, numTM1, numTM2,
                seqlen1, seqlen2, seqidt, record['pfamid-inter'])) 
            if g_params['isPrintCountPairInFam']:
                for pfamid in record['pfamid-inter']:
                    if pfamid in g_params['countPairInFam'][idxGroup]:
                        g_params['countPairInFam'][idxGroup][pfamid] += 1
                    else:
                        g_params['countPairInFam'][idxGroup][pfamid] = 1

            if g_params['isPrintDIFFPair']:
                if record['cmpclass'] == "DIFF" and seqidt >= 50 :
                    td = (record['seqidt1'], record['id1'], record['id2'],
                            record['pfamid1'], record['pfamid2'],
                            record['seqdef1'],record['seqdef2'])
#                 print record['id1'], record['seqidt1'], record['seqdef1']
#                 print record['id2'], record['seqidt2'], record['seqdef2']
                    g_params['DIFFPairList'].append(td)
#}}}
def AnaPiarCmpResultCmpClass_ps_bin(recordList, dataTable, #{{{
        classList, psAllGroup, rltyDict):
    # ps (rlty) of the pair is calculated in three ways, avg, max, min
    numGroup = len(psAllGroup)/2
    numClass = len(classList)
    seqidttype = g_params['seqidttype']
    binsListSeqIDT = SEQIDT_GROUP_4
    numBinsSeqIDT = len(binsListSeqIDT)/2

    for record in recordList:
        if record == {}:
            continue
        id1 = record['id1']
        id2 = record['id2']
        # ps is using the minimum of the pair
        try:
            ps1 = rltyDict[id1]
        except KeyError:
            ps1 = -100.0
        try:
            ps2 = rltyDict[id2]
        except KeyError:
            ps2 = -100.0

        cmpclass = record['cmpclass']
        if cmpclass.find('|') != -1:
            cmpclass = cmpclass.split('|')[0]
        seqidt = lcmp.GetSeqIDT(record, seqidttype)
        idxClass = GetClassIndex(cmpclass, classList)
        if idxClass == -1:
            continue
        idxSeqidt = GetIndexOfBins(seqidt, binsListSeqIDT, numBinsSeqIDT)
        for key in ['min_ps', 'avg_ps', 'max_ps']:
            if key == 'min_ps':
                ps = min(ps1,ps2)
            elif key == 'max_ps':
                ps = max(ps1,ps2)
            elif key == 'avg_ps':
                ps = (ps1+ps2)/2
            idxGroup = GetIndexOfBins(ps , psAllGroup, numGroup)
            if g_params['isPrintFileRltyCmpclass']:
                try:
                    g_params['fpout_rlty_cmpclass_'+ key].write("%.1f %s\n" %(ps, cmpclass))
                except IOError:
                    pass
            if idxGroup != numGroup and idxClass != numClass:
                try:
                    dataTable[key]['freq'][idxGroup][idxClass] += 1
                    dataTable[key]['freq_by_seqidt'][idxGroup][idxSeqidt] += 1
                    dataTable[key]['subsum'][idxGroup] += 1
                except IndexError:
                    msg =  "Error! %s - %s idxGroup=%d  ps=%g  idxSeqidt=%d seqidt=%g"\
                            " idxClass=%d SIZE1=%d SIZE2=%d"
                    print >> sys.stderr, msg%(id1, id2, idxGroup, ps,
                            idxSeqidt, seqidt, idxClass, len(dataTable[key]['freq_by_seqidt']),
                            len(dataTable[key]['freq_by_seqidt'][idxGroup]))
                    raise
#}}}
def AnaHelixPairCmpResultCmpClass(recordList, dataTable, classList, seqIDTGroupList):#{{{
# ChangeLog 2015-03-23
# analyze also the DG scores of aligned TM helix region for x-y plotting
# e.g.
#Q7YJX5    0 miiafqlavfaliatssillisvpvvfaspdgwsnnknvVFSGTSLWIGLVFLVAILNSLis   61
#            | ||||||||||||||||||||||||||| ||||.|||||||||||||||||||||||||||
#Q5IHA9    0 mtiafqLAVFALIATSSILLISVPVVFassdgwssnknvVFSGTSLWIGLVFLVAILNSLis   61
# Get DG score of 
# lavfaliatssillisvpvv
# LAVFALIATSSILLISVPVV
# 
    freq = dataTable['freq']
    subsum = dataTable['subsum']
    numGroup = len(seqIDTGroupList)/2
    numClass = len(classList)
    seqidttype = g_params['seqidttype']
    for record in recordList:
        if record == {}:
            continue

        seqidt = lcmp.GetSeqIDT(record, seqidttype)
        idxGroup = GetSeqIDTGroupIndex(seqidt, seqIDTGroupList)
        if "mapArray" in record:
            for li in record["mapArray"]:
                for mp in li:
                    idxClass = GetClassIndex(mp, classList)
                    if idxClass == -1:
                        continue
                    if idxGroup != numGroup and idxClass != numClass:
                        freq[idxGroup][idxClass] += 1
                        subsum[idxGroup] += 1
#}}}
def AnaPairCmpTM2GAP_add_term(recordList, dataTable, seqIDTGroupList):#{{{
    freq = dataTable['freq']
    subsum = dataTable['subsum']
    numGroup = len(seqIDTGroupList)/2
    seqidttype = g_params['seqidttype']
    for record in recordList:
        if record == {}:
            continue
        id1 = record['id1']
        id2 = record['id2']
        cmpclass = record['cmpclass']

        if record['alignrange'] != g_params['alignrange']:
            continue
        if cmpclass.find("TM2GAP|nonDUP") == -1:
            continue

        seqidt = lcmp.GetSeqIDT(record, seqidttype)
        idxGroup = GetSeqIDTGroupIndex(seqidt, seqIDTGroupList)
        if "mapArray" in record:
# first get alignedMapArray. e.g.
# TMMap B3S493, numTM1 =  7:  1  1  1  0  0  0  0
# TMMap A8DVP9, numTM2 =  4:           0  0  0  0
# agree_str = 1110000
# then get the segment of 1, which is for indels
            alignedMapArray = GetAlignedMapArray_mp1(record['mapArray'])
            agree_list = []
            lengthAlignedMapArray = len(alignedMapArray[0])
            for j in xrange(len(alignedMapArray[0])):
                if alignedMapArray[0][j][0] == alignedMapArray[1][j][0]:
                    agree_list.append("0")
                else:
                    agree_list.append("1")
            agree_str = "".join(agree_list)
            segList = myfunc.GetSegPos(agree_str, "1")
            numTMSegList = [x[1]-x[0] for x in segList]
            cls = 3
# numClass = 4 (diff_by_term_1, diff_by_term_larger_than_1, diff_by_inter, other)
#                  0                  1                          2          3
            if len(segList) == 1:
                if numTMSegList[0] == 1:
                    if segList[0][0] == 0 or segList[0][1] == lengthAlignedMapArray:
                        cls = 0
                    else:
                        cls = 2
                else: #numTM>=2
                    if segList[0][0] == 0 or segList[0][1] == lengthAlignedMapArray:
                        cls = 1
                    else:
                        cls = 2
            elif len(segList) == 2:
                if max(numTMSegList) == 1:
                    if segList[0][0] == 0 and segList[1][1] == lengthAlignedMapArray:
                        cls = 0
                    elif segList[0][0] > 0 and segList[1][1] < lengthAlignedMapArray:
                        cls = 2
                    else:
                        cls = 3
                else: # at least one indels has >=2 TM
                    if segList[0][0] == 0 and segList[1][1] == lengthAlignedMapArray:
                        cls = 1
                    elif segList[0][0] > 0 and segList[1][1] < lengthAlignedMapArray:
                        cls = 2
                    else:
                        cls = 3
            else: # len(segList)>2
                isAllInter = True
                for seg in segList:
                    if seg[0] == 0 or seg[1] == lengthAlignedMapArray:
                        isAllInter = False
                        break
                if isAllInter:
                    cls = 2
                else:
                    cls = 3
            # print to std
            print "class = ", cls
            print record['mapTMline'][0]
            print record['mapTMline'][1]
            if idxGroup != numGroup:
                freq[idxGroup][cls] += 1
                subsum[idxGroup] += 1
#}}}

def AnaHelixPairCmpResultCmpClass_ps_bin(recordList,  #{{{
        dataTable,  classList, psAllGroup, rltyDict):
    numGroup = len(psAllGroup)/2
    numClass = len(classList)
    seqidttype = g_params['seqidttype']
    binsListSeqIDT = SEQIDT_GROUP_4
    numBinsSeqIDT = len(binsListSeqIDT)/2
    for record in recordList:
        if record == {}:
            continue
        id1 = record['id1']
        id2 = record['id2']
        # ps is using the minimum of the pair
        try:
            ps1 = rltyDict[id1]
        except KeyError:
            ps1 = -100.0
        try:
            ps2 = rltyDict[id2]
        except KeyError:
            ps2 = -100.0

        seqidt = lcmp.GetSeqIDT(record, seqidttype)
        idxSeqidt = GetIndexOfBins(seqidt, binsListSeqIDT, numBinsSeqIDT)
        try: 
            for key in ['min_ps', 'avg_ps', 'max_ps']:
                if key == 'min_ps':
                    ps = min(ps1,ps2)
                elif key == 'max_ps':
                    ps = max(ps1,ps2)
                elif key == 'avg_ps':
                    ps = (ps1+ps2)/2

                fpout_rlty_helixcmpclass = None
                if g_params['isPrintFileRltyHelixCmpclass']:
                    try:
                        fpout_rlty_helixcmpclass = g_params['fpout_rlty_helixcmpclass_'+ key]
                    except KeyError:
                        pass
                idxGroup = GetIndexOfBins(ps , psAllGroup, numGroup)
                freq = dataTable[key]['freq']
                freq_by_seqidt = dataTable[key]['freq_by_seqidt']
                subsum = dataTable[key]['subsum']
                for li in record["mapArray"]:
                    for mp in li:
                        idxClass = GetClassIndex(mp, classList)
                        if idxClass == -1:
                            continue
                        if fpout_rlty_helixcmpclass != None:
                            fpout_rlty_helixcmpclass.write("%.1f %s\n"%(ps, mp))
                        if idxGroup != numGroup and idxClass != numClass:
                            freq[idxGroup][idxClass] += 1
                            freq_by_seqidt[idxGroup][idxSeqidt] += 1
                            subsum[idxGroup] += 1
        except KeyError, IndexError:
            pass
#}}}
def AnaHelixPairCmpResultDGScore(recordList, dataTable, classList, seqIDTGroupList):#{{{
    table = dataTable['data']
    table_noshift = dataTable['data_noshift']
    table_shift = dataTable['data_shift']
    numGroup = len(seqIDTGroupList)/2
    numClass = len(classList)
    seqidttype = g_params['seqidttype']
    for record in recordList:
        if record == {}:
            continue
        seqidt = lcmp.GetSeqIDT(record, seqidttype)
        idxGroup = GetSeqIDTGroupIndex(seqidt, seqIDTGroupList)
        if "mapArray" in record and "dgscore" in record:
            alignedMapArray = GetAlignedMapArray_mp1(record['mapArray'])
#debug
#             if record['cmpclass'].find( "TM2SEQ") != -1:
#                 print record['id1'], record['id2']
#                 print alignedMapArray[0]
#                 print alignedMapArray[1]

            for i in xrange(2):
                idx1 = i
                idx2 = (i+1)%2
                numTM1 = len(record['mapArray'][idx1])
                numTM2 = len(record['mapArray'][idx2])
                array = record['mapArray'][i]
                dgscores = record['dgscore'][i]
#                 print dgscores
                if len(array) == len(dgscores) and len(array) > 0:
                    for j in xrange(len(array)):
                        mp = array[j]
                        dg = dgscores[j]
                        idxClass = GetClassIndex(mp, classList)
                        if idxClass == -1:
                            continue
                        if idxGroup != numGroup and idxClass != numClass:
                            table[idxClass][idxGroup].append(dg)
                    for j in xrange(len(alignedMapArray[i])):
                        tup = alignedMapArray[i][j]
                        mp = tup[0]
                        if mp in [0, 1, 2]:
                            idxClass = GetClassIndex(mp, classList)
                            dg = dgscores[tup[1]]
                            if mp in [0, 1]:
                                table_noshift[idxClass][idxGroup].append(dg)
                                table_shift[idxClass][idxGroup].append(dg)
                            elif mp == 2:
                                if (((j-1)>=0 and alignedMapArray[(i+1)%2][j-1][0] in [2]) 
                                    or ((j+1)<len(alignedMapArray[i]) and
                                        alignedMapArray[(i+1)%2][j+1][0] in [2])
                                    ):
                                    isShiftedTM = True
                                else:
                                    isShiftedTM = False
                                #print record['id1'], record['id2'], "isShiftedTM = ", isShiftedTM #debug
                                if not isShiftedTM:
                                    table_noshift[idxClass][idxGroup].append(dg)
                                else:
                                    table_shift[idxClass][idxGroup].append(dg)
                else:
                    print >> sys.stderr, "Bad array"
                    print >> sys.stderr, "  Map array", array
                    print >> sys.stderr, "  dgscore", dgscores
#}}}

def AnaPairCmpResultUnmappedTMPosition(recordList, dataTable, item,  #{{{
        topoalnDict, rltyDict, dupPairDict):
# analyse TM helix mapping, show frequency of non-TM2TM mapped helices
# must distinguish single-spanning and multi-spanning
# single spanning TM protein is not counted in (2012-10-11)
# item == 0:
#   position calculation using gapless sequence
#       pp = (beginTM + endTM)/2/seqLength 
# item == 1:
#   position calculation using aligned sequence
#   anything before the aligned region
#        pp = 0
#   anything after the aligned region
#        pp = 1
#   For the rest:
#        pp = (beginAlnTM + endAlnTM - 2*bAlnRegion)/2/lengthAlnRegion
# item == 2:
#   position calculation using gapless sequence
#   For the first TM:
#        pp = 0
#   For the last TM:
#        pp = 1
#   For the rest:
#        pp = (beginTM + endTM)/2/seqLength 
# item == 3:
#   position calculation using aligned sequence
#   determine the internal region
#   beginInternalRegion = max(beginAlignedRegion, endFirstAlnTM)
#   endInternalRegion = min(endAlignedRegion, beginLastAlnTM)
#   Before beginInternalRegion:
#       pp = 0
#   After endInternalRegion:
#       pp = 1
#   else:
#       pp = (beginAlnTM + endAlnTM - 2*beginInternalRegion)/2/lengthInternalRegion
# item == 4:(default method for determining N- & C- terminal) 
#   position calculation using aligned sequence
#   posFirstTMEnd = max(posFirstTMEnd1, posFirstTMEnd2)
#   posLastTMBeg = min(posLastTMBeg1, posLastTMBeg2)
#   interRegion = (posFirstTMEnd, posLastTMBeg)
#   lengthInterRegion = posLastTMBeg - posFirstTMEnd
#   anything after the aligned region
#   if pos < posFirstTMEnd:
#        pp = 0 (N-terminal)
#   elif pos > posLastTMBeg:
#        pp = 1 (C-terminal)
#   else:
#        pp = (beginAlnTM + endAlnTM - 2*posFirstTMEnd)/2/lengthInterRegion
# item == 5 (created 2013-06-07)
#   anything before the first TM2TM is N-terminal
#   anything after the last TM2TM is C-terminal
#    e.g. for the alignment
#           0   0   0   0  1  1 
#    2 1    0   0   0   0
#    N-Term                C-Term
#

    tablePosition = dataTable['dataPosition']
    tableNumContinuousTM  = dataTable['dataNumContinuousTM']
    seqidttype = g_params['seqidttype']
    for record in recordList:
        if record == {}:
            continue
        seqidt = lcmp.GetSeqIDT(record, seqidttype)
        rlty = GetRlty(record, rltyDict, "min")
        cmpclass = record['cmpclass']


        if not "mapArray" in  record:
            continue

        if (cmpclass.find('UNALIGNED') == 0 
                or cmpclass.find('AMBIGUOUS') == 0):
            continue

        cmpclass = record['cmpclass']
# This just works for pairwise_comparison_method 1 and 3. Ana5 is not working
        if cmpclass.find("TM2GAP|DUP") == 0:
            cmpclass = "DUP"
        else:
            cmpclass = cmpclass.split('|')[0]

        id1 = record['id1']
        id2 = record['id2']
        idList = [id1, id2]
        key = "%s-%s"%(id1, id2)

        try:
            pairaln = topoalnDict[key]
        except KeyError:
            print >> sys.stderr, "topology alignment does not find for %s"%(key)
            continue

        if len(pairaln[0]) != len(pairaln[1]):
            print >> sys.stderr, "Error! unequal length %d (%s) != %d (%s)" %(
                    id1, len(pairaln[0]), id2, len(pairaln[1]))
            continue
        posTMList_gapless = []
        for i in range(2): # get posTM for gapless seq
            posTMList_gapless.append(myfunc.GetTMPosition(pairaln[i].replace(GAP,"")))

        duppair = [(id1, len(pairaln[0].replace(GAP,""))),
                (id2, len(pairaln[1].replace(GAP,"")))]
# sort dupair as shorter - longer
        duppair = sorted(duppair, key=lambda x:x[1], reverse=False)
        key = (duppair[0][0], duppair[1][0])
        idList_duppair = [duppair[0][0], duppair[1][0]]
        try:
            duphit = dupPairDict[key]['hit']
#             print duphit
        except KeyError:
            duphit = []


        if item == 5:#(the default method, created 2013-06-07, updated 2013-06-07)
# first get the (ib1, ie1), (ib2, ie2) of alignned TM index
            mapArrayList = record['mapArray']
            numTM2TMList = [x.count(0) for x in mapArrayList]
            if numTM2TMList[0] != numTM2TMList[1]:
                msg = "numTM2TM not equal! numTM2TM of %s (%d) != numTM2TM2 of"\
                " %s (%id). Ignore"
                print >> sys.stderr, msg%(id1, numTM2TMList[0], id2,
                        numTM2TMList[1])
                continue
            firstLastTM2TMIndexList = []
            for i in xrange(2):
                firstLastTM2TMIndexList.append(GetFirstLastTM2TM(mapArrayList[i]))

# analyze the relationship of N-terminal i/o status and the number of
# difference in N terminal helices
            for side1 in ["Nterm", "Cterm"]:
                for side2 in ["Nterm", "Cterm"]:
                    if side1 == "Nterm":
                        diff_numTM = abs(firstLastTM2TMIndexList[0][0] -
                                firstLastTM2TMIndexList[1][0])
                    else:
                        diff_numTM = abs(
                                (len(mapArrayList[0])-firstLastTM2TMIndexList[0][1]-1) -
                                (len(mapArrayList[1])-firstLastTM2TMIndexList[1][1]-1))

                    if side2 == "Nterm":
                        state1 =  lcmp.GetNtermState(pairaln[0]) 
                        state2 =  lcmp.GetNtermState(pairaln[1]) 
                    else:
                        state1 =  lcmp.GetCtermState(pairaln[0]) 
                        state2 =  lcmp.GetCtermState(pairaln[1]) 

                    if state1 == state2:
                        diff_state = 1
                    else:
                        diff_state = -1
                    dataTable['dataDiffNumTM%s_%sStatus'%(side1,
                        side2)].append((diff_numTM, diff_state))

            for i in xrange(2):  # iterate over two proteins in the pair
# analyze locate for each unmapped TM helix
                # get the index of this seqid in the duplication pair
                idx_in_dupair = idList_duppair.index(idList[i])
                for j in xrange(len(mapArrayList[i])):
                    if j < firstLastTM2TMIndexList[i][0]:
                        pp = 0.0
                    elif j > firstLastTM2TMIndexList[i][1]:
                        pp = 1.0
                    else:
                        pp = 0.5
                    mp = mapArrayList[i][j]
                    if mp == _TM2GAP:
                        tablePosition[0].append((pp, seqidt, rlty))
                    elif mp == _TM2SEQ:
                        tablePosition[1].append((pp, seqidt, rlty))
                    idx = GetUnmappedTMPositionIndex(mp, cmpclass)
                    if idx >= 0:
                        tablePosition[idx].append((pp, seqidt, rlty))
# analyze distribution of continious TM2GAP or TM2SEQ
                str_maparray_list = []
                for x in mapArrayList[i]:
                    if x == -1:
                        str_maparray_list.append('X')
                    else:
                        str_maparray_list.append("%d"%x)
                str_maparray = "".join(str_maparray_list)
                for st in [1, 2]:  # for TM2GAP or TM2SEQ helix mapping
                    posContList = myfunc.GetSegPos(str_maparray, "%d"%st)
                    if len(posContList) > 0:
                        for (b,e) in posContList:
                            # check if this TM to gap indel is a duplication
                            isDup = IsIndelTMDuplicated(posTMList_gapless[i],
                                    range(b,e), duphit, idx_in_dupair)
                            if isDup:
                                ss_isDup = 'y'
                            else:
                                ss_isDup = 'n'
                            numTMSeg = e-b
                            if e-1 < firstLastTM2TMIndexList[i][0]:
                                pp = 0.0
                            elif b > firstLastTM2TMIndexList[i][1]:
                                pp = 1.0
                            else:
                                if (myfunc.coverage(b, e,
                                    firstLastTM2TMIndexList[i][0],
                                    firstLastTM2TMIndexList[i][1]+1) < e-b):
                                    msg = "Bad NCTerm for %s - %s, b,e (%d, %d)"\
                                            " but FstLst (%d, %d)"
                                    print >> sys.stderr, msg%(id1, id2, b, e-1,
                                            firstLastTM2TMIndexList[i][0],
                                            firstLastTM2TMIndexList[i][1])
                                pp = 0.5

                            if st == _TM2GAP:
                                tableNumContinuousTM[0].append((numTMSeg, pp,
                                    seqidt, rlty, id1, id2, ss_isDup))
                            elif st == _TM2SEQ:
                                tableNumContinuousTM[1].append((numTMSeg, pp,
                                    seqidt, rlty, id1, id2, ss_isDup))
                            idx = GetUnmappedTMPositionIndex(st, cmpclass)
                            if idx >= 0:
                                tableNumContinuousTM[idx].append((numTMSeg, pp,
                                    seqidt, rlty, id1, id2, ss_isDup))
        else:
            if item in [1, 3,4]:
                (interB, interE) = GetInterRegion_pairwise(pairaln, item)
            if DEBUG_UNMAPPED_TM_POSITION:
                print "Unmapped TM position for %s - %s" %(id1, id2)
                print "mapArray1:", record['mapArray'][0]
                print "mapArray2:", record['mapArray'][1]
                print "TopoAln1:", pairaln[0]
                print "TopoAln2:", pairaln[1]
            for i in xrange(2):
                idx_in_dupair = idList_duppair.index(idList[i])
                mapArray = record['mapArray'][i]
                seqid = record['id%d'%(i+1)]
                numTM  = len(mapArray)
                if numTM <= 1: #ignore single spanning TM protein in location analysis
                    continue
                if item in [0, 2]: # gapless
                    alnTopo = pairaln[i].replace("-","")
                    posTM = myfunc.GetTMPosition_gapless(alnTopo)
                    seqLength = record['seqLength%s'%(i+1)]
                elif item in [1,3,4]: # aligned
                    alnTopo = pairaln[i]
                    posTM = myfunc.GetTMPosition(alnTopo)
                    lengthInterRegion = (interE - interB)


                str_maparray_list = ["%d"%x for x in mapArray]
                str_maparray = "".join(str_maparray_list)
                neighbour_char = "0"
                for st in [1, 2]:
                    posContList = myfunc.GetSegPos(str_maparray, "%d"%st)
                    posContList = FilterSegPos(posContList, str_maparray, neighbour_char)
                    if len(posContList) > 0:
                        for (b,e) in posContList:
                            # check if this TM to gap indel is a duplication
                            isDup = IsIndelTMDuplicated(posTMList_gapless[i],
                                    range(b,e), duphit, idx_in_dupair)
                            if isDup:
                                ss_isDup = 'y'
                            else:
                                ss_isDup = 'n'
                            numTMSeg = e-b
                            locTM =  (posTM[b][0] + posTM[e-1][1])/2.0
                            if item == 0:
                                pp = locTM/seqLength
                            elif item == 2:
                                if b == 0:
                                    pp = 0.0
                                elif e == numTM:
                                    pp = 1.0
                                else:
                                    pp = locTM/seqLength
                            elif item in [1,3,4]:
                                if locTM < interB:
                                    pp = 0.0
                                elif locTM > interE:
                                    pp = 1.0
                                else:
                                    pp = float(locTM - interB)/lengthInterRegion
                            if st == _TM2GAP:
                                tableNumContinuousTM[0].append((numTMSeg, pp,
                                    seqidt, rlty, id1, id2, ss_isDup))
                            elif st == _TM2SEQ:
                                tableNumContinuousTM[1].append((numTMSeg, pp,
                                    seqidt, rlty, id1, id2, ss_isDup))
                            idx = GetUnmappedTMPositionIndex(st, cmpclass)
                            if idx >= 0:
                                tableNumContinuousTM[idx].append((numTMSeg, pp,
                                    seqidt, rlty, id1, id2, ss_isDup))
                for j in xrange(len(mapArray)):
                    mp = mapArray[j]
                    locTM = (posTM[j][0]+posTM[j][1])/2.0
                    if item == 0:
                        pp = locTM / seqLength
                    elif item == 2:
                        if j == 0:
                            pp = 0.0
                        elif j == numTM-1:
                            pp = 1.0
                        else:
                            pp = locTM/seqLength
                    elif item in [1,3,4]:
                        if locTM < interB:
                            pp = 0.0
                        elif locTM > interE:
                            pp = 1.0
                        else:
                            pp = float(locTM - interB)/lengthInterRegion
                    if mp == _TM2GAP:
                        tablePosition[0].append((pp, seqidt, rlty))
                    elif mp == _TM2SEQ:
                        tablePosition[1].append((pp, seqidt, rlty))

                    idx = GetUnmappedTMPositionIndex(mp, cmpclass)
                    if idx >= 0:
                        tablePosition[idx].append((pp, seqidt, rlty))
                    if DEBUG_UNMAPPED_TM_POSITION and item == 1:
                        print "position (%s) = %6.3f" %(seqid, pp)
                    elif mp != 0 and DEBUG_UNMAPPED_TM_POSITION :
                        print >> sys.stderr, seqid, ": index < 0. TMMapclass = ", \
                                mp, "cmpclass=", cmpclass, "index=", idx
    return 0
#}}}

def CalHist_obsolete(li, binsize):#{{{
    if len(li) <= 0:
        return ([], [])
    else:
        min_v = min(li)
        max_v = max(li)
        numbin = max(1, int(ceil((max_v - min_v) / binsize)))
        return numpy.histogram(li, bins=numbin, density=True)
        #if numpy.__version__ >= "1.6.0":
        #else:
        #    return numpy.histogram(li, bins=numbin, normed=True)
#}}}
def CalHist_grouping(li):#{{{
    """
    """
    hist = {}
    for x in li:
        if not x in hist:
            hist[x] = 0
        hist[x] += 1
    return hist

#}}}
def GetBinMax(x, binsize):
    halfbin = binsize/2.0
    return int(x/halfbin)*halfbin+halfbin

def GetBinMin(x, binsize):
    halfbin = binsize/2.0
    return int(x/halfbin)*halfbin-halfbin

def CalHist(li, binsize, binList):#{{{
    if len(li) <= 0:
        return ([], [])
    else:
        if binList == []:
            if binsize > 0:
                min_v = min(li)
                max_v = max(li)
                histmin = GetBinMin(min_v, binsize)
                histmax = GetBinMax(max_v, binsize)
                numBin = max(1, int(ceil(histmax - histmin)/binsize))
                for i in xrange(numBin):
                    binList.append(histmin + binsize * i)
                    binList.append(histmin + binsize * (i+1))
            else:
                return ([],[])

        numBin = len(binList)/2
        freq = [0] * numBin
        edges = []
        edges.append(binList[0])
        for i in xrange(numBin):
            edges.append(binList[2*i+1])
        for x in li:
            idx = GetIndexOfBins(x, binList, numBin)
            if idx < numBin and idx >= 0:
                freq[idx] += 1
        return (freq, edges)
#}}}
def CalHistDGScore(table2d, seqIDTGroupList):#{{{
    isEvodist = g_params['isEvodist'] 
    if isEvodist:
        thHigh = g_params['thHigh_evodist']
    else:
        thHigh = g_params['thHigh_seqidt']
    numGroup = len(seqIDTGroupList)/2
    listLow = []
    listHigh = []
    listAll = []
    histDict = {}
    binsize = 1.0
    for i in xrange(numGroup):
        item = "%g-%g"%(seqIDTGroupList[2*i], seqIDTGroupList[2*i+1])
        th = seqIDTGroupList[2*i]
        listAll += table2d[i]
        if th < thHigh:
            listLow += table2d[i]
        else:
            listHigh += table2d[i]
        histDict[item] = CalHist(table2d[i], binsize, [])
    histDict['low'] = CalHist(listLow, binsize, [])
    histDict['high'] = CalHist(listHigh, binsize, [])
    histDict['all'] = CalHist(listAll, binsize, [])
    return histDict
#}}}
def InitXY(N):#{{{
    lx = range(N+1)
    ly = [0]*(N+1)
    return (lx, ly)
    #}}}
def CalHistNumTM(data, positionBins, isCmpDup):#{{{
    """
    data is a list of tuples, each tuple contains 
    (numTM, pos, seqidt, rlty, id1, id2, ss_isDup)
    positionBins 
    the format of the reture data
    dt[0] dt[1..3] according to bins
    1 3
    2 50
    2 10
    """
# the total number of classes of histograms is (ncls+1)
    MAX_NUMTM = 100
    dt = {}
    numBin = len(positionBins)/2
    if isCmpDup:
        ncls = numBin + 1
    else:
        ncls = numBin

    for i in xrange(ncls):
        dt[i] = InitXY(MAX_NUMTM)
    dt['all'] = InitXY(MAX_NUMTM)

    for x in data:
        if x[0] < MAX_NUMTM:
            if x[6] == 'y' and isCmpDup:
                idx = numBin
                dt[idx][1][x[0]] += 1
                dt['all'][1][x[0]] += 1
            else:
                idx = GetIndexOfBins(x[1], positionBins, numBin)
                if idx < numBin and idx >= 0:
                    dt[idx][1][x[0]] += 1
                    dt['all'][1][x[0]] += 1
    return dt
#}}}
def CalHistNumContinuousTM(table2d, termmode, seqidtmode,#{{{
        th_rlty, isCmpDup):
    """
    Calculate the histogram of the frequency of number of TM helices in the indels
    when isCmpDup is True, there will be four columns for this histogram
    N-term, Internal, C-term and Duplicated
    """
    histDict = {}
    thHigh = g_params['thHigh_seqidt']
    positionBins = [
            0, 1e-9,
            1e-9, 1.0,
            1.0, 1.0+1e-9]
    numGroup = len(table2d)
    numBin = len(positionBins)/2
    isRLTYSupplied = g_params['isRLTYSupplied']
    for i in xrange(numGroup):
        dt = []
        for j in xrange(len(table2d[i])):
            pp =     table2d[i][j][1] 
            seqidt = table2d[i][j][2]
            rlty   = table2d[i][j][3]
            if ((seqidtmode == "low" and seqidt >= thHigh) or 
                    (seqidtmode == "high" and seqidt < thHigh)):
                continue
            if isRLTYSupplied and rlty < th_rlty:
                continue
            if termmode == "withoutterm":
                if not (pp == 0.0 or pp == 1.0):
                    dt.append(table2d[i][j])
            else:
                dt.append(table2d[i][j])
        histDict[i] = CalHistNumTM(dt, positionBins, isCmpDup)
    return histDict
#}}}
def CalHistUnmappedTMPosition(table2d, termmode, seqidtmode, th_rlty, binList):#{{{
    histDict = {}
    binsize = 0.05
    numGroup = len(table2d)
    isRLTYSupplied = g_params['isRLTYSupplied']
    thHigh = g_params['thHigh_seqidt']
    for i in xrange(numGroup):
        dt = []
        for j in xrange(len(table2d[i])):
            seqidt = table2d[i][j][1]
            rlty = table2d[i][j][2]
            if ((seqidtmode == "low" and seqidt >= thHigh) or 
                    (seqidtmode == "high" and seqidt < thHigh)):
                continue
            if isRLTYSupplied and rlty < th_rlty:
                continue
            if termmode == "withoutterm":
                if not (table2d[i][j][0] == 0.0 or table2d[i][j][0] == 1.0):
                    dt.append(table2d[i][j][0])
            else:
                dt.append(table2d[i][j][0])
        histDict[i] = CalHist(dt, binsize, binList)
    return histDict
#}}}

def CountNCtermInter(ana, idxGroup, classList, numClass, freq, subsum):#{{{
    if ana == {}:
        return
    for cls in ['Nterm', 'Cterm']:
        if cls in ana and ana[cls]['numTMunmapped'] >0:
            idxClass = GetClassIndex(cls, classList)
            if idxClass == -1:
                continue
            if idxClass != numClass:
                freq[idxGroup][idxClass] += 1
                subsum[idxGroup] += 1
    if 'internal' in ana:
        idxClass = GetClassIndex('internal', classList)
        if idxClass != -1 and idxClass != numClass:
            for m in range(len(ana['internal'])):
                freq[idxGroup][idxClass] += 1
                subsum[idxGroup] += 1
#}}}
def AnaPairCmpResultNCTermInter(recordList, dataTable, classList, seqIDTGroupList):#{{{
# freq1 is counting for all occurrences of unmapped TMs at N, C-terminals and
# internal regions.
#
# freq2 is counting for pairs with unmapped TMs at N,C-terminals and internal
# regions. In this case,  for each pair, occurrence of unmapped TMs at internal
# regions is counted once even there are several.
    numGroup = len(seqIDTGroupList)/2
    numClass = len(classList)
    freq1 = dataTable['freq1']
    subsum1 = dataTable['subsum1']
    freq2 = dataTable['freq2']
    subsum2 = dataTable['subsum2']
    seqidttype = g_params['seqidttype']
    for record in recordList:
        if record == {}:
            continue
        if record['cmpclass'] != 'DIFF':
            continue
        seqidt = lcmp.GetSeqIDT(record, seqidttype)
        idxGroup = GetSeqIDTGroupIndex(seqidt, seqIDTGroupList)
        if idxGroup == numGroup:
            continue
        # count for this record
        localFreq1 = []
        localSubsum1 = []
        for i in xrange(numGroup):
            localFreq1.append([0]*numClass)
            localSubsum1.append(0)

        CountNCtermInter(record['ana1'], idxGroup, classList, numClass,
                localFreq1, localSubsum1)
        CountNCtermInter(record['ana2'], idxGroup, classList, numClass,
                localFreq1, localSubsum1)

        # subsum1 and subsum2 are the number of pairs in each group.
        for i in xrange(numGroup):
            #subsum1[i] += localSubsum1[i]
            #subsum1[i] += bool(sum(localSubsum1))
            subsum1[i] += bool(localSubsum1[i])
            subsum2[i] += bool(localSubsum1[i])
            for j in xrange(numClass):
                freq1[i][j] += localFreq1[i][j]
                freq2[i][j] += bool(localFreq1[i][j])

#}}}
def CountNumTMDiffRegion(ana, idxGroup, classList, numClass,freq, subsum, MAX_NUMTM):#{{{
    if ana == {}:
        return
    for cls in ['Nterm', 'Cterm']:
        if cls in ana and ana[cls]['numTMunmapped'] >0:
            idxClass = GetClassIndex(cls, classList)
            if idxClass != numClass:
                numTMunmapped = ana[cls]['numTMunmapped']
                if numTMunmapped < MAX_NUMTM:
                    freq[idxGroup][idxClass][numTMunmapped] += 1
                    subsum[idxGroup][idxClass] += 1
    if 'internal' in ana:
        idxClass = GetClassIndex('internal', classList)
        if idxClass != numClass:
            for m in range(len(ana['internal'])):
                numTMunmapped = ana['internal'][m]['numTMunmapped']
                if numTMunmapped < MAX_NUMTM:
                    freq[idxGroup][idxClass][numTMunmapped] += 1
                    subsum[idxGroup][idxClass] += 1
#}}}

def AnaPairCmpResultNumTMDistribution(recordList, dataTable, classList, #{{{
        seqIDTGroupList, MAX_NUMTM):
    numGroup = len(seqIDTGroupList)/2
    numClass = len(classList)
    freq = dataTable['freq']
    subsum = dataTable['subsum']
    seqidttype = g_params['seqidttype']
    for record in recordList:
        if record == {}:
            continue
        if record['cmpclass'] != 'DIFF':
            continue

        seqidt = lcmp.GetSeqIDT(record, seqidttype)
        idxGroup = GetSeqIDTGroupIndex(seqidt, seqIDTGroupList)
        if idxGroup == numGroup:
            continue
        # count for this record
        CountNumTMDiffRegion(record['ana1'], idxGroup, classList, numClass,
                freq, subsum, MAX_NUMTM)
        CountNumTMDiffRegion(record['ana2'], idxGroup, classList, numClass,
                freq, subsum, MAX_NUMTM)

#}}}
def AnaPairCmpResultNumTMHeatMap(recordList, dataTable, classList):#{{{
    for record in recordList:
        if record == {}:
            continue
        try:
            cmpclass = record['cmpclass']
            numTM1 = record['numTM1']
            numTM2 = record['numTM2']
        except KeyError:
            print >> sys.stderr, "KeyError for the record %s - %s" %(
                    record['id1'], record['id2'])
            continue


        if cmpclass.find('|') != -1:
            cmpclass = cmpclass.split('|')[0]

        for cls in classList:
            dt = dataTable[cls]
            data = dt['data']
            if ((cls == "ALL") or 
                    (cls == "ALL_VARY" and cmpclass in ["TM2GAP", "TM2SEQ",
                        "TM2GAP_AND_TM2SEQ"]) or 
                    (cls == cmpclass) or
                    (cls == "ONLY_DUP" and record['isDup'] == True)):
                data[min(numTM1,numTM2)][max(numTM1, numTM2)] += 1
                if max(numTM1, numTM2) > dt['maxNumTM']:
                    dt['maxNumTM'] = max(numTM1, numTM2)
                dt['numPair'] += 1
    return 0
#}}}

def WriteTM2SEQDGScoreSegment(dataTable, outfile):#{{{
    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    for rd in dataTable:
# item description
# seqid1 seqid2 cmpclass seqidt begin_TM end_TM segTM_nogap segNonTM_nogap segTM segNonTM seqidt_seg 
        fpout.write("%s %s %19s %6.1f %4d %4d %21s %21s %s %s %6.1f\n"%(
            rd[0], rd[1],
            rd[2], rd[3], rd[4], rd[5], 
            rd[6].replace("-",""), rd[7].replace("-",""),
            rd[6],rd[7],rd[8]))
    myfunc.myclose(fpout)
#}}}
def WriteTM2GAPDGScoreSegment(dataTable, outfile):#{{{
    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    for rd in dataTable:
# item description
# seqid1 seqid2 cmpclass seqidt begin_TM end_TM segTM_nogap segNonTM_nogap segTM segNonTM seqidt_seg 
        segNonTM_nogap = "A"*21
        fpout.write("%s %s %19s %6.1f %4d %4d %21s %21s %s %s %6.1f\n"%(
            rd[0], rd[1],
            rd[2], rd[3], rd[4], rd[5], 
            rd[6].replace("-",""), segNonTM_nogap,
            rd[6],rd[7],rd[8]))
    myfunc.myclose(fpout)
#}}}

def WriteTM2SEQDGScoreSegment2(dataTable, outfile):#{{{
    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    cnt = 0
    for rd in dataTable:
        fpout.write("%4d : %s %s %19s %6.1f %6.1f %4d %4d %21s %21s\n"%(cnt+1, rd[0], rd[1],
            rd[2], rd[3], rd[8], rd[4], rd[5], rd[6].replace("-",""),
            rd[7].replace("-","")))
        fpout.write("%s\n"%rd[9])
        fpout.write("%s\n"%rd[6])
        fpout.write("%s\n"%rd[7])
        fpout.write("%s\n"%rd[10])
        fpout.write("\n")
        cnt += 1
    myfunc.myclose(fpout)
#}}}

def WritePairInfo(pairinfoList, outfile):#{{{
    try:
        fpout = open(outfile, "w")
        for tup in pairinfoList:
#seqid1 seqid2 NtermState1 NtermState2 numTM1 numTM2 len1 len2 seqidt1 pfamid1 ...
            print >> fpout, "%s %s %2s %2s %2d %2d %4d %4d %6.1f %s"%(tup[0], tup[1], tup[2],
                    tup[3], tup[4], tup[5], tup[6], tup[7], tup[8],
                    " ".join(tup[9]))
        fpout.close()
    except IOError:
        print >> sys.stderr, "Failed to write to file %s"%outfile
#}}}
def WriteTable2D(freq, subsum, classList, seqIDTGroupList, outfile):#{{{
    #works with dataCmpClass and dataNCtermInter
    try:
        fpout = open(outfile, "wb")
    except IOError:
        print >> sys.stderr, "Failed to write to file %s"%outfile
        return 1

    numGroup = len(seqIDTGroupList)/2
    numClass = len(classList)
    fpout.write("%4s %7s"%("#Idx","SeqIDT"))
    for cls in classList:
        fpout.write(" %9s"%(cls))
    fpout.write(" %9s"%"Maximum")
    fpout.write(" %10s"%"Occurrence")
    fpout.write("\n")
    for i in xrange(numGroup):
        stridtrange="%g-%g"%(seqIDTGroupList[i*2], seqIDTGroupList[i*2+1])
        fpout.write("%-4d %7s"%(i,stridtrange))
        for j in xrange(numClass):
            fpout.write(" %9.3f"%(myfunc.FloatDivision(freq[i][j],subsum[i])*100))
        fpout.write(" %9.3f"%max( [ myfunc.FloatDivision(freq[i][j],subsum[i])*100 for j in range(numClass)] ))
        fpout.write(" %10d"%subsum[i])
        fpout.write("\n")

    fpout.write("%-12s"%("#sum"))
    totalOccur = [0]*numClass
    totalSum = sum(subsum)
    for j in xrange(numClass):
        for i in xrange(numGroup):
            totalOccur[j] += freq[i][j]
        fpout.write(" %9.3f"%(myfunc.FloatDivision(totalOccur[j],totalSum)*100))
    fpout.write(" %9.3f"%max( [ myfunc.FloatDivision(totalOccur[j],totalSum)*100 for j
        in range(numClass)] ))
    fpout.write(" %10d"%totalSum)
    fpout.write("\n")

    fpout.close()
    return 0
#}}}
def WriteHistogramData(freq, subsum, tag, classList, binsList, outfile):#{{{
    #works with dataCmpClass and dataNCtermInter
    try:
        fpout = open(outfile, "wb")
    except IOError:
        print >> sys.stderr, "Failed to write to file %s"%outfile
        return 1

    numGroup = len(binsList)/2
    numClass = len(classList)
    fpout.write("%4s %7s"%("#Idx",tag))
    for cls in classList:
        fpout.write(" %9s"%(cls))
    fpout.write(" %9s"%"Maximum")
    fpout.write(" %10s"%"Occurrence")
    fpout.write("\n")
    for i in xrange(numGroup):
        stridtrange="%g"%((binsList[i*2]+binsList[i*2+1])/2)
        fpout.write("%-4d %7s"%(i,stridtrange))
        for j in xrange(numClass):
            fpout.write(" %9.3f"%(myfunc.FloatDivision(freq[i][j],subsum[i])*100))
        fpout.write(" %9.3f"%max( [ myfunc.FloatDivision(freq[i][j],subsum[i])*100 for j in range(numClass)] ))
        fpout.write(" %10d"%subsum[i])
        fpout.write("\n")

    fpout.write("%-12s"%("#sum"))
    totalOccur = [0]*numClass
    totalSum = sum(subsum)
    for j in xrange(numClass):
        for i in xrange(numGroup):
            totalOccur[j] += freq[i][j]
        fpout.write(" %9.3f"%(myfunc.FloatDivision(totalOccur[j],totalSum)*100))
    fpout.write(" %9.3f"%max( [ myfunc.FloatDivision(totalOccur[j],totalSum)*100 for j
        in range(numClass)] ))
    fpout.write(" %10d"%totalSum)
    fpout.write("\n")

    fpout.close()
    return 0
#}}}
def WriteHistogram_grouping(hist, outfile):#{{{
    try:
        fpout = open(outfile, "w")
        totalcount = 0
        for key in hist:
            totalcount += hist[key]
        print >> fpout, "#%s %s %s %s"%("Idx", "NumDiffTM", "%#", "Count")
        cnt = 0
        for key in hist:
            print >> fpout, "%2d %3d %6.3f %5d"%(
                    cnt, key,
                    myfunc.FloatDivision(hist[key], totalcount),
                    hist[key])
            cnt += 1
        fpout.close()
    except IOError:
        print >> sys.stderr, "Failed to write to file %s"%(outfile)
#}}}
def WriteNumPair2D(freq, binsListX, binsListY, outfile):#{{{
    try:
        fpout = open(outfile, "wb")
    except IOError:
        print >> sys.stderr, "Failed to write to file %s"%outfile
        return 1

    numGroupX = len(binsListX)/2
    numGroupY = len(binsListY)/2
    fpout.write("%4s %7s"%("#Idx","RLTY"))
    for j in  xrange(numGroupY):
        strx="%g"%((binsListY[j*2]+binsListY[j*2+1])/2)
        fpout.write(" %9s"%(strx))
    fpout.write(" %10s"%"Subsum")
    fpout.write("\n")
    for i in xrange(numGroupX):
        strx="%g"%((binsListX[i*2]+binsListX[i*2+1])/2)
        fpout.write("%-4d %7s"%(i,strx))
        for j in xrange(numGroupY):
            fpout.write(" %9d"%(freq[i][j]))
        fpout.write(" %10d"%(sum(freq[i])))
        fpout.write("\n")
# last line
    fpout.write("%-12s"%("#sum"))
    totalOccur = [0]*numGroupY
    for j in xrange(numGroupY):
        for i in xrange(numGroupX):
            totalOccur[j] += freq[i][j]
        fpout.write(" %9d"%(totalOccur[j]))
    totalSum = sum(totalOccur)
    fpout.write(" %10d"%totalSum)
    fpout.write("\n")
    
#}}}
def WriteHelixDGScore(table2d, seqIDTGroupAll, outfile):#{{{
    try:
        fpout = open(outfile, "w")
        numSeqIDTGroupAll = len(seqIDTGroupAll)/2
        for i in xrange(numSeqIDTGroupAll):
            fpout.write("%d-%d"%(seqIDTGroupAll[2*i], seqIDTGroupAll[2*i+1]))
            fpout.write(" %5d"%len(table2d[i]))
            for j in xrange(len(table2d[i])):
                fpout.write(" %6.3f"%table2d[i][j])
            fpout.write("\n")
        fpout.close()
        return 0
    except IOError:
        print >> sys.stderr, "Failed to write to file %s"%outfile
        return 1
#}}}
def WriteUnmappedTMPosition(table2d,  outfile, termmode, seqidtmode, th_rlty):#{{{
    try:
        isRLTYSupplied = g_params['isRLTYSupplied']
        fpout = open(outfile, "w")
        numRow = len(table2d)
        thHigh = g_params['thHigh_seqidt']
        for i in xrange(numRow):
            dt = []
            for j in xrange(len(table2d[i])):
                seqidt = table2d[i][j][1]
                rlty = table2d[i][j][2]
                if ((seqidtmode == "low" and seqidt >= thHigh) or 
                        (seqidtmode == "high" and seqidt < thHigh)):
                    continue
                if isRLTYSupplied and rlty < th_rlty:
                    continue
                if termmode == "withoutterm":
                    if not (table2d[i][j][0] == 0.0 or table2d[i][j][0] == 1.0):
                        dt.append(table2d[i][j][0])
                else:
                    dt.append(table2d[i][j][0])

            fpout.write("%-3d"%(i))
            fpout.write(" %6d"%len(dt))
            for j in xrange(len(dt)):
                fpout.write(" %6.3f"%dt[j])
            fpout.write("\n")
        fpout.close()
        return 0
    except IOError:
        print >> sys.stderr, "Failed to write to %s in function %s"%(outfile,
                sys._getframe().f_code.co_name)
        return 1
#}}}
def WriteNumContinuousTM_obsolete(table2d, outfile, termmode, seqidtmode, th_rlty):#{{{
    try:
        fpout = open(outfile, "w")
        thHigh = g_params['thHigh_seqidt']
        numRow = len(table2d)
        isRLTYSupplied = g_params['isRLTYSupplied']
        for i in xrange(numRow):
            dt = []
            for j in xrange(len(table2d[i])):
                pp =  table2d[i][j][1] 
                seqidt = table2d[i][j][2]
                rlty = table2d[i][j][3]
                if ((seqidtmode == "low" and seqidt >= thHigh) or 
                        (seqidtmode == "high" and seqidt < thHigh)):
                    continue
                if isRLTYSupplied and rlty < th_rlty:
                    continue
                if termmode == "withoutterm":
                    if not (pp == 0.0 or pp == 1.0):
                        dt.append(table2d[i][j])
                else:
                    dt.append(table2d[i][j])

            fpout.write("%-3d"%(i))
            fpout.write(" %6d"%len(dt))
            for j in xrange(len(dt)):
                fpout.write(" %3.0f %6.3f"%(dt[j][0], dt[j][1]))
            fpout.write("\n")
        fpout.close()
        return 0
    except IOError:
        print >> sys.stderr, "Failed to write to %s in function %s"%(outfile,
                sys._getframe().f_code.co_name)
        return 1
#}}}
def WriteNumContinuousTM(table2d, outfile, termmode, seqidtmode, th_rlty, #{{{
        seqid2pfamidDict, numTMDict):
    try:
        fpout = open(outfile, "w")
        thHigh = g_params['thHigh_seqidt']
        numRow = len(table2d)
        isRLTYSupplied = g_params['isRLTYSupplied']
        for i in xrange(numRow):
            dt = []
            for j in xrange(len(table2d[i])):
                pp =  table2d[i][j][1] 
                seqidt = table2d[i][j][2]
                rlty = table2d[i][j][3]
                if ((seqidtmode == "low" and seqidt >= thHigh) or 
                        (seqidtmode == "high" and seqidt < thHigh)):
                    continue
                if isRLTYSupplied and rlty < th_rlty:
                    continue
                if termmode == "withoutterm":
                    if not (pp == 0.0 or pp == 1.0):
                        dt.append(table2d[i][j])
                else:
                    dt.append(table2d[i][j])


            fpout.write("#%-3d %6d\n"%(i, len(dt)))
            fpout.write("#id1  id2   Type   nTM  Pos  nTM1  nTM2  SeqIDT isDup PfamID\n")
            for j in xrange(len(dt)):
                id1 = dt[j][4]
                id2 = dt[j][5]
                try:
                    pfamidlist1 = seqid2pfamidDict[id1]
                    pfamidlist2 = seqid2pfamidDict[id2]
                    common_pfamidlist = list(set(pfamidlist1) & set(pfamidlist2))
                except KeyError:
                    common_pfamidlist = []
                fpout.write("%s %s %2d %3.0f %6.3f %2d %2d %6.2f %6s "%(
                    dt[j][4], dt[j][5], i, dt[j][0], dt[j][1], numTMDict[id1],
                    numTMDict[id2], dt[j][2], dt[j][6]))
                for pfamid in common_pfamidlist:
                    fpout.write("%s "%(pfamid))
                fpout.write("\n")
            fpout.write("\n")
        fpout.close()
        return 0
    except IOError:
        print >> sys.stderr, "Failed to write to %s in function %s"%(outfile,
                sys._getframe().f_code.co_name)
        return 1
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

        for i in xrange(0, maxNumTM+1):
            if mode == "norm_diag":
                scale = scale_norm_diag
            elif mode == "norm_row":
                scale = myfunc.FloatDivision(count,sum(data[i]))
            for j in xrange(0, maxNumTM+1):
                if mode == "norm_col":
                    scale = scale_norm_col_list[j]
                fpout.write(" %6.3g"%(myfunc.FloatDivision(data[i][j],count)*scale*100))
            fpout.write("\n")
        fpout.close()
        return 0
    except IOError:
        print >> sys.stderr, "Failed to write to file %s"%outfile
        return 1
#}}}
def WriteHelixDGScoreHistogram(histTableDict, seqIDTGroupAll, outfile):#{{{
    try:
        fpout = open(outfile, "w")
        fpout.write("#")
        numSeqIDTGroupAll = len(seqIDTGroupAll)/2
        for i in xrange(numSeqIDTGroupAll):
            ss = "%g-%g"%(seqIDTGroupAll[2*i], seqIDTGroupAll[2*i+1])
            if i == 0:
                fpout.write("#%-7s %7s"%(ss+"X", ss+"Y"))
            else:
                fpout.write(" %7s %7s"%(ss+"X", ss+"Y"))
        fpout.write(" %7s %7s"%("LowX", "LowY"))
        fpout.write(" %7s %7s"%("HighX", "HighY"))
        fpout.write(" %7s %7s"%("AllX", "AllY"))
        fpout.write("\n")

        binsList = []
        for i in xrange(numSeqIDTGroupAll):
            ss = "%g-%g"%(seqIDTGroupAll[2*i], seqIDTGroupAll[2*i+1])
            binsList.append(len(histTableDict[ss][0]))
        binsList.append(len(histTableDict['low'][0]))
        binsList.append(len(histTableDict['high'][0]))
        binsList.append(len(histTableDict['all'][0]))
        maxBin = max(binsList)
        #maxBin = 50

        for i in xrange(maxBin):
            for j in xrange(numSeqIDTGroupAll):
                ss = "%g-%g"%(seqIDTGroupAll[2*j], seqIDTGroupAll[2*j+1])
                (hist, edges) = histTableDict[ss]
                total = sum(hist)
                numbin = len(hist)
                if i < numbin:
                    fpout.write(" %7.2f %7.3f"%((edges[i]+edges[i+1])/2.0,
                        hist[i]/float(total)))
                else:
                    fpout.write(" %7s %7s"%("?", "?"))
            for ss in ['low', 'high', 'all']:
                (hist, edges) = histTableDict[ss]
                total = sum(hist)
                numbin = len(hist)
                fpout.write(" ")
                if i < numbin:
                    fpout.write(" %7.2f %7.3f"%((edges[i]+edges[i+1])/2.0,
                        hist[i]/float(total)))
                else:
                    fpout.write(" %7s %7s"%("?", "?"))
            fpout.write("\n")
        fpout.close()
        return 0
    except IOError:
        print >> sys.stderr, "Failed to write to file %s"%outfile
        return 1
#}}}
def WriteNumContinuousTMHistogram(histTableDict, outfile, isNorm):#{{{
    try:
        fpout = open(outfile, "w")
        numItem1 = len(histTableDict)
#         print
#         print histTableDict
#         print 
        for i in xrange(numItem1):
            numItem2 = len(histTableDict[i])
            for j in xrange(numItem2-1):
                ss = "%d-%d"%(i,j)
                if i == 0 and j == 0:
                    fpout.write("#%-7s %7s"%(ss+"X", ss+"Y"))
                else:
                    fpout.write(" %7s %7s"%(ss+"X", ss+"Y"))
            for j in ['all']:
                ss = "%d-%s"%(i,j)
                fpout.write(" %7s %7s"%(ss+"X", ss+"Y"))
        fpout.write("\n")

        maxBin = len(histTableDict[0][0][0])

        for irow in xrange(1, maxBin):
            for i in xrange(numItem1):
                numItem2 = len(histTableDict[i])
                for j in xrange(numItem2-1):
                    (lx, ly) = histTableDict[i][j]
                    if isNorm:
                        count = sum(ly)
                        fpout.write(" %7.0f %7.3f"%(lx[irow],
                            myfunc.FloatDivision(ly[irow],count)*100))
                    else:
                        fpout.write(" %7.0f %7.0f"%(lx[irow],ly[irow]))
                for j in ['all']:
                    (lx, ly) = histTableDict[i][j]
                    if isNorm:
                        count = sum(ly)
                        fpout.write(" %7.0f %7.3f"%(lx[irow],
                            myfunc.FloatDivision(ly[irow],count)*100))
                    else:
                        fpout.write(" %7.0f %7.0f"%(lx[irow],ly[irow]))
            fpout.write("\n")

        for i in xrange(numItem1):
            numItem2 = len(histTableDict[i])
            for j in xrange(numItem2-1):
                ss = "Total"
                (lx, ly) = histTableDict[i][j]
                count = float(sum(ly))
                if i == 0 and j == 0:
                    fpout.write("#%-7s %7.0f"%(ss, count))
                else:
                    fpout.write(" %7s %7.0f"%(ss, count))
            for j in ['all']:
                ss = "Total"
                (lx, ly) = histTableDict[i][j]
                count = float(sum(ly))
                fpout.write(" %7s %7.0f"%(ss, count))
        fpout.write("\n")
        fpout.close()
        return 0
    except IOError:
        print >> sys.stderr, "Failed to write to file %s"%outfile
        return 1
#}}}
def WriteNumContinuousTMHistogram_grouped(histTableDict, outfile, isNorm):#{{{
    try:
        typeList=[
        "N-terminus",
        "Internal",
        "C-terminus",
        "Duplicated"]

        fpout = open(outfile, "w")
        numItem1 = len(histTableDict) #number of groups, e.g. TM2GAP in all, TM2SEQ in all
#         print
#         print histTableDict
#         print 
        # print header line
        fpout.write("#%-10s"%("Type"))
        for i in xrange(numItem1):
            fpout.write(" %6s %6s %6s"%("1TM", "2TM", "+2TM"))
        fpout.write("\n")

        # print content
        numItem2 = len(histTableDict[0])
        count_total_list = [0]*numItem1 # total count for each group
        for i in xrange(numItem1):
            for j in xrange(numItem2-1):
                (lx, ly) = histTableDict[i][j]
                count_total_list[i] += sum(ly)

        for j in xrange(numItem2-1):


            fpout.write("%-11s"%(typeList[j]))
            for i in xrange(numItem1):
                (lx, ly) = histTableDict[i][j]
                count_1TM = ly[1]
                count_2TM = ly[2]
                count_more = sum(ly[3:])
                if not isNorm:
                    fpout.write(" %6d %6d %6d"%(count_1TM, count_2TM, count_more))
                else:
                    freq_1TM = myfunc.FloatDivision(count_1TM, count_total_list[i])
                    freq_2TM = myfunc.FloatDivision(count_2TM, count_total_list[i])
                    freq_more = myfunc.FloatDivision(count_more, count_total_list[i])
                    fpout.write(" %6.3f %6.3f %6.3f"%(freq_1TM, freq_2TM, freq_more))
            fpout.write("\n")
        fpout.close()
        return 0
    except IOError:
        print >> sys.stderr, "Failed to write to file %s"%outfile
        return 1
#}}}
def WriteUnmappedTMPositionHistogram(histTableDict, outfile):#{{{
    try:
        fpout = open(outfile, "w")

        if outfile.find("3bin.") != -1:
            is3bin = True #3bin: N-terminal, Internal, C-terminal
        else:
            is3bin = False

        numItem = len(histTableDict)
        for i in xrange(numItem):
            ss = "%d"%(i)
            if i == 0:
                fpout.write("#%-7s %7s"%(ss+"X", ss+"Y"))
            else:
                fpout.write(" %7s %7s"%(ss+"X", ss+"Y"))
        fpout.write("\n")

        binsList = []
        for i in xrange(numItem):
            binsList.append(len(histTableDict[i][0]))
        maxBin = max(binsList)

        for i in xrange(maxBin):
            for j in xrange(numItem):
                (hist, edges) = histTableDict[j]
                total = sum(hist)
                numbin = len(hist)
                if i < numbin:
                    if not is3bin:
                        xstr = "%7.3f"%((edges[i]+edges[i+1])/2.0)
                    else:
                        if i==0:
                            xstr = "N-terminus"
                        elif i==1:
                            xstr = "Internal"
                        elif i==2:
                            xstr = "C-terminus"
                        else:
                            xstr = "?"
                    fpout.write(" %s %7.3f"%(xstr, hist[i]/float(total)))
                else:
                    fpout.write(" %7s %7s"%("?", "?"))
            fpout.write("\n")
        fpout.close()
        return 0
    except IOError:
        print >> sys.stderr, "Failed to write to %s in function %s"%(outfile,
                sys._getframe().f_code.co_name)
        return 1
#}}}
def WriteDiffNumTM(data, outfile):#{{{
    try:
        fpout = open(outfile, "w")
        for tup in data:
            fpout.write("%2d %2d\n"%(tup[0], tup[1]))
        fpout.close()
    except IOError:
        print >> sys.stderr, "Failed to write to file %s"%outfile
        return 1
#}}}
def WriteDiffNumTMHistogram(histSameDirection, histInverted, outfile):#{{{
    try:
        fpout = open(outfile, "w")
        keys = list(set(histSameDirection.keys()+histInverted.keys()))
        fpout.write("%-10s %14s %14s\n"%("#DiffNumTM", "SameDireciton", "DiffDirection"))
        for key in keys:
            fpout.write("%-10d"%(key))
            try:
                fpout.write("% 14d"%(histSameDirection[key]))
            except KeyError:
                fpout.write("% 14s"%("?"))
            try:
                fpout.write("% 14d"%(histInverted[key]))
            except KeyError:
                fpout.write("% 14s"%("?"))
            fpout.write("\n")
        fpout.close()
    except IOError:
        print >> sys.stderr, "Failed to write to file %s"%outfile
        return 1
#}}}


def WriteTable2D_1(freq, subsum, classList, MAX_NUMTM, outfile):#{{{
    #works with dataTableNumTMDistribution
    try:
        fpout = open(outfile, "wb")
    except IOError:
        print >> sys.stderr, "Failed to write to file %s"%outfile
        return 1

    numClass = len(classList)
    fpout.write("%4s %5s"%("#Idx","numTM"))
    for cls in classList:
        fpout.write(" %9s"%(cls))
    fpout.write("\n")
    for i in range(1, MAX_NUMTM):
        fpout.write("%-4d %5d"%(i-1,i))
        for j in xrange(numClass):
            fpout.write(" %9.3f"%(myfunc.FloatDivision(freq[j][i], subsum[j])*100))
        fpout.write("\n")

    fpout.write("%-10s"%"#Total")
    for j in xrange(numClass):
        fpout.write("%10d"%subsum[j])
    fpout.write("\n")

    fpout.close()
    return 0
#}}}
def WriteDIFFPair():#{{{
    fpout = open(g_params['outDIFFPairFile'], "w")
    for td in g_params['DIFFPairList']:
        pfamid_union = list(set(td[3]) | set(td[4]))
        if len(pfamid_union) == 1:
            famtype = 'singlefam'
        else:
            famtype = 'multifam'

        fpout.write("%-10s %6.1f %s" %(td[1], td[0], td[5]))
        print >> fpout, ", ", famtype, ", ", td[3]
        fpout.write("%-10s %6.1f %s" %(td[2], td[0], td[6]))
        print >> fpout, ", ", famtype, ", ", td[4]
        fpout.write("\n")
    fpout.close()
    print "diffpair output to %s" %g_params['outDIFFPairFile']
#}}}
def WriteCountPairInFam(seqIDTGroupAll):#{{{
    numSeqIDTGroupAll = len(seqIDTGroupAll)/2
    fpout = open(g_params['outCountPairInFam'], "w")
    fpout.write("#%-10s %8s\n"%("seqIDT", "CountFam")) 
    for i in range(numSeqIDTGroupAll):
        fpout.write("%-11s %8d\n" % ("%d-%d"%(seqIDTGroupAll[2*i],
            seqIDTGroupAll[2*i+1]), len(g_params['countPairInFam'][i]))) 

    for i in range(numSeqIDTGroupAll):
        fpout.write("\n") 
        fpout.write("%-11s %8d\n" % ("%d-%d"%(seqIDTGroupAll[2*i],
            seqIDTGroupAll[2*i+1]), len(g_params['countPairInFam'][i]))) 
        d = g_params['countPairInFam'][i]
        sortedTupList = sorted(d.items(), key=itemgetter(1), reverse=True)
        for j in range(len(sortedTupList)):
            fpout.write("%-4d %-10s %6d\n"%(j+1, sortedTupList[j][0],
                sortedTupList[j][1]))
    fpout.close()
    print "countPairInFam file output to %s" %g_params['outCountPairInFam']
#}}}
def OutputPairInfoFile(dataTableCmpClass, cmpClassList, rootname, addname, outpath, isCmpdup):
    if "pairinfo" in dataTableCmpClass:
        newAddname = addname
        if isCmpdup:
            newAddname = ".cmpdup%s"%(newAddname)
        for i in range(len(cmpClassList)):
            outfilepairinfo = (outpath + os.sep + rootname + '_' + newAddname +
                    ".%s.pairinfo.txt"%cmpClassList[i])
            WritePairInfo(dataTableCmpClass['pairinfo'][i], outfilepairinfo)
            if cmpClassList[i] == "INV":
                outfilepairinfowithkrbias =  (outpath + os.sep + rootname + '_'
                        + newAddname + ".%s.pairinfo.withkrbias.txt"%cmpClassList[i])
                cmd = "%s/pairinfoAppendKRbias.sh %s -o %s"%(binpath,
                        outfilepairinfo, outfilepairinfowithkrbias)
#                 if os.path.exists(outfilepairinfowithkrbias):
#                     cmd = "%s"
def main(g_params):#{{{
    argv = sys.argv
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    infile = ""
    outpath = "./"
    isQuiet=False
    diffseqidtgroup = "0"
    tableinfoFile = ""
    rltyinfoFile = ""
    seqDefFile = ""
    seq2famMapfile = ""

    gomapfile = "%s/wk/MPTopo/pfamAna_refpro/pfammap_from_uniprot/Pfam-A-full.seqfrompfamfasta.percentTMpro_scampi.perTM75_nseq20.Family.nr100.filter.fragmented.uniq.pfam.goinfowithancestor.txt"%(DATADIR3)
    gotermfile = "%s/wk/MPTopo/pfamAna_refpro/GO_analysis/GO_term.txt"%(DATADIR3)
    

    #NCInterDefList = [0, 1, 2, 3, 4]
    #th_rlty_list = [0.0, 50,0, 80.0, 90.0]
    #seqidtmode_list = ["high", "low", "all"]

    NCInterDefList = [5]
    th_rlty_list = [0.0]
    #seqidtmode_list = ["all"]
    seqidtmode_list = ["high", "low", "all"]

# pfamdef use the file  /data3/data/pfam/pfam26.0/Pfam-A.clans.tsv
#    seqDefFile = '/data3/wk/MPTopo/pfamAna/pfam2-selTM-giid-refseqid-pfamid-description.txt'
    #seqDefFile = '/data3/wk/MPTopo/pfamAna_uniref/result_pfamscan/Pfam-A-full.perTM75_nseq20.seqdeffile'
#    seq2famMapfile = '/data3/wk/MPTopo/pfamAna/pfamfullseq.selTM_uniq.idmap2'
    #seq2famMapfile = '/data3/wk/MPTopo/pfamAna_uniref/result_pfamscan/uniref90_refpro_20120521.random.seqid2pfamid'
    seqid2pfamidMapFile = ""

    dgscoreFile = ""
    topoalnFile = ""
    seqalnFile = ""
    signalpFile = ""
    dupFile = ""
    selectIDListFile = ""
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
            elif argv[i] in [ "-outpath", "--outpath"]:
                (outpath, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in [ "-signalp", "--signalp"]:
                (signalpFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in [ "-dup", "--dup", "-dupfile", "--dupfile"]:
                (dupFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in [ "-rmsp", "--rmsp"]:
                g_params['isRemoveSignalP']  = True; i+=1
            elif argv[i] in [ "-rmdup", "--rmdup"]:
                g_params['isRemoveDup']  = True; i+=1
            elif argv[i] in ["-seq2fammap", "--seq2fammap"]:
                (seq2famMapfile, i) = myfunc.my_getopt_str(argv, i)
            elif (sys.argv[i] in ["-printdiffseq", "--printdiffseq"]):
                g_params['isPrintDIFFPair'] = True
                i += 1
            elif (argv[i] in ["-printcountpair", "--printcountpair"]):
                g_params['isPrintCountPairInFam'] = True
                i += 1
            elif (argv[i] in ["-diffseqidtgroup", "--diffseqidtgroup"]):
                (diffseqidtgroup, i) = myfunc.my_getopt_str(argv,i)
            elif argv[i] in ["-type", "--type"]:
                (g_params['selecttype'], i) = myfunc.my_getopt_str(argv,i)
            elif argv[i] in ["-mp", "--mp"]:
                g_params['pairwise_comparison_method'], i = myfunc.my_getopt_int(argv, i)
            elif argv[i] in ["-seqidttype", "--seqidttype"]:
                g_params['seqidttype'], i = myfunc.my_getopt_int(argv,i)
            elif argv[i] in ["-tableinfo", "--tableinfo"]:
                tableinfoFile, i = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-rltyinfo", "--rltyinfo"]:
                rltyinfoFile, i = myfunc.my_getopt_str(argv,i)
            elif argv[i] in ["-seqdef", "--seqdef"]:
                seqDefFile, i = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-seq2fammap", "--seq2fammap"]:
                seq2famMapfile, i = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-gap",  "--gap"]:
                g_params['minGapFraction'], i  = myfunc.my_getopt_float(argv, i)
            elif argv[i] == "-ps" or argv[i] == "--ps":
                g_params['minRLTY'], i  = myfunc.my_getopt_float(argv, i)
            elif argv[i] in [ "-seqidtmode", "--seqidtmode"]:
                tmpstr, i = myfunc.my_getopt_str(argv, i)
                seqidtmode_list = tmpstr.split()
            elif argv[i] in [ "-ncinterdef", "--ncinterdef"]:
                tmpstr, i = myfunc.my_getopt_str(argv, i)
                NCInterDefList = [int(x) for x in tmpstr.split()]
            elif argv[i] in [ "-thrltylist", "--thrltylist"]:
                tmpstr, i = myfunc.my_getopt_str(argv, i)
                th_rlty_list = [float(x) for x in tmpstr.split()]
            elif argv[i] in ["-dg", "--dg"]:
                g_params['maxDGvalue'], i  = myfunc.my_getopt_float(argv, i)
            elif argv[i] in ["-min-seqidt", "--min-seqidt"]:
                g_params['minSeqIDT'], i  = myfunc.my_getopt_float(argv, i)
            elif argv[i] in ["-max-seqidt", "--max-seqidt"]:
                g_params['maxSeqIDT'], i  = myfunc.my_getopt_float(argv, i)
            elif argv[i] in ["-evodist", "--evodist"]:
                g_params['isEvodist'] = True
                i += 1
            elif argv[i] in ["-dgscore", "--dgscore"]:
                dgscoreFile,i  = myfunc.my_getopt_str(argv,i)
            elif argv[i] in ["-heatmap", "--heatmap"]:
                g_params['numTMHeatMapMode'],i  = myfunc.my_getopt_str(argv,i)
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
            elif argv[i] in ["-topoaln", "--topoaln"]:
                topoalnFile, i  = myfunc.my_getopt_str(argv,i )
            elif argv[i] in ["-selidlist", "--selidlist"]:
                selectIDListFile, i  = myfunc.my_getopt_str(argv,i )
            elif argv[i] in ["-seqid2pfamid", "--seqid2pfamid"]:
                seqid2pfamidMapFile, i  = myfunc.my_getopt_str(argv,i )
            elif argv[i] in ["-seqaln", "--seqlan"]:
                seqalnFile, i  = myfunc.my_getopt_str(argv,i )
            elif argv[i] in ["-print-rlty-cmpclass", "--print-rlty-cmpclass"]:
                g_params['isPrintFileRltyCmpclass'] = True
                i += 1
            elif argv[i] in ["-print-rlty-helixcmpclass", "--print-rlty-helixcmpclass"]:
                g_params['isPrintFileRltyHelixCmpclass'] = True
                i += 1
            elif argv[i] in [ "-testlevel", "--testlevel"]:
                g_params['testlevel'], i  = myfunc.my_getopt_int(argv,i )
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

    selectIDListSet = []
    if selectIDListFile != "":
        selectIDSet = set(myfunc.ReadIDList(selectIDListFile))

    # try to obtain Pfam family tagPfamType
    tagPfamType = ""
    if infile.find(".Family.") != -1:
        tagPfamType = ".Family"
    elif infile.find(".Domain.") != -1:
        tagPfamType = ".Domain"
    elif infile.find(".Repeat.") != -1:
        tagPfamType = ".Repeat"
    elif infile.find(".Repeat.") != -1:
        tagPfamType = ".Repeat"
    elif infile.find(".Motif.") != -1:
        tagPfamType = ".Motif"
    else:
        tagPfamType = ""
    if seqid2pfamidMapFile == "":
        seqid2pfamidMapFile = "%s/wk/MPTopo/pfamAna_refpro/pfammap_from_uniprot/Pfam-A-full.seqfrompfamfasta.percentTMpro_scampi.perTM75_nseq20%s.nr100.filter.fragmented.seqid2pfamid"%(DATADIR3, tagPfamType)

    seqid2pfamidDict = myfunc.ReadFam2SeqidMap(seqid2pfamidMapFile)

    pairalnStat = {}
    if g_params['seqidttype'] != 0:
        if tableinfoFile == "" or not os.path.exists(tableinfoFile):
            print >> sys.stderr, "tableinfoFile must be set when seqidttype is set to 1 or 2"
            print >> sys.stderr, "but seqidttype = %d is set. Exit."%g_params['seqidttype']
            return -1
        pairalnStat = lcmp.ReadPairAlnTableInfo(tableinfoFile)

    cmpClassList=["OK","SHIFT","INV","INV_SHIFT","DUP", "SIGNALP", "DIFF1", "DIFF2"]
    cmpClassList_method1 = ["IDT","INV","TM2GAP", "TM2SEQ", "TM2GAP_AND_TM2SEQ"]
    cmpClassList_method3 = ["IDT","INV","TM2GAP", "TM2SEQ", "TM2SP", "Mixed"]
    cmpClassList_mp3_cmpdup = ["IDT","INV","DUP", "TM2GAP", "TM2SEQ", "TM2SP", "Mixed"]

    cmpClassList_mp3_cmpdup_SADI = ["IDT","INV","DUP"]
    for ss1 in ["TM2GAP", "TM2SEQ", "TM2SP", "Mixed"]: # analyze also SAME and DIFF
        for ss2 in ["SAME","DIFF"]:
            cmpClassList_mp3_cmpdup_SADI.append("%s|%s"%(ss1,ss2))

    cmpClassList_cmpsp = []
    for ss1 in cmpClassList_method1:
        for ss2 in ["noSP", "SP2TM", "SP2GAP", "SP2SEQ", "unalignedSP"]:
            cmpClassList_cmpsp.append("%s|%s"%(ss1, ss2))

    cmpClassList_cmpsp_rmUnalignedSP = []
    for ss1 in cmpClassList_method1:
        for ss2 in ["noSP", "SP2TM", "SP2GAP", "SP2SEQ"]:
            cmpClassList_cmpsp_rmUnalignedSP.append("%s|%s"%(ss1, ss2))
    cmpClassList_cmpsp_5 = ["IDT", "INV", "TM2GAP", "TM2SEQ", "SP2TM", "Mixed"]

    if g_params['pairwise_comparison_method'] == 0:
        helixCmpClassList = [0, 1, 9]
    elif g_params['pairwise_comparison_method'] == 1:
        helixCmpClassList = [_TM2TM, _TM2GAP, _TM2SEQ]
    elif g_params['pairwise_comparison_method'] == 3:
        helixCmpClassList = [_TM2TM, _TM2GAP, _TM2SEQ, _TM2SP]

# DIFF1: different topology, but with the same number of TM helices
# DIFF2: different topology, and with different number of TM helices
# used to check whether the different topology is mainly caused by sequence
# alignemnt errors

    diffClassList=["Nterm", "Cterm", "internal"]

    if g_params['pairwise_comparison_method'] == 1:
        cmpClassList = cmpClassList_method1
    elif g_params['pairwise_comparison_method'] == 3:
        cmpClassList = cmpClassList_method3

    numCmpClass = len(cmpClassList)
    numCmpClass_mp3_cmpdup = len(cmpClassList_mp3_cmpdup)
    numCmpClass_mp3_cmpdup_SADI = len(cmpClassList_mp3_cmpdup_SADI)
    numCmpClass_cmpsp = len(cmpClassList_cmpsp)
    numCmpClass_cmpsp_rmUnalignedSP = len(cmpClassList_cmpsp_rmUnalignedSP)
    numCmpClass_cmpsp_5 = len(cmpClassList_cmpsp_5)
    numDiffClass = len(diffClassList)
    numHelixCmpClass = len(helixCmpClassList)
    seqIDTGroupAll = [
            0 ,10,
           10 ,20,
           20 ,30,
           30 ,40,
           40 ,50,
           50 ,60,
           60 ,70,
           70 ,80,
           80 ,90,
           90 ,100 
            ]
    evodistGroupAll = [
            0.0,  0.5,
            0.5,  1.0,
            1.0,  1.5,
            1.5,  2.0,
            2.0,  2.5,
            2.5,  3.0,
            3.0,  3.5,
            3.5,  4.0
            ]
    if g_params['isEvodist']:
        seqIDTGroupAll = evodistGroupAll

    if diffseqidtgroup == "0":
        seqIDTGroupDIFF = SEQIDT_GROUP_ALL
    elif diffseqidtgroup == "1":
        seqIDTGroupDIFF = SEQIDT_GROUP_1
    elif diffseqidtgroup == "2":
        seqIDTGroupDIFF = SEQIDT_GROUP_2
    elif diffseqidtgroup == "3":
        seqIDTGroupDIFF = SEQIDT_GROUP_3
    else:
        seqIDTGroupDIFF = SEQIDT_GROUP_ALL

    numSeqIDTGroupDIFF = len(seqIDTGroupDIFF)/2
    numSeqIDTGroupAll = len(seqIDTGroupAll)/2
    MAX_NUMTM = 20

    dataTableUnaligned = {}
    dataTableCmpClass = {}
    dataTableCmpClass_mp3_cmpdup = {}
    dataTableCmpClass_mp3_cmpdup_SADI = {}
    dataTableCmpClass_cmpsp = {}
    dataTableCmpClass_cmpsp_rmUnalignedSP = {}
    dataTableCmpClass_cmpsp_5 = {}
    dataTableNCTermInter = {}
    dataTableNumTMDistribution = {}
    dataTableHelixCmpClass = {} #for counting helix comparison, state 0 1 2
    dataTableHelixDGScore = {}
    dataTableNumTMHeatMap = {}
    dataTableUnmappedTMPosition = {} #0:position calculation using unaligned sequence
                                      #1:position calculation using aligned sequence
    dataTableCmpClassPSBin = {}
    dataTableHelixCmpClassPSBin = {}

    dataTableTM2SEQ = {}
    dataTableTM2TM = {} 
    dataTableTM2SP = {}
    dataTableTM2GAP = {}

    dataTableTM2GAP_add_terminal_numTM = {} # calculate the fraction of different topology 
                                            # differed by adding/removing 1 TM
                                            # at the terminal

    # numGroup = 4 (seqidt 0-20, 20-30, 30-40, 40-100)
    # numClass = 4 (diff_by_term_1, diff_by_term_larger_than_1, diff_by_inter, other)
    numGroup_TM2GAP_add_term = len(SEQIDT_GROUP_5)/2
    numClass_TM2GAP_add_term = 4
    InitTableTM2GAP_add_terminal_numTM(dataTableTM2GAP_add_terminal_numTM,
            numGroup_TM2GAP_add_term, numClass_TM2GAP_add_term)

    InitTableTM2SEQ(dataTableTM2SEQ)
    InitTableTM2TM(dataTableTM2TM)
    InitTableTM2GAP(dataTableTM2GAP)
    InitTableTM2SP(dataTableTM2SP)

    InitTableUnaligned(dataTableUnaligned)
    InitTableCmpClass(dataTableCmpClass, numSeqIDTGroupAll, numCmpClass)

    InitTableCmpClass(dataTableCmpClass_mp3_cmpdup, numSeqIDTGroupAll, numCmpClass_mp3_cmpdup)
    InitTableCmpClass(dataTableCmpClass_mp3_cmpdup_SADI, numSeqIDTGroupAll, numCmpClass_mp3_cmpdup_SADI)

    InitTableCmpClass(dataTableCmpClass_cmpsp, numSeqIDTGroupAll, numCmpClass_cmpsp)
    InitTableCmpClass(dataTableCmpClass_cmpsp_rmUnalignedSP, numSeqIDTGroupAll, numCmpClass_cmpsp_rmUnalignedSP)
    InitTableCmpClass(dataTableCmpClass_cmpsp_5, numSeqIDTGroupAll, numCmpClass_cmpsp_5)
    InitTableNCTermInter(dataTableNCTermInter, numSeqIDTGroupDIFF, numDiffClass)
    InitTableNumTMDistribution(dataTableNumTMDistribution, numSeqIDTGroupDIFF, numDiffClass, MAX_NUMTM)
    InitTableHelixCmpClass(dataTableHelixCmpClass, numSeqIDTGroupAll, numHelixCmpClass)
    InitTableHelixDGScore(dataTableHelixDGScore, numSeqIDTGroupAll, numHelixCmpClass)

    classList_TableNumTMHeatMap = ["ALL", "ONLY_DUP", "ALL_VARY", "TM2GAP", "TM2SEQ", "TM2GAP_AND_TM2SEQ"] 
#    "ALL": everything
#    "ALL_VARY": everything except those categories with identity numTM

    InitTableNumTMHeatMap(dataTableNumTMHeatMap, classList_TableNumTMHeatMap, 100)
    InitTableUnmappedTMPosition(dataTableUnmappedTMPosition, NCInterDefList)

    ps_bin_size = 5.0
    psAllGroup = []
    xv = 0
    while xv <= 100.0:
        psAllGroup += [xv - ps_bin_size/2, xv + ps_bin_size/2]
        xv += ps_bin_size
    numPSGroup = len(psAllGroup)/2
    InitTableCmpClass_ps_bin(dataTableCmpClassPSBin, numPSGroup, numCmpClass)
    InitTableCmpClass_ps_bin(dataTableHelixCmpClassPSBin, numPSGroup, numHelixCmpClass)

    rootname=os.path.basename(os.path.splitext(infile)[0])
    os.system("mkdir -p %s"%outpath)

    #print "binpath=%s"%binpath
# Read in rlty file
    rltyDict = {}
    if rltyinfoFile != "" and os.path.exists(rltyinfoFile):
        rltyDict = lcmp.ReadRLTYInfo(rltyinfoFile)
        #rltyDict = ReadInRLTYInfo_old(rltyinfoFile)
    if rltyDict != {}:
        g_params['isRLTYSupplied'] = True

    seqInfoDict = {}
    if seqDefFile != "" and os.path.exists(seqDefFile):
        if seqDefFile.find('idwithanno') != -1:
            seqInfoDict = ReadSeqDefInfo_idwithanno(seqDefFile)
        elif seqDefFile.find('giid') != -1 or seqDefFile.find('refseq') != -1:
            seqInfoDict = ReadSeqDefInfo_refseq(seqDefFile)
        elif seqDefFile.find('uniref') != -1 or seqDefFile.find('uniprot') != -1:
            seqInfoDict = ReadSeqDefInfo_uniref(seqDefFile)
        else:
            print >> sys.stderr, "un-recognized format for seqDefFile ", seqDefFile

    seq2famDict = {}
    if seq2famMapfile != "" and os.path.exists(seq2famMapfile):
        seq2famDict = ReadIDMap2(seq2famMapfile)

    dgScoreDict = {}
    if dgscoreFile != "":
        dgScoreDict = lcmp.ReadDGScore(dgscoreFile)
        #dgScoreDict = ReadDGScore_old(dgscoreFile)

    topoalnDict = {}
    if topoalnFile != "":
        topoalnDict = ReadTopoAln(topoalnFile)

    seqalnDict = {}
    if seqalnFile != "":
        seqalnDict = ReadTopoAln(seqalnFile)
    signalpDict = {}
    if signalpFile != "":
        signalpDict = lcmp.ReadSignalPDict(signalpFile)
    if signalpDict != {}:
        g_params['isSignalPSet'] = True

    dupPairList = []
    dupPairDict = {} #dictoinary storing duplications
    if dupFile != "":
        dupPairList = lcmp.ReadDupPairList(dupFile)
        dupPairDict = lcmp.ReadDupPairDict(dupFile)
    if len(dupPairList) > 0:
        g_params['isDupSet'] = True
    dupPairSet = set(dupPairList)

    numTMDict = {}

    if g_params['isPrintCountPairInFam']:
        for i in range(numSeqIDTGroupAll):
            g_params['countPairInFam'].append({})

    if g_params['isPrintFileRltyCmpclass']:
        OpenFileRltyCmpclass(g_params, outpath, rootname)
    if g_params['isPrintFileRltyHelixCmpclass']:
        OpenFileRltyHelixCmpclass(g_params, outpath, rootname)

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

        AddSeqDefInfo(pairCmpRecordList, seqInfoDict)
        AddIDMap2Info(pairCmpRecordList, seq2famDict)
        AddTableInfo(pairCmpRecordList, pairalnStat)
        AddDGScore(pairCmpRecordList, dgScoreDict)
        AddSignalPInfo(pairCmpRecordList, signalpDict)
        AddDupInfo(pairCmpRecordList, dupPairSet)
        pairCmpRecordList = FilterPairCmpResult(pairCmpRecordList, rltyDict, selectIDListSet)

#         for rd in pairCmpRecordList:
#             if rd['cmpclass'] == "DIFF":
#                 print rd
        if len(pairCmpRecordList) > 0: 
#            (status, cntTotalOutputRecord) = lcmp.WritePairCmpRecord(pairCmpRecordList, cntTotalOutputRecord, sys.stdout)

            AnaPairCmpTM2GAP_add_term(pairCmpRecordList,
                    dataTableTM2GAP_add_terminal_numTM, SEQIDT_GROUP_5)

            isCmpSP = False
            isRmUnalignedSP = False
            isAna5 = False
            isCmpDup = False
            isAnaSADI = False
            AnaPairCmpResultCmpClass(pairCmpRecordList, dataTableCmpClass,
                    cmpClassList, seqIDTGroupAll, topoalnDict, 
                    isCmpSP, isCmpDup, isRmUnalignedSP, isAna5, isAnaSADI)

            if g_params['pairwise_comparison_method'] == 3:
                isCmpSP = False
                isRmUnalignedSP = False
                isAna5 = False
                isCmpDup = True
                isAnaSADI = False
                AnaPairCmpResultCmpClass(pairCmpRecordList, 
                        dataTableCmpClass_mp3_cmpdup, cmpClassList_mp3_cmpdup, 
                        seqIDTGroupAll, topoalnDict,
                        isCmpSP, isCmpDup, isRmUnalignedSP, isAna5, isAnaSADI)

                isAnaSADI = True
                AnaPairCmpResultCmpClass(pairCmpRecordList, 
                        dataTableCmpClass_mp3_cmpdup_SADI, cmpClassList_mp3_cmpdup_SADI, 
                        seqIDTGroupAll, topoalnDict,
                        isCmpSP, isCmpDup, isRmUnalignedSP, isAna5, isAnaSADI)

            if g_params['pairwise_comparison_method'] != 3:
                isCmpSP = True
                isRmUnalignedSP = False
                isAna5 = False
                isCmpDup = False
                AnaPairCmpResultCmpClass(pairCmpRecordList,
                        dataTableCmpClass_cmpsp, cmpClassList_cmpsp,
                        seqIDTGroupAll, topoalnDict, 
                        isCmpSP, isCmpDup, isRmUnalignedSP, isAna5, isAnaSADI)

                isCmpSP = True
                isRmUnalignedSP = True
                isAna5 = False
                isCmpDup = False
                AnaPairCmpResultCmpClass(pairCmpRecordList,
                        dataTableCmpClass_cmpsp_rmUnalignedSP, cmpClassList_cmpsp_rmUnalignedSP, 
                        seqIDTGroupAll, topoalnDict, 
                        isCmpSP, isCmpDup, isRmUnalignedSP, isAna5, isAnaSADI)

                isCmpSP = True
                isRmUnalignedSP = True
                isAna5 = True
                isCmpDup = False
                AnaPairCmpResultCmpClass(pairCmpRecordList,
                        dataTableCmpClass_cmpsp_5, cmpClassList_cmpsp_5, 
                        seqIDTGroupAll, topoalnDict, 
                        isCmpSP, isCmpDup, isRmUnalignedSP, isAna5, isAnaSADI)

            AnaHelixPairCmpResultCmpClass(pairCmpRecordList,
                    dataTableHelixCmpClass,  helixCmpClassList,
                    seqIDTGroupAll)
            for rd in pairCmpRecordList:
                numTMDict[rd['id1']] = rd['numTM1']
                numTMDict[rd['id2']] = rd['numTM2']

            if g_params['minRLTY'] == 0.0 and g_params['maxRLTY'] == 100.0:
                AnaPiarCmpResultCmpClass_ps_bin(pairCmpRecordList,
                        dataTableCmpClassPSBin, cmpClassList, psAllGroup,
                        rltyDict)
                AnaHelixPairCmpResultCmpClass_ps_bin(pairCmpRecordList,
                        dataTableHelixCmpClassPSBin,  helixCmpClassList, psAllGroup, rltyDict)

            AnaHelixPairCmpResultDGScore(pairCmpRecordList,
                    dataTableHelixDGScore,  helixCmpClassList,
                    seqIDTGroupAll)
            AnaPairCmpResultNCTermInter(pairCmpRecordList,
                    dataTableNCTermInter, diffClassList, seqIDTGroupDIFF)
            AnaPairCmpResultNumTMDistribution(pairCmpRecordList,
                    dataTableNumTMDistribution, diffClassList, seqIDTGroupDIFF,
                    MAX_NUMTM)
            AnaPairCmpResultNumTMHeatMap(pairCmpRecordList,
                    dataTableNumTMHeatMap, classList_TableNumTMHeatMap)
            AnaTM2SEQSegment(pairCmpRecordList,
                    dataTableTM2SEQ, topoalnDict, seqalnDict)
            AnaTM2TMSegment(pairCmpRecordList,
                    dataTableTM2TM, topoalnDict, seqalnDict)
            AnaTM2GAPSegment(pairCmpRecordList,
                    dataTableTM2GAP, topoalnDict, seqalnDict)
            AnaTM2SPSegment(pairCmpRecordList,
                    dataTableTM2SP, topoalnDict, seqalnDict)
            AnaUnalignedTerminal(pairCmpRecordList, dataTableUnaligned)
            if (g_params['pairwise_comparison_method'] in [1,3]
                    and topoalnDict != {}):
                for item in NCInterDefList:
                    AnaPairCmpResultUnmappedTMPosition(pairCmpRecordList,
                            dataTableUnmappedTMPosition[item], item,
                            topoalnDict, rltyDict, dupPairDict)
            cntTotalReadInRecord += len(pairCmpRecordList)
        if isEOFreached == True:
            break

    fpin.close()

    #return 0

#     if g_params['numTMHeatMapMode'] == "full":
#         FillSymmetricDataTableNumTMHeatMap(dataTableNumTMHeatMap, classList_TableNumTMHeatMap)

    print "cntTotalReadInRecord =", cntTotalReadInRecord

    #print "cntTotalOutputRecord =", cntTotalOutputRecord

#   write out the result
#=================================================
    if g_params['pairwise_comparison_method'] == 0:
        addname = 'type_%s_dg%.1f_gap%.1f_ps%.1f'%(g_params['selecttype'],
                g_params['maxDGvalue'], g_params['minGapFraction'],
                g_params['minRLTY'])
    else:
        addname = ""
    if g_params['isRemoveSignalP']:
        addname += ".RMSP"
    if g_params['alignrange'] != 'all':
        addname += ".%s"%(g_params['alignrange'])
    if g_params['isRemoveDup']:
        addname += ".RMDUP"
    if g_params['isEvodist']:
        addname += ".evodist"
    outFileCmpClass = outpath + os.sep + rootname + '_' + addname + "-cmpclass.table.txt"; 
    outFileCmpClass_mp3_cmpdup = outpath + os.sep + rootname + '_' + addname + "-cmpclass.cmpdup.table.txt"; 
    outFileCmpClass_mp3_cmpdup_SADI = outpath + os.sep + rootname + '_' + addname + "-cmpclass.cmpdup.SADI.table.txt"; 
    outFileCmpClass_cmpsp = outpath + os.sep + rootname + '_' + addname + "-cmpclass.cmpsp.table.txt"; 
    outFileCmpClass_cmpsp_rmUnalignedSP = outpath + os.sep + rootname + '_' + addname + "-cmpclass.cmpsp.rmUnAlnSP.table.txt"; 
    outFileCmpClass_cmpsp_5 = outpath + os.sep + rootname + '_' + addname + "-cmpclass.cmpsp.ana5.table.txt"; 
    outFileCmpClassWithoutIDT =  (outpath + os.sep + rootname + '_' + addname + "-cmpclass.table.withoutIDT.txt")
    outFileNCtermInter1 = (outpath + os.sep + rootname +  '_' + addname + 
            "-NCtermInter.table1.c%s.txt"%diffseqidtgroup); 
    outFileNCtermInter2 = (outpath + os.sep + rootname +   '_' + addname +  
            "-NCtermInter.table2.c%s.txt"%diffseqidtgroup); 

    outFileTM2GAP_addterm = outpath + os.sep + rootname + '_' + addname + ".TM2GAP_Add_Term.table.txt"; 


#----------------------------------------------
# TM2GAP add term analysis
    classList_TM2GAP_add = ['TERM_1', 'TERM_2','Inter','Other']
    if (WriteTable2D(dataTableTM2GAP_add_terminal_numTM['freq'], 
            dataTableTM2GAP_add_terminal_numTM['subsum'], 
            classList_TM2GAP_add, SEQIDT_GROUP_5, outFileTM2GAP_addterm) == 0):
        print "outFileTM2GAP_addterm %s output" %(outFileTM2GAP_addterm)
#----------------------------------------------
# TM2SEQ seqidt analysis
#     print "dataTableTM2SEQ:", dataTableTM2SEQ
    for item in ["", "_noshift", "_shift"]:
        outfileTM2SEQSegment = (outpath + os.sep + rootname + "-" + addname +
                "-TM2SEQ-Segment%s.txt"%(item))
        WriteTM2SEQDGScoreSegment(dataTableTM2SEQ['data%s'%item], outfileTM2SEQSegment)
        print "outfileTM2SEQSegment %s output"%(outfileTM2SEQSegment)

        # append DG score
        outfileTM2SEQSegment_withdgscore = (outpath + os.sep + rootname + "-" + addname +
                "-TM2SEQ-Segment%s.withdgscore.txt"%(item))
        cmd = ["%s/TM2SEQSegmentAppendDGvalue.sh"%(rundir), outfileTM2SEQSegment, "-o", outfileTM2SEQSegment_withdgscore]
        try:
            subprocess.check_output(cmd)
        except subprocess.CalledProcessError, e:
            print e

        outfileTM2SEQSegment2 = (outpath + os.sep + rootname + "-" + addname +
                "-TM2SEQ-Segment%s.2.txt"%item)
        WriteTM2SEQDGScoreSegment2(dataTableTM2SEQ['data%s'%item], outfileTM2SEQSegment2)
        print "outfileTM2SEQSegment2 %s output"%(outfileTM2SEQSegment2)
#----------------------------------------------
#----------------------------------------------
# TM2GAP seqidt analysis
    for item in [""]:
        outfileTM2GAPSegment = (outpath + os.sep + rootname + "-" + addname +
                "-TM2GAP-Segment%s.txt"%(item))
        WriteTM2GAPDGScoreSegment(dataTableTM2GAP['data%s'%item], outfileTM2GAPSegment)
        print "outfileTM2GAPSegment %s output"%(outfileTM2GAPSegment)

        # append DG score
        outfileTM2GAPSegment_withdgscore = (outpath + os.sep + rootname + "-" + addname +
                "-TM2GAP-Segment%s.withdgscore.txt"%(item))
        cmd = ["%s/TM2SEQSegmentAppendDGvalue.sh"%(rundir), outfileTM2GAPSegment, "-o", outfileTM2GAPSegment_withdgscore]
        try:
            subprocess.check_output(cmd)
        except subprocess.CalledProcessError, e:
            print e

        outfileTM2GAPSegment2 = (outpath + os.sep + rootname + "-" + addname +
                "-TM2GAP-Segment%s.2.txt"%item)
        WriteTM2SEQDGScoreSegment2(dataTableTM2GAP['data%s'%item], outfileTM2GAPSegment2)
        print "outfileTM2GAPSegment2 %s output"%(outfileTM2GAPSegment2)
#----------------------------------------------
#----------------------------------------------
# TM2TM seqidt analysis
    for item in [""]:
        outfileTM2TMSegment = (outpath + os.sep + rootname + "-" + addname +
                "-TM2TM-Segment%s.txt"%(item))
        WriteTM2SEQDGScoreSegment(dataTableTM2TM['data%s'%item], outfileTM2TMSegment)
        print "outfileTM2TMSegment %s output"%(outfileTM2TMSegment)

        # append DG score
        outfileTM2TMSegment_withdgscore = (outpath + os.sep + rootname + "-" + addname +
                "-TM2TM-Segment%s.withdgscore.txt"%(item))
        cmd = ["%s/TM2SEQSegmentAppendDGvalue.sh"%(rundir), outfileTM2TMSegment, "-o", outfileTM2TMSegment_withdgscore]
        try:
            subprocess.check_output(cmd)
        except subprocess.CalledProcessError, e:
            print e

        outfileTM2TMSegment2 = (outpath + os.sep + rootname + "-" + addname +
                "-TM2TM-Segment%s.2.txt"%item)
        WriteTM2SEQDGScoreSegment2(dataTableTM2TM['data%s'%item], outfileTM2TMSegment2)
        print "outfileTM2TMSegment2 %s output"%(outfileTM2TMSegment2)
#----------------------------------------------
#----------------------------------------------
# TM2TM seqidt analysis
    for item in [""]:
        outfileTM2SPSegment = (outpath + os.sep + rootname + "-" + addname +
                "-TM2SP-Segment%s.txt"%(item))
        WriteTM2SEQDGScoreSegment(dataTableTM2SP['data%s'%item], outfileTM2SPSegment)
        print "outfileTM2SPSegment %s output"%(outfileTM2SPSegment)

        # append DG score
        outfileTM2SPSegment_withdgscore = (outpath + os.sep + rootname + "-" + addname +
                "-TM2SP-Segment%s.withdgscore.txt"%(item))
        cmd = ["%s/TM2SEQSegmentAppendDGvalue.sh"%(rundir), outfileTM2SPSegment, "-o", outfileTM2SPSegment_withdgscore]
        try:
            subprocess.check_output(cmd)
        except subprocess.CalledProcessError, e:
            print e

        outfileTM2SPSegment2 = (outpath + os.sep + rootname + "-" + addname +
                "-TM2SP-Segment%s.2.txt"%item)
        WriteTM2SEQDGScoreSegment2(dataTableTM2SP['data%s'%item], outfileTM2SPSegment2)
        print "outfileTM2SPSegment2 %s output"%(outfileTM2SPSegment2)
#----------------------------------------------

# debug 2015-04-20
    if g_params['testlevel'] == 2:
        sys.exit(0)
#========================

    if g_params['isPrintDIFFPair'] and g_params['pairwise_comparison_method'] == 0:
        g_params['outDIFFPairFile'] = (outpath + os.sep + rootname + '_' +
                addname + '.diffpair.txt')
        g_params['DIFFPairList'] = sorted(g_params['DIFFPairList'], key=lambda
                tup:tup[0], reverse=True)
        WriteDIFFPair()

    if (g_params['isPrintCountPairInFam'] and
            g_params['pairwise_comparison_method'] in [0,3]):
        g_params['outCountPairInFam'] = (outpath + os.sep + rootname + '_'
                + addname + '.pairInFam.txt')
        WriteCountPairInFam(seqIDTGroupAll)

    # write cmpclass with RLTY bins, added 2012-11-19
    if g_params['minRLTY'] == 0.0 and g_params['maxRLTY'] == 100.0:
        for key in ['min_ps', 'avg_ps', 'max_ps']:
            outfile_psbin_cmpclass = (outpath + os.sep + rootname + '_' 
                    + addname + key +'bin.cmpclass.table.txt')
            outfile_psbin_helixcmpclass = (outpath + os.sep + rootname + '_'
                    + addname + key + 'bin.helixcmpclass.table.txt')
            outfile_numpair_cmpclass = (outpath+os.sep+rootname + '_'
                    + addname + key + 'bin.numpair.cmpclass.table.txt')
            outfile_numpair_helixcmpclass = (outpath+os.sep+rootname + '_'
                    + addname + key + 'bin.numpair.helixcmpclass.table.txt')
            if WriteHistogramData(dataTableCmpClassPSBin[key]['freq'],
                    dataTableCmpClassPSBin[key]['subsum'], "RLTY",
                    cmpClassList, psAllGroup, outfile_psbin_cmpclass) == 0:
                print "table file %s output" %(outfile_psbin_cmpclass)
            if WriteHistogramData(dataTableHelixCmpClassPSBin[key]['freq'],
                    dataTableHelixCmpClassPSBin[key]['subsum'], "RLTY",
                    helixCmpClassList, psAllGroup, 
                    outfile_psbin_helixcmpclass) == 0:
                print "table file %s output" %(outfile_psbin_helixcmpclass)
            if WriteNumPair2D(dataTableCmpClassPSBin[key]['freq_by_seqidt'],
                        psAllGroup, SEQIDT_GROUP_4,
                        outfile_numpair_cmpclass) == 0:
                print "table file %s output" %(outfile_numpair_cmpclass)
            if WriteNumPair2D(dataTableHelixCmpClassPSBin[key]['freq_by_seqidt'],
                        psAllGroup, SEQIDT_GROUP_4,
                        outfile_numpair_helixcmpclass) == 0:
                print "table file %s output" %(outfile_numpair_helixcmpclass)

    if dgScoreDict != {}:
        outfile = ""#{{{
        if g_params['isEvodist']:
            addtxt = ".evodist"
        else:
            addtxt = ""
        for item in ["", "_noshift", "_shift"]:
            for i in xrange(numHelixCmpClass):
                outFileHelixDGScore = (outpath + os.sep + rootname + "_" + addname +
                        ".TMMap%s"%(helixCmpClassList[i]) + addtxt + ".dgscore%s.txt"%item )
                WriteHelixDGScore(dataTableHelixDGScore['data%s'%item][i],
                        seqIDTGroupAll, outFileHelixDGScore)
                print "table file %s output" %(outFileHelixDGScore)
                histTableDict = CalHistDGScore(dataTableHelixDGScore['data%s'%item][i],
                        seqIDTGroupAll)
                outFileHelixDGScoreHistogram = (outpath + os.sep + rootname + "_" + addname +
                        ".TMMap%s"%(helixCmpClassList[i]) + addtxt + ".dgscore%s.hist"%item )
                outfile = outFileHelixDGScoreHistogram
                WriteHelixDGScoreHistogram(histTableDict, seqIDTGroupAll,
                        outFileHelixDGScoreHistogram)
                print "table file %s output" %(outFileHelixDGScoreHistogram)
            if os.path.exists(outfile):
                cmd = "%s/plotHelixDGScoreHistMerged_mp1.sh %s" %(binpath, outfile)
                os.system(cmd)
#}}}
# write pairinfo list
    OutputPairInfoFile(dataTableCmpClass, cmpClassList, rootname, addname, outpath, False)
    OutputPairInfoFile(dataTableCmpClass_mp3_cmpdup, cmpClassList_mp3_cmpdup, rootname, addname, outpath, True)  #output also pairinfo list with cmpdup

    if WriteTable2D(dataTableCmpClass['freq'], dataTableCmpClass['subsum'],
            cmpClassList, seqIDTGroupAll, outFileCmpClass) == 0:
        if g_params['isEvodist']:
            xlabel = "Evolutionary distance"
        else:
            xlabel = "Sequence identity"

        if g_params['pairwise_comparison_method'] == 0:
            os.system("%s/plotCmpClass.sh %s -xlabel \"%s\" -outstyle eps  -outpath %s" %(
                binpath, outFileCmpClass, xlabel, outpath) )
            os.system("awk '/^[^#]/{sum1=$4+$5+$6+$7+$8+$9+$10; print $1,$2, $4/sum1*100, $5/sum1*100, $6/sum1*100, $7/sum1*100, $8/sum1*100, $9/sum1*100, $10/sum1*100, int($12*sum1/100+0.5) }' %s > %s" %(outFileCmpClass, outFileCmpClassWithoutIDT))
            os.system("%s/plotCmpClass_noOK.sh %s -xlabel \"%s\" -outstyle eps  -outpath %s" %(
                binpath, outFileCmpClassWithoutIDT, xlabel, outpath) )
            print "table file %s output" %outFileCmpClass
            print "table file %s output" %outFileCmpClassWithoutIDT
        elif  g_params['pairwise_comparison_method'] == 1:
            os.system("%s/plotCmpClass_mp1.sh %s -xlabel \"%s\" -outstyle eps  -outpath %s -plot1" %(
                binpath, outFileCmpClass, xlabel, outpath) )
            os.system("%s/plotCmpClass_line_mp1.sh %s -xlabel \"%s\" -outstyle eps  -outpath %s" %(
                binpath, outFileCmpClass, xlabel, outpath) )
            print "table file %s output" %outFileCmpClass
        elif  g_params['pairwise_comparison_method'] == 3:
            cmd = "%s/plotCmpClass_mp3.sh %s -xlabel \"%s\" -outstyle eps"\
            " -outpath %s -plot1 -multiplot" %( binpath, outFileCmpClass, xlabel, outpath)
            os.system(cmd)
            GOAnalysis_mp3(dataTableCmpClass, cmpClassList, outpath, rootname,
                    addname, gomapfile, gotermfile)

    if g_params['pairwise_comparison_method'] == 3:
        if WriteTable2D(dataTableCmpClass_mp3_cmpdup['freq'],
                dataTableCmpClass_mp3_cmpdup['subsum'],
                cmpClassList_mp3_cmpdup, seqIDTGroupAll,
                outFileCmpClass_mp3_cmpdup) == 0:
            cmd = "%s/plotCmpClass_mp3_cmpdup.sh %s -xlabel \"%s\" -outstyle eps"\
            " -outpath %s -plot1 -multiplot" %(binpath, outFileCmpClass_mp3_cmpdup, xlabel, outpath)
            os.system(cmd)
            GOAnalysis_mp3(dataTableCmpClass_mp3_cmpdup,
                    cmpClassList_mp3_cmpdup, outpath, rootname,
                    addname, gomapfile, gotermfile)

        if WriteTable2D(dataTableCmpClass_mp3_cmpdup_SADI['freq'],
                dataTableCmpClass_mp3_cmpdup_SADI['subsum'],
                cmpClassList_mp3_cmpdup_SADI, seqIDTGroupAll,
                outFileCmpClass_mp3_cmpdup_SADI) == 0:
            cmd = "%s/plotCmpClass_mp3_cmpdup_SADI.sh %s -xlabel \"%s\" -outstyle eps"\
            " -outpath %s -plot1 -multiplot -color" %(binpath, outFileCmpClass_mp3_cmpdup_SADI, xlabel, outpath)
            print cmd
            os.system(cmd)

    if g_params['testlevel'] == 1:
        return 1 # DEBUG 2014-06-23

#================================= BEGIN CMPSP for mp 1 =============
    if g_params['pairwise_comparison_method'] == 1:
        if WriteTable2D(dataTableCmpClass_cmpsp['freq'], dataTableCmpClass_cmpsp['subsum'],
                cmpClassList_cmpsp, seqIDTGroupAll, outFileCmpClass_cmpsp) == 0:
            if g_params['isEvodist']:
                xlabel = "Evolutionary distance"
            else:
                xlabel = "Sequence identity"
            os.system("%s/plotCmpClass_mp1_cmpsp.sh %s -xlabel \"%s\" -outstyle eps  -outpath %s -plot1" %(
                binpath, outFileCmpClass_cmpsp, xlabel, outpath) )

        if WriteTable2D(dataTableCmpClass_cmpsp_rmUnalignedSP['freq'],
                dataTableCmpClass_cmpsp_rmUnalignedSP['subsum'],
                cmpClassList_cmpsp_rmUnalignedSP, seqIDTGroupAll,
                outFileCmpClass_cmpsp_rmUnalignedSP) == 0:
            if g_params['isEvodist']:
                xlabel = "Evolutionary distance"
            else:
                xlabel = "Sequence identity"
            os.system("%s/plotCmpClass_mp1_cmpsp_rmUnalignedSP.sh %s -xlabel \"%s\" -outstyle eps  -outpath %s -plot1" %(
                binpath, outFileCmpClass_cmpsp_rmUnalignedSP, xlabel, outpath) )

#==============CMPSP ana5 BEGIN==============================================
        if WriteTable2D(dataTableCmpClass_cmpsp_5['freq'],
                dataTableCmpClass_cmpsp_5['subsum'],
                cmpClassList_cmpsp_5, seqIDTGroupAll,
                outFileCmpClass_cmpsp_5) == 0:
            if g_params['isEvodist']:
                xlabel = "Evolutionary distance"
            else:
                xlabel = "Sequence identity"
            cmd =  "%s/plotCmpClass_mp1_cmpsp_5.sh %s -xlabel \"%s\""\
                    " -outstyle eps  -outpath %s -plot1" %( binpath,
                            outFileCmpClass_cmpsp_5, xlabel, outpath) 
            os.system(cmd)
            GOAnalysis_mp3(dataTableCmpClass_cmpsp_5, cmpClassList_cmpsp_5,
                    outpath, rootname, addname+".cmpsp.ana5", gomapfile,
                    gotermfile)
#==============CMPSP ana5 END==============================================


    outFileHelixCmpClass = (outpath + os.sep + rootname + '_' + addname +
            "-helixcmpclass.table.txt")
    if WriteTable2D(dataTableHelixCmpClass['freq'], dataTableHelixCmpClass['subsum'],
            helixCmpClassList, seqIDTGroupAll, outFileHelixCmpClass) == 0:
        print "%s output"%outFileHelixCmpClass
        if g_params['isEvodist']:
            xlabel = "Evolutionary distance"
        else:
            xlabel = "Sequence identity"
        if g_params['pairwise_comparison_method'] == 1:
            cmd = "%s/plotHelixCmpClass_mp1.sh %s -xlabel \"%s\" -plot1"\
                    " -multiplot -outpath %s -outstyle eps"%(binpath, 
                    outFileHelixCmpClass, xlabel, outpath)
        elif g_params['pairwise_comparison_method'] == 3:
            cmd = "%s/plotHelixCmpClass_mp3.sh %s -xlabel \"%s\" -plot1"\
                    " -multiplot -outpath %s -outstyle eps"%(binpath, 
                    outFileHelixCmpClass, xlabel, outpath)
        os.system(cmd)

    if g_params['pairwise_comparison_method'] == 0:
        if WriteTable2D(dataTableNCTermInter['freq1'],
                dataTableNCTermInter['subsum1'], diffClassList, seqIDTGroupDIFF,
                outFileNCtermInter1) == 0:
            os.system("%s/plotNCtermInterOccu.sh %s -outstyle eps  -outpath %s" %(
                binpath, outFileNCtermInter1, outpath) )
            print "table file %s output" %outFileNCtermInter1

        if WriteTable2D(dataTableNCTermInter['freq2'],
                dataTableNCTermInter['subsum2'], diffClassList, seqIDTGroupDIFF,
                outFileNCtermInter2) == 0:
            os.system("%s/plotNCtermInterOccu.sh %s -outstyle eps  -outpath %s" %(
                binpath, outFileNCtermInter2, outpath) )
            print "table file %s output" %outFileNCtermInter2

        for i in range(numSeqIDTGroupDIFF):
            seqidtrange = "%d-%d" % (seqIDTGroupDIFF[2*i], seqIDTGroupDIFF[2*i+1])
            outFileNumTMDistri = (outpath + os.sep + rootname +   '_' + addname +  
                    "-NumTMDistri.seqidt%s.txt"%(seqidtrange))
            if WriteTable2D_1(dataTableNumTMDistribution['freq'][i],
                    dataTableNumTMDistribution['subsum'][i], diffClassList,
                    MAX_NUMTM, outFileNumTMDistri)  == 0:
                os.system("%s/plotNumTMDistri.sh %s -outstyle eps  -outpath %s" %(
                    binpath, outFileNumTMDistri, outpath) )
                print "table file %s output"%outFileNumTMDistri
    # Count numTM - numTM heat map
    if g_params['isPrintNumTMHeatMap']:
        for cls in classList_TableNumTMHeatMap: 
            for mode_norm in ["norm_row", "norm_col", "norm_diag"]:
                if mode_norm == "norm_diag":
                    heatmapmode = 'half'
                else:
                    heatmapmode = 'full'
                outFileNumTMHeatMap = (outpath + os.sep + rootname + '_' +
                        addname + ".numTMHeatMap.%s.%s.%s.txt" %
                        (heatmapmode, cls, mode_norm))
                if heatmapmode == 'full':
                    mtx = myfunc.FillSymmetricMatrix(
                            dataTableNumTMHeatMap[cls]['data'],
                            dataTableNumTMHeatMap[cls]['maxNumTM'])
                else:
                    mtx = dataTableNumTMHeatMap[cls]['data']
                if WriteNumTMHeatMap(mtx,
                        dataTableNumTMHeatMap[cls]['maxNumTM'],
                        dataTableNumTMHeatMap[cls]['numPair'], mode_norm,
                        outFileNumTMHeatMap) == 0:
                    print "heatmap %s output"%(outFileNumTMHeatMap) 
                    cmd = "%s/plotNumTMHeatMap.sh %s" %(binpath,
                            outFileNumTMHeatMap)
                    os.system(cmd)

    if g_params['pairwise_comparison_method'] in [1,3] and topoalnDict != {}:
        isNorm = False
        for seqidtmode in seqidtmode_list:
            if seqidtmode == "low":
                str_seqidtmode = "lt%.0f"%(g_params['thHigh_seqidt'])
            elif seqidtmode == "high":
                str_seqidtmode = "ge%.0f"%(g_params['thHigh_seqidt'])
            else:
                str_seqidtmode = seqidtmode
            for item in NCInterDefList:
                for th_rlty in th_rlty_list:
                    for termmode in ["withterm"]:
                        outFileUnmappedTMPosition =\
                        "%s%s%s_%s.UnmappedTMPosition_%d_%s_%s_ps%.0f.txt"%(outpath,
                                os.sep, rootname, addname, item, termmode,
                                str_seqidtmode, th_rlty)
                        WriteUnmappedTMPosition(
                                dataTableUnmappedTMPosition[item]['dataPosition'],
                                outFileUnmappedTMPosition, termmode,
                                seqidtmode, th_rlty)
                        print "table file %s output" %(outFileUnmappedTMPosition)

#----------------------------------------------------------------------------------
#output histogram, first use a more detailed binList
                        binList = [
                                0.0, 1e-9,
                                1e-9, 0.1,
                                0.1, 0.2,
                                0.2, 0.3,
                                0.3, 0.4,
                                0.4, 0.5,
                                0.5, 0.6,
                                0.6, 0.7,
                                0.7, 0.8,
                                0.8, 0.9,
                                0.9, 1.0,
                                1.0, 1.0+1e-9
                                ]
                        histTableDict = CalHistUnmappedTMPosition(
                                dataTableUnmappedTMPosition[item]['dataPosition'],
                                termmode, seqidtmode, th_rlty, binList)
                        outFileUnmappedTMPositionHistogram =\
                        "%s%s%s_%s.UnmappedTMPosition_%d_%s_%s_ps%.0f.hist"%(outpath,
                                os.sep, rootname, addname, item, termmode,
                                str_seqidtmode, th_rlty)
                        WriteUnmappedTMPositionHistogram(histTableDict,
                                outFileUnmappedTMPositionHistogram)
                        print "table file %s output" %(outFileUnmappedTMPositionHistogram)
                        if os.path.exists(outFileUnmappedTMPositionHistogram):
                            cmd = "%s/plotHelixUnmappedTMPosition.sh %s" %(binpath,
                                    outFileUnmappedTMPositionHistogram)
                            os.system(cmd)
#then use 3 bin, Nterm, Cterm, Inside
                        binList = [
                                0.0, 1e-9,
                                1e-9, 1.0,
                                1.0, 1.0+1e-9
                                ]
                        histTableDict = CalHistUnmappedTMPosition(
                                dataTableUnmappedTMPosition[item]['dataPosition'],
                                termmode, seqidtmode, th_rlty, binList)
                        outFileUnmappedTMPositionHistogram =\
                                "%s%s%s_%s.UnmappedTMPosition_%d_%s_%s_ps%.0f.3bin.hist"%(
                                        outpath, os.sep, rootname, addname,
                                        item, termmode, str_seqidtmode, th_rlty)
                        WriteUnmappedTMPositionHistogram(histTableDict,
                                outFileUnmappedTMPositionHistogram)
                        print "table file %s output" %(outFileUnmappedTMPositionHistogram)
                        if os.path.exists(outFileUnmappedTMPositionHistogram):
                            cmd = "%s/plotHelixUnmappedTMPosition.sh %s" %(binpath,
                                    outFileUnmappedTMPositionHistogram)
                            os.system(cmd)
#-------------------------------------------------------------------------------
                        outFileNumContinuous =\
                                "%s%s%s_%s.NumContUnmappedTM_%d_%s_%s_ps%.0f.txt"%(
                                        outpath, os.sep, rootname, addname,
                                        item, termmode, str_seqidtmode, th_rlty)
                        WriteNumContinuousTM(
                                dataTableUnmappedTMPosition[item]['dataNumContinuousTM'],
                                outFileNumContinuous, termmode, seqidtmode,
                                th_rlty, seqid2pfamidDict, numTMDict)
                        print "data file %s output" %(outFileNumContinuous)
#-------------------------------------------------------------------------------
                        # without cmpdup
                        isCmpDup = False
                        histNumContinuousTableDict = CalHistNumContinuousTM(
                                dataTableUnmappedTMPosition[item]['dataNumContinuousTM'],
                                termmode, seqidtmode, th_rlty, isCmpDup)
                        outFileNumContinuousHistogram =\
                                "%s%s%s_%s.NumContUnmappedTM_%d_%s_%s_ps%.0f.hist"%(
                                        outpath, os.sep, rootname, addname,
                                        item, termmode, str_seqidtmode, th_rlty)
                        WriteNumContinuousTMHistogram(histNumContinuousTableDict,
                                outFileNumContinuousHistogram, isNorm)
                        print "table file %s output" %(outFileNumContinuousHistogram)
                        if os.path.exists(outFileNumContinuousHistogram):
                            cmd = "%s/plotHelixNumContUnmappedTM.sh %s" %(binpath,
                                    outFileNumContinuousHistogram)
                            os.system(cmd)

                            cmd = "%s/plotHelixNumContUnmappedTM.withHistogram.sh %s %s" %(binpath,
                                    outFileNumContinuousHistogram,
                                    outFileUnmappedTMPositionHistogram)
                            os.system(cmd)
#-------------------------------------------------------------------------------
                        # with cmpdup, include classificaiton of duplicated
                        # pairs
                        isCmpDup = True
                        histNumContinuousTableDict = CalHistNumContinuousTM(
                                dataTableUnmappedTMPosition[item]['dataNumContinuousTM'],
                                termmode, seqidtmode, th_rlty, isCmpDup)

                        outFileNumContinuousHistogram =\
                                "%s%s%s_%s.NumContUnmappedTM.cmpdup_%d_%s_%s_ps%.0f.hist"%(
                                        outpath, os.sep, rootname, addname,
                                        item, termmode, str_seqidtmode, th_rlty)
                        WriteNumContinuousTMHistogram(histNumContinuousTableDict,
                                outFileNumContinuousHistogram, isNorm)
                        print "table file %s output" %(outFileNumContinuousHistogram)

                        # not normalized
                        outFileNumContinuousHistogram2 =\
                                "%s%s%s_%s.NumContUnmappedTM.cmpdup_%d_%s_%s_ps%.0f.non_norm.hist"%(
                                        outpath, os.sep, rootname, addname,
                                        item, termmode, str_seqidtmode, th_rlty)
                        WriteNumContinuousTMHistogram_grouped(histNumContinuousTableDict,
                                outFileNumContinuousHistogram2, False)
                        # normalized
                        outFileNumContinuousHistogram3 =\
                                "%s%s%s_%s.NumContUnmappedTM.cmpdup_%d_%s_%s_ps%.0f.norm.hist"%(
                                        outpath, os.sep, rootname, addname,
                                        item, termmode, str_seqidtmode, th_rlty)
                        WriteNumContinuousTMHistogram_grouped(histNumContinuousTableDict,
                                outFileNumContinuousHistogram3, True)

                        if os.path.exists(outFileNumContinuousHistogram):
                            cmd = "%s/plotHelixNumContUnmappedTM_cmpdup.sh %s" %(binpath,
                                    outFileNumContinuousHistogram)
                            os.system(cmd)

                            cmd = "%s/plotHelixNumContUnmappedTM_cmpdup.withHistogram.sh %s %s" %(
                                    binpath,
                                    outFileNumContinuousHistogram,
                                    outFileNumContinuousHistogram2)
                            os.system(cmd)

                            cmd = "%s/plotHelixNumContUnmappedTM_cmpdup.withHistogram.sh %s %s" %(
                                    binpath,
                                    outFileNumContinuousHistogram,
                                    outFileNumContinuousHistogram3)
                            os.system(cmd)
#-------------------------------------------------------------------------------
                        for side1 in ["Nterm", "Cterm"]:
                            for side2 in ["Nterm", "Cterm"]:
                                tmpdataTable = dataTableUnmappedTMPosition[
                                        item]['dataDiffNumTM%s_%sStatus'%(side1, side2)]
                                if len(tmpdataTable) > 0:
                                    outFileDiffNumTM =\
                                            "%s%s%s_%s.DiffNumTM%s_%sStatus_%d_%s_%s_ps%.0f.txt"%(
                                                    outpath, os.sep, rootname,
                                                    addname, side1, side2,
                                                    item, termmode,
                                                    str_seqidtmode, th_rlty)
                                    WriteDiffNumTM(tmpdataTable, outFileDiffNumTM)
                                    print "data file %s output" %(outFileDiffNumTM)

                                    liSameDirection = []
                                    liInverted = []
                                    for tup in tmpdataTable:
                                        if tup[1] == 1:
                                            liSameDirection.append(tup[0])
                                        else:
                                            liInverted.append(tup[0])
                                    histSameDirection = CalHist_grouping(liSameDirection)
                                    histInverted = CalHist_grouping(liInverted)
                                    outFileDiffNumTMHistogram =\
                                            "%s%s%s_%s.DiffNumTM%s_%sStatus_%d_%s_%s_ps%.0f.hist"%(
                                                    outpath, os.sep, rootname,
                                                    addname,side1, side2,
                                                    item, termmode,
                                                    str_seqidtmode, th_rlty)
                                    WriteDiffNumTMHistogram(histSameDirection, histInverted,
                                            outFileDiffNumTMHistogram)
                                    print "table file %s output" %(outFileDiffNumTMHistogram)
                                    if os.path.exists(outFileDiffNumTMHistogram):
                                        cmd = "%s/plotDiffNumTMNterm.sh %s" %(binpath,
                                                outFileDiffNumTMHistogram)
                                        os.system(cmd)

    # output unaligned terminal analysis
    if dataTableUnaligned['numPair'] > 0:
        dataTableUnaligned['diff_numTM_Nterm_hist'] =\
                CalHist_grouping(dataTableUnaligned['diff_numTM_Nterm'])
        dataTableUnaligned['diff_numTM_Cterm_hist'] =\
                CalHist_grouping(dataTableUnaligned['diff_numTM_Cterm'])
        tmpfilelist = []
        for item in ['Nterm', 'Cterm']:
            outFileUnalignedTerm = "%s%s%s_%s.diffNumTM_Unaligned_%s.hist"%(
                            outpath, os.sep, rootname, addname, item)
            WriteHistogram_grouping(dataTableUnaligned['diff_numTM_%s_hist'%(item)], 
                    outFileUnalignedTerm)
            if os.path.exists(outFileUnalignedTerm):
                cmd = "%s/plotDiffUnalignedNumTM.sh %s"%(binpath, outFileUnalignedTerm)
                os.system(cmd)
                tmpfilelist.append(outFileUnalignedTerm)
        cmd = "%s/plotDiffUnalignedNumTM.sh %s"%(binpath, " ".join(tmpfilelist))
        os.system(cmd)

    if g_params['isPrintFileRltyCmpclass']:
        CloseFileRltyCmpclass()
    if g_params['isPrintFileRltyHelixCmpclass']:
        CloseFileRltyHelixCmpclass()

    return 0

#}}}
def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isDEBUG'] = False
    g_params['selecttype'] = 'all'
    g_params['outpath'] = ""
    g_params['testlevel'] = 0

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
    g_params['alignrange'] = 'all'
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

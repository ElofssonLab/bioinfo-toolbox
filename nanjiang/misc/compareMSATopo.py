#!/usr/bin/env python

# Filename: compareMSATopo.py
# Description:
#   Finding topology variations and further classify topology variations by
#   compareing aligned topologies, either multiple alignment or pairwise
#   alignment
#
# Author:
#   Nanjiang Shu  nanjiang.shu@scilifelab.se

# 2011-10-25
# Note that length of TM regions do not vary a lot.
# it is about 22.3(3.3) residues and SCAMPI or PRODIV usually predict with a 
# fixed length of about 21(1.1)
# therefore, when mapping TM regions, one can consider that the length of TM
# regions are the identical.

# comapre topologies in the context of multiple sequence alignment
# 1st. Select the largest group with the identical topology by all to all
#      pairwise topology comparison
# 2nd. Get the consensus topology based on topologies in the largest IDT
#      (identical) group.
# 3rd. Compare the resting topologies to the consensus topology and analyze 
#      the N/C terminal insertion/deletions and internal i/ds
#

# Three methods for pairwise topology comparison are used.
# 1. Gapless
# 2. Local 
# 3. Global

# ChangeLog#{{{
# ChangeLog 2011-10-21 
#   major change: when method_getIDTgroup == 1 and method_comparison == 1
#   (which are default values), original topologies instead of trimmed topology
#   are enough, since for method_getIDTgroup = 1, topology comparison is based
#   on posTM, not the topology itself. And the locaiton of TM regions can be
#   obtained from the original topology (with gaps). Note that trimTopoMSA can
#   take very long time in some cases. 
#
# ChangeLog 2011-10-24
#   GetTMPosition updated, the new algorithm is much faster than using finditer
#
# ChangeLog 2011-10-25
#   A bug found for MappingTM, topology with different number of TM regions are
#   classifed as OK or SHIFT, which is wrong. E.g. for PF01944
#   This bug is fixed. MapAlignedTMregion is rewritten.
#   Also, the program has been speed up by 2x in general
#
# ChangeLog 2011-10-26 
#   dgscore will be supplied in the fasta format topology alignment file,
#   therefore -tdg option is not necessary. If dg values are supplied, it will
#   be used, otherwise, listed as NA (not available)
#
# ChangeLog 2011-10-27 
#   Support for pairwise comparison, pairwise alignments are placed
#   sequentially in the input Fasta file
#   option -mode INT is added
#
# ChangeLog 2011-10-30 
#   1. Sequence identity also read out from topology alignment file
#   2. The return value of ReadTopoWithDGScore() is merged into recordList
#
# ChangeLog 2011-11-03
#   For pairwise comparison, the input file size can be very large when there
#   are lots of pairwise aligments. Read in all alinged  topoligies at once is
#   not necessary and might cause memory overflow. This has been modified.
#
# ChangeLog 2011-11-07 
#   For pairwise comparison, write //Begin //End for each record
#   WriteAna(ana1,fpout, "1"), write 1 or 2 instead of ID
#
# ChangeLog 2011-11-15
#   1. A bug in GetConsensusTopo() is fixed, the correct i,o sate are assigned at
#      the beginnings and ends of TM helices to avoid the confusion of TM
#      helices when all resideus between two TM helices are predicted as gaps.
#      In addition, proof reading is added to make sure that the i, o state
#      between two TM helices are consistent.
#   2. IsIdenticalTopology is updated, a minimal overlap of TM helices
#      is set. 
#   3. posTMcons are not taken as the average of posTMs from all identical
#   topologies, but from 70% of middle posTMs sorted by position. Similar as
#   average of (remove minimum and remove maximum)
#
# ChnageLog 2011-11-16 
#   1. option -origtopomsa removed, all input is considered as original
#      topology multiple sequence alignment, trimmed topology is not recommened
#      as input, but if supplied, it will be considered as a MSA with no gaps.
#
# ChangeLog 2011-11-20
#   grouped functionality added with the options
#   -og FILE
#   -wog FILE
#
# ChangeLog 2011-11-21
#   1. Remove unnecessary gaps in the resorted alignment output
#   2. Write "group of ... "only when there are more than 2 group members
#
# ChangeLog 2012-03-16
#   New option added
#   -woc    FILE   Output the clustered topomsa to file, only for mode 1
#   output also alignment in clusters (clusters of identical topology)
#
# ChangeLog 2012-05-15
#   For pairwise comparison, a simple method added for duplications detection
#
# ChangeLog 2012-07-05
#   For duplication, using hhsearch
#   For variations happens only at N-terminal, check whether it is signal
#   peptide
#
# ChangeLog 2012-08-30
#   For pairwise comparison, method 1 is created.
#
# ChangeLog 2012-10-01
#   Bug fixed
#       if (min(mapArray1) == _TM2TM and min(mapArray2) == _TM2TM and
#   changed to
#       if (max(mapArray1) == _TM2TM and max(mapArray2) == _TM2TM and 
#   since _TM2TM is the smallest, 
#
# ChangeLog 2012-10-22
#   compare_method 2 is added for multiple comparison, classes are classified
#   as IDT, INV, Only TM2GAP, Only TM2SEQ, Both TM2GAP and TM2SEQ
#
# ChangeLog 2013-02-06
#   1. Fixed keyError 'repIndex', groupList[0] is the largest identical group,
#   do not need to output, it is already the consensus
#
# ChangeLog 2013-03-20
#   option -rmsp is added, when enabled, signal peptide will be masked, e.g.
#   for protien with topology
#       iiiiiMMMMMMMMMMMMMoooooooooooooMMMMMMMMMiiiiiii
#   if 1-24 is signal peptide, then the new topology will be 
#       oooooooooooooooooooooooooooooooMMMMMMMMMiiiiiii
# ChangeLog 2013-05-15
#   option -localali FILE
#   added, this option is at the most works only for pairwise comparison
#   when local alignment is added, aligned resides in lower characters are
#   considered as unaligned.
# ChangeLog 2013-05-28
#   for pairwise comparison method 1, threshold_TM2TM is disabled, a consistent
#   definition of matched TM helice is used, that two TM helices should be
#   considered as matched only when they are sharing min_TM_overlap residues.
# ChangeLog 2013-07-11
#   option "-woinv FILE" added
# ChangeLog 2013-09-04
#   option "-cmpsp" is added. This enables the comparison of Signal Peptide,
#   that is, there are intotal 4 types: SP, TM, GAP, SEP
#   Then there will be 7 types of comparison: 
#   TM2TM, TM2GAP, TM2SEQ, SP2SP, SP2TM, SP2GAP, SP2SEQ
#   When doing signal peptide comparison, the topology itself should have
#   signal peptide filtered, and then provide the signal peptide location to
#   determine the signal peptide compairson category in one of the SP2SP,
#   SP2TM, SP2GAP, SP2SEQ
# ChangeLog 2013-10-04
#   method 2
#   Helix mapping in the groups of TM2TM, TM2GAP, TM2SEQ, TM2SP
#   protein level classified as, IDT, INV, TM2GAP, TM2SEQ, TM2SP, Mixed 
# ChangeLog 2013-10-22
#   option "-cmpdup" is added, This will add DUP or nonDUP information to
#   cmpclass
# ChangeLog 2013-11-26
#   option "-tableinfo" and "-seqidttype" is added, this is used to obtain
#   sequence identities for each compared pair in pairwise mode
# #}}}

#import argparse
import os
import sys
import re
import comptopo as ct
import myfunc
import libtopologycmp as lcmp
import math
import tempfile
from operator import itemgetter
import copy
#import cProfile
import subprocess


BLOCK_SIZE = myfunc.BLOCK_SIZE

#define TM map value, method 0
_UN_MAPPED = 0
_DIRECT_MAPPED = 9
_SHIFT_MAPPED = 1

#
GAP = "-"


# 
TH_RATIO = 1.0

#define TM map value, method 1
_UNALIGNED = -1
_TM2TM = 0
_TM2GAP = 1
_TM2SEQ = 2
_TM2SP = 3

_SP2SP = 4
_SP2TM = 5
_SP2GAP = 6
_SP2SEQ = 7

#init_dgvalue, meaning that the value is not set
INIT_DGVALUE = -99999.0
INIT_SEQUENCE_IDENTITY = -99.0; # Init value for sequence identity, 
                                # added 2011-10-30

comparisonClassNameListAll = {}
# for method_comparison 0, 1,2, 11 ...
comparisonClassNameListAll[0] = ["OK","SHIFT","INV","DIFF"];  
comparisonClassNameListAll[1] = ["OK","SHIFT","INV","INV_SHIFT", "DIFF"]
comparisonClassNameListAll[2] = ["IDT","INV","TM2GAP","TM2SEQ", "TM2GAP_AND_TM2SEQ"]
comparisonClassNameListAll[3] = ["IDT","INV","TM2GAP","TM2SEQ", "TM2SP", "Mixed"]
comparisonClassNameListAll[11] = ["SAME_NTERM","DIFF_NTERM"]; 

progname = os.path.basename(sys.argv[0])
usage="""
Usage: %s topology-alignment-file

Description:
    Compare aligned topology, either pairwise alignment or multiple alignment.
    Note that topology-alignment-file should be in Fasta format
    When the topo-MSA is untrimmed, origmsa is not needed

Options:
  -mode       INT   Comparison mode, (default: 1)
                    0: pairwise. Pairwise alignments are placed sequentially
                    1: multiple. The given file is a multiple sequence alignment
                       of topologies.  get the consensus from the largest identical
                       group and then compare the rest to the consensus
  -mp         INT   Pairwise comparison method, 0 or 1 (default: 0)
                    3: helix level, TM2TM, TM2GAP, TM2SEQ, TM2SP
                       protein level, IDT, INV, TM2GAP, TM2SEQ, TM2SP, Mixed
  -mm         INT   Method for multiple comparison, (default: 0)
                        0: comparing topologies based on alignment
                        9: comparing topologies based on number of TM helices
  -cmpsp            Compare also signal peptide, when this is enabled, signalp
                    file must be supplied. Currently, this works only for pairwise
                    comparison -mp=2. SP-TM, SP-GAP, SP-SEQ are listed
  -cmpdup           Show also duplication or non duplication in cmpclass
  -dupfile   FILE   duplication result file
  -signalp   FILE   Set signal peptide prediction file
  -localali  FILE   Local alignment file
  -tableinfo FILE   Set pairwise alignment table info, get more pairwise statistics
  -seqidttype INT   Type of sequence identity, (default: 1)
                      0: seqIDT = numIDTRes /alnLength
                      1: seqIDT = numIDTRes / min(len1, len2)
                      2: seqIDT = numIDTRes / (alnLength - NumGAP)
  -o         FILE   Output the result to file
  -obad      FILE   Output bad mapped pairs
  -og        FILE   Output the grouped result file. Identical topologies are
                    grouped
  -wo        FILE   Output the re-sorted topomsa to file, only for mode 1
  -wog       FILE   Output the grouped re-sorted topomsa to file, only for mode 1
  -woc       FILE   Output the clustered topomsa to file, only for mode 1
  -woinv     FILE   Output inverted topology
  -log       FILE   Output the logfile 
  -outpath    DIR   Output the result to outpath, (default: the same folder as the input file)
  -v          INT   Set verbose level, (default: 1)
  -h, --help        Print this help message and exit

Advanced options
see documentation for details
  -mcmp   INT    method when comparing to the consensus. (default: 1)
                 This is for multiple comparison
                 0: Four class comparison
                 1: Six class comparison, TM regions are mapped by direct
                    mapping followed by shifted mapping
                 2: Topology comparisons are classified as IDT, INV, Only
                 TM2GAP, Only TM2SEQ, Both TM2SEQ
                 when this method is enabled, midt == 1

  -mcons  INT    method for getting consensus topology, (default: 2)
                 0: determined by the state with the highest percentage at
                    each column
                 1: first find the core TM and then extend at both ends
                 2: find the average location of all TM regions at each
                    aligned position

  -midt   INT    method for getting the largest IDT group,(default: 1)
                 0: get the largest identical group by trimmed topology, slow.
                 1: get the largest identical group by locations of TM
                    regions, fast.
  -mino, -minoverlap INT
                 Topologies are considered identical only when all TM helices
                 are overlapping with at least minimal residues.

  -wt, -write-sorted-trimmedmsa FILE
                 Output the re-sorted topomsa to the specified file
  -wohtml  FILE  Output the re-sorted topomsa to the specified file in html
                 format
  -wthtml  FILE  Output the re-sorted topomsa to the specified file in html
                 format
  -pi, --print-index-idt
                 Print the index of sequences in the identity group
  -pdgcons       Print the dg values of the consensusTopo
  -maxdgdiff     Set the maximum allowed DG value difference for a pair of
                 topologies with similar DG values, (default: 0.5)
  -c-tm2tm       Threshold_TM2TM, (default: 1/3, i.e. 33.33%%)
  -c-tm2gap      Threshold_TM2GAP, (default: 1/2, i.e. 50%%) 
  -pdbtosp FILE  PDB to swissprot id (uniprot ac) maplist 
  -sprot   FILE  Supply swissprot aclist

Debug options:
  -debug-consensus
  -debug-tmmapping
  -debug-grouping

Created 2011-08-01, updated 2013-10-07, Nanjiang Shu

Examples:
    # compare topologies in a multiple sequence alignment
      compareMSATopo.py test.topomsa.topowithdgscore -o tmp1/t1.rst \\
              -wo tmp1/t1.fa

    # compare topologies of a number of pairwise alignments
      compareMSATopo.py test.pairaln -mode 0 -o tmp1/t1.rst
"""%(progname)
def PrintHelp():#{{{
    print usage
#}}}
def Get_IOState_upstream(topo, begin_TM):#{{{
    i = begin_TM
    while i >= 0:
        if topo[i] in ['i','o']:
            return topo[i]
        i -= 1
    return ''
    #}}}
def StatIOMFreq(topo):#{{{
    length = len(topo)
    cntGap = topo.count(GAP)
    cntTM = topo.count('M')
    cntSeq = length - cntGap - cntTM
    freqTM = myfunc.FloatDivision(cntTM , length)
    freqGap = myfunc.FloatDivision(cntGap , length)
    freqSeq = myfunc.FloatDivision(cntSeq , length)
    return (cntTM, cntGap, cntSeq, freqTM, freqGap, freqSeq)
#}}}

def WriteDetailedDIFFTopo_old(anaDIFFConsList, anaDIFFQueryList, idList,  #{{{
        numTMList, seqLenList, indexDIFF, fpout):
    if fpout != None:
        for i in range(len(anaDIFFConsList)):
            if  anaDIFFConsList[i] != {} or anaDIFFQueryList[i] != {}:
                numLineToPrint=0
                if  anaDIFFConsList[i] != {} :
                    numLineToPrint += (
                            bool(anaDIFFConsList[i]['Nterm']['numTMunmapped'])
                            + bool(anaDIFFConsList[i]['Cterm']['numTMunmapped'])
                            + len(anaDIFFConsList[i]['internal']))
                if anaDIFFQueryList[i] != {}:
                    numLineToPrint += (
                            bool(anaDIFFQueryList[i]['Nterm']['numTMunmapped'])
                            +
                            bool(anaDIFFQueryList[i]['Cterm']['numTMunmapped'])
                            + len(anaDIFFQueryList[i]['internal']))
                fpout.write("SeqID %-15s | numTM  %2d | " 
                        % (idList[indexDIFF[i]], numTMList[indexDIFF[i]]) 
                        + "length  %5d | numLineRecord %d\n" 
                        % (seqLenList[indexDIFF[i]], numLineToPrint))
            if  anaDIFFConsList[i] != {} :
                lcmp.WriteAna(anaDIFFConsList[i],fpout , "Cons") 
            if  anaDIFFQueryList[i] != {}:
                lcmp.WriteAna(anaDIFFQueryList[i],fpout , "Query") 
#}}}
def WriteDetailedDIFFTopo(recordList, fpout):  #{{{
    if fpout != None:
        for i in xrange(len(recordList)):
            record = recordList[i]
            print >> fpout, "//Begin record", i+1
            lcmp.WriteOverallInfo_pairwise(record['id1'], record['id2'],
                    INIT_SEQUENCE_IDENTITY, record['cmpclass'],
                    record['numTM1'], record['numTM2'], record['seqLength1'],
                    record['seqLength2'], fpout,
                    g_params['uniprot2pdbMap'], g_params['swissprotAcSet'])
            PrintMappedArray(record['mapArray1'], record['mapArray2'],
                record['id1'], record['id2'], fpout)
            print >> fpout, "NtermTopo1", record['NtermTopo1']
            print >> fpout, "NtermTopo2", record['NtermTopo2']
            if 'ana1' in record and record['ana1'] != {}:
                lcmp.WriteAna(record['ana1'], fpout, "1") 
            if 'ana2' in record and record['ana2'] != {}:
                lcmp.WriteAna(record['ana2'], fpout, "2") 
            print >> fpout, "//End record", i+1
#}}}
def WritePairwiseRecord_method1(recordList, fpout):  #{{{
# here
    if fpout != None:
        cnt = 0
        seqidt_null = 0.0
        for rd in recordList:
            if rd['cmpclass'] == "":
                continue
            cnt += 1
# Write overall information
            print >> fpout, "//Begin record", cnt
            lcmp.WriteOverallInfo_pairwise(rd['id1'], rd['id2'],
                    seqidt_null,  rd['cmpclass'], rd['numTM1'],
                    rd['numTM2'], rd['seqLength1'], rd['seqLength2'],
                    fpout, g_params['uniprot2pdbMap'], g_params['swissprotAcSet'])
            PrintMappedArray_method1_1(rd['mapArray1'], rd['mapArray2'],
                    rd['id1'], rd['id2'], fpout)
            print >> fpout, "//End record", cnt
#}}}
def WriteGroupedDetailedDIFFTopo(recordList, groupList, dgScoreList, #{{{
        idxFull2OtherClass, fpout): 
    if fpout != None:
        cntRecord = 0
        for grp in groupList[1:]: #groupList[0] is for the IDT group
            record = recordList[idxFull2OtherClass[grp['repIndex']]]
            print >> fpout, "//Begin record", cntRecord + 1
            lcmp.WriteOverallInfo_pairwise(record['id1'], record['id2'],
                    INIT_SEQUENCE_IDENTITY, record['cmpclass'],
                    record['numTM1'], record['numTM2'], record['seqLength1'],
                    record['seqLength2'], fpout, g_params['uniprot2pdbMap'],
                    g_params['swissprotAcSet'])
            if len(grp['index-members-without-rep']) >= 0:
                cntMember = 0; # member 0 is the representative
                sizeIDList = [len(recordList[idxFull2OtherClass[j]]['id2']) for
                        j in [grp['repIndex']] +
                        grp['index-members-without-rep']]
                maxSizeID = max(sizeIDList)
                for j in [grp['repIndex']] + grp['index-members-without-rep']:
                    rd = recordList[idxFull2OtherClass[j]]
                    fpout.write("Member-%d: " % cntMember + 
                            "%*s %4d " % (maxSizeID, rd['id2'], rd['seqLength2']))
                    fpout.write(" dgscore (%d): "  % rd['numTM2'])
                    for k in range(rd['numTM2']):
                        fpout.write("%6.3f " % dgScoreList[j][k])
                    cntMember += 1
                    fpout.write('\n')

            PrintMappedArray(record['mapArray1'], record['mapArray2'],
                record['id1'], record['id2'], fpout)
            print >> fpout, "NtermTopo1", record['NtermTopo1']
            print >> fpout, "NtermTopo2", record['NtermTopo2']
            if 'ana1' in record and record['ana1'] != {}:
                lcmp.WriteAna(record['ana1'], fpout, "1") 
            if 'ana2' in record and record['ana2'] != {}:
                lcmp.WriteAna(record['ana2'], fpout, "2") 
            print >> fpout, "//End record", cntRecord+1
            cntRecord += 1
        if g_params['verbose'] > 0:
            print ("Grouped result file output to %s" %
                    g_params['outGroupedResultFile'])
#}}}


def WriteHTMLHeader(fpout, title):#{{{
    print >> fpout, "<HTML>"
    print >> fpout, "<title>%s</title>"%(title)
    print >> fpout, "<style type=\"text/css\">"
    print >> fpout, "<!--"
    print >> fpout, ("td {font-family: \"SansSerif\", \"SansSerif\", " + 
            "mono; font-size: 10px; text-align: center; }")
    print >> fpout, "-->"
    print >> fpout, "</style>"
    print >> fpout, "<BODY>"
    print >> fpout, "<table border=\"1\"><tr><td>"
    print >> fpout, ""
    print >> fpout, ("<table border=\"0\"  cellpadding=\"0\"" +
            " cellspacing=\"0\">  ")
#}}}
def WriteHTMLScaleBar(fpout, lengthMSA):#{{{
    print >> fpout,"<tr><td colspan=\"6\"></td>"; 
    i=10
    while i < lengthMSA:
        print >> fpout,"<td colspan=\"9\">%d<br>|</td><td></td> "%(i)
        i+=10
    print >> fpout,"</tr>"; 
#}}}
def WriteHTMLTail(fpout):#{{{
    print >> fpout, "</table>"
    print >> fpout, "</td></tr></table>"
    print >> fpout, "</body>"
    print >> fpout, "</html>"
#}}}
def WriteHTMLColorChar(fpout, ch):#{{{
    if ch == "M":
        print >> fpout, "<td bgcolor=\"#FF0000\">%s</td>" %(ch);# red
    elif ch == "i":
        print >> fpout, "<td bgcolor=\"#FFFFCC\">%s</td>" %(ch);# faded yellow
    elif ch == "o":
        print >> fpout, "<td bgcolor=\"#CCFFFF\">%s</td>" %(ch);# faded blue
    else:
        print >> fpout, "<td>%s</td>" %(" ");   # faded blue
#}}}
def WriteHTMLID(fpout, seqid):#{{{
    print >> fpout, "<tr><td nowrap>%s&nbsp;&nbsp;</td> " % (seqid)
#}}}

def WriteSortedOrigTopoMSA(outSortedOrigTopoMSAFile, idList,   #{{{
        topoSeqList,posTMList, numIDTTopo, consensusTopo, numTM_IDT,
        indexIDTTopo, indexClass,comparisonClassNameList, infoDIFF):
# @params:
#   topoSeqList  Untrimmed topology sequences
    numSeq = len(topoSeqList)
    newAnnoList = []
    newTopoSeqList = []
    # 1. consensus
    newAnnoList.append("%s of %d, nTM=%d"%("Consensus", numIDTTopo,
        myfunc.CountTM(consensusTopo)))
    newTopoSeqList.append(consensusTopo)
    # 2. IDTgroup
    for i in range(len(indexIDTTopo)):
        newAnnoList.append("%s IDT %d nTM=%d"%(idList[indexIDTTopo[i]],
            i+1, numTM_IDT))
        newTopoSeqList.append(topoSeqList[indexIDTTopo[i]])
    # 3. cmpClassOtherTopoList
    if "DIFF" in comparisonClassNameList:
        rng =  len(comparisonClassNameList)-1
    else:
        rng = len(comparisonClassNameList)
    for icls in range(rng):
        for i in range(len(indexClass[icls])):
            newAnnoList.append("%s %s %d nTM=%d" % (
                idList[indexClass[icls][i]], comparisonClassNameList[icls],
                i+1, len(posTMList[indexClass[icls][i]])))
            newTopoSeqList.append(topoSeqList[indexClass[icls][i]])
    # 3.2 for DIFF
    if "DIFF" in comparisonClassNameList:
        idx=comparisonClassNameList.index("DIFF")
        indexDIFF=indexClass[idx]
        for i in range(len(indexDIFF)):
            tmpanno = ("%s DIFF %d nTM=%d "%(idList[indexDIFF[i]], i+1,
                len(posTMList[indexDIFF[i]])))
            if len(infoDIFF) > 0:
                info=infoDIFF[i]
                tmpanno += ("N_ins=%d N_del=%d C_ins=%d C_del=%d" % (
                    info["N_ins"], info["N_del"], info["C_ins"],
                    info["C_del"]))
                tmpanno += ("  | I_ins:")
                for ii in range(numTM_IDT-1):
                    tmpanno += (" %2d"%info["i_ins"][ii])
                tmpanno += (" | I_del: ")
                for ii in range(numTM_IDT-2):
                    tmpanno += (" %2d"%info["i_del"][ii])

            newAnnoList.append(tmpanno)
            newTopoSeqList.append(topoSeqList[indexDIFF[i]])
    newTopoSeqList = lcmp.RemoveUnnecessaryGap(newTopoSeqList)
#data checking
    numIncludedSeq= len(indexIDTTopo)+sum(len(x) for x in indexClass)
    if numIncludedSeq != numSeq:
        sys.stderr.write(
                "Error! For file %s, numIncludedSeq (%d) != numSeq (%d)\n",
                inFile, numIncludedSeq, numSeq) 
    try:
        fpout = open (outSortedOrigTopoMSAFile, "w")
        for i in xrange(len(newAnnoList)):
            fpout.write(">%s\n" % newAnnoList[i])
            fpout.write("%s\n" % newTopoSeqList[i])
        fpout.close()
        if g_params['verbose'] > 0:
            print ("Sorted original Topology MSA file output to %s" %
                    (outSortedOrigTopoMSAFile))
        return 0
    except IOError:
        print >> sys.stderr, "Failed to write to file %s"%outSortedOrigTopoMSAFile
        return -1

#}}}
def WriteSortedClusteredTopoMSA(filename, #{{{
        idList, topoSeqList, posTMList, Mcmp):
    numSeq = len(idList)

#     cntIDTTopoTupleList = []
#     for i in xrange(numSeq):
#         cntIDTTopoTupleList.append((i, sum(Mcmp[i])))
#     sortedCntIDTTopoTupleList = sorted(cntIDTTopoTupleList, key=lambda x: x[1], reverse=True)
#     for i in xrange(numSeq):
#         sys.stdout.write("%-4d: "%sortedCntIDTTopoTupleList[i][1])
#         for j in xrange(numSeq):
#             sys.stdout.write(" %d"%Mcmp[sortedCntIDTTopoTupleList[i][0]][j])
#         sys.stdout.write("\n")

    tmpMcmp = copy.deepcopy(Mcmp)
    markArray = [1]*numSeq
    clusterList = []
    cntCluster = 0
    while 1:
        cntIDTTopo = [sum(ii) for ii in tmpMcmp]
        maxCnt = max(cntIDTTopo)
        maxIndex = cntIDTTopo.index(maxCnt)
        clusterList.append([])
        for j in xrange(numSeq):
            if tmpMcmp[maxIndex][j] == 1:
                clusterList[cntCluster].append(j)
                markArray[j] = 0
                for i in xrange(numSeq):
                    tmpMcmp[i][j] = 0
        cntCluster += 1
        if sum(markArray) == 0:
            break
#     for i in xrange(len(clusterList)):
#         print len(clusterList[i]), clusterList[i]

    # print sequences

    newAnnoList = []
    newTopoSeqList = []
    for i in xrange(len(clusterList)):
        for idx in clusterList[i]:
            newAnnoList.append("%s, nTM=%d ClusterNo=%d numSeqInCluster=%d"%(idList[idx],
                len(posTMList[idx]), i+1, len(clusterList[i])))
            newTopoSeqList.append(topoSeqList[idx])
    try:
        fpout = open (filename, "w")
        for i in xrange(len(newAnnoList)):
            fpout.write(">%s\n" % newAnnoList[i])
            fpout.write("%s\n" % newTopoSeqList[i])
        fpout.close()
        if g_params['verbose'] > 0:
            sys.stdout.write(("Clustered sorted original Topology MSA file " +
                "output to %s\n" % (filename)))
    except IOError:
        print >> sys.stderr, ("Failed to write to file %s" % filename)


#}}}
def WriteInvertedTopologyGroup_MSA(outfile, #{{{
        idList, topoSeqList, posTMList, Mcmp):
    numSeq = len(idList)
# calculate inverted topology
# cluster proteins with identical topology, find the largest two groups with
# inverted topology
# 2013-07-11 : not finished
    tmpMcmp = copy.deepcopy(Mcmp)
    markArray = [1]*numSeq
    clusterList = []
    cntCluster = 0
    while 1:
        cntIDTTopo = [sum(ii) for ii in tmpMcmp]
        maxCnt = max(cntIDTTopo)
        maxIndex = cntIDTTopo.index(maxCnt)
        clusterList.append([])
        tmpCluster = []
        for j in xrange(numSeq):
            if tmpMcmp[maxIndex][j] == 1:
                tmpCluster.append(j)
                markArray[j] = 0
                for i in xrange(numSeq):
                    tmpMcmp[i][j] = 0
        cntCluster += 1
        if sum(markArray) == 0:
            break
#     for i in xrange(len(clusterList)):
#         print len(clusterList[i]), clusterList[i]

    # print sequences

    newAnnoList = []
    newTopoSeqList = []
    for i in xrange(len(clusterList)):
        for idx in clusterList[i]:
            newAnnoList.append("%s, nTM=%d ClusterNo=%d numSeqInCluster=%d"%(idList[idx],
                len(posTMList[idx]), i+1, len(clusterList[i])))
            newTopoSeqList.append(topoSeqList[idx])
    try:
        fpout = open (outfile, "w")
        for i in xrange(len(newAnnoList)):
            fpout.write(">%s\n" % newAnnoList[i])
            fpout.write("%s\n" % newTopoSeqList[i])
        fpout.close()
        if g_params['verbose'] > 0:
            sys.stdout.write(("Clustered sorted original Topology MSA file " +
                "output to %s\n" % (filename)))
    except IOError:
        print >> sys.stderr, ("Failed to write to file %s" % filename)


#}}}
def WriteInvertedTopologyPairwise_MSA(outfile, #{{{
        idList, topoSeqList, posTMList, NtermStateList):
    numSeq = len(idList)
# given all-to-all pairwise comparison, calculate the fraction of pairwise
# comparison with inverted topology
# 2013-07-11
    Mcmp = GetInvertedTopologyMatrix(topoSeqList, NtermStateList, posTMList)
    seqLenList=[ len(tp.replace(GAP,'')) for tp in topoSeqList]; 
    numINVpair = sum(a.count(1) for a in Mcmp)
    numAllpair = numSeq*(numSeq-1)/2
    pairIdxList = []
    uniq_inv_idset = set([])
    for i in xrange (numSeq):
        for j in xrange(i+1, numSeq):
            if Mcmp[i][j] == 1:
                pairIdxList.append((i,j))
                uniq_inv_idset.add(idList[i])
                uniq_inv_idset.add(idList[j])
    try:
        fpout = open(outfile, "w")
        fpout.write("General: %4d %6d %6.3f\n"%(numINVpair, numAllpair,
                myfunc.FloatDivision(numINVpair, numAllpair)))
        fpout.write("UniqIDList: %4d "%(len(uniq_inv_idset)))
        for idd in uniq_inv_idset:
            fpout.write(" %s"%(idd))
        fpout.write("\n")
        for tup in pairIdxList:
            fpout.write("Pair: %s %s %s %s %2d %4d %4d\n"%(
                idList[tup[0]], idList[tup[1]], 
                NtermStateList[tup[0]], NtermStateList[tup[1]],
                len(posTMList[tup[0]]),
                seqLenList[tup[0]], seqLenList[tup[1]]))
        fpout.close()
    except IOError:
        print >> sys.stderr, "Failed to write to file %s"%outfile
        pass
    outfile2 = outfile + ".INVseq.idlist"
    try:
        fpout = open(outfile2, "w")
        for idd in uniq_inv_idset:
            fpout.write("%s\n"%(idd))
        fpout.close()
    except IOError:
        print >> sys.stderr, "Failed to write to file %s"%outfile2
        pass
    
#}}}
def WriteGroupedSortedOrigTopoMSA(filename, groupList, idList,  #{{{ 
        topoSeqList, numIDTTopo, consensusTopo, indexIDTTopo,
        comparisonClassNameList):
# @params:
#   topoSeqList  Untrimmed topology sequences
    numSeq = len(topoSeqList)
    newAnnoList = []
    newTopoSeqList = []
    # 1. consensus
    newAnnoList.append("%s of %d, nTM=%d"%("Consensus", numIDTTopo,
        myfunc.CountTM(consensusTopo)))
    newTopoSeqList.append(consensusTopo)
    # 3. cmpClassOtherTopoList
    for cmpclass in comparisonClassNameList:
        cnt = 0
        numGroup = len(groupList)
        for i in xrange(1, numGroup):
            grp = groupList[i]
            if grp['cmpclass'] == cmpclass:
                tmpanno = (">%s %s %d nTM=%d " % (idList[grp['repIndex']],
                    cmpclass, cnt+1, grp['repNumTM']))
                if len(grp['index-members-without-rep']) > 0:
                    tmpanno += ("group of %d topologies, Also " %
                            len(grp['index-members']))
                    for j in grp['index-members-without-rep']:
                        tmpanno += ("%s " % idList[j])
                newAnnoList.append(tmpanno)
                newTopoSeqList.append(topoSeqList[grp['repIndex']])
                cnt += 1
    newTopoSeqList = lcmp.RemoveUnnecessaryGap(newTopoSeqList)
#     print newAnnoList
    try:
        fpout = open (filename, "w")
        for i in xrange(len(newAnnoList)):
            fpout.write(">%s\n" % newAnnoList[i])
            fpout.write("%s\n" % newTopoSeqList[i])
        fpout.close()
        if g_params['verbose'] > 0:
            sys.stdout.write(("Grouped sorted original Topology MSA file " +
                "output to %s\n" % (filename)))
    except IOError:
        print >> sys.stderr, ("Failed to write to file %s" % filename)
#}}}
def WriteGroupedTopoAnaResult(numSeq, cmpToConsRecordList, groupList, #{{{
        numTM_IDT, indexClass, indexIDTTopo, dgScoreList,
        comparisonClassNameList, fpout):
# Write overall information
    fulllist=[x for x in range(numSeq)]
    indexOtherTopo=list(set(fulllist)-set(indexIDTTopo)); # set subtraction
    idxFull2OtherClass = {}; # a dictionary
    for i in range(len(indexOtherTopo)):
        idxFull2OtherClass[indexOtherTopo[i]] = i
    numIDTTopo = len(indexIDTTopo)

    numTMcons = groupList[0]['repNumTM']
    fpout.write("//Begin CMPMSA\n")
    lcmp.WriteOverallInfo_msa(comparisonClassNameList, g_params['rootname'],
            numSeq, numTM_IDT, numTMcons, numIDTTopo, indexClass, fpout)
# Write detailed topology variation info
    WriteGroupedDetailedDIFFTopo(cmpToConsRecordList, groupList, dgScoreList, 
            idxFull2OtherClass, fpout)
    fpout.write("//End CMPMSA\n")

#}}}


def WriteSortedTrimmedTopoMSA(outSortedTrimmedTopoMSAFile, idList,   #{{{
        topoSeqList, posTMList, numIDTTopo, consensusTopo, numTM_IDT,
        indexIDTTopo, indexClass,comparisonClassNameList, infoDIFF):
    try: 
        numSeq = len(topoSeqList)
        fpout = open (outSortedTrimmedTopoMSAFile, "w")
        # 1. consensus
        fpout.write(">%s of %d, nTM=%d\n" % ("Consensus", numIDTTopo,
            myfunc.CountTM(consensusTopo)))
        fpout.write("%s\n"%consensusTopo)
        # 2. IDTgroup
        for i in range(len(indexIDTTopo)):
            fpout.write(">%s IDT %d nTM=%d\n"%(idList[indexIDTTopo[i]], i+1,
                numTM_IDT))
            fpout.write("%s\n"%topoSeqList[indexIDTTopo[i]])

        # 3. cmpClassOtherTopoList
        for icls in range(len(comparisonClassNameList)-1):
            for i in range(len(indexClass[icls])):
                fpout.write(">%s %s %d nTM=%d\n" % (idList[indexClass[icls][i]],
                    comparisonClassNameList[icls], i+1,
                    len(posTMList[indexClass[icls][i]])))
                fpout.write("%s\n"%topoSeqList[indexClass[icls][i]])
        # 3.2 for DIFF
        idx = comparisonClassNameList.index("DIFF")
        indexDIFF = indexClass[idx]
        for i in range(len(indexDIFF)):
            fpout.write(">%s DIFF %d nTM=%d " % (idList[indexDIFF[i]], i+1,
                len(posTMList[indexDIFF[i]])))
            if len(infoDIFF) > 0:
                info=infoDIFF[i]
                fpout.write("N_ins=%d N_del=%d C_ins=%d C_del=%d"
                        % (info["N_ins"], info["N_del"], info["C_ins"],
                            info["C_del"]))
                fpout.write("  | I_ins:")
                for ii in range(numTM_IDT-1):
                    fpout.write(" %2d"%info["i_ins"][ii])
                fpout.write(" | I_del: ")
                for ii in range(numTM_IDT-2):
                    fpout.write(" %2d"%info["i_del"][ii])

            fpout.write("\n")
            fpout.write("%s\n"%topoSeqList[indexDIFF[i]])
#data checking
        numIncludedSeq= len(indexIDTTopo)+sum(len(x) for x in indexClass)
        if numIncludedSeq != numSeq:
            sys.stderr.write("Error! For file %s"%inFile 
                    + ", numIncludedSeq (%d) != numSeq (%d)\n"
                    % (numIncludedSeq, numSeq)) 
        fpout.close()
        if g_params['verbose'] > 0:
            print ("Sorted Trimmed Topology MSA file output to \"%s\""
                    %(outSortedTrimmedTopoMSAFile))
    except IOError:
        print >> sys.stderr, ("Failed to write file %s" %
                outSortedTrimmedTopoMSAFile)
#}}}
def WriteSortedTrimmedTopoMSAHTML(outSortedTrimmedTopoMSAFileHTML, #{{{
        lengthAlignment, idList, topoSeqList, indexIDTTopo):
    try: 
        fpout = open (outSortedTrimmedTopoMSAFileHTML, "w")
        title=""
        WriteHTMLHeader(fpout, title)
        WriteHTMLScaleBar(fpout, lengthAlignment)
        fulllist=[x for x in range(len(idList))]
        restlist=list(set(fulllist) - set(indexIDTTopo))
        for i in range(len(indexIDTTopo)):
            seqid=idList[indexIDTTopo[i]]
            seq=topoSeqList[indexIDTTopo[i]]
            WriteHTMLID(fpout,seqid)
            for i in range (lengthAlignment):
                WriteHTMLColorChar(fpout, seq[i])

        for i in range(len(restlist)):
            seqid=idList[restlist[i]]
            seq=topoSeqList[restlist[i]]
            WriteHTMLID(fpout,seqid)
            for i in range (lengthAlignment):
                WriteHTMLColorChar(fpout, seq[i])
        WriteHTMLTail(fpout)
        fpout.close()
        if verbose > 0:
            print ("Sorted Trimmed Topology MSA HTML file output to \"%s\""  %
                    (outSortedTrimmedTopoMSAFileHTML))
    except IOError:
        print >> sys.stderr, ("Failed to write to file %s" %
                outSortedTrimmedTopoMSAFileHTML)
#}}}
def WriteSortedOrigTopoMSAHTML(outSortedTrimmedTopoMSAFileHTML,  #{{{
        origTopoMSASeqList, lengthAlignment,idList, topoSeqList,
        indexIDTTopo):
    try:
        fpout = open (outSortedOrigTopoMSAFileHTML, "w")
        title=""
        WriteHTMLHeader(fpout, title)
        WriteHTMLScaleBar(fpout, lengthAlignment)
        fulllist=[x for x in range(len(idList))]
        restlist=list(set(fulllist) - set(indexIDTTopo))
        for i in range(len(indexIDTTopo)):
            seqid=idList[indexIDTTopo[i]]
            seq=seqList[indexIDTTopo[i]]
            WriteHTMLID(fpout,seqid)
            for i in range (lengthAlignment):
                WriteHTMLColorChar(fpout, seq[i])

        for i in range(len(restlist)):
            seqid=idList[restlist[i]]
            seq=seqList[restlist[i]]
            WriteHTMLID(fpout,seqid)
            for i in range (lengthAlignment):
                WriteHTMLColorChar(fpout, seq[i])
        WriteHTMLTail(fpout)
        fpout.close()
        if verbose > 0:
            print ("Sorted original Topology MSA HTML file output to \"%s\""
                    % (outSortedOrigTopoMSAFileHTML))
    except IOError:
        print >> sys.stderr, ("Failed to write to file %s" %
                outSortedOrigTopoMSAFileHTML)
#}}}

def PrintMappedArray(mapArray1, mapArray2, seqID1, seqID2, fpout):#{{{
# Print the mapped arry in aligned form
    if fpout:
        numTM1=len(mapArray1)
        numTM2=len(mapArray2)
        alignList1=[]
        alignList2=[]
        cnt1=0
        cnt2=0
        while cnt1 < numTM1 and cnt2 < numTM2:
            if mapArray1[cnt1] == _UN_MAPPED and mapArray2[cnt2]== _UN_MAPPED:
                alignList1.append("%d"%mapArray1[cnt1])
                alignList2.append("%d"%mapArray2[cnt2])
                cnt1+=1
                cnt2+=1
            else:
                if  (mapArray1[cnt1] != _UN_MAPPED 
                        and mapArray2[cnt2]  != _UN_MAPPED):
                    alignList1.append("%d"%mapArray1[cnt1])
                    alignList2.append("%d"%mapArray2[cnt2])
                    cnt1+=1
                    cnt2+=1
                elif (mapArray1[cnt1] != _UN_MAPPED 
                        and mapArray2[cnt2]  == _UN_MAPPED): 
                    alignList1.append("-")
                    alignList2.append("%d"%mapArray2[cnt2])
                    cnt2+=1
                else:
                    alignList1.append("%d"%mapArray1[cnt1])
                    alignList2.append("-")
                    cnt1+=1
        #tail
        while cnt1 < numTM1:
            alignList1.append("%d" % mapArray1[cnt1])
            alignList2.append("-")
            cnt1 +=1
        while cnt2 < numTM2:
            alignList1.append("-")
            alignList2.append("%d" % mapArray2[cnt2])
            cnt2 +=1
    #print
        maxSizeID = max(len(seqID1), len(seqID2))
        if len(alignList1) == len(alignList2):
            fpout.write("TMMap %*s, numTM1 = %2d:"%(maxSizeID, seqID1, numTM1))
            for i in alignList1:
                fpout.write(" %3s"%i)
            fpout.write("\n")

            fpout.write("TMMap %*s, numTM2 = %2d:"%(maxSizeID, seqID2, numTM2))
            for i in alignList2:
                fpout.write(" %3s"%i)
            fpout.write("\n")
        else:
            sys.stderr.write("Error! alignList1 and alignList2 " 
                    + "not the same size for %s and %s\n" % (seqID1, seqID2))
            print >> sys.stderr, "alignList1:", alignList1
            print >> sys.stderr, "alignList2:", alignList2
            print >> sys.stderr, "mapArray1:",mapArray1 
            print >> sys.stderr, "mapArray2:", mapArray2
#}}}
def PrintMappedArray_method1(mapArray1, mapArray2, seqID1, seqID2, fpout):#{{{
# Print the mapped arry in unaligned form
    if fpout:
        numTM1=len(mapArray1)
        numTM2=len(mapArray2)
        maxSizeID = max(len(seqID1), len(seqID2))
        fpout.write("TMMap %*s, numTM1 = %2d:"%(maxSizeID, seqID1, numTM1))
        for i in mapArray1:
            fpout.write(" %3s"%i)
        fpout.write("\n")

        fpout.write("TMMap %*s, numTM2 = %2d:"%(maxSizeID, seqID2, numTM2))
        for i in mapArray2:
            fpout.write(" %3s"%i)
        fpout.write("\n")
#}}}
def PrintMappedArray_method1_1(mapArray1, mapArray2, seqID1, seqID2, fpout):#{{{
# Print the mapped arry in aligned form
    if fpout:
        numTM1=len(mapArray1)
        numTM2=len(mapArray2)
        maxSizeID = max(len(seqID1), len(seqID2))
        mapstrlist1 = []
        mapstrlist2 = []

        cnt1=0
        cnt2=0
        isHead1 = True
        isHead2 = True
        while cnt1 < numTM1 and cnt2 < numTM2:
            if mapArray1[cnt1] == _UNALIGNED or mapArray2[cnt2] == _UNALIGNED:
                if  mapArray2[cnt2] != _UNALIGNED:
                    if isHead2:
                        mapstrlist1.append("%2d"%mapArray1[cnt1])
                        mapstrlist2.append("%2s"%"")
                        cnt1 += 1
                    else:
                        mapstrlist2.append("%2d"%mapArray2[cnt2])
                        mapstrlist1.append("%2s"%"")
                        cnt2 += 1
                elif mapArray1[cnt1] != _UNALIGNED:
                    if isHead1:
                        mapstrlist2.append("%2d"%mapArray2[cnt2])
                        mapstrlist1.append("%2s"%"")
                        cnt2 += 1
                    else:
                        mapstrlist1.append("%2d"%mapArray1[cnt1])
                        mapstrlist2.append("%2s"%"")
                        cnt1 += 1
                else:
                    mapstrlist1.append("%2d"%mapArray1[cnt1])
                    mapstrlist2.append("%2d"%mapArray2[cnt2])
                    cnt1 += 1
                    cnt2 += 1
            else: # none of them are unaligned
                if mapArray1[cnt1] != _TM2TM:
                    mapstrlist1.append("%2d"%mapArray1[cnt1])
                    mapstrlist2.append("%2s"%"")
                    cnt1 += 1
                elif mapArray2[cnt2] != _TM2TM:
                    mapstrlist2.append("%2d"%mapArray2[cnt2])
                    mapstrlist1.append("%2s"%"")
                    cnt2 += 1
                else:
                    mapstrlist1.append("%2d"%mapArray1[cnt1])
                    mapstrlist2.append("%2d"%mapArray2[cnt2])
                    cnt1 += 1
                    cnt2 += 1
                isHead1 = False
                isHead2 = False

        #tail
        while cnt1 < numTM1:
            mapstrlist1.append("%2d"%mapArray1[cnt1])
            mapstrlist2.append("%2s"%"")
            cnt1 +=1
        while cnt2 < numTM2:
            mapstrlist1.append("%2s"%"")
            mapstrlist2.append("%2d" % mapArray2[cnt2])
            cnt2 +=1



        fpout.write("TMMap %*s, numTM1 = %2d: %s\n"%(maxSizeID, seqID1, numTM1,
            " ".join(mapstrlist1)))
        fpout.write("TMMap %*s, numTM2 = %2d: %s\n"%(maxSizeID, seqID2, numTM2,
            " ".join(mapstrlist2)))
#}}}
def PrintMappedArray_method1_2(mapArray1, mapArray2, posTM1, posTM2, seqID1, seqID2, fpout):#{{{
# Print the mapped arry in aligned form
# using relative positions of TM helices
# not finished
    if fpout:
        numTM1=len(mapArray1)
        numTM2=len(mapArray2)
        maxSizeID = max(len(seqID1), len(seqID2))
        mapstrlist1 = []
        mapstrlist2 = []

        cnt1=0
        cnt2=0
        isHead1 = True
        isHead2 = True
        while cnt1 < numTM1 and cnt2 < numTM2:
            if mapArray1[cnt1] == _UNALIGNED or mapArray2[cnt2] == _UNALIGNED:
                if  mapArray2[cnt2] != _UNALIGNED:
                    if isHead2:
                        mapstrlist1.append("%2d"%mapArray1[cnt1])
                        mapstrlist2.append("%2s"%"")
                        cnt1 += 1
                    else:
                        mapstrlist2.append("%2d"%mapArray2[cnt2])
                        mapstrlist1.append("%2s"%"")
                        cnt2 += 1
                elif mapArray1[cnt1] != _UNALIGNED:
                    if isHead1:
                        mapstrlist2.append("%2d"%mapArray2[cnt2])
                        mapstrlist1.append("%2s"%"")
                        cnt2 += 1
                    else:
                        mapstrlist1.append("%2d"%mapArray1[cnt1])
                        mapstrlist2.append("%2s"%"")
                        cnt1 += 1
                else:
                    mapstrlist1.append("%2d"%mapArray1[cnt1])
                    mapstrlist2.append("%2d"%mapArray2[cnt2])
                    cnt1 += 1
                    cnt2 += 1
            else: # none of them are unaligned
                if mapArray1[cnt1] != _TM2TM:
                    mapstrlist1.append("%2d"%mapArray1[cnt1])
                    mapstrlist2.append("%2s"%"")
                    cnt1 += 1
                elif mapArray2[cnt2] != _TM2TM:
                    mapstrlist2.append("%2d"%mapArray2[cnt2])
                    mapstrlist1.append("%2s"%"")
                    cnt2 += 1
                else:
                    mapstrlist1.append("%2d"%mapArray1[cnt1])
                    mapstrlist2.append("%2d"%mapArray2[cnt2])
                    cnt1 += 1
                    cnt2 += 1
                isHead1 = False
                isHead2 = False

        #tail
        while cnt1 < numTM1:
            mapstrlist1.append("%2d"%mapArray1[cnt1])
            mapstrlist2.append("%2s"%"")
            cnt1 +=1
        while cnt2 < numTM2:
            mapstrlist1.append("%2s"%"")
            mapstrlist2.append("%2d" % mapArray2[cnt2])
            cnt2 +=1



        fpout.write("TMMap %*s, numTM1 = %2d: %s\n"%(maxSizeID, seqID1, numTM1,
            " ".join(mapstrlist1)))
        fpout.write("TMMap %*s, numTM2 = %2d: %s\n"%(maxSizeID, seqID2, numTM2,
            " ".join(mapstrlist2)))
#}}}
def ExtractFromSeqWithAnno(seqWithAnno):#{{{
    """
    Extract information from the record seqWithAnno
    Return (seqID, anno, seq, seqIdentity, dgscore)
    ==updated 2013-03-20
    """
    posAnnoEnd = seqWithAnno.find('\n')
    anno = seqWithAnno[1:posAnnoEnd]
    anno = anno.lstrip('>')
    seqID = myfunc.GetSeqIDFromAnnotation(anno)

    # extract seqIDT from annotation line, this is for pairwise aligned
    # topology file 
    seqIdentity = INIT_SEQUENCE_IDENTITY
    m = re.search('seqIDT=[^\s]* ', anno)
    if m != None:
        seqIdentity = float(m.group(0).split('=')[1])

    dgscore = []; # a list of DG scores of the TM regions of the topology
    m = re.search('{[\s]*dgscore.*}', seqWithAnno[posAnnoEnd:].replace('\n', ' '))
    if m != None:
        strdgscore = m.group(0).lstrip('{').rstrip('}')
        strdgscore = strdgscore.split('dgscore')[1].lstrip(':')
        strs = strdgscore.split()
        for s in strs:
            try:
                dgscore.append(float(s))
            except ValueError, TypeError:
                pass

    seq = seqWithAnno[posAnnoEnd+1:]
    seq = seq.replace('\n','').replace(' ','')
    seq = re.sub("{.*}", '',seq)
    return (seqID, anno, seq, seqIdentity, dgscore)
#}}}
def ReadTopoWithDGScore(infile, BLOCK_SIZE=100000): #{{{
    """
    Read enhanced fasta with dg values at the end of sequences and enclosed by
    {dgvalue }
    Return recordList
    ==updated 2011-10-30
    """
    recordList = []; #recordList is a list of n-tuples
    fpin = None
    try:
        fpin=open(infile,"rb")
    except IOError:
        msg = "Failed to read fastaWithDGScore file %s"
        print >> sys.stderr, msg%(infile)
        return []
    buff = fpin.read(BLOCK_SIZE)
    brokenSeqWithAnnoLine=""; ##for the annotation line broken by BLOCK read
    while buff:
        beg=0
        end=0
        while 1:
            if brokenSeqWithAnnoLine:
                if brokenSeqWithAnnoLine[len(brokenSeqWithAnnoLine)-1] == "\n":
                    end=buff.find(">")
                else:
                    end=buff.find("\n>")
                if end >= 0:
                    seqWithAnno = brokenSeqWithAnnoLine + buff[0:end]
                    recordList.append(ExtractFromSeqWithAnno(seqWithAnno))
                    brokenSeqWithAnnoLine = ""
                    beg=end
                else:
                    brokenSeqWithAnnoLine += buff
                    break

            beg=buff.find(">",beg)
            end=buff.find("\n>",beg+1)
            if beg >= 0:
                if end >=0:
                    seqWithAnno=buff[beg:end]
                    recordList.append(ExtractFromSeqWithAnno(seqWithAnno))
                    beg=end
                else:
                    brokenSeqWithAnnoLine=buff[beg:]
                    break
            else:
                break
        buff = fpin.read(BLOCK_SIZE)

    if brokenSeqWithAnnoLine:
        seqWithAnno=brokenSeqWithAnnoLine
        recordList.append( ExtractFromSeqWithAnno(seqWithAnno))
    fpin.close()
    return recordList
#}}}
def RemoveSignalPeptide(topoRecordList, signalpDict):#{{{
    """
    remove signal peptide in topology, content of topoRecordList will be changed
    @params
    topoRecordList  A list of topology records of n-tuple
                    (seqID, anno, seq, seqIdentity, dgscore)
    """
    newRecordList = []
    for rd in topoRecordList:
        isToRemove = False
        newtopo = ""
        newDGList = []
        seqid = rd[0]
        try:
            sp_pos = signalpDict[seqid]
        except KeyError:
            sp_pos = -1
            pass
        if sp_pos != -1:
            newtopo = lcmp.FilterSignalPeptideInTopology(rd[2], sp_pos)
            if newtopo != rd[2]:
                newDGList = rd[4][1:]
                isToRemove = True
        if isToRemove:
            newRecordList.append((rd[0], rd[1], newtopo, rd[3], newDGList))
        else:
            newRecordList.append(rd)
    return newRecordList
#}}}

def ReadSeqPathMapDict(infile):#{{{
    try:
        seqPathMapDict = {}
        fpin = open(infile, "r")
        line = fpin.readline()
        while line:
            if line[0] != "#":
                strs = line.split("\t")
                if len(strs) == 2:
                    seqPathMapDict[strs[0].strip()] = strs[1].strip()
            line = fpin.readline()
        fpin.close()
        return seqPathMapDict
    except IOError:
        print >> sys.stderr, "Failed to read infile %s"%infile
        return {}
#}}}


def ReadTopoWithDGScoreFromBuffer(buff,recordList, isEOFreached):#{{{
    """ 
    Return (unprocessedBuffer)"""
    if not buff:
        return ""
    brokenSeqWithAnnoLine=""
    beg=0
    end=0
    while 1:
        beg=buff.find(">",beg)
        end=buff.find("\n>",beg+1)
        if beg >= 0:
            if end >=0:
                seqWithAnno=buff[beg:end]
                recordList.append(ExtractFromSeqWithAnno(seqWithAnno))
                beg=end
            else:
                brokenSeqWithAnnoLine=buff[beg:]
                break
        else:
            break
    if isEOFreached and brokenSeqWithAnnoLine:
        seqWithAnno=brokenSeqWithAnnoLine
        recordList.append( ExtractFromSeqWithAnno(seqWithAnno))
        brokenSeqWithAnnoLine=""
    return brokenSeqWithAnnoLine
#}}}

def IsTrimmedMSA(topoSeqList):#{{{
    for i in range (len(topoSeqList)):
        if topoSeqList[i].find(GAP) != -1:
            return False
    return True
#}}}

def GetTMPosition_slow(topo):#{{{
# Get the position of TM helices given the topology
# The return value is a list of 2-tuples
# [ (beg, end), (beg, end)...]
# this version works for topologys with and without gaps
    posTM=[]
    m=re.finditer("([M-]*M[M-]*)",topo)
    #m=pTM.finditer(topo)
    for i in m:
        posTM.append((i.start(0), i.end(0)))
    return posTM
#}}}
def GetUnAlignedString(alnseq1, alnseq2):#{{{
    # for local alignment
    sslist = []
    if len(alnseq1) != len(alnseq2):
        print >> sys.stderr, "alnseq1 and alnseq2 length not matched"
        return ""
    else:
        for j in xrange(len(alnseq1)):
            if alnseq1[j].islower() or alnseq2[j].islower():
                sslist.append('0')
            else:
                sslist.append('1')
        return "".join(sslist)
#}}}
def VerifyUnalignedRegion_obsolete(topo1, topo2, unaligned_str):#{{{
# verify unaligned region, for those unaligned region satisfying 
# 1. at most one topology has TM region
# 2. the one without TM region is very short <=10
# 3. or the one without TM regions is shorter than the distance to the first
# unaligned TM of the other TM 
    alignedPosList = myfunc.GetSegPos(unaligned_str, "1")
    lengthAln = len(topo1)
    if len(alignedPosList) != 1:
        print >> sys.stderr, "aligned region not equal 1 for %s %s"%(id1,id2)
        return unaligned_str
    else:
        alignedPos = alignedPosList[0]
    if alignedPos[0] == 0 and alignedPos[1] == lengthAln:
        return unaligned_str
    else:
        topo_Nterm1 = topo1[:alignedPos[0]]
        topo_Cterm1 = topo1[alignedPos[1]:]
        topo_Nterm2 = topo2[:alignedPos[0]]
        topo_Cterm2 = topo2[alignedPos[1]:]
        posTM_Nterm1 = myfunc.GetTMPosition(topo_Nterm1)
        posTM_Cterm1 = myfunc.GetTMPosition(topo_Cterm1)
        posTM_Nterm2 = myfunc.GetTMPosition(topo_Nterm2)
        posTM_Cterm2 = myfunc.GetTMPosition(topo_Cterm2)

        if len(posTM_Nterm1) >0 and len(posTM_Nterm2) > 0:
            s_Nterm = 2
        elif len(posTM_Nterm1) >0 or len(posTM_Nterm2) > 0:
            s_Nterm = 1
        else:
            s_Nterm = 0

        if len(posTM_Cterm1) >0 and len(posTM_Cterm2) > 0:
            s_Cterm = 2
        elif len(posTM_Cterm1) >0 or len(posTM_Cterm2) > 0:
            s_Cterm = 1
        else:
            s_Cterm = 0

        if s_Nterm == 2 or s_Cterm == 2:
            return unaligned_str
        else:
            if len(posTM_Nterm1) > 0:
                num_res_unaligned_Nterm = len(topo_Nterm2.replace("-",""))
                numTM_unaligned_Nterm = len(posTM_Nterm1)
                num_res_to_TM_Nterm = len(topo_Nterm1) - posTM_Nterm1[len(posTM_Nterm1)-1][1]
            elif len(posTM_Nterm2) > 0:
                num_res_unaligned_Nterm = len(topo_Nterm1.replace("-",""))
                numTM_unaligned_Nterm = len(posTM_Nterm2)
                num_res_to_TM_Nterm = len(topo_Nterm2) - posTM_Nterm2[len(posTM_Nterm2)-1][1]
            else:
                num_res_unaligned_Nterm = 0
                numTM_unaligned_Nterm = 0
                num_res_to_TM_Nterm = 0
            if (num_res_unaligned_Nterm <= 10 or
                    myfunc.FloatDivision(num_res_unaligned_Nterm,
                        num_res_to_TM_Nterm)<=0.5):
                new_unaligned_Nterm = "1"*len(topo_Nterm1)
            else:
                new_unaligned_Nterm = "0"*len(topo_Nterm1)


            if len(posTM_Cterm1) > 0:
                num_res_unaligned_Cterm = len(topo_Cterm2.replace("-",""))
                numTM_unaligned_Cterm = len(posTM_Cterm1)
                num_res_to_TM_Cterm = posTM_Cterm1[0][0]
            elif len(posTM_Cterm2) > 0:
                num_res_unaligned_Cterm = len(topo_Cterm1.replace("-",""))
                numTM_unaligned_Cterm = len(posTM_Cterm2)
                num_res_to_TM_Cterm = posTM_Cterm2[0][0]
            else:
                num_res_unaligned_Cterm = 0
                numTM_unaligned_Cterm = 0
                num_res_to_TM_Cterm = 0

            if (num_res_unaligned_Cterm <= 10 or
                    myfunc.FloatDivision(num_res_unaligned_Cterm,
                        num_res_to_TM_Cterm)<=0.5):
                new_unaligned_Cterm = "1"*len(topo_Cterm1)
            else:
                new_unaligned_Cterm = "0"*len(topo_Cterm1)  

            new_unaligned_str = (new_unaligned_Nterm +
                    unaligned_str[alignedPos[0]:alignedPos[1]] +
                    new_unaligned_Cterm)
            return new_unaligned_str
#}}}
def CountCommonM_Aligned(alignedPos, topo1, topo2, helix1, helix2):#{{{
    common_b = max(alignedPos[0], helix1[0], helix2[0])
    common_e = min(alignedPos[1], helix1[1], helix2[1])
    cnt_TM_common = 0
    for i in xrange(common_b, common_e):
        if topo1[i] == "M" and topo2[i] == "M":
            cnt_TM_common += 1
    return cnt_TM_common
#}}}
def CountSEQ_Aligned(alignedPos, topo1, topo2, helix):#{{{
# helix1 is TM 
# helix2 is segment
    b = max(alignedPos[0], helix[0])
    e = min(alignedPos[1], helix[1])
    cntSEQ = 0
    for i in xrange(b, e):
        if topo1[i] == "M" and topo2[i] in ["i", "o"]:
            cntSEQ += 1
    return cntSEQ
#}}}
def CountTM2SEQRes_Aligned(alignedPos, sp_pos2, topo1, topo2, helix):#{{{
# helix1 is TM 
# helix2 is segment
    b = max(alignedPos[0], helix[0])
    e = min(alignedPos[1], helix[1])
    cnt = 0
    for i in xrange(b, e):
        if topo1[i] == "M" and topo2[i] in ["i", "o"] and i < sp_pos2:
            cnt += 1
    return cnt
#}}}
def CountNumberOf_TM2SEQ_residues_at_borders(helix, alignedPos, alingedTopo):#{{{
    b = max(alignedPos[0], helix[0])
    e = min(alignedPos[1], helix[1])
    cnt = 0
    for i in xrange(b,e):
        if alingedTopo[i] not in ["-", "M"]:
            cnt += 1
    return cnt
#}}}
def IsTMAlignDetermined_Nterm(topo1, topo2, posTM1, posTM2, topo_Nterm1,#{{{
        topo_Nterm2, posTM_Nterm1, posTM_Nterm2, alignedPos):
    isStatusDetermiend = False
    cnt_common_M_aln = CountCommonM_Aligned(alignedPos, topo1, topo2,
            posTM1[max(0, len(posTM_Nterm1)-1)], posTM2[max(0,
                len(posTM_Nterm2)-1)])
    if len(posTM_Nterm1) == 0:
        if len(posTM_Nterm2) == 0:
            isStatusDetermiend = True
        else:
            # overlapping of TM region at aligned region
            cnt_seqres_aln = CountSEQ_Aligned(alignedPos, topo1, topo2,
                    posTM2[0])
            num_res_unaligned_Nterm = len(topo_Nterm1.replace("-",""))
            num_res_to_TM_Nterm = len(topo_Nterm2) - posTM_Nterm2[len(posTM_Nterm2)-1][1]
            numTMRes = alignedPos[0] - posTM_Nterm2[len(posTM_Nterm2)-1][0]
            if len(posTM_Nterm2) == 1:
                if (cnt_common_M_aln >=5 
                    or (  numTMRes <= 21 
                        and myfunc.FloatDivision(cnt_seqres_aln +
                            min(num_res_unaligned_Nterm, alignedPos[0] -
                                posTM_Nterm2[len(posTM_Nterm2)-1][0]),
                            21-cnt_common_M_aln) < 0.5) 
                    or ( numTMRes > 21 and
                        (num_res_unaligned_Nterm <= 10 or
                            myfunc.FloatDivision(num_res_unaligned_Nterm,
                                num_res_to_TM_Nterm) <= TH_RATIO))
                ):
                    isStatusDetermiend = True
            else:
                if (( numTMRes <= 21 and myfunc.FloatDivision(cnt_seqres_aln +
                            min(num_res_unaligned_Nterm, alignedPos[0] -
                                posTM_Nterm2[len(posTM_Nterm2)-1][0]),
                            21-cnt_common_M_aln) < 0.5 and num_res_unaligned_Nterm <= 10) 
                    or ( numTMRes > 21 and
                        (num_res_unaligned_Nterm <= 10 or
                            myfunc.FloatDivision(num_res_unaligned_Nterm,
                                num_res_to_TM_Nterm) <= TH_RATIO))
                ):
                    isStatusDetermiend = True
    elif len(posTM_Nterm1) == 1:
        if len(posTM_Nterm2) == 0:
            cnt_seqres_aln = CountSEQ_Aligned(alignedPos, topo1, topo2, posTM1[0])
            num_res_unaligned_Nterm = len(topo_Nterm2.replace("-",""))
            num_res_to_TM_Nterm = len(topo_Nterm1) - posTM_Nterm1[len(posTM_Nterm1)-1][1]
            numTMRes = alignedPos[0] - posTM_Nterm1[len(posTM_Nterm1)-1][0]
            if (cnt_common_M_aln >= 5
                or (  numTMRes <= 21
                    and myfunc.FloatDivision(cnt_seqres_aln +
                        min(num_res_unaligned_Nterm, alignedPos[0] -
                            posTM_Nterm1[len(posTM_Nterm1)-1][0]),
                        21-cnt_common_M_aln) < 0.5) 

                or ( numTMRes > 21 and
                    (num_res_unaligned_Nterm <= 10 or
                        myfunc.FloatDivision(num_res_unaligned_Nterm,
                            num_res_to_TM_Nterm) <= TH_RATIO))
               ):
                isStatusDetermiend = True
        elif len(posTM_Nterm2) == 1:
            if cnt_common_M_aln >= 5:
                isStatusDetermiend = True
        else:
            num_res_unaligned_Nterm = len(topo_Nterm1.replace("-",""))
            numTMRes1 = alignedPos[0] - posTM_Nterm1[len(posTM_Nterm1)-1][0] 
            numTMRes2 = alignedPos[0] - posTM_Nterm2[len(posTM_Nterm2)-1][0] 
            num_res_to_TM_Nterm = len(topo_Nterm2) - posTM_Nterm2[len(posTM_Nterm2)-2][1]
            if (numTMRes1 < 21 and numTMRes2 < 21 and 
                    (num_res_unaligned_Nterm <= 10
                        or  myfunc.FloatDivision(num_res_unaligned_Nterm,
                            num_res_to_TM_Nterm) < TH_RATIO)):
                isStatusDetermiend = True
    else:
        if len(posTM_Nterm2) == 0:
            num_res_unaligned_Nterm = len(topo_Nterm2.replace("-",""))
            num_res_to_TM_Nterm = len(topo_Nterm1) - posTM_Nterm1[len(posTM_Nterm1)-1][1]
            numTMRes = alignedPos[0] - posTM_Nterm1[len(posTM_Nterm1)-1][0] 
            cnt_seqres_aln = CountSEQ_Aligned(alignedPos, topo1, topo2,
                    posTM1[max(0,len(posTM_Nterm1)-1)])
            if ((numTMRes <= 21 and myfunc.FloatDivision(cnt_seqres_aln +
                            min(num_res_unaligned_Nterm, alignedPos[0] -
                                posTM_Nterm1[len(posTM_Nterm1)-1][0]),
                            21-cnt_common_M_aln) < 0.5 and num_res_unaligned_Nterm <= 10) 
                or ( numTMRes > 21 and
                    (num_res_unaligned_Nterm <= 10 or
                        myfunc.FloatDivision(num_res_unaligned_Nterm,
                            num_res_to_TM_Nterm) <= TH_RATIO))
               ):
                isStatusDetermiend = True
        elif len(posTM_Nterm2) == 1:
            num_res_unaligned_Nterm = len(topo_Nterm2.replace("-",""))
            numTMRes1 = alignedPos[0] - posTM_Nterm1[len(posTM_Nterm1)-1][0] 
            numTMRes2 = alignedPos[0] - posTM_Nterm2[len(posTM_Nterm2)-1][0] 
            num_res_to_TM_Nterm = len(topo_Nterm1) - posTM_Nterm1[len(posTM_Nterm1)-2][1]
            if (numTMRes2 < 21 and numTMRes1 < 21 and 
                    (num_res_unaligned_Nterm <= 10
                        or  myfunc.FloatDivision(num_res_unaligned_Nterm,
                            num_res_to_TM_Nterm) < TH_RATIO)):
                isStatusDetermiend = True
    return isStatusDetermiend
#}}}
def IsTMAlignDetermined_Nterm3(topo1, topo2, posTM1, posTM2, topo_Nterm1,#{{{
        topo_Nterm2, posTM_Nterm1, posTM_Nterm2, sp_pos1, sp_pos2, alignedPos):
    """
    This verifycation of unaligned region is only for when mapping topology 1
    to topology 2, therefore, the TM regions of topology 2 is not need to 
    be verified
    """
    isStatusDetermiend = False
    cnt_common_M_aln = CountCommonM_Aligned(alignedPos, topo1, topo2,
            posTM1[max(0, len(posTM_Nterm1)-1)], posTM2[max(0,
                len(posTM_Nterm2)-1)])
    if len(posTM_Nterm1) == 0:
        isStatusDetermiend = True
    elif len(posTM_Nterm1) == 1:
        if len(posTM_Nterm2) == 0:
            cnt_seqres_aln = CountSEQ_Aligned(alignedPos, topo1, topo2, posTM1[0])
            cnt_TM2SPres_aln = CountTM2SEQRes_Aligned(alignedPos, sp_pos2,
                    topo1, topo2, posTM1[0])
            num_res_unaligned_Nterm = len(topo_Nterm2.replace("-",""))
            num_res_to_TM_Nterm = len(topo_Nterm1) - posTM_Nterm1[len(posTM_Nterm1)-1][1]
            numTMRes = alignedPos[0] - posTM_Nterm1[len(posTM_Nterm1)-1][0]
            num_res_to_SP = alignedPos[0] - sp_pos2
            if (cnt_common_M_aln >= 5
                    or ( cnt_TM2SPres_aln + numTMRes - num_res_to_SP > 10
                        )
               ):
                isStatusDetermiend = True
        elif len(posTM_Nterm2) == 1:
            if cnt_common_M_aln >= 5:
                isStatusDetermiend = True
        else: #numTM_Nterm2 >= 2
            num_res_unaligned_Nterm = len(topo_Nterm1.replace("-",""))
            numTMRes1 = alignedPos[0] - posTM_Nterm1[len(posTM_Nterm1)-1][0] 
            numTMRes2 = alignedPos[0] - posTM_Nterm2[len(posTM_Nterm2)-1][0] 
            num_res_to_TM_Nterm = len(topo_Nterm2) - posTM_Nterm2[len(posTM_Nterm2)-2][1]
            if (numTMRes1 < 21 and numTMRes2 < 21 and 
                    (num_res_unaligned_Nterm <= 10
                        or  myfunc.FloatDivision(num_res_unaligned_Nterm,
                            num_res_to_TM_Nterm) < TH_RATIO)):
                isStatusDetermiend = True
    else: #numTM_Nterm1 >= 2
        if len(posTM_Nterm2) == 0:
            num_res_unaligned_Nterm = len(topo_Nterm2.replace("-",""))
            num_res_to_TM_Nterm = len(topo_Nterm1) - posTM_Nterm1[len(posTM_Nterm1)-1][1]
            numTMRes = alignedPos[0] - posTM_Nterm1[len(posTM_Nterm1)-1][0] 
            cnt_seqres_aln = CountSEQ_Aligned(alignedPos, topo1, topo2,
                    posTM1[max(0,len(posTM_Nterm1)-1)])
            num_res_to_SP = alignedPos[0] - sp_pos2
            if (num_res_to_SP <= 10):  # must be TM2SP
                isStatusDetermiend = True
    return isStatusDetermiend
#}}}
def IsTMAlignDetermined_Cterm(topo1, topo2, posTM1, posTM2, topo_Cterm1,#{{{
        topo_Cterm2, posTM_Cterm1, posTM_Cterm2, alignedPos):
    isStatusDetermiend = False
    cnt_common_M_aln = CountCommonM_Aligned(alignedPos, topo1,
            topo2, posTM1[len(posTM1)-max(1, len(posTM_Cterm1))],
            posTM2[len(posTM2)-max(1, len(posTM_Cterm2))])
    if len(posTM_Cterm1) == 0:
        if len(posTM_Cterm2) == 0:
            isStatusDetermiend = True
        else:
            # overlapping of TM region at aligned region
            cnt_seqres_aln = CountSEQ_Aligned(alignedPos, topo1, topo2,
                    posTM2[len(posTM2)-1])
            num_res_unaligned_Cterm = len(topo_Cterm1.replace("-",""))
            num_res_to_TM_Cterm = posTM_Cterm2[0][0]-alignedPos[1]
            numTMRes =  posTM_Cterm2[0][1] - alignedPos[1]
            if len(posTM_Cterm2) == 1:
                if (cnt_common_M_aln >=5 
                    or ( numTMRes <= 21 and myfunc.FloatDivision(cnt_seqres_aln +
                            min(num_res_unaligned_Cterm, numTMRes),
                            21-cnt_common_M_aln) < 0.5) 
                    or ( numTMRes > 21 and
                        (num_res_unaligned_Cterm <= 10 or
                            myfunc.FloatDivision(num_res_unaligned_Cterm,
                                num_res_to_TM_Cterm) <= TH_RATIO))
                ):
                    isStatusDetermiend = True
            else: #numTM2 > 1
                if ((numTMRes <= 21 and myfunc.FloatDivision(cnt_seqres_aln +
                            min(num_res_unaligned_Cterm, numTMRes),
                            21-cnt_common_M_aln) < 0.5 and num_res_unaligned_Cterm <= 10)
                    or ( numTMRes > 21 and
                        (num_res_unaligned_Cterm <= 10 or
                            myfunc.FloatDivision(num_res_unaligned_Cterm,
                                num_res_to_TM_Cterm) <= TH_RATIO))
                ):
                    isStatusDetermiend = True
    elif len(posTM_Cterm1) == 1:
        if len(posTM_Cterm2) == 0:
            cnt_seqres_aln = CountSEQ_Aligned(alignedPos, topo1, topo2,
                    posTM1[len(posTM1)-max(1, len(posTM_Cterm1))])
            num_res_unaligned_Cterm = len(topo_Cterm2.replace("-",""))
            num_res_to_TM_Cterm = posTM_Cterm1[0][0] - alignedPos[1]
            numTMRes = posTM_Cterm1[0][1] - alignedPos[1]
            if (cnt_common_M_aln >=5 
                or ( numTMRes <= 21 
                    and myfunc.FloatDivision(cnt_seqres_aln +
                        min(num_res_unaligned_Cterm, numTMRes),
                        21-cnt_common_M_aln) < 0.5) 
                or ( numTMRes > 21 and
                    (num_res_unaligned_Cterm <= 10 or
                        myfunc.FloatDivision(num_res_unaligned_Cterm,
                            num_res_to_TM_Cterm) <= TH_RATIO))
               ):
                isStatusDetermiend = True
        elif len(posTM_Cterm2) == 1:
            if cnt_common_M_aln >= 5:
                isStatusDetermiend = True
        else:
            num_res_unaligned_Cterm = len(topo_Cterm1.replace("-",""))
            numTMRes1 = posTM_Cterm1[0][1] - alignedPos[1] 
            numTMRes2 = posTM_Cterm2[0][1] - alignedPos[1] 
            num_res_to_TM_Cterm = posTM_Cterm2[1][0] - alignedPos[1]
            if (numTMRes1 < 21 and numTMRes2 < 21 and 
                    (num_res_unaligned_Cterm <= 10
                        or  myfunc.FloatDivision(num_res_unaligned_Cterm,
                            num_res_to_TM_Cterm )<TH_RATIO)):
                isStatusDetermiend = True
    else:
        if len(posTM_Cterm2) == 0:
            num_res_unaligned_Cterm = len(topo_Cterm2.replace("-",""))
            num_res_to_TM_Cterm = posTM_Cterm1[0][0] - alignedPos[1]
            numTMRes = posTM_Cterm1[0][1] - alignedPos[1]
            cnt_seqres_aln = CountSEQ_Aligned(alignedPos, topo1, topo2,
                    posTM1[len(posTM1)-max(1, len(posTM_Cterm1))])
            if ((numTMRes <= 21 and myfunc.FloatDivision(cnt_seqres_aln +
                            min(num_res_unaligned_Cterm, numTMRes),
                            21-cnt_common_M_aln) < 0.5 and num_res_unaligned_Cterm <= 10) 
                or ( numTMRes > 21 and
                    (num_res_unaligned_Cterm <= 10 or
                        myfunc.FloatDivision(num_res_unaligned_Cterm,
                            num_res_to_TM_Cterm) <= TH_RATIO))
               ):
                isStatusDetermiend = True
        elif len(posTM_Cterm2) == 1:
            num_res_unaligned_Cterm = len(topo_Cterm2.replace("-",""))
            numTMRes1 = posTM_Cterm1[0][1] - alignedPos[1] 
            numTMRes2 = posTM_Cterm2[0][1] - alignedPos[1] 
            num_res_to_TM_Cterm = posTM_Cterm1[1][0] - alignedPos[1]
            if (numTMRes2 < 21 and numTMRes1 < 21 and 
                    (num_res_unaligned_Cterm <= 10
                        or myfunc.FloatDivision(num_res_unaligned_Cterm,
                            num_res_to_TM_Cterm)<TH_RATIO)):
                isStatusDetermiend = True
    return isStatusDetermiend
#}}}



def AdjustPosTM(posTM, adjust):#{{{
    newPosTM = []
    for item in posTM:
        newPosTM.append((item[0]+adjust, item[1]+adjust))
    return newPosTM
#}}}
def VerifyUnalignedRegion(topo1, topo2, unaligned_str, id1, id2):#{{{
# updated 2013-07-01
    alignedPosList = myfunc.GetSegPos(unaligned_str, "1")
    lengthAln = len(topo1)
    if len(alignedPosList) != 1:
        print >> sys.stderr, "aligned region not equal 1 for %s %s, unaligned_str=\"%s\""%(id1,id2, unaligned_str)
        return unaligned_str
    else:
        alignedPos = alignedPosList[0]
    if alignedPos[0] == 0 and alignedPos[1] == lengthAln:
        return unaligned_str
    else:
        topo_Nterm1 = topo1[:alignedPos[0]]
        topo_Cterm1 = topo1[alignedPos[1]:]
        topo_Nterm2 = topo2[:alignedPos[0]]
        topo_Cterm2 = topo2[alignedPos[1]:]
        posTM_Nterm1 = myfunc.GetTMPosition(topo_Nterm1)
        posTM_Cterm1 = myfunc.GetTMPosition(topo_Cterm1)
        posTM_Nterm2 = myfunc.GetTMPosition(topo_Nterm2)
        posTM_Cterm2 = myfunc.GetTMPosition(topo_Cterm2)

        posTM_Cterm1 = AdjustPosTM(posTM_Cterm1, alignedPos[1])
        posTM_Cterm2 = AdjustPosTM(posTM_Cterm2, alignedPos[1])

        posTM1 = myfunc.GetTMPosition(topo1)
        posTM2 = myfunc.GetTMPosition(topo2)

        if IsTMAlignDetermined_Nterm(topo1, topo2, posTM1, posTM2, topo_Nterm1,
                topo_Nterm2, posTM_Nterm1, posTM_Nterm2, alignedPos):
            new_unaligned_Nterm = "1"*len(topo_Nterm1)
        else:
            new_unaligned_Nterm = "0"*len(topo_Nterm1)

        if IsTMAlignDetermined_Cterm(topo1, topo2, posTM1, posTM2, topo_Cterm1,
                topo_Cterm2, posTM_Cterm1, posTM_Cterm2, alignedPos):
            new_unaligned_Cterm = "1"*len(topo_Cterm1)
        else:
            new_unaligned_Cterm = "0"*len(topo_Cterm1)

        return  (new_unaligned_Nterm +
                unaligned_str[alignedPos[0]:alignedPos[1]] +
                new_unaligned_Cterm)


#}}}
def VerifyUnalignedRegion3(topo1, topo2, sp_pos1, sp_pos2,  #{{{
        unaligned_str, id1,  id2):
    """
    Change unaligned region to be aligned if the status of TM helix mapping can
    be determined.
    created 2013-10-07, updated 2013-10-07
    """
    alignedPosList = myfunc.GetSegPos(unaligned_str, "1")
    lengthAln = len(topo1)
    if len(alignedPosList) != 1:
        print >> sys.stderr, "aligned region not equal 1 for %s %s, unaligned_str=\"%s\""%(id1,id2, unaligned_str)
        return unaligned_str
    else:
        alignedPos = alignedPosList[0]
    if alignedPos[0] == 0 and alignedPos[1] == lengthAln:
        return unaligned_str
    else:
        topo_Nterm1 = topo1[:alignedPos[0]]
        topo_Cterm1 = topo1[alignedPos[1]:]
        topo_Nterm2 = topo2[:alignedPos[0]]
        topo_Cterm2 = topo2[alignedPos[1]:]
        posTM_Nterm1 = myfunc.GetTMPosition(topo_Nterm1)
        posTM_Cterm1 = myfunc.GetTMPosition(topo_Cterm1)
        posTM_Nterm2 = myfunc.GetTMPosition(topo_Nterm2)
        posTM_Cterm2 = myfunc.GetTMPosition(topo_Cterm2)

        posTM_Cterm1 = AdjustPosTM(posTM_Cterm1, alignedPos[1])
        posTM_Cterm2 = AdjustPosTM(posTM_Cterm2, alignedPos[1])

        posTM1 = myfunc.GetTMPosition(topo1)
        posTM2 = myfunc.GetTMPosition(topo2)


        if ((sp_pos2 == -1 and IsTMAlignDetermined_Nterm(topo1, topo2, posTM1,
            posTM2, topo_Nterm1, topo_Nterm2, posTM_Nterm1, posTM_Nterm2,
            alignedPos)) 
            or (sp_pos2 >= 0 and IsTMAlignDetermined_Nterm3(topo1,
                topo2, posTM1, posTM2, topo_Nterm1, topo_Nterm2, posTM_Nterm1,
                posTM_Nterm2, sp_pos1, sp_pos2, alignedPos))):
            new_unaligned_Nterm = "1"*len(topo_Nterm1)
        else:
            new_unaligned_Nterm = "0"*len(topo_Nterm1)

        if IsTMAlignDetermined_Cterm(topo1, topo2, posTM1, posTM2, topo_Cterm1,
                topo_Cterm2, posTM_Cterm1, posTM_Cterm2, alignedPos):
            new_unaligned_Cterm = "1"*len(topo_Cterm1)
        else:
            new_unaligned_Cterm = "0"*len(topo_Cterm1)

        return  (new_unaligned_Nterm +
                unaligned_str[alignedPos[0]:alignedPos[1]] +
                new_unaligned_Cterm)

#}}}


def MapAlignedTMregionForward( posTM1, posTM2, commonMarray, seqID1 = "", seqID2 = ""):#{{{
    fpLog = g_params['fpLog']
    DEBUG_TMMAPPING = g_params['DEBUG_TMMAPPING']
    numTM1 = len(posTM1)
    numTM2 = len(posTM2)
    ## Annotation: mapToIndexList is index of TMs in topo2 that are selected by
    ## mapping topo1 to topo2
    mapToIndexList=[-1]*numTM1;  
    #lengthOrigTMList=[]
    #for (b2,e2) in posTM2:
        ## Note: lengthOrigTM is always 21 for SCAMPI
        #lengthOrigTMList.append(len(topo2[b2:e2].replace(GAP,''))); 
    for i in xrange (numTM1):
        (b1,e1) = posTM1[i]
        overlapTMList=[]; # list of tuples [(index: fraction),...]
        for j in range (numTM2):
            (b2,e2) = posTM2[j]
            covb = max(b1,b2)
            cove = min(e1,e2)
            numCommomM = 0
            if cove>covb: # coverage is not zero
                numCommomM=sum(commonMarray[covb:cove])
#                 overlapTMList.append((j,
#                     myfunc.FloatDivision(float(numCommomM),
#                         lengthOrigTMList[j])))
                try:
                    freq =  numCommomM/float(21)
                except ZeroDivisionError:
                    freq = 0.0
                overlapTMList.append((j, freq))

            if DEBUG_TMMAPPING and fpLog :
                print >> fpLog, ("%s - %s "%(seqID1, seqID2), "TM1", (b1,e1),
                        "TM2", (b2,e2), "coverage", (covb,cove), "numCommomM=",
                        numCommomM)

        if len(overlapTMList) > 0:
             ## sort in descending order
            overlapTMList.sort(key=lambda tup:tup[1], reverse=True); 
            (idxTM, fraction) = overlapTMList[0]
            mapToIndexList[i] = idxTM
            if DEBUG_TMMAPPING and fpLog != None:
                print >> fpLog, ("%s - %s , TM %d -> %d ," % 
                        (seqID1, seqID2, i, idxTM) + 
                        " overlap fraction = %.3f"% (fraction)) 
    return mapToIndexList

#}}}
def MapAlignedTMregion1(posTM1, posTM2, commonMarray, mapArray1, mapArray2, #{{{
        seqID1 = "", seqID2 = ""):
    numTM1 = len(posTM1)
    numTM2 = len(posTM2)
    mapToIndexList1 = MapAlignedTMregionForward(posTM1,posTM2,
            commonMarray,seqID1, seqID2)
    mapToIndexList2 = MapAlignedTMregionForward(posTM2,posTM1,
            commonMarray,seqID1, seqID2)
    forwardMap1 = set([ (i, mapToIndexList1[i]) for i in range(numTM1) ])
    backMap2 = set([ (mapToIndexList2[i], i) for i in range(numTM2) ])
    commonMap = forwardMap1.intersection(backMap2)
    for (i,j) in commonMap:
        if i >=0 and j >=0:
            mapArray1[i] = _DIRECT_MAPPED
            mapArray2[j] = _DIRECT_MAPPED
    return (mapArray1, mapArray2)
#}}}
def MapAlignedTMregion(posTM1, posTM2, topo2, mapArray1, mapArray2, #{{{
        seqID1="", seqID2=""):
    """
    One-to-One mapping of TM regions of two topologies. The mapped TM regions
    must have overlaps 
    mapArray1 and mapArray2 are initialized to _UN_MAPPED
    Generic algorithm description:      
        For each TM in topo1
            Get TMs in topo2 that overlap with TM
            if numTMWithOverlap <= 0
                status of this TM = UNMAPPED
            elif numTMWithOverlap == 1
                if TM on topo2 is MAPPED
                    status of this TM = UNMAPPED
                else 
                    status of this TM = DIRECT_MAPPED
            else (>1)
                fl1 = Get fraction of overlap region for the TM in topo1
                fl2 = Get fraction of overlap region for TMs in topo2
                avgfl = average(fl1, fl2)
                Sort (avgfl)
                for i in range(numTMWithOverlap)
                    if state(indexsortedTM) == UNMAPPED
                        map this 
                        break
                    i +=1

    Since the length of predicted TM regions by SCAMPI is fixed to 21 (and
    actually the true length of TM regions varies slightly 22.3 (3.3). A
    simpler algorithm can be used. fraction on topo1 is the same as fraction on topo2
    Algorithm 2:
        For each TM in topo1
            Get TMs in topo2 that overlap with TM
            if numTMWithOverlap >= 1
                fl2 = Get fraction of overlap region for TMs in topo2
                Sort (fl2)
                for i in range(numTMWithOverlap)
                    if state(indexsortedTM) == UNMAPPED
                        map this 
                        break
                    i +=1
    """
# @params
# topo2   untrimmed topology of sequence 2
    fpLog = g_params['fpLog']
    DEBUG_TMMAPPING = g_params['DEBUG_TMMAPPING']
    numTM1 = len(posTM1)
    numTM2 = len(posTM2)

    lengthOrigTMList=[]
    for (b2,e2) in posTM2:
        ## this is always 21 for SCAMPI
        lengthOrigTMList.append(len(topo2[b2:e2].replace(GAP,''))); 
    for i in range (numTM1):
        (b1,e1) = posTM1[i]
        overlapTMList=[]; # list of tuples [(index: fraction),...]
        for j in range (numTM2):
            (b2,e2) = posTM2[j]
            covb = max(b1,b2)
            cove = min(e1,e2)
            numCommomM = 0
            if cove>covb: # coverage is not zero
                numCommomM=topo2[covb:cove].count("M")
                overlapTMList.append((j,myfunc.FloatDivision(numCommomM, lengthOrigTMList[j])))
                #overlapTMList.append((j, myfunc.FloatDivision(numCommomM, 21)))

            if DEBUG_TMMAPPING and fpLog != None:
                print >> fpLog, ("%s - %s "%(seqID1, seqID2), "TM1", (b1,e1),
                        "TM2", (b2,e2), "coverage", (covb,cove), "numCommomM=",
                        numCommomM)

        if len(overlapTMList) > 0:
            ## sort in descending order
            overlapTMList.sort(key=lambda tup:tup[1], reverse=True); 
            for (idxTM, fraction) in overlapTMList:
                if mapArray2[idxTM] == _UN_MAPPED:
                    mapArray1[i] = _DIRECT_MAPPED
                    mapArray2[idxTM] = _DIRECT_MAPPED
                    if DEBUG_TMMAPPING and fpLog != None:
                        print >> fpLog, ("%s - %s TM %d of topo1 "%(seqID1,
                            seqID2, i), "indexMaxFrac=", idxTM)
                        print >> fpLog, ("%s - %s TM %d of topo1 "%(seqID1,
                            seqID2, i), "maxFrac %4.3f"%(fraction))
                    break
    return (mapArray1, mapArray2)
#}}}
def MapShiftedTMregion(numTM1, numTM2, posTM1, posTM2, mapArray1, mapArray2,#{{{
        seqID1="", seqID2=""):
    """ map shifted (closest) TM regions"""
# find the closest TM region
# if UN_MAPPED
#     map this 
# else 
#     find the closest TM region on the other side
#     if UN_MAPPED
#         map this 
    fpLog = g_params['fpLog']
    DEBUG_TMMAPPING = g_params['DEBUG_TMMAPPING']
    for i in range (numTM1):
        if mapArray1[i] != _UN_MAPPED: continue
        (b1,e1)=posTM1[i]
        dist=[999999]*numTM2
        for j in range (numTM2):
            #if mapArray2[j] != _UN_MAPPED: continue
            (b2,e2)=posTM2[j]
            dist[j] = -(min(e1,e2)-max(b1,b2))
        minDist=min(dist)
        indexMinDist=dist.index(minDist)
        if mapArray2[indexMinDist] == _UN_MAPPED:
            mapArray1[i]  = _SHIFT_MAPPED
            mapArray2[indexMinDist] = _SHIFT_MAPPED
            if DEBUG_TMMAPPING and fpLog != None:
                fpLog.write("%s - %s shiftMapped %d -> %d"%(seqID1, seqID2, i,
                    indexMinDist))
        else:
            (b2,e2) = posTM2[indexMinDist]
            indexMinDistOnTheOtherSide = -1
# Note that overlap is allowed even in the shifted mapping. One TM helix can
# have overlap with multiple TM helices to the opposite topology sequence. When
# doing the direct map, only the one with the largest overlaping region are
# mapped. Therefore, the rest umapped TM helices may still be encountered in
# the shifted mapping 
            if (b2 < b1 and e2 < e1) or (b1 < b2 and e1 < e2 ):
                if b2<b1 and e2<e1: #TM2 is to the left of TM1
                    indexMinDistOnTheOtherSide=indexMinDist +1; 
                elif b1 <b2 and e1 < e2: #TM2 is to the right of TM1
                    indexMinDistOnTheOtherSide=indexMinDist -1
                if (indexMinDistOnTheOtherSide >=0 
                        and indexMinDistOnTheOtherSide < numTM2 
                        and mapArray2[indexMinDistOnTheOtherSide] == _UN_MAPPED):
                    mapArray1[i]  = _SHIFT_MAPPED
                    mapArray2[indexMinDistOnTheOtherSide] = _SHIFT_MAPPED
                    if DEBUG_TMMAPPING and fpLog != None:
                        fpLog.write("%s - %s shiftMapped %d -> %d"%(seqID1,
                            seqID2, i, indexMinDist))
            elif b2 <= b1 and e2 >= e1: #search in both sides
                if DEBUG_TMMAPPING and fpLog != None:
                    fpLog.write("MapShiftedTMregion %s: %d(%d, %d) is contained by %s: %d(%d, %d), but TM2 is mapped by other TMs, this one will remain unmapped\n" %(seqID1, i, b1,e1, seqID2, indexMinDist, b2,e2))
                    PrintMappedArray(mapArray1, mapArray2, seqID1, seqID2, fpLog)
            else:
                sys.stderr.write("Error in MapShiftedTMregion %s: %d(%d, %d) covers entirely %s: %d(%d, %d), but the second TM is mapped by other TMs. This is a bug\n" %(seqID1, i, b1,e1, seqID2, indexMinDist, b2,e2))
                PrintMappedArray(mapArray1, mapArray2, seqID1, seqID2, sys.stderr)
    return (mapArray1, mapArray2)
#}}}
def MapAlignedTopo_method1_obsolete1(topo1, topo2, posTM1, posTM2,   #{{{
        unaligned_str, seqID1, seqID2):
    # unaligned_str is a string of 0 and 1s for unaligned regions
    # e.g. 0000011100000111000, indicates those 0s as unaligned regions
# updated 2013-05-27 
# when the aligned segment contains multiple TM helices, e.g. 
# MMMMMMMMMMMMMM
# MMMooooMMMMMMM
# take the one with the same insertion direction
    mapArray = []
    threshold_TM2TM = g_params['threshold_TM2TM']
    threshold_TM2GAP = g_params['threshold_TM2GAP']
    isLocalAlignment = False
    if unaligned_str != "":
        isLocalAlignment = True

    if len(posTM1) > 0:
        for (b, e) in posTM1:
            segList2 = []
            cntTM = 0
            cntGap = 0
            cntSeq = 0
            cntUnaligned = 0
            cntAligned = 0
            for j in xrange(b,e):
                if topo1[j] == 'M':
                    segList2.append(topo2[j])
                    if isLocalAlignment and unaligned_str[j] == '0':
                        cntUnaligned += 1
                    else:
                        cntAligned += 1
                        if topo2[j] == 'M':
                            cntTM += 1
                        elif topo2[j] == GAP:
                            cntGap += 1
                        else:
                            cntSeq += 1
            sizeSeg = len(segList2)
            try:
                freqTM = cntTM / float(cntAligned)
                freqGap = cntGap / float(cntAligned)
                freqSeq = cntSeq / float(cntAligned)
                freqUnaligned = cntUnaligned/float(sizeSeg)
            except ZeroDivisionError:
                freqTM = 0.0
                freqGap = 0.0
                freqSeq = 0.0
                freqUnaligned = 0.0

            if cntUnaligned > 5:
                print "cntUnaligned = %2d / %2d = %5.1f %21s #M=%d #GAP=%d #SEQ=%d"%(
                        cntUnaligned,
                        sizeSeg, freqUnaligned*100, "".join(segList2),
                        segList2.count('M'), segList2.count('-'),
                        sizeSeg-(segList2.count('M')+segList2.count('-')))
            #if cntUnaligned / float(sizeSeg) >= 0.5:
            if freqUnaligned >= 1- threshold_TM2TM:
                cls = _UNALIGNED
            elif freqTM >= threshold_TM2TM:
                cls = _TM2TM
            elif freqGap >= threshold_TM2GAP:
                cls = _TM2GAP
            else:
                cls = _TM2SEQ
            mapArray.append(cls)
    return mapArray
#}}}
def MapAlignedTopo_method1_obsolete2(topo1, topo2, posTM1, posTM2,   #{{{
        unaligned_str, seqID1, seqID2):
    # unaligned_str is a string of 0 and 1s for unaligned regions
    # e.g. 0000011100000111000, indicates those 0s as unaligned regions
# updated 2013-05-27 
# when the aligned segment contains multiple TM helices, e.g. 
# MMMMMMMMMMMMMM
# MMMooooMMMMMMM
# take the one with the same insertion direction
    mapArray = []
    threshold_TM2TM = g_params['threshold_TM2TM']
    threshold_TM2GAP = g_params['threshold_TM2GAP']
    isLocalAlignment = False
    if unaligned_str != "":
        isLocalAlignment = True

    if len(posTM1) > 0:
        for (b, e) in posTM1:
            segList2 = []
            aligned_segList2 = []
            cntTM = 0
            cntGap = 0
            cntSeq = 0
            cntUnaligned = 0
            cntAligned = 0
            begin_aligned_segment = -1
            for j in xrange(b,e):
                if topo1[j] == 'M':
                    segList2.append(topo2[j])
                    if isLocalAlignment and unaligned_str[j] == '0':
                        cntUnaligned += 1
                    else:
                        cntAligned += 1
                        aligned_segList2.append(topo2[j])
                        if begin_aligned_segment == -1:
                            begin_aligned_segment = j
                        if topo2[j] == 'M':
                            cntTM += 1
                        elif topo2[j] == GAP:
                            cntGap += 1
                        else:
                            cntSeq += 1
            sizeSeg = len(segList2)
            freqTM = myfunc.FloatDivision(cntTM, cntAligned)
            freqGap = myfunc.FloatDivision(cntGap, cntAligned)
            freqSeq = myfunc.FloatDivision(cntSeq , cntAligned)
            freqUnaligned = myfunc.FloatDivision(cntUnaligned, sizeSeg)
            if cntUnaligned > 5 and 0:
                print "cntUnaligned = %2d / %2d = %5.1f %21s #M=%d #GAP=%d #SEQ=%d"%(
                        cntUnaligned,
                        sizeSeg, freqUnaligned*100, "".join(segList2),
                        segList2.count('M'), segList2.count('-'),
                        sizeSeg-(segList2.count('M')+segList2.count('-')))
            #if cntUnaligned / float(sizeSeg) >= 0.5:
            str_segment = "".join(aligned_segList2)
            seg_posTM = myfunc.GetTMPosition(str_segment)
            method_multi_helix_overlap = 1
# 
# more than one TM helices in the segment, 
# method_multi_helix_overlap = 0 take the one with the same insertion direction 
# method_multi_helix_overlap = 1 take the one with the longest overlap
            if len(seg_posTM) > 1:
                # more than one TM helices in the segment, take the one with
                # the same insertion direction 
                if method_multi_helix_overlap == 0:
                    Nterm_thisTM = Get_IOState_upstream(topo1, b)
                    tmp_Nterm = Get_IOState_upstream(topo2,
                            seg_posTM[0][0]+begin_aligned_segment)
                    if Nterm_thisTM == tmp_Nterm:
                        new_str_segment = str_segment[0:(seg_posTM[1][0])]
                    else:
                        new_str_segment = str_segment[seg_posTM[0][1]:]

                    print "SAME_DIRECTION: %s %s %21s %21s %5d"%(seqID1,
                            seqID2, str_segment, new_str_segment, cntUnaligned)
                else:
                    cntM1 = str_segment[seg_posTM[0][0]:seg_posTM[0][1]].count("M")
                    cntM2 = str_segment[seg_posTM[1][0]:seg_posTM[1][1]].count("M")
                    if cntM1 > cntM2:
                        new_str_segment = str_segment[0:(seg_posTM[1][0])]
                    elif cntM2 > cntM1:
                        new_str_segment = str_segment[seg_posTM[0][1]:]
                    else:
                        Nterm_thisTM = Get_IOState_upstream(topo1, b)
                        tmp_Nterm = Get_IOState_upstream(topo2,
                                seg_posTM[0][0]+begin_aligned_segment)
                        if Nterm_thisTM == tmp_Nterm:
                            new_str_segment = str_segment[0:(seg_posTM[1][0])]
                        else:
                            new_str_segment = str_segment[seg_posTM[0][1]:]
                    print "MAX: %s %s %21s %21s %5d"%(seqID1, seqID2,
                            str_segment, new_str_segment, cntUnaligned)

                (cntTM, cntGap, cntSeq, freqTM, freqGap, freqSeq) = StatIOMFreq(
                        new_str_segment)

            if freqUnaligned >= 1- threshold_TM2TM:
                cls = _UNALIGNED
            elif freqTM >= threshold_TM2TM:
                cls = _TM2TM
            elif freqGap >= threshold_TM2GAP:
                cls = _TM2GAP
            else:
                cls = _TM2SEQ
            mapArray.append(cls)
    return mapArray
#}}}
def MapAlignedTopo_method1(topo1, topo2, posTM1, posTM2,   #{{{
        unaligned_str, seqID1, seqID2):
    # unaligned_str is a string of 0 and 1s for unaligned regions
    # e.g. 0000011100000111000, indicates those 0s as unaligned regions
# updated 2013-05-27 
# when the aligned segment contains multiple TM helices, e.g. 
# MMMMMMMMMMMMMM
# MMMooooMMMMMMM
# take the one with the same insertion direction
    mapArray = []
#     threshold_TM2TM = g_params['threshold_TM2TM']
#     threshold_TM2GAP = g_params['threshold_TM2GAP']
    min_TM_overlap = g_params['min_TM_overlap']
    isLocalAlignment = False
    if unaligned_str != "":
        isLocalAlignment = True
        unaligned_str = VerifyUnalignedRegion(topo1, topo2, unaligned_str, seqID1,
            seqID2)
#       unaligned_str = VerifyUnalignedRegion_obsolete(topo1, topo2, unaligned_str)

    if len(posTM1) > 0:
        for (b, e) in posTM1:
            segList2 = []
            aligned_segList2 = []
            cntTM = 0
            cntGap = 0
            cntSeq = 0
            cntUnaligned = 0
            cntAligned = 0
            begin_aligned_segment = -1
            for j in xrange(b,e):
                if topo1[j] == 'M':
                    segList2.append(topo2[j])
                    if isLocalAlignment and unaligned_str[j] == '0':
                        cntUnaligned += 1
                    else:
                        cntAligned += 1
                        aligned_segList2.append(topo2[j])
                        if begin_aligned_segment == -1:
                            begin_aligned_segment = j
                        if topo2[j] == 'M':
                            cntTM += 1
                        elif topo2[j] == GAP:
                            cntGap += 1
                        else:
                            cntSeq += 1
            sizeSeg = len(segList2)
            freqTM = myfunc.FloatDivision(cntTM, cntAligned)
            freqGap = myfunc.FloatDivision(cntGap, cntAligned)
            freqSeq = myfunc.FloatDivision(cntSeq , cntAligned)
            freqUnaligned = myfunc.FloatDivision(cntUnaligned, sizeSeg)
            if cntUnaligned > 5 and 0:
                print "cntUnaligned = %2d / %2d = %5.1f %21s #M=%d #GAP=%d #SEQ=%d"%(
                        cntUnaligned,
                        sizeSeg, freqUnaligned*100, "".join(segList2),
                        segList2.count('M'), segList2.count('-'),
                        sizeSeg-(segList2.count('M')+segList2.count('-')))
            #if cntUnaligned / float(sizeSeg) >= 0.5:
            str_segment = "".join(aligned_segList2)
            seg_posTM = myfunc.GetTMPosition(str_segment)
            #method_multi_helix_overlap = 1
            method_multi_helix_overlap = 0
# 
# more than one TM helices in the segment, 
# method_multi_helix_overlap = 0 take the one with the same insertion direction 
# method_multi_helix_overlap = 1 take the one with the longest overlap
            if len(seg_posTM) > 1:
                # more than one TM helices in the segment, take the one with
                # the same insertion direction 
                if method_multi_helix_overlap == 0:
                    Nterm_thisTM = Get_IOState_upstream(topo1, b)
                    tmp_Nterm = Get_IOState_upstream(topo2,
                            seg_posTM[0][0]+begin_aligned_segment)
                    if Nterm_thisTM == tmp_Nterm:
                        new_str_segment = str_segment[0:(seg_posTM[1][0])]
                    else:
                        new_str_segment = str_segment[seg_posTM[0][1]:]

                    print "SAME_DIRECTION: %s %s %21s %21s %5d"%(seqID1,
                            seqID2, str_segment, new_str_segment, cntUnaligned)
                else:
                    cntM1 = str_segment[seg_posTM[0][0]:seg_posTM[0][1]].count("M")
                    cntM2 = str_segment[seg_posTM[1][0]:seg_posTM[1][1]].count("M")
                    if cntM1 > cntM2:
                        new_str_segment = str_segment[0:(seg_posTM[1][0])]
                    elif cntM2 > cntM1:
                        new_str_segment = str_segment[seg_posTM[0][1]:]
                    else:
                        Nterm_thisTM = Get_IOState_upstream(topo1, b)
                        tmp_Nterm = Get_IOState_upstream(topo2,
                                seg_posTM[0][0]+begin_aligned_segment)
                        if Nterm_thisTM == tmp_Nterm:
                            new_str_segment = str_segment[0:(seg_posTM[1][0])]
                        else:
                            new_str_segment = str_segment[seg_posTM[0][1]:]
                    print "MAX: %s %s %21s %21s %5d"%(seqID1, seqID2,
                            str_segment, new_str_segment, cntUnaligned)
                (cntTM, cntGap, cntSeq, freqTM, freqGap, freqSeq) = StatIOMFreq(
                        new_str_segment)
            if cntUnaligned > sizeSeg - min_TM_overlap:
            #if cntUnaligned > 11:
            #if freqUnaligned > 0.8:
                cls = _UNALIGNED
            elif cntTM >= min_TM_overlap:
                cls = _TM2TM
            else:
                if freqGap >= freqSeq:
                    cls = _TM2GAP
                else:
                    cls = _TM2SEQ
            mapArray.append(cls)
    return mapArray
#}}}
def MapAlignedTopo_method3(topo1, topo2, posTM1, posTM2,   #{{{
        sp_pos1, sp_pos2, unaligned_str, seqID1, seqID2):
    """
    unaligned_str is a string of 0 and 1s for unaligned regions
    e.g. 0000011100000111000, indicates those 0s as unaligned regions
    when the aligned segment contains multiple TM helices, e.g. 
    MMMMMMMMMMMMMM
    MMMooooMMMMMMM
    take the one with the same insertion direction
    created 2013-10-07
    """
    mapArray = []
    min_TM_overlap = g_params['min_TM_overlap']
    isLocalAlignment = False
    if unaligned_str != "":
        isLocalAlignment = True
        unaligned_str = VerifyUnalignedRegion3(topo1, topo2, sp_pos1, sp_pos2,
                unaligned_str, seqID1, seqID2)
#       unaligned_str = VerifyUnalignedRegion_obsolete(topo1, topo2, unaligned_str)

    for iTM in xrange(len(posTM1)):
        (b,e) = posTM1[iTM]
        segList2 = []
        aligned_segList2 = []
        cntTM = 0
        cntSP = 0
        cntGap = 0
        cntSeq = 0
        cntUnaligned = 0
        cntAligned = 0
        begin_aligned_segment = -1
        for j in xrange(b,e):
            if topo1[j] == 'M':
                segList2.append(topo2[j])
                if isLocalAlignment and unaligned_str[j] == '0':
                    cntUnaligned += 1
                else:
                    cntAligned += 1
                    aligned_segList2.append(topo2[j])
                    if begin_aligned_segment == -1:
                        begin_aligned_segment = j
                    if j <= sp_pos2:
                        if topo2[j] != GAP:
                            cntSP += 1
                    elif topo2[j] == 'M':
                        cntTM += 1
                    elif topo2[j] == GAP:
                        cntGap += 1
                    else:
                        cntSeq += 1
        sizeSeg = len(segList2)
        freqTM = myfunc.FloatDivision(cntTM, cntAligned)
        freqSP = myfunc.FloatDivision(cntSP, cntAligned)
        freqGap = myfunc.FloatDivision(cntGap, cntAligned)
        freqSeq = myfunc.FloatDivision(cntSeq , cntAligned)
        freqUnaligned = myfunc.FloatDivision(cntUnaligned, sizeSeg)
        if cntUnaligned > 5 and 0:
            print "cntUnaligned = %2d / %2d = %5.1f %21s #M=%d #GAP=%d #SEQ=%d"%(
                    cntUnaligned,
                    sizeSeg, freqUnaligned*100, "".join(segList2),
                    segList2.count('M'), segList2.count('-'),
                    sizeSeg-(segList2.count('M')+segList2.count('-')))
        #if cntUnaligned / float(sizeSeg) >= 0.5:
        str_segment = "".join(aligned_segList2)
        seg_posTM = myfunc.GetTMPosition(str_segment)
        #method_multi_helix_overlap = 1
        method_multi_helix_overlap = 0
# 
# more than one TM helices in the segment, 
# method_multi_helix_overlap = 0 take the one with the same insertion direction 
# method_multi_helix_overlap = 1 take the one with the longest overlap
        if len(seg_posTM) > 1:
            # more than one TM helices in the segment, take the one with
            # the same insertion direction 
            if method_multi_helix_overlap == 0:
                Nterm_thisTM = Get_IOState_upstream(topo1, b)
                tmp_Nterm = Get_IOState_upstream(topo2,
                        seg_posTM[0][0]+begin_aligned_segment)
                if Nterm_thisTM == tmp_Nterm:
                    new_str_segment = str_segment[0:(seg_posTM[1][0])]
                else:
                    new_str_segment = str_segment[seg_posTM[0][1]:]

                print "SAME_DIRECTION: %s %s %21s %21s %5d"%(seqID1,
                        seqID2, str_segment, new_str_segment, cntUnaligned)
            else:
                cntM1 = str_segment[seg_posTM[0][0]:seg_posTM[0][1]].count("M")
                cntM2 = str_segment[seg_posTM[1][0]:seg_posTM[1][1]].count("M")
                if cntM1 > cntM2:
                    new_str_segment = str_segment[0:(seg_posTM[1][0])]
                elif cntM2 > cntM1:
                    new_str_segment = str_segment[seg_posTM[0][1]:]
                else:
                    Nterm_thisTM = Get_IOState_upstream(topo1, b)
                    tmp_Nterm = Get_IOState_upstream(topo2,
                            seg_posTM[0][0]+begin_aligned_segment)
                    if Nterm_thisTM == tmp_Nterm:
                        new_str_segment = str_segment[0:(seg_posTM[1][0])]
                    else:
                        new_str_segment = str_segment[seg_posTM[0][1]:]
                print "MAX: %s %s %21s %21s %5d"%(seqID1, seqID2,
                        str_segment, new_str_segment, cntUnaligned)
            (cntTM, cntGap, cntSeq, freqTM, freqGap, freqSeq) = StatIOMFreq(
                    new_str_segment)
        #if cntUnaligned > 11:
        #if freqUnaligned > 0.8:
        if cntUnaligned > sizeSeg - min_TM_overlap:
            cls = _UNALIGNED
        elif cntTM >= min_TM_overlap:
            cls = _TM2TM
        #elif iTM == 0 and cntSP > 10:
        elif iTM == 0 and freqSP > freqGap and freqSP > freqSeq:
            cls = _TM2SP
        elif freqGap >= freqSeq:
            cls = _TM2GAP
        else:
            cls = _TM2SEQ
# debuging
        if sp_pos2 != -1 and b < sp_pos2:
            print "TM - SP (%s - %s), cls=%d, cntSP=%d"%(seqID1, seqID2, cls, cntSP)
            tmpli = []
            for jj in range(b,e):
                if jj < sp_pos2 and topo2[jj] != GAP:
                    tmpli.append("P")
                else:
                    tmpli.append(topo2[jj])
            print "%s: %s"%(seqID1, topo1[b:e])
            print "%s: %s"%(seqID2, "".join(tmpli))
# debuging
        mapArray.append(cls)
    return mapArray
#}}}
def MapAlignedSeqIndex(sp, aligned_topo):
# Created 2013-09-04, updated 2013-09-04, Nanjiang Shu  
# given a sequence position, get the position in the aligned seq
    j = 0
    cnt = 0
    while cnt < sp:
        if aligned_topo[j] != GAP:
            cnt += 1
        j += 1
    return j

def MapAlignedSP(topo1, topo2, posTM1, posTM2,   #{{{
        mapArray1, mapArray2, sp1, sp2, unaligned_str, seqID1, seqID2):
# Created 2013-09-05, updated 2013-09-04, Nanjiang Shu 
# unaligned_str is a string of 0 and 1s for unaligned regions
# e.g. 0000011100000111000, indicates those 0s as unaligned regions
# map signal peptide to aligned topology, classified as 
# SP2TM, SP2GAP, SP2SEQ and unalignedSP

    if unaligned_str != "":
        isLocalAlignment = True
        alignedPosList = myfunc.GetSegPos(unaligned_str, "1")
        alignedPos = alignedPosList[0]
    else:
        isLocalAlignment = False

    aligned_sp1 = MapAlignedSeqIndex(sp1, topo1)
    aligned_sp2 = MapAlignedSeqIndex(sp2, topo2)

    segList2 = []
    aligned_segList2 = []
    cntSP = 0
    cntTM = 0
    cntGap = 0
    cntSeq = 0
    cntUnaligned = 0
    cntAligned = 0
    begin_aligned_segment = -1
    for j in xrange(0, aligned_sp1):
        segList2.append(topo2[j])
        if isLocalAlignment and unaligned_str[j] == '0':
            cntUnaligned += 1
        else:
            cntAligned += 1
            aligned_segList2.append(topo2[j])
            if begin_aligned_segment == -1:
                begin_aligned_segment = j
            if j < aligned_sp2 and topo2[j] != GAP:
                cntSP += 1

            if topo2[j] == 'M':
                cntTM += 1
            elif topo2[j] == GAP:
                cntGap += 1
            else:
                cntSeq += 1
    sizeSeg = len(segList2)
    freqSP = myfunc.FloatDivision(cntSP, cntAligned)
    freqTM = myfunc.FloatDivision(cntTM, cntAligned)
    freqGap = myfunc.FloatDivision(cntGap, cntAligned)
    freqSeq = myfunc.FloatDivision(cntSeq , cntAligned)
    freqUnaligned = myfunc.FloatDivision(cntUnaligned, sizeSeg)
    if cntAligned < 1:
        unalignedNterm1 = topo1[:alignedPos[0]]
        unalignedNterm2 = topo2[:alignedPos[0]]
        if unalignedNterm2.count('M') >= 5 and mapArray2[0] != _TM2TM:
            cls = "SP2TM"
        elif sp2 > 0:
            cls = "SP2SP"
        elif (len(unalignedNterm2) - unalignedNterm2.count(GAP)) < sp1/2.0:
            cls = "SP2GAP"
        elif unalignedNterm2.count(GAP) < sp1/2.0:
            cls = "SP2SEQ"
        else:
            cls = "unalignedSP"

    elif cntTM >= 5 and mapArray2[0] != _TM2TM:
        cls = "SP2TM"
    elif cntSP > 0:
        cls = "SP2SP"
    else:
        if freqGap >= freqSeq:
            cls = "SP2GAP"
        else:
            cls = "SP2SEQ"
    return cls
#}}}

def GetGapSequenceConsensus(origTopoMSASeqList, indexIDTTopo):#{{{
    """Get a sequence of gap fractions of the consensus topology """
    lengthAlignment = len(origTopoMSASeqList[0]); 
    cntGapList=[0] * lengthAlignment
    numIDTTopo = len(indexIDTTopo)
    for i in xrange(0, numIDTTopo):
        topo=origTopoMSASeqList[indexIDTTopo[i]]
        j=0
        for s in topo:
            if s == '-':
                cntGapList[j] += 1
            j+=1
    gapSeqCons=[]
    for j in xrange(0,lengthAlignment):
        try:
            freq = cntGapList[j]/float(numIDTTopo)
        except ZeroDivisionError:
            freq = 0.0
        gapSeqCons.append(freq)
    return gapSeqCons
#}}}
def GetGapSequenceConsensus_old(origTopoMSASeqList, indexIDTTopo):#{{{
    """Get a sequence of gap fractions of the consensus topology """
    gapSeqCons =[]; 
    apd=gapSeqCons.append
    lengthAlignment = len(origTopoMSASeqList[0]); 
    numIDTTopo = len(indexIDTTopo)
    for i in xrange(0, lengthAlignment):
        cntGap=0
        for j in xrange(0, numIDTTopo):
            if origTopoMSASeqList[indexIDTTopo[j]][i] == '-':
                cntGap += 1
        try:
            freq =  cntGap/float(numIDTTopo)
        except ZeroDivisionError:
            freq = 0
        apd(freq)
    return gapSeqCons
#}}}
def GetDGvalueTMconsensus(dgScoreList, indexIDTTopo, numTM_IDT):#{{{
    """Get the average DGvalue of each aligned TM at the IDTgroup"""
    DGvalueTMcons=[INIT_DGVALUE]*numTM_IDT; 
    numIDTTopo = len(indexIDTTopo)
    sumDGList = [0.0] * numTM_IDT
    cntTopoWithValidDGScore = 0
    for i in xrange(numIDTTopo):
        dglist = dgScoreList[indexIDTTopo[i]]
        if dglist and len(dglist) == numTM_IDT:
            cntTopoWithValidDGScore += 1
            for j in xrange(numTM_IDT):
                sumDGList[j] += dgScoreList[indexIDTTopo[i]][j]
    if cntTopoWithValidDGScore > 0:
        DGvalueTMcons = [x/float(cntTopoWithValidDGScore) for x in sumDGList]
    return DGvalueTMcons
#}}}

def MappingTM(posTMcons, posTMquery, commonMarray, seqID1="", seqID2=""):#{{{
    """
    Map TM helices of the consensus topology to the query topology
    Return (indexMappedTMcons, indexMappedTMquery)
    @param
        posTMcons  Location of TM helices of the consensusTopo
        posTMquery Location of TM helices of the trimmed query topology
        topoquery  Untrimmed (with gaps) aligned topology of the query
    """
# Mapping algorithm description
# 1. find the mapping of TM regions with the biggest overlaps 
# for the rest 
# if n1 < n2 
#    start from topo1
#    find the mapping to the closest TM helices, but if both the closest
#    right and left TM regions of the targeting topology have been
#    assigned, then this TM helix is left unmapped
# else 
#    start from topo2
#    do the same thing as for the case n1 < n2
    fpLog = g_params['fpLog']
    DEBUG_TMMAPPING = g_params['DEBUG_TMMAPPING']
    numTMcons=len(posTMcons)
    numTMquery=len(posTMquery)
    #init map array 
    mapArraycons=[_UN_MAPPED]*numTMcons
    mapArrayquery=[_UN_MAPPED]*numTMquery

# Step 1: direct map
#   map TMs that have common region, if there are more than one have common
#   regions, take the one with the largest common region. The fraction is
#   according to the orignal msa, not the trimmed msa
    (mapArraycons, mapArrayquery) = MapAlignedTMregion1(posTMcons,
            posTMquery,commonMarray, mapArraycons, mapArrayquery, seqID1,
            seqID2)
    if DEBUG_TMMAPPING and fpLog :
        fpLog.write("After Directmap:\n")
        PrintMappedArray(mapArraycons, mapArrayquery, seqID1, seqID2, fpLog)

# Step 2, map shifted (closest) TM regions
#
# find the closest TM region
# if UN_MAPPED
#     map this 
# else 
#     find the closest TM region on the other side
#     if UN_MAPPED
#         map this 
# 
# this means if the closest TM regions on both sides are already mapped, this
# TM region of the query will be remained unmapped start mapping the one with
# smaller number of TM regions and then map the rest in the reversed way
    if numTMcons <= numTMquery:
        (mapArraycons, mapArrayquery) = MapShiftedTMregion(numTMcons,
                numTMquery, posTMcons, posTMquery, mapArraycons, mapArrayquery,
                seqID1, seqID2)
        (mapArrayquery, mapArraycons) = MapShiftedTMregion(numTMquery,
                numTMcons, posTMquery, posTMcons, mapArrayquery,
                mapArraycons,seqID2, seqID1)
    else:
        (mapArrayquery, mapArraycons) = MapShiftedTMregion(numTMquery,
                numTMcons, posTMquery, posTMcons, mapArrayquery, mapArraycons,
                seqID2, seqID1)
        (mapArraycons, mapArrayquery) = MapShiftedTMregion(numTMcons,
                numTMquery, posTMcons, posTMquery, mapArraycons, mapArrayquery,
                seqID1, seqID2)
    return (mapArraycons, mapArrayquery)
#}}}
def MappingTM_method1(posTM1, posTM2, topo1, topo2, seqID1, seqID2):#{{{
    """
    Map TM helices of the topo1 to topo2 
    Return (mapArray1, mapArray2)
    @param
        topo1  Topology of seq1
        topo2  Topology of seq2
        posTM1 Location of TM helices of topo1
        posTM2 Location of TM helices of topo2
    Note that topology is untrimmed (with gaps) aligned topology
    """
# Mapping algorithm description
# For TM in TMList1:
#   classify its mapping status by the composition of the gapless segment of
#   in the aligned seq, into
#   if freqTM >=  1/3 (one third):
#       class = TM2TM
#   elif freqGap >= 1/2:
#       class = TM2Gap
#   else: 
#       class = TM2Seq
# updated 2013-05-27 
# when the aligned segment contains multiple TM helices, e.g. 
# MMMMMMMMMMMMMM
# MMMooooMMMMMMM
# take the one with the same insertion direction

    fpLog = g_params['fpLog']
    DEBUG_TMMAPPING = g_params['DEBUG_TMMAPPING']
    unaligned_str = ""
    if g_params['localseqpairDict'] != {}:
        try:
            rd = g_params['localseqpairDict'][(seqID1, seqID2)]
            unaligned_str = rd[2]
        except KeyError:
            msg = "Failed to find local alignment for (%s,%s)"
            print >> sys.stderr, msg%(seqID1, seqID2)
            unaligned_str = ""

    mapArray1 = MapAlignedTopo_method1(topo1, topo2, posTM1, posTM2,
            unaligned_str, seqID1, seqID2)
    mapArray2 = MapAlignedTopo_method1(topo2, topo1, posTM2, posTM1,
            unaligned_str, seqID1, seqID2)
    return (mapArray1, mapArray2)
#}}}
def MappingTM_method3(posTM1, posTM2, topo1, topo2, sp_pos1, sp_pos2, #{{{
        seqID1, seqID2):
    """
    Map TM helices of the topo1 to topo2 
    Return (mapArray1, mapArray2)
    @param
        topo1    Topology of seq1
        topo2    Topology of seq2
        posTM1   Location of TM helices of topo1
        posTM2   Location of TM helices of topo2
        sp_pos1  position of the signal peptide for protein 1
        sp_pos2  position of the signal peptide for protein 2
    created 2013-10-07
    """
# Mapping algorithm description
# For TM in TMList1:
#   classify its mapping status by the composition of the gapless segment of
#   in the aligned seq, into
#   if commonM >=  5 
#       class = TM2TM
#   elif isFirstTM and shareWithSignalP > 10:
#       class = TM2SP   #print all TM2SP to check
#   elif freqGAP >= freqSEQ:
#       class = TM2GAP
#   else: 
#       class = TM2SEQ
# updated 2013-05-27 
# when the aligned segment contains multiple TM helices, e.g. 
# MMMMMMMMMMMMMM
# MMMooooMMMMMMM
# take the one with the same insertion direction

    fpLog = g_params['fpLog']
    DEBUG_TMMAPPING = g_params['DEBUG_TMMAPPING']
    unaligned_str = ""
    if g_params['localseqpairDict'] != {}:
        try:
            rd = g_params['localseqpairDict'][(seqID1, seqID2)]
            unaligned_str = rd[2]
        except KeyError:
            msg = "Failed to find local alignment for (%s,%s)"
            print >> sys.stderr, msg%(seqID1, seqID2)
            unaligned_str = ""

    mapArray1 = MapAlignedTopo_method3(topo1, topo2, posTM1, posTM2,
            sp_pos1, sp_pos2, unaligned_str, seqID1, seqID2)
    mapArray2 = MapAlignedTopo_method3(topo2, topo1, posTM2, posTM1,
            sp_pos2, sp_pos1, unaligned_str, seqID1, seqID2)
    return (mapArray1, mapArray2)
#}}}

def IsDGScoreSimilar(dgScore1, dgScore2, maxDGdifference): #{{{
    if len(dgScore1) != len(dgScore2):
        return False
    else:
        for i in range(len(dgScore1)):
            if abs(dgScore1[i] - dgScore2[i]) > maxDGdifference:
                return False
    return True
#}}}

def IsIdenticalTopology_old(Nterm1, Nterm2, numTM1, numTM2, posTM1, posTM2):#{{{
    """check whether topo1 and topo2 are identical"""
# if has overlap 
    #numTM1=len(posTM1)
    #numTM2=len(posTM2)
    if numTM1 != numTM2:
        return False
    else:
        if Nterm1 != Nterm2:
            return False
        else:
            for i in range (numTM1): 
                (b1,e1)=posTM1[i]
                (b2,e2)=posTM2[i]
                if e2<b1 or e1<b2:
                    return False
    return True
#}}}

def CompareTrimmedToposGlobally(strTop1, strTop2, strProtein1, strProtein2):#{{{
    if len(strTop1) <= 0 and len(strTop2) <=0:
        return ("DIFF", 0,0 )
    elif len(strTop1)*len(strTop2) == 0 and len(strTop1)+len(strTop2) > 0:
        print >> sys.stderr, ("%s %s global length does not match" %
                (strProtein1, strProtein2))
        sys.exit(1)
    (intNumMem1,Nterm1)=ct.counttopo(strTop1)
    (intNumMem2,Nterm2)=ct.counttopo(strTop2)
    return ct.compareTopos(intNumMem1,intNumMem2,strTop1,strTop2, Nterm1,
            Nterm2)
#}}}
def ClassifyTopoComparison(mapArray1, mapArray2, Nterm1, Nterm2):#{{{
# classify the comparsion by 
# if all not _UN_MAPPED
#     if Nterm1 == Nterm2
#         if all _DIRECT_MAPPED
#              1. OK        # Nterm1 == Nterm2 and all _DIRECT_MAPPED
#         else
#              2. SHIFT     # Nterm1 == Nterm2 and all not _UN_MAPPED
#     else
#         if all _DIRECT_MAPPED
#              3. INV       # Nterm1 != Nterm2 and all _DIRECT_MAPPED
#         else
#              4. INV_SHIFT # Nterm1 != Nterm2 and all not _UN_MAPPED
# else
#     5. DIFF      # having _UN_MAPPED
    if min(mapArray1) != _UN_MAPPED and min(mapArray2) != _UN_MAPPED:
        if Nterm1 == Nterm2:
            if min(mapArray1) == _DIRECT_MAPPED:
                return "OK"
            else:
                return "SHIFT"
        else:
            if min(mapArray1) == _DIRECT_MAPPED:
                return "INV"
            else:
                return "INV_SHIFT"
    else:
        return "DIFF"
#}}}
def IsDuplicated0(mapArray1, mapArray2, posTM1, posTM2, seqLength1,#{{{
                seqLength2):
    numTM1 = len(posTM1)
    numTM2 = len(posTM2)
    if numTM1 > numTM2:
        mapArrayLong = mapArray1
        mapArrayShort = mapArray2
        posTMLong = posTM1
        posTMShort = posTM2
        seqLengthLong = seqLength1
        seqLengthShort = seqLength2
    else:
        mapArrayLong = mapArray2
        mapArrayShort = mapArray1
        posTMLong = posTM2
        posTMShort = posTM1
        seqLengthLong = seqLength2
        seqLengthShort = seqLength1
    numTMLong = len(posTMLong)
    numTMShort = len(posTMShort)

    if numTMLong < numTMShort * 2:
        return False
    else:
        ratioSeqlen = float(seqLengthLong)/float(seqLengthShort)
        ratioNumTM = float (numTMLong)/float(numTMShort)
        diffRatio = math.fabs(ratioSeqlen - ratioNumTM)
#        print "diffRatio=", diffRatio, "ratioNumTM=", ratioNumTM, "ratioSeqlen", ratioSeqlen
#        print seqLength1, seqLength2
        if (min(mapArrayShort) != _UN_MAPPED and diffRatio <= 0.4):
            return True
        else:
            return False
#}}}
def IsDuplicated1(mapArray1, mapArray2, posTM1, posTM2, seqLength1,#{{{
                seqLength2, seqid1, seqid2, dupPairSet):
# updated 2012-07-05, using also hhsearch to determine whether it is duplicated
# algorithm:
# if lengthLong/lengthShort > 1.7 and numTMLong%numTMShort <= 1
# check if it is duplicated by hhsearch
    if (seqid1, seqid2) in dupPairSet or (seqid2, seqid1) in dupPairSet:
        isDuplicated = True
    else:
        isDuplicated = False
    return isDuplicated

#}}}
def IsVariedAtSignalPeptide(mapArray1, mapArray2, posTM1, posTM2, #{{{
        seqid1, seqid2, signalpDict):

    isSigPep = False
    if (len(mapArray1) >= 1 and len(mapArray2) >= 1 and 
            (mapArray1[0] == _UN_MAPPED or mapArray2[0] == _UN_MAPPED) and 
            (mapArray1[0] != mapArray2[0])): 
# at least of one of the N-term TM should be unmapped
# note that both N-term TM can not be unmapped
        if mapArray1[0] == _UN_MAPPED:
            if min(mapArray2) > _UN_MAPPED:
                j = 1
                while j < len(mapArray1) and mapArray1[j] == _UN_MAPPED:
                    j += 1
                numNtermUnmappedTM = j
                if min(mapArray1[numNtermUnmappedTM:]) > _UN_MAPPED:
                    if seqid1 in signalpDict:
                        pos_sigp = signalpDict[seqid1]
                        if posTM1[numNtermUnmappedTM-1][1] <= (pos_sigp + 5):
                            isSigPep = True
        elif mapArray2[0] == _UN_MAPPED:
            if min(mapArray1) > _UN_MAPPED:
                j = 1
                while j < len(mapArray2) and mapArray2[j] == _UN_MAPPED:
                    j += 1
                numNtermUnmappedTM = j
                if min(mapArray2[numNtermUnmappedTM:]) > _UN_MAPPED:
                    if seqid1 in signalpDict:
                        pos_sigp = signalpDict[seqid1]
                        if posTM2[numNtermUnmappedTM-1][1] <= (pos_sigp + 5):
                            isSigPep = True
    return isSigPep
#}}}


def ClassifyTopoComparison_pairwise(mapArray1, mapArray2, Nterm1, Nterm2, #{{{
        posTM1, posTM2, seqLength1, seqLength2, seqid1, seqid2, dupPairSet,
        signalpDict):
# classify the comparsion by 
# if all not _UN_MAPPED
#     if Nterm1 == Nterm2
#         if all _DIRECT_MAPPED
#              1. OK        # Nterm1 == Nterm2 and all _DIRECT_MAPPED
#         else
#              2. SHIFT     # Nterm1 == Nterm2 and all not _UN_MAPPED
#     else
#         if all _DIRECT_MAPPED
#              3. INV       # Nterm1 != Nterm2 and all _DIRECT_MAPPED
#         else
#              4. INV_SHIFT # Nterm1 != Nterm2 and all not _UN_MAPPED
# else
#     5. DIFF      # having _UN_MAPPED
    if min(mapArray1) != _UN_MAPPED and min(mapArray2) != _UN_MAPPED:
        if Nterm1 == Nterm2:
            if min(mapArray1) == _DIRECT_MAPPED:
                return "OK"
            else:
                return "SHIFT"
        else:
            if min(mapArray1) == _DIRECT_MAPPED:
                return "INV"
            else:
                return "INV_SHIFT"
    else:
        if IsDuplicated1(mapArray1, mapArray2, posTM1, posTM2, seqLength1,
                seqLength2, seqid1, seqid2, dupPairSet):
            return "DUP"
        elif IsVariedAtSignalPeptide(mapArray1, mapArray2, posTM1, posTM2,
                seqid1, seqid2, signalpDict):
            return "SIGNALP"
        else:
            return "DIFF"
#}}}
def GetNumTMUnaligned_NCterm(mapArray):#{{{
    numTM = len(mapArray)
    i = 0
    while i < numTM:
        if mapArray[i] != _UNALIGNED:
            break
        i += 1
    numNterm = i

    i = numTM - 1
    while i >= 0:
        if mapArray[i] != _UNALIGNED:
            break
        i -= 1
    numCterm = numTM - i - 1
    return (numNterm, numCterm)
#}}}
def ClassifyTopoComparison_pairwise_method1(mapArray1, mapArray2, Nterm1,#{{{
        Nterm2, seqID1, seqID2, cmpclassSP):
# return string for global classification:
# e.g.
# format:

# CLS_OF_ALIGNED_REGION;
# FULL_OR_PART_ALIGNED;
# NUMTM_ALIGNED_1;
# NUMTM_ALIGNED_2;
# numTM_UNALIGNED_NTERM1
# numTM_UNALIGNED_NTERM2
# numTM_UNALIGNED_CTERM1
# numTM_UNALIGNED_CTERM2

# TM2GAP;FULL_ALIGNED;4;6;0;0;0;0
# TM2GAP;PART_ALIGNED;4;6;2;1;1;1

    numTM1 = len(mapArray1)
    numTM2 = len(mapArray2)
#     DEBUG_CLS=1
#     if DEBUG_CLS:
#         print seqID1, mapArray1, "min=", min(mapArray1)
#         print seqID2, mapArray2, "min=", min(mapArray2)
    if min(mapArray1 + mapArray2) == _UNALIGNED:
        isFullAligned = False
    else:
        isFullAligned = True

    if isFullAligned:
        numTM_unaligned_Nterm1 = 0
        numTM_unaligned_Nterm2 = 0
        numTM_unaligned_Cterm1 = 0
        numTM_unaligned_Cterm2 = 0
        aligned_Nterm1 = Nterm1
        aligned_Nterm2 = Nterm2
        aligned_numTM1 = numTM1
        aligned_numTM2 = numTM2
        aligned_mapArray1 = mapArray1
        aligned_mapArray2 = mapArray2
    else:
        aligned_mapArray1 = filter(lambda a:  a!= _UNALIGNED, mapArray1)
        aligned_mapArray2 = filter(lambda a:  a!= _UNALIGNED, mapArray2)
        aligned_numTM1 = len(aligned_mapArray1)
        aligned_numTM2 = len(aligned_mapArray2)
        (numTM_unaligned_Nterm1, numTM_unaligned_Cterm1
                ) = GetNumTMUnaligned_NCterm(mapArray1)
        (numTM_unaligned_Nterm2, numTM_unaligned_Cterm2
                ) = GetNumTMUnaligned_NCterm(mapArray2)
        li_io = ['i','o']
        aligned_Nterm1 =  li_io[(li_io.index(Nterm1) + numTM_unaligned_Nterm1)%2] 
        aligned_Nterm2 =  li_io[(li_io.index(Nterm2) + numTM_unaligned_Nterm2)%2] 



    cls_aligned_region = ""
    if (aligned_numTM1 <= 0 or aligned_numTM2 <= 0 or
            aligned_mapArray1.count(_TM2TM) !=
            aligned_mapArray2.count(_TM2TM)):
        print >> g_params['fpout_badmap'], "Bad TM mapping %s - %s"%(seqID1, seqID2)
        PrintMappedArray_method1_1(mapArray1, mapArray2, seqID1, seqID2,
                g_params['fpout_badmap'])
        if(aligned_numTM1 <= 0 or aligned_numTM2 <= 0):
            cls_aligned_region = "UNALIGNED"
        elif (aligned_mapArray1.count(_TM2TM) !=
                aligned_mapArray2.count(_TM2TM)):
            cls_aligned_region = "AMBIGUOUS"
    else:
        if (max(aligned_mapArray1) == _TM2TM and max(aligned_mapArray2) ==
                _TM2TM and aligned_numTM1 == aligned_numTM2):
            if aligned_Nterm1 == aligned_Nterm2:
                cls_aligned_region =  "IDT"
            else:
                cls_aligned_region =  "INV"
        else:
            uniqMergedSet = set(aligned_mapArray1 + aligned_mapArray2)
            if _TM2GAP in uniqMergedSet and _TM2SEQ in uniqMergedSet:
                cls_aligned_region =  "TM2GAP_AND_TM2SEQ"
            elif _TM2GAP in uniqMergedSet:
                cls_aligned_region =  "TM2GAP"
            elif _TM2SEQ in uniqMergedSet:
                cls_aligned_region =  "TM2SEQ"
            else:
                print >> g_params['fpout_badmap'], "Bad TM mapping %s - %s"%(seqID1, seqID2)
                PrintMappedArray_method1_1(mapArray1, mapArray2, seqID1, seqID2,
                        g_params['fpout_badmap'])
                cls_aligned_region = "UNKNOWN"
    if g_params['isCompareSP']:
        cls_aligned_region = "%s|%s"%(cls_aligned_region, cmpclassSP)

    if cls_aligned_region != "":
        if isFullAligned:
            str_align = "FULL_ALIGNED"
        else:
            str_align = "PART_ALIGNED"

        cls = "%s;%s;%d;%d;%d;%d;%d;%d"%(cls_aligned_region, str_align, aligned_numTM1,
                aligned_numTM2, numTM_unaligned_Nterm1, numTM_unaligned_Nterm2,
                numTM_unaligned_Cterm1, numTM_unaligned_Cterm2)
        return cls
    else:
        return ""
#}}}
def ClassifyTopoComparison_pairwise_method3(mapArray1, mapArray2, Nterm1,#{{{
        Nterm2, seqID1, seqID2):
    """
    return string for global classification:
    e.g.
    format:

    CLS_OF_ALIGNED_REGION;
    FULL_OR_PART_ALIGNED;
    NUMTM_ALIGNED_1;
    NUMTM_ALIGNED_2;
    numTM_UNALIGNED_NTERM1
    numTM_UNALIGNED_NTERM2
    numTM_UNALIGNED_CTERM1
    numTM_UNALIGNED_CTERM2

    TM2GAP;FULL_ALIGNED;4;6;0;0;0;0
    TM2GAP;PART_ALIGNED;4;6;2;1;1;1
    Created 2013-10-07, updated 2013-10-07, Nanjiang Shu 
    """

    numTM1 = len(mapArray1)
    numTM2 = len(mapArray2)
#     DEBUG_CLS=1
#     if DEBUG_CLS:
#         print seqID1, mapArray1, "min=", min(mapArray1)
#         print seqID2, mapArray2, "min=", min(mapArray2)
    if min(mapArray1 + mapArray2) == _UNALIGNED:
        isFullAligned = False
    else:
        isFullAligned = True

    if isFullAligned:
        numTM_unaligned_Nterm1 = 0
        numTM_unaligned_Nterm2 = 0
        numTM_unaligned_Cterm1 = 0
        numTM_unaligned_Cterm2 = 0
        aligned_Nterm1 = Nterm1
        aligned_Nterm2 = Nterm2
        aligned_numTM1 = numTM1
        aligned_numTM2 = numTM2
        aligned_mapArray1 = mapArray1
        aligned_mapArray2 = mapArray2
    else:
        aligned_mapArray1 = filter(lambda a:  a!= _UNALIGNED, mapArray1)
        aligned_mapArray2 = filter(lambda a:  a!= _UNALIGNED, mapArray2)
        aligned_numTM1 = len(aligned_mapArray1)
        aligned_numTM2 = len(aligned_mapArray2)
        (numTM_unaligned_Nterm1, numTM_unaligned_Cterm1
                ) = GetNumTMUnaligned_NCterm(mapArray1)
        (numTM_unaligned_Nterm2, numTM_unaligned_Cterm2
                ) = GetNumTMUnaligned_NCterm(mapArray2)
        li_io = ['i','o']
        aligned_Nterm1 =  li_io[(li_io.index(Nterm1) + numTM_unaligned_Nterm1)%2] 
        aligned_Nterm2 =  li_io[(li_io.index(Nterm2) + numTM_unaligned_Nterm2)%2] 



    cls_aligned_region = ""
    if (aligned_numTM1 <= 0 or aligned_numTM2 <= 0 or
            aligned_mapArray1.count(_TM2TM) !=
            aligned_mapArray2.count(_TM2TM)):
        print >> g_params['fpout_badmap'], "Bad TM mapping %s - %s"%(seqID1, seqID2)
        PrintMappedArray_method1_1(mapArray1, mapArray2, seqID1, seqID2,
                g_params['fpout_badmap'])
        if(aligned_numTM1 <= 0 or aligned_numTM2 <= 0):
            cls_aligned_region = "UNALIGNED"
        elif (aligned_mapArray1.count(_TM2TM) !=
                aligned_mapArray2.count(_TM2TM)):
            cls_aligned_region = "AMBIGUOUS"
    else:
        if (max(aligned_mapArray1) == _TM2TM and max(aligned_mapArray2) ==
                _TM2TM and aligned_numTM1 == aligned_numTM2):
            if aligned_Nterm1 == aligned_Nterm2:
                cls_aligned_region =  "IDT"
            else:
                cls_aligned_region =  "INV"
        else:
            uniqMergedSet = set(aligned_mapArray1 + aligned_mapArray2)
            if (uniqMergedSet == set([_TM2TM, _TM2GAP]) 
                    or uniqMergedSet == set([_TM2GAP])):
                cls_aligned_region = "TM2GAP"
            elif (uniqMergedSet == set([_TM2TM, _TM2SEQ]) 
                    or uniqMergedSet == set([_TM2SEQ])):
                cls_aligned_region = "TM2SEQ"
            elif (uniqMergedSet == set([_TM2TM, _TM2SP]) 
                    or uniqMergedSet == set([_TM2SP])):
                cls_aligned_region = "TM2SP"
            else:
                cls_aligned_region =  "Mixed"

    dupPairSet = g_params['dupPairSet']
    if g_params['isCompareDup'] and dupPairSet != {}:
        if (seqID1, seqID2) in dupPairSet or (seqID2, seqID1) in dupPairSet:
            cls_aligned_region = "%s|%s"%(cls_aligned_region, "DUP")
        else:
            cls_aligned_region = "%s|%s"%(cls_aligned_region, "nonDUP")
    if cls_aligned_region != "":
        if isFullAligned:
            str_align = "FULL_ALIGNED"
        else:
            str_align = "PART_ALIGNED"

        cls = "%s;%s;%d;%d;%d;%d;%d;%d"%(cls_aligned_region, str_align, aligned_numTM1,
                aligned_numTM2, numTM_unaligned_Nterm1, numTM_unaligned_Nterm2,
                numTM_unaligned_Cterm1, numTM_unaligned_Cterm2)
        return cls
    else:
        return ""
#}}}


def AnaDIFFTopology(topoquery,numTMquery, consensusTopo, numTMcons):#{{{
    N_ins=0
    N_del=0
    C_ins=0
    C_del=0
    i_ins=[0]*numTMcons
    i_del=[0]*numTMcons
#  find the location of TMs in the consensusTopo
    m=re.finditer("(M+)", consensusTopo)
    startCons=[]
    endCons=[]
    for jr in m:
        startCons.append(jr.start(0))
        endCons.append(jr.end(0))
    if len(startCons) != numTMcons:
        sys.stderr.write("Error! numTM of the consensusTopo is inconsistent\n")
# find the location of TMs in the query topo
    m=re.finditer("(M+)", topoquery)
    startTopo=[]
    endTopo=[]
    for jr in m:
        startTopo.append(jr.start(0))
        endTopo.append(jr.end(0))
    if len(startTopo) != numTMquery:
        sys.stderr.write("Error! numTM of the query topo is inconsistent\n")

    
#1. determine the condition at N-terminal 
    if endTopo[0] < startCons[0]:
        N_ins = 1
    if startTopo[0] > endCons[0]:
        N_del = 1
#2. determine the condition at C-terminal 
    if startTopo[numTMtopo-1] > endCons[numTMcons-1]:
        C_ins = 1
    if endTopo[numTMtopo-1] < startCons[numTMcons-1]:
        C_del = 1

#3. determine the condition for internal deletion
    for i in range(numTMcons-2):
        substr=topoquery[endCons[i]:startCons[i+2]]
        if substr.count("M") == 0:
            i_del[i] = 1
#4. determine the condition for internal insertion
    for i in range(numTMcons-1):
        substr=topoquery[startCons[i]:endCons[i+1]]
        if myfunc.CountTM(substr) >= 3:
            i_ins[i] = 1

    return (N_ins,N_del,C_ins,C_del, i_ins, i_del)
    #}}}
def AnaDIFFTopology1(mapArray ):#{{{
# Since it is hard to say insertions or deletions without knowing the
# phylogenetic tree, we just count numbers here
# At terminals
#   numUnmappedTMatNterm, indices of TM helices
#   numUnmappedTMatCterm, indices of TM helices
# for internal regions
#   numUnmapedTMatInternal_1, indices of TM helices
#   numUnmapedTMatInternal_2, indices of TM helices
#   numUnmapedTMatInternal_3, indices of TM helices 
#   ...
# Note that we can also do pairwise comparison
# e.g. sequence identity vs number of variations by means of the above numbers
    numTM = len(mapArray)
    ana = {}
    ana["Nterm"] ={}
    ana["Nterm"]["numTMunmapped"] = 0
    ana["Nterm"]["index"] = []
    ana["Cterm"] = {}
    ana["Cterm"]["numTMunmapped"] = 0
    ana["Cterm"]["index"] = []
    ana["internal"] = []
    string=""
    for i in mapArray:
        if i == _UN_MAPPED: 
            string += "0"
        else:
            string += "1"
    posUnmappedTM=[]
    m=re.finditer("(0+)", string)
    cntinter=0
    for i in m:
        posUnmappedTM.append((i.start(0), i.end(0)))
    for i in range(len(posUnmappedTM)):
        (b,e)=posUnmappedTM[i]
        if b == 0:
            ana["Nterm"]["numTMunmapped"] = e-b
            ana["Nterm"]["index"]=range(b,e)
        elif e == numTM  :
            ana["Cterm"]["numTMunmapped"] = e-b
            ana["Cterm"]["index"]=range(b,e)
        else:
            ana["internal"].append({})
            ana["internal"][cntinter]["numTMunmapped"] = e-b
            ana["internal"][cntinter]["index"] = range(b,e)
            cntinter +=1
    return ana
#}}}
def CheckGapOfMSA(ana, posTM,gapSeqOpposite):#{{{
#@params
# ana           : result for pairwise topology comparison, structured in
#               : dictionary
# posTM         : list of position of TM regions in 2-tuples (beg, end)
# gapSeqOpposite: sequence of gap percentages of the opposite topology
#               : (compared topology), for consensus topology, it is the
#               : percentages of each column, for a real topology, it is a
#               : series of 0 and 1
#Check the gap percentage of the opposite topology at the unmapped TM region
    if ana["Nterm"]["numTMunmapped"] > 0:
        ana["Nterm"]["gapFraction"] =[]
        for i in range(ana["Nterm"]["numTMunmapped"] ):
            (b,e) = posTM[ana["Nterm"]["index"][i]]
            try: 
                gapFrac=sum(gapSeqOpposite[b:e])/float(e-b)
            except ZeroDivisionError:
                gapFrac = 0.0
            ana["Nterm"]["gapFraction"].append(gapFrac)
    if ana["Cterm"]["numTMunmapped"] > 0:
        ana["Cterm"]["gapFraction"] =[]
        for i in range(ana["Cterm"]["numTMunmapped"] ):
            (b,e) = posTM[ana["Cterm"]["index"][i]]
            try:
                gapFrac=sum(gapSeqOpposite[b:e])/float(e-b)
            except ZeroDivisionError:
                gapFrac = 0.0
            ana["Cterm"]["gapFraction"].append(gapFrac)
    for m in range ( len(ana["internal"]) ):
        if ana["internal"][m]["numTMunmapped"] > 0:
            ana["internal"][m]["gapFraction"] =[]
            for i in range(ana["internal"][m]["numTMunmapped"] ):
                (b,e) = posTM[ana["internal"][m]["index"][i]]
                try:
                    gapFrac=sum(gapSeqOpposite[b:e])/float(e-b)
                except ZeroDivisionError:
                    gapFrac = 0.0
                ana["internal"][m]["gapFraction"].append(gapFrac)
    return ana
#}}}
def CheckDGOfMSA(ana, DGvalueTM):#{{{
#@params
# ana           : result for pairwise topology comparison, structured in
#               : dictionary
# DGvalueTM     : list of DG values for all TM helices 
#               : for consensus topology, it is the average of each column, 
#               : for a real topology, it is the DG value of each TM region
#Check the gap percentage of the opposite topology at the unmapped TM region
    if ana["Nterm"]["numTMunmapped"] > 0:
        ana["Nterm"]["DGvalue"] =[]
        for i in range(ana["Nterm"]["numTMunmapped"] ):
            ana["Nterm"]["DGvalue"].append(DGvalueTM[ana["Nterm"]["index"][i]])
    if ana["Cterm"]["numTMunmapped"] > 0:
        ana["Cterm"]["DGvalue"] =[]
        for i in range(ana["Cterm"]["numTMunmapped"] ):
            ana["Cterm"]["DGvalue"].append(DGvalueTM[ana["Cterm"]["index"][i]])
    for m in range ( len(ana["internal"]) ):
        if ana["internal"][m]["numTMunmapped"] > 0:
            ana["internal"][m]["DGvalue"] =[]
            for i in range(ana["internal"][m]["numTMunmapped"] ):
                ana["internal"][m]["DGvalue"].append(
                        DGvalueTM[ana["internal"][m]["index"][i]])
    return ana
#}}}

def GetTopoStateFraction_slow(idtTopoSeqList):#{{{
    "this version is slow"
    lengthAlignment=len(idtTopoSeqList[0])
    numIDTTopo = len(idtTopoSeqList)
    numIDTTopo_float = float(numIDTTopo)
    statestr = 'oiM-'
    perStateList=[]
    for kk in range(len(statestr)):
        binaryMatrix=[]
        apd=binaryMatrix.append
        state=statestr[kk]
        for i in xrange(numIDTTopo):
            apd ([(s==state) for s in idtTopoSeqList[i]]); 
        cntState=[sum( [binaryMatrix[i][j] for i in xrange(numIDTTopo)]) for j
                in xrange(lengthAlignment)]; 
        perStateList.append(  [x/numIDTTopo_float for x in cntState ] ); 

    return (perStateList[0], perStateList[1], perStateList[2], perStateList[3])

#}}}

def GetConsensusTopo(idtTopoSeqList, idtPosTMList, method_consensus):#{{{
# @params:
#   idtTopoSeqList  topology sequences of the largest identical group
#   posTMList       location of TM regions for all sequences
# method_consensus:
#   2   get consensus based on the agreement of posTM
    numIDTTopo = len(idtTopoSeqList)
    if numIDTTopo <= 0:
        return (None, None)

    lengthAlignment=len(idtTopoSeqList[0])
    (cnt_i, cnt_o, cnt_M, cnt_GAP, 
            per_i, per_o, per_M, per_GAP) = lcmp.GetTopoStateFraction(
                    idtTopoSeqList)

    consensusTopo=""
    if method_consensus == 0:
        state="";#{{{
        for i in range(lengthAlignment):
            if per_GAP[i] > 0.5:
                state=GAP
            else:
                if per_M[i] >= per_o[i] and per_M[i] >= per_i[i]:
                    state="M"
                elif per_i[i] >= per_M[i] and per_i[i] >= per_o[i]:
                    state="i"
                else:
                    state="o"
            consensusTopo+=state;#}}}
    elif method_consensus == 1:
        TM_extend_threshold=0.5;#{{{
        strlist=[GAP]*lengthAlignment
        for i in range(lengthAlignment):
            if per_M[i] >= 1.0:
                strlist[i] = "M"
        newstr="".join(strlist)

        posCoreTM = myfunc.GetTMPosition(newstr)
        numTM = len(posCoreTM)
#         print "newseq=",newstr ; #debug
#         print "posCoreTM=",posCoreTM; #debug
        for i in range(numTM):
            if i - 1 < 0:
                leftborder=0
            else:
                leftborder=posCoreTM[i-1][1]+1
            if i +1 > numTM-1:
                rightborder=lengthAlignment-1
            else:
                rightborder=posCoreTM[i+1][0]-1

            j = posCoreTM[i][0]
            while j>leftborder:
                if per_GAP[j] <= 0.5:
                    if per_M[j] > (per_i[j]+per_o[j]):
                        strlist[j] = "M"
                    else:
                        break
                j-=1
            j = posCoreTM[i][1]
            while j<rightborder:
                if per_GAP[j] <= 0.5:
                    if per_M[j] > (per_i[j]+per_o[j]):
                        strlist[j] = "M"
                    else:
                        break
                j+=1
        for j in range(lengthAlignment):
            if strlist[j] != 'M' and per_GAP[j] <= 0.5:
                if per_i[j] >= per_o[j]:
                    strlist[j]="i"
                else:
                    strlist[j]="o"
        consensusTopo="".join(strlist)
#}}}
    elif method_consensus == 2:
        strlist=[GAP]*lengthAlignment
        numTM = len(idtPosTMList[0]);#{{{
        posTMCons=[]; # 2-tuple list
        cntTMpos=0
        for i in xrange(numTM):
            posBegs = [idtPosTMList[j][i][0] for j in xrange(numIDTTopo)]
            posEnds = [idtPosTMList[j][i][1] for j in xrange(numIDTTopo)]
            posBegs.sort()
            posEnds.sort()
            idxTMPosToBeUsed_begin = 0
            idxTMPosToBeUsed_end = 1
            if numIDTTopo <= 2 :
                idxTMPosToBeUsed_begin = 0
                idxTMPosToBeUsed_end = numIDTTopo
            else:
                idxTMPosToBeUsed_begin = max(1,int(round(numIDTTopo*0.15)));  #take only 70% of the middle positions
                idxTMPosToBeUsed_end = numIDTTopo - idxTMPosToBeUsed_begin
            #print [idxTMPosToBeUsed_begin, idxTMPosToBeUsed_end], numIDTTopo

            avgBeg = sum(posBegs[idxTMPosToBeUsed_begin: idxTMPosToBeUsed_end])/float(idxTMPosToBeUsed_end-idxTMPosToBeUsed_begin)
            avgEnd = sum(posEnds[idxTMPosToBeUsed_begin: idxTMPosToBeUsed_end])/float(idxTMPosToBeUsed_end-idxTMPosToBeUsed_begin)
            avgBeg = int(round(avgBeg))
            avgEnd = int(round(avgEnd))
            if len(posTMCons) > 0 and avgBeg <= posTMCons[cntTMpos-1][1]:
                avgBeg = posTMCons[i-1][1] + 1
            if avgEnd <= avgBeg:
                print >> sys.stderr, ("TM %d is overlapping and neglected." %
                        (i+1))
            else:
                posTMCons.append((avgBeg, avgEnd))
                cntTMpos += 1
        NtermStateCons = lcmp.GetNtermState(idtTopoSeqList[0])
        if not NtermStateCons in 'io':
            print >> sys.stderr, ("GetNtermState failed when "
                    + "making consensusTopo.")
            return (None, None)
# assign 'M' state and 'i' 'o' sate at the beginnings and ends of TM
# regions for the consensus 
        stateiostr = "io"
        state_index= stateiostr.index(NtermStateCons)
        for pos in posTMCons:
            for j in range(pos[0],pos[1]):
                if per_GAP[j] <= 0.5:
                    strlist[j] = 'M'
            # set the beginnings and ends
            if pos[0] == 0:
                strlist[pos[0]] = stateiostr[state_index]
            else:
                strlist[pos[0]-1] = stateiostr[state_index]
            state_index = (state_index+1)%2
            if pos[1] == lengthAlignment:
                strlist[pos[1]-1] = stateiostr[state_index]
            else:
                strlist[pos[1]] = stateiostr[state_index]

        for j in range(lengthAlignment):
            if per_GAP[j] <= 0.5 and strlist[j] == GAP:
                if per_i[j] >= per_o[j]:
                    strlist[j]="i"
                else:
                    strlist[j]="o"
        # proof reading the consensus topology, the i,o sate between two TM
        # helices should be consistent.
        state_index= stateiostr.index(NtermStateCons)
        for i in xrange(numTM+1):
            beg = 0
            end = 0
            if i > 0:
                beg=posTMCons[i-1][1]
            else:
                beg = 0
            if i < numTM:
                end = posTMCons[i][0]
            else:
                end = lengthAlignment -1
            correctState =  stateiostr[state_index]; 
            for j in xrange(beg, end):
                if strlist[j] != GAP and strlist[j] != correctState:
                    strlist[j] == correctState
            state_index = (state_index+1)%2

        consensusTopo="".join(strlist)
        #}}}
    return (per_GAP, consensusTopo) 
#}}}
def GetIDTTopoGroup(idList,topoSeqList):#{{{
    #return (consensusTopo, indexIDTTopo,  numTM_IDT, Mcmp)
# 0 for different
# 1 for identical
    numSeq=len(idList)
    # Initialize a unit matrix of the size (numSeq x numSeq)
    Mcmp= [[0]*numSeq for x in xrange(numSeq)]
    for i in range(numSeq): 
        Mcmp[i][i]=1
    # All-to-all comparison
    for i in range (numSeq):
        for j in range(i+1, numSeq):
            strTop1 = topoSeqList[i]
            strTop2 = topoSeqList[j]
            strProtein1 = idList[i]
            strProtein2 = idList[j]
            (class_global, num1_global, num2_global) = CompareTrimmedToposGlobally(strTop1, strTop2, strProtein1, strProtein2)
            if class_global == "OK":
                Mcmp[i][j] = 1
            else:
                Mcmp[i][j] = 0
    #copy the symmetric matrix
    for i in range (numSeq):
        for j in range(i+1, numSeq):
            Mcmp[j][i]=Mcmp[i][j]

    #get the largest group with the identical topology
    cntIDTTopo=[sum(ii) for ii in Mcmp]
    maxCnt=max(cntIDTTopo)
    maxIndex=cntIDTTopo.index(maxCnt)
    numIDTTopo=maxCnt
    indexIDTTopo=[]
    for j in range(numSeq):
        if Mcmp[maxIndex][j]==1:
            indexIDTTopo.append(j)

    numIDTTopo=len(indexIDTTopo)
# To get the consensusTopo
    idtTopoSeqList=[]
    idtPosTMList=[]
    for idx in indexIDTTopo:
        idtTopoSeqList.append(topoSeqList[idx])
        idtPosTMList.append(posTMList[idx])
    (gapSeqCons, consensusTopo) = GetConsensusTopo(idtTopoSeqList,
            idtPosTMList, method_consensus)

# check whether the numTM in consensusTopo is the same as the average
# number of numTM of the largest identical topology group
    numTMList = []
    for i in range(numIDTTopo):
        numTMList.append(myfunc.CountTM(topoSeqList[indexIDTTopo[i]]))
    # get the most common element in the list
    counter={}
    for i in numTMList: counter[i] = counter.get(i, 0) + 1
    (maxcnt, numTM_IDT) = sorted([(cnt, num) for num,cnt in counter.items()],
            reverse=True)[0]
    numTMcons = myfunc.CountTM(consensusTopo)
    if numTMcons != numTM_IDT:
        sys.stderr.write("%s: numTMcons (%d) != numTM_largestgroup (%d)\n" %
                (inFile, numTMcons, numTM_IDT))

#debug code
    fpLog = g_params['fpLog']
    DEBUG_CONSENSUS = g_params['DEBUG_CONSENSUS']
    if fpLog != None and DEBUG_CONSENSUS:
        for i in range (numSeq):
            for j in range(numSeq):
                fpLog.write("%d "%Mcmp[i][j])
            fpLog.write("\n")
        print >> fpLog, "numIDTTopo=",numIDTTopo
        print >> fpLog, "indexIDTTopo:\n",indexIDTTopo
        print >> fpLog
        for i in range (numIDTTopo):
            print >> fpLog, topoSeqList[indexIDTTopo[i]]
        print >> fpLog, "ConsensusTopo:"
        print >> fpLog, consensusTopo

        conf_o=[max(0,int(x*10)-1) for x in per_o]; 
        conf_M=[max(0,int(x*10)-1) for x in per_M]; 
        conf_i=[max(0,int(x*10)-1) for x in per_i]; 
        print >> fpLog, ''.join(["%d"%x for x in conf_M])
        print >> fpLog, ''.join(["%d"%x for x in conf_o])
        print >> fpLog, ''.join(["%d"%x for x in conf_i])

        print >> fpLog, conf_o
        print >> fpLog, conf_M
        print >> fpLog, conf_i
        print >> fpLog, "numTMcons=", numTMcons
        print >> fpLog, "numTM_largestgroup=", numTM_IDT
#debug code
    return (consensusTopo, indexIDTTopo,  numTM_IDT, Mcmp)
#}}}
def GetIDTTopoGroup1(topoSeqList, NtermStateList, posTMList):#{{{
    #return (consensusTopo, indexIDTTopo,  gapSeqCons, numTM_IDT, Mcmp)
# return ()
# 0 for different
# 1 for identical
    numSeq=len(posTMList)
    numTMList=[ len (posTM) for posTM in posTMList]

    # Initialize a unit matrix of the size (numSeq x numSeq)
#   Mcmp = [ [(i==j) for i in xrange(numSeq)] for j in xrange(numSeq) ]  
    Mcmp= [[0]*numSeq for x in xrange(numSeq)]
    for i in xrange(numSeq): 
        Mcmp[i][i]=1
    # All-to-all comparison
    if g_params['method_comparison'] == 2: # modified 2012-11-22 
        min_TM_overlap = int(g_params['threshold_TM2TM']*21+0.5)
    else:
        min_TM_overlap = g_params['min_TM_overlap']
    for i in xrange (numSeq):
        for j in xrange(i+1, numSeq):
            #if lcmp.IsIdenticalTopology(NtermStateList[i], NtermStateList[j], numTMList[i], numTMList[j], posTMList[i], posTMList[j]):
            if lcmp.IsIdenticalTopology(NtermStateList[i], NtermStateList[j],
                    numTMList[i], numTMList[j], posTMList[i], posTMList[j],
                    topoSeqList[i], topoSeqList[j], min_TM_overlap):
                Mcmp[i][j] = 1
            else:
                Mcmp[i][j] = 0
    #copy the symmetric matrix
    for i in xrange (numSeq):
        for j in xrange(i+1, numSeq):
            Mcmp[j][i]=Mcmp[i][j]

    #get the largest group with the identical topology, bug fixed 2011-09-19 22:26:53 Monday Week 38
    cntIDTTopo=[sum(ii) for ii in Mcmp]
    maxCnt=max(cntIDTTopo)
    maxIndex=cntIDTTopo.index(maxCnt)
    numIDTTopo=maxCnt
    indexIDTTopo=[]
    for j in range(numSeq):
        if Mcmp[maxIndex][j]==1:
            indexIDTTopo.append(j)

    numIDTTopo = len(indexIDTTopo)
# To get the consensusTopo
    idtTopoSeqList=[]
    idtPosTMList=[]
    for idx in indexIDTTopo:
        idtTopoSeqList.append(topoSeqList[idx])
        idtPosTMList.append(posTMList[idx])
    (gapSeqCons, consensusTopo) = GetConsensusTopo(idtTopoSeqList, idtPosTMList, g_params['method_consensus'])

    # check whether the numTM in consensusTopo is the same as the average
    # number of numTM of the largest identical topology group
    numTM_IDT = len(posTMList[indexIDTTopo[0]]) 
    numTMcons = myfunc.CountTM(consensusTopo)
    if numTMcons != numTM_IDT:
        sys.stderr.write("%s: numTMcons (%d) != numTM_largestgroup (%d)\n"%(inFile, numTMcons, numTM_IDT))

#debug code begin#{{{
    fpLog = g_params['fpLog']
    DEBUG_CONSENSUS = g_params['DEBUG_CONSENSUS']
    if fpLog != None and DEBUG_CONSENSUS:
        for i in range (numSeq):
            for j in range(numSeq):
                fpLog.write("%d "%Mcmp[i][j])
            fpLog.write("\n")
        print >> fpLog, "numIDTTopo=",numIDTTopo
        print >> fpLog, "indexIDTTopo:\n",indexIDTTopo
        print >> fpLog
        for i in range (numIDTTopo):
            print >> fpLog, topoSeqList[indexIDTTopo[i]]
        print >> fpLog, "ConsensusTopo:"
        print >> fpLog, consensusTopo

        print >> fpLog, "numTMcons=", numTMcons
        print >> fpLog, "numTM_largestgroup=", numTM_IDT
#debug code ends#}}}

    return (consensusTopo, indexIDTTopo,  gapSeqCons, numTM_IDT, Mcmp)
#}}}
def GetInvertedTopologyMatrix(topoSeqList, NtermStateList, posTMList):#{{{
    #return (Mcmp)
# return ()
# 0 for not-inverted
# 1 for inverted
    numSeq=len(posTMList)
    numTMList=[ len (posTM) for posTM in posTMList]
    # Initialize a unit matrix of the size (numSeq x numSeq)
#   Mcmp = [ [(i==j) for i in xrange(numSeq)] for j in xrange(numSeq) ]  
    Mcmp= [[0]*numSeq for x in xrange(numSeq)]
    # All-to-all comparison
    min_TM_overlap = 5
    for i in xrange (numSeq):
        for j in xrange(i+1, numSeq):
            if lcmp.IsInvertedTopology(NtermStateList[i], NtermStateList[j],
                    numTMList[i], numTMList[j], posTMList[i], posTMList[j],
                    topoSeqList[i], topoSeqList[j], min_TM_overlap):
                Mcmp[i][j] = 1
            else:
                Mcmp[i][j] = 0
# it is a half matrix
    return Mcmp
#}}}
def CompareToConsensus(topoSeqList, indexIDTTopo,consensusTopo, idList, #{{{
        numTM_IDT, comparisonClassNameList):
    """compare each topology not in the largest identical group to the
    consensusTopo"""
    numSeq=len(idList)
    numIDTTopo=len(indexIDTTopo)
    fulllist=[x for x in range(numSeq)]
    indexOtherTopo=list(set(fulllist)-set(indexIDTTopo)); # set subtraction
    numOtherTopo=len(indexOtherTopo)
    cmpClassOtherTopoList=[];# a list of classes for all topologies not in the IDTgroup
    numTMcons = myfunc.CountTM(consensusTopo)
    numCmpClass=len(comparisonClassNameList); #number of defined classes
    indexClass=[]
    for i in range(numCmpClass):
        indexClass.append([])
    infoDIFF=[]

    cntN_ins=0
    cntN_del=0
    cntC_ins=0
    cntC_del=0
    cnti_ins=[0]*numTMcons
    cnti_del=[0]*numTMcons

    cntdiff = 0; # counter for the proteins with DIFF topology
    for i in range(numOtherTopo):
        topoquery=topoSeqList[indexOtherTopo[i]]
        numTMquery = myfunc.CountTM(topoquery)
        (class_global, num1_global, num2_global) = CompareTrimmedToposGlobally(
                topoquery, consensusTopo, idList[indexOtherTopo[i]],
                "Consensus", fpLog); # OK, DIFF, INV, SHIFT analysis
        cmpClassOtherTopoList.append(class_global)

        idx = comparisonClassNameList.index(class_global)
        indexClass[idx].append(indexOtherTopo[i])

        # for DIFF topology, do further analysis. Analyzing whether it has
        # N-terminal indels, C-terminal indels and the number of internal
        # indels
        if class_global=="DIFF":
            (N_ins,N_del,C_ins,C_del, i_ins, i_del) = AnaDIFFTopology(
                    topoquery , numTMquery, consensusTopo, numTMcons)
            cntN_ins+=N_ins
            cntN_del+=N_del
            cntC_ins+=C_ins
            cntC_del+=C_del
            for ii in range(numTMcons):
                cnti_ins[ii] += i_ins[ii]
                cnti_del[ii] += i_del[ii]

            infoDIFF.append({})
            infoDIFF[cntdiff]["N_ins"] = N_ins
            infoDIFF[cntdiff]["N_del"] = N_del
            infoDIFF[cntdiff]["C_ins"] = C_ins
            infoDIFF[cntdiff]["C_del"] = C_del
            infoDIFF[cntdiff]["i_ins"] = i_ins
            infoDIFF[cntdiff]["i_del"] = i_del
            cntdiff += 1

    # note that maximum internal insertions is numTMcons - 1
    # while the maximum internal deletions is numTMcons - 2
    normI_ins = sum(x for x in cnti_ins)/float(numTMcons-1+1e-6); 
    normI_del = sum(x for x in cnti_del)/float(numTMcons-2+1e-6)
    # get the occurrence of each element in the list
    counter={}
    for i in cmpClassOtherTopoList: counter[i] = counter.get(i, 0) + 1
    sortedOccurrence1=sorted([(status, cnt) for status,cnt in counter.items()], reverse=True)
    cntOccu=[0]*numCmpClass
    for i in range(numCmpClass):
        for j in range(len(sortedOccurrence1)):
            if comparisonClassNameList[i] == sortedOccurrence1[j][0]:
                cntOccu[i]=sortedOccurrence1[j][1]

    # header line
    fpout = g_params['fpout']
    fpout.write("%-7s %5s %6s %9s %8s" %("PfamID","nSeq", "numIDT", "numIDT/N%", "numTMIDT"))
    for i in range(numCmpClass):
        fpout.write(" %5s"% comparisonClassNameList[i])
    for i in range(numCmpClass):
        fpout.write(" %6s"% (comparisonClassNameList[i]+"%"))
    fpout.write(" %4s %4s %4s %4s %8s %8s"%("Nins", "Ndel", "Cins", "Cdel", "normIins", "normIdel"))
    fpout.write(" %s"%("persite_info"))
    fpout.write("\n")
        

    # content
    fpout.write("%-7s %5d %6d %9.1f %8s" %(rootname, numSeq, numIDTTopo, numIDTTopo/float(numSeq)*100, numTM_IDT))
    for i in range(numCmpClass):
        fpout.write(" %5d"% cntOccu[i])
    for i in range(numCmpClass):
        fpout.write(" %6.1f"% (cntOccu[i]/float(numSeq)*100))
    fpout.write(" %4d %4d %4d %4d %8.1f %8.1f"%(cntN_ins, cntN_del, cntC_ins, cntC_del, normI_ins, normI_del))
    fpout.write(" |")
    for i in range(numTMcons-1):
        fpout.write(" %2d"%cnti_ins[i])
    fpout.write(" |")
    for i in range(numTMcons-2):
        fpout.write(" %2d"%cnti_del[i])
    
    fpout.write("\n")

    return (indexClass, infoDIFF)
#    print sortedOccurrence1

#}}}
def CompareToConsensus1(topoSeqList, NtermStateList, #{{{
        indexIDTTopo,consensusTopo, gapSeqCons, idList, posTMList, dgScoreList,
        comparisonClassNameList):
    """ @param
        topoSeqList    Untrimmed (with gaps) aligned topology sequences
        NtermStateList N-terminal topology state of all sequences
        indexIDTTopo   Index of topology sequence in the largest identity group
        consensusTopo  The consensus topology built from the largest identity
                       group
        idList         Identities of all sequences
        posTMList      Location of TM helices in 2-tuple (beg, end) for all
                       untrimmed topology sequences
        dgScoreList    Delta G values of all TM regions of all sequences 
        ComparisonClassNameList 
                       Name list of comparison classes. e.g. 
                       ["OK","SHIFT","INV","DIFF" ]
    """

# Compare each topology not in the largest identical group to the consensusTopo
# brief description, 
# 1. Map every TM region to the consensus topology
#
# 2. For the unmapped TM regions of the query topology, analyze its condition
#   2.1. If it is an insertion, count the number of inserted TM regions and
#        check
#      2.1.1. Whether it is mostly gaps at the aligned regions of the IDTgroup
#      2.1.2. DG value distribution
#      2.1.3. topcons_single prediction
#      2.1.4. Check gene duplication
#
# 3. For the unmapped TM regions of the consensus topology 
#    (deletion should be analyzed in the other way)
#   3.1. If it is a deletion, count the number of deleted TM regions and check
#     3.1.1. DG value distribution of the topologies in the IDTgroup which
#            generate this TM helix
#     3.1.2. DG value of the query topology, to see if it is a mis-prediction
#     3.1.3. topcons_single prediction 
#     3.1.4. Check gene duplication
#     3.1.5. If the above is check qualified, check if it is a gap at the query
#            topology

    fpLog = g_params['fpLog']
    fpout = g_params['fpout']
    DEBUG_CONSENSUS = g_params['DEBUG_CONSENSUS']

    numSeq = len(idList)
    numIDTTopo = len(indexIDTTopo)
    fulllist = [x for x in range(numSeq)]
    indexOtherTopo = list(set(fulllist)-set(indexIDTTopo)); # set subtraction
    numOtherTopo = len(indexOtherTopo)
    cmpClassOtherTopoList = [];# comparison class list

    posTMcons = myfunc.GetTMPosition(consensusTopo)
    numTMcons = len(posTMcons)

    # Get number of TM helices for the largest identity group
    numTM_IDT = len (posTMList[indexIDTTopo[0]]); 

    DGvalueTMcons = GetDGvalueTMconsensus(dgScoreList, indexIDTTopo, numTM_IDT)
    if g_params['isPrintDGCons']:
        print "DGvalueTMcons (%d) =" % len(DGvalueTMcons), ["%.3f"%x for x in
                DGvalueTMcons]
    numCmpClass = len(comparisonClassNameList); #number of defined classes

    # Get length of un-aligned sequences
    seqLenList = [ len(tp.replace(GAP,'')) for tp in topoSeqList]; 
    seqLenConsensus = len(consensusTopo.replace(GAP,''))

    indexClass=[]
    for i in range(numCmpClass):
        indexClass.append([])
    infoDIFF=[]

#     anaDIFFConsList=[]
#     anaDIFFQueryList=[]
# recordList for the pairwise comparison to the consensus
# each record is a dictionary
# record['cmpclass'] = 'DIFF'
# record['ana1'] are only for when cmpclass = DIFF, ana1 is for the consensus
# record['ana2'] are only for when cmpclass = DIFF, ana2 is for the query
    cmpToConsRecordList = []

    Ntermcons = lcmp.GetNtermState(consensusTopo)

    for i in xrange(numOtherTopo):

        cmpToConsRecordList.append({})
        cmprecord = cmpToConsRecordList[i]

        seqID = idList[indexOtherTopo[i]]
        posTMquery = posTMList[indexOtherTopo[i]]
        numTMquery = len(posTMquery)
        topoquery = topoSeqList[indexOtherTopo[i]]
        lengthQuery = len(topoquery)
        seqLenQuery = seqLenList[i]
        DGvalueTMquery = [INIT_DGVALUE] * numTMquery
        dglist =  dgScoreList[indexOtherTopo[i]]; 
        if dglist and len(dglist) == numTMquery:
            DGvalueTMquery = dglist
        gapSeqQuery = [ (s==GAP) for s in topoquery]; # avoid loop, it is faster
        commonMarray = [ (topoquery[j] == 'M' and consensusTopo[j] == 'M') for
                j in range(lengthQuery)]; 

        (mapArraycons, mapArrayquery) = MappingTM(posTMcons,
                posTMquery,commonMarray, "Consensus",
                idList[indexOtherTopo[i]])
        Ntermquery = NtermStateList[indexOtherTopo[i]] 
        class_global = ClassifyTopoComparison(mapArraycons, mapArrayquery,
                Ntermcons, Ntermquery)
        if len(mapArraycons) != len(mapArrayquery) and class_global == "OK":
            print seqID, mapArrayquery
            print "cons", mapArraycons

        cmpClassOtherTopoList.append(class_global)
        idx = comparisonClassNameList.index(class_global)

        if DEBUG_CONSENSUS and fpLog != None:
            print >> fpLog, "SeqID", seqID, "class=", class_global, "idx=", idx
            print >> fpLog, ("SeqID %s: %s"%(idList[indexOtherTopo[i]],
                    class_global)) 
            PrintMappedArray(mapArraycons, mapArrayquery,
                    "Consensus",idList[indexOtherTopo[i]], fpLog)

        indexClass[idx].append(indexOtherTopo[i])

        # for DIFF topology, do further analysis. Analyzing whether it has
        # N-terminal indels, C-terminal indels and the number of internal
        # indels
        cmprecord['cmpclass'] = class_global
        cmprecord['numTM1'] = numTMcons
        cmprecord['numTM2'] = numTMquery
        cmprecord['NtermTopo1'] = Ntermcons
        cmprecord['NtermTopo2'] = Ntermquery
        cmprecord['mapArray1'] = mapArraycons
        cmprecord['mapArray2'] = mapArrayquery
        cmprecord['id1'] = 'Cons'
        cmprecord['id2'] = seqID
        cmprecord['seqLength1'] = seqLenConsensus
        cmprecord['seqLength2'] = seqLenQuery
        if class_global == "DIFF":
            anaCons={}
            anaQuery={}
            anaCons=AnaDIFFTopology1(mapArraycons)
            anaQuery=AnaDIFFTopology1(mapArrayquery); 
#check with DG, topcons_single and gap 
# 1. check if the compared topology are mostly gaps at the aligned region
            anaCons  = CheckGapOfMSA(anaCons, posTMcons,gapSeqQuery)
            anaQuery = CheckGapOfMSA(anaQuery, posTMquery,gapSeqCons)
        
# 2. check the DG value of the TM region to see if the prediction is reliable
# we should probably also check the DG value of the opposite topology, but it
# is not implemented at this stage
            anaCons  = CheckDGOfMSA(anaCons, DGvalueTMcons)
            anaQuery = CheckDGOfMSA(anaQuery, DGvalueTMquery)
# 3. check with the predictions from topocons_single 
# not implemented yet
            cmprecord['ana1'] = anaCons
            cmprecord['ana2'] = anaQuery

    if DEBUG_CONSENSUS and fpLog != None:
        for icls in range(len(comparisonClassNameList)):
            print >> fpLog, ("%-10s" % (comparisonClassNameList[icls]),
                    indexClass[icls])

# Write overall information
    fpout.write("//Begin CMPMSA\n")
    lcmp.WriteOverallInfo_msa(comparisonClassNameList, g_params['rootname'],
            numSeq, numTM_IDT, numTMcons, numIDTTopo, indexClass, fpout)
# Write detailed topology variation info
    WriteDetailedDIFFTopo(cmpToConsRecordList, fpout)
    fpout.write("//End CMPMSA\n")

    return (cmpToConsRecordList, indexClass, infoDIFF)
#}}}
def CompareToConsensus2(topoSeqList, NtermStateList, #{{{
        indexIDTTopo,consensusTopo, gapSeqCons, idList, posTMList, dgScoreList,
        comparisonClassNameList):
    """ @param
        topoSeqList    Untrimmed (with gaps) aligned topology sequences
        NtermStateList N-terminal topology state of all sequences
        indexIDTTopo   Index of topology sequence in the largest identity group
        consensusTopo  The consensus topology built from the largest identity
                       group
        idList         Identities of all sequences
        posTMList      Location of TM helices in 2-tuple (beg, end) for all
                       untrimmed topology sequences
        dgScoreList    Delta G values of all TM regions of all sequences 
        ComparisonClassNameList 
                       Name list of comparison classes. e.g. 
                       ["OK","SHIFT","INV","DIFF" ]
    """

    fpLog = g_params['fpLog']
    fpout = g_params['fpout']
    DEBUG_CONSENSUS = g_params['DEBUG_CONSENSUS']

    numSeq = len(idList)
    numIDTTopo = len(indexIDTTopo)
    fulllist = [x for x in range(numSeq)]
    indexOtherTopo = list(set(fulllist)-set(indexIDTTopo)); # set subtraction
    numOtherTopo = len(indexOtherTopo)
    cmpClassOtherTopoList = [];# comparison class list

    posTMcons = myfunc.GetTMPosition(consensusTopo)
    numTMcons = len(posTMcons)

    # Get number of TM helices for the largest identity group
    numTM_IDT = len (posTMList[indexIDTTopo[0]]); 

    DGvalueTMcons = GetDGvalueTMconsensus(dgScoreList, indexIDTTopo, numTM_IDT)
    if g_params['isPrintDGCons']:
        print "DGvalueTMcons (%d) =" % len(DGvalueTMcons), ["%.3f"%x for x in
                DGvalueTMcons]
    numCmpClass = len(comparisonClassNameList); #number of defined classes

    # Get length of un-aligned sequences
    seqLenList = [ len(tp.replace(GAP,'')) for tp in topoSeqList]; 
    seqLenConsensus = len(consensusTopo.replace(GAP,''))

    indexClass=[]
    for i in range(numCmpClass):
        indexClass.append([])
    infoDIFF=[]

    cmpToConsRecordList = []

    Ntermcons = lcmp.GetNtermState(consensusTopo)

    for i in xrange(numOtherTopo):

        cmpToConsRecordList.append({})
        cmprecord = cmpToConsRecordList[i]

        seqID = idList[indexOtherTopo[i]]
        posTMquery = posTMList[indexOtherTopo[i]]
        numTMquery = len(posTMquery)
        topoquery = topoSeqList[indexOtherTopo[i]]
        lengthQuery = len(topoquery)
        seqLenQuery = seqLenList[i]
        DGvalueTMquery = [INIT_DGVALUE] * numTMquery
        dglist =  dgScoreList[indexOtherTopo[i]]; 
        if dglist and len(dglist) == numTMquery:
            DGvalueTMquery = dglist
        gapSeqQuery = [ (s==GAP) for s in topoquery]; # avoid loop, it is faster
        commonMarray = [ (topoquery[j] == 'M' and consensusTopo[j] == 'M') for
                j in range(lengthQuery)]; 

        idCons = "Consensus"
        idQuery =  idList[indexOtherTopo[i]]

        (mapArraycons, mapArrayquery) = MappingTM_method1(posTMcons,
                posTMquery, consensusTopo, topoquery, idCons, idQuery)
        Ntermquery = NtermStateList[indexOtherTopo[i]] 
        cmpclassSP = "" # comparison of signal peptide is not fulfilled for multiple comparison
        class_global = ClassifyTopoComparison_pairwise_method1(mapArraycons,
                mapArrayquery, Ntermcons, Ntermquery, idCons, idQuery,
                cmpclassSP)

        #print comparisonClassNameList
        class_global = class_global.split(";")[0]
        cmpClassOtherTopoList.append(class_global)
        idx = comparisonClassNameList.index(class_global)

        if DEBUG_CONSENSUS and fpLog != None:
            print >> fpLog, "SeqID", seqID, "class=", class_global, "idx=", idx
            print >> fpLog, ("SeqID %s: %s"%(idList[indexOtherTopo[i]],
                    class_global)) 
            PrintMappedArray(mapArraycons, mapArrayquery,
                    "Consensus",idList[indexOtherTopo[i]], fpLog)

        indexClass[idx].append(indexOtherTopo[i])

        cmprecord['cmpclass'] = class_global
        cmprecord['numTM1'] = numTMcons
        cmprecord['numTM2'] = numTMquery
        cmprecord['NtermTopo1'] = Ntermcons
        cmprecord['NtermTopo2'] = Ntermquery
        cmprecord['mapArray1'] = mapArraycons
        cmprecord['mapArray2'] = mapArrayquery
        cmprecord['id1'] = 'Cons'
        cmprecord['id2'] = idQuery
        cmprecord['seqLength1'] = seqLenConsensus
        cmprecord['seqLength2'] = seqLenQuery

    if DEBUG_CONSENSUS and fpLog != None:
        for icls in range(len(comparisonClassNameList)):
            print >> fpLog, ("%-10s" % (comparisonClassNameList[icls]),
                    indexClass[icls])

# Write overall information
    fpout.write("//Begin CMPMSA\n")
    lcmp.WriteOverallInfo_msa(comparisonClassNameList, g_params['rootname'],
            numSeq, numTM_IDT, numTMcons, numIDTTopo, indexClass, fpout)
# Write detailed topology variation info
    WritePairwiseRecord_method1(cmpToConsRecordList, fpout)
    fpout.write("//End CMPMSA\n")

    return (cmpToConsRecordList, indexClass, infoDIFF)
#}}}
def CompareToConsensus11(topoSeqList, indexIDTTopo,consensusTopo, idList, #{{{
        numTM_IDT, comparisonClassNameList):
    numSeq=len(idList)
    numIDTTopo=len(indexIDTTopo)
    fulllist=[x for x in range(numSeq)]
    indexOtherTopo=list(set(fulllist)-set(indexIDTTopo)); # set subtraction
    numOtherTopo=len(indexOtherTopo)
    numCmpClass=len(comparisonClassNameList); #number of defined classes
    indexClass=[]
    for i in range(numCmpClass):
        indexClass.append([])

    for i in range(numOtherTopo):
        topo=topoSeqList[indexOtherTopo[i]]
        class_global=None
        if topo[0] == consensusTopo[0]: 
            class_global="SAME_NTERM"
        else:
            class_global="DIFF_NTERM"

        idx = comparisonClassNameList.index(class_global)
        indexClass[idx].append(indexOtherTopo[i])
    return indexClass
#}}}
def CompareSignalPeptide(seqID1, seqID2, topo1, topo2, posTM1, posTM2,
        mapArray1, mapArray2, localseqpairDict, signalpDict):
# compare signal peptide
# Created 2013-09-04, updated 2013-09-04, Nanjiang Shu
# 
    unaligned_str = ""
    if localseqpairDict != {}:
        try:
            rd = localseqpairDict[(seqID1, seqID2)]
            unaligned_str = rd[2]
        except KeyError:
            msg = "Failed to find local alignment for (%s,%s)"
            print >> sys.stderr, msg%(seqID1, seqID2)
            unaligned_str = ""
    sp1 = -1
    sp2 = -1
    try:
        sp1 = signalpDict[seqID1]
    except KeyError:
        sp1 = -1
    try:
        sp2 = signalpDict[seqID2]
    except KeyError:
        sp2 = -1

    cls = ""
    cls1 = ""
    cls2 = ""
    if sp1 == -1 and sp2 == -1:
        cls = "noSP"
    else:
        if sp1 >= 0:
            cls1 = MapAlignedSP(topo1, topo2, posTM1, posTM2, mapArray1,
                    mapArray2, sp1, sp2, unaligned_str, seqID1, seqID2)
        if sp2 >= 0:
            cls2 = MapAlignedSP(topo2, topo1, posTM2, posTM1, mapArray2,
                    mapArray1, sp2, sp1, unaligned_str, seqID2, seqID1)

        li = [cls1, cls2]
        li = myfunc.uniquelist(li)
        li = filter(None, li)
        cls = "|".join(li)
    return cls



def GetGroupedTopoMSA(numSeq, idList, topoSeqList, dgScoreList, #{{{
        cmpToConsRecordList, indexIDTTopo, consensusTopo, Mcmp, indexClass,
        comparisonClassNameList):

    fpLog = g_params['fpLog']
    DEBUG_GROUPING = g_params['DEBUG_GROUPING']

    numSeq = len(idList)
    numIDTTopo = len(indexIDTTopo)
    fulllist = [x for x in range(numSeq)]
    indexOtherTopo = list(set(fulllist)-set(indexIDTTopo)); # set subtraction
    numOtherTopo = len(indexOtherTopo)
#check
    if len(cmpToConsRecordList) != numOtherTopo:
        msg = "Error! len(cmpToConsRecordList) ({}) != numOtherTopo ({})",
        print >> sys.stderr, msg.format(len(cmpToConsRecordList), numOtherTopo)
    groupList = []
    cnt = 0
#groupList[0] is the largest identical group
    groupList.append({})
    grp = groupList[cnt]
    grp['cmpclass'] = "IDT"; 
    grp['reptopo'] = consensusTopo; # representative topology
    grp['repID'] = "Cons"; # representative ID
    grp['repNumTM'] = myfunc.CountTM(consensusTopo)
    grp['index-members'] = indexIDTTopo; # index of the members of the group
    cnt += 1

# map the index in the original all sequences to that in cmpToConsRecordList
    idxFull2OtherClass = {}; # a dictionary
    for i in range(len(indexOtherTopo)):
        idxFull2OtherClass[indexOtherTopo[i]] = i
# now for the topologies in the other classes
#     print "idxFull2OtherClass=", idxFull2OtherClass
#     print "indexOtherTopo=", indexOtherTopo

    maxDGdifference = g_params['maxDGdifference']
    if fpLog and DEBUG_GROUPING:
        print >> fpLog, "maxDGdifference = ", maxDGdifference
    for icls in range(len(comparisonClassNameList)):
        indexThisClass = indexClass[icls]
        if len(indexThisClass) < 1: 
            continue
        indexThisClassSet = set(indexThisClass)
# 1. For each topology in this class, find its similar topologies
#    if they are Identical topology according to Mcmp and 
#    if their DG values differ by <= maxDGdifference
# 2. Add the largest similar group to groupList 
#    do 
#        add the next largest similar group to groupList
#        until all items in this class are grouped:

        similarTopoForEach = []
        for pivot in indexThisClass:
            lst = []
            lst.append(pivot)
            for j in indexThisClass:
                if j == pivot:
                    continue
                if (Mcmp[pivot][j] == 1 and
                        IsDGScoreSimilar(dgScoreList[pivot], dgScoreList[j],
                            maxDGdifference)):
                    lst.append(j)
            similarTopoForEach.append(lst)

        groupedIndexSet = set([]); 
        while 1:
            if len(groupedIndexSet) == len(indexThisClass):
                break
            numSimilarTopoForEachList = [len(l) for l in similarTopoForEach]
            idxLargestGroup = numSimilarTopoForEachList.index(
                    max(numSimilarTopoForEachList))
            largestGroup = similarTopoForEach[idxLargestGroup]
            subDGScoreList = [sum(dgScoreList[j]) for j in largestGroup]
            idxRep = largestGroup[subDGScoreList.index(min(subDGScoreList))]
            groupList.append({})
            grp = groupList[cnt]
#         print "indexGrpMemberList=",indexGrpMemberList
#         print "idxRep=", idxRep
            cmprecord = cmpToConsRecordList[idxFull2OtherClass[idxRep]] 
            grp['cmpclass'] = comparisonClassNameList[icls]
            grp['reptopo'] = topoSeqList[idxRep]; # representative topology
            grp['repID'] = idList[idxRep]; # representative ID
            grp['repIndex'] = idxRep; # index of the representative
            grp['repNumTM'] = myfunc.CountTM(topoSeqList[idxRep])
            # index of the members of the group, including the representative
            grp['index-members'] = largestGroup
            if len(largestGroup) > 1:
                grp['index-members-without-rep'] = [x for x in largestGroup
                    if x != idxRep]
            else:
                grp['index-members-without-rep'] = []
            cnt += 1
            groupedIndexSet.update(largestGroup)

            if fpLog and DEBUG_GROUPING:
                print >> fpLog, "cmpclass=", comparisonClassNameList[icls]
                for i in range(len(similarTopoForEach)):
                    print >> fpLog, "similarTopoForEach[%d]"%i, similarTopoForEach[i]
                print >> fpLog, "largestGroup=", largestGroup
                print >> fpLog

# remove the items that are already grouped
            largestGroupSet = set(largestGroup)
            for i in range(len(similarTopoForEach)):
                similarTopoForEach[i] = list(set(similarTopoForEach[i]) -
                        largestGroupSet)

    return groupList
#}}}
def PairwiseTopologyComparison(topoRecordList, g_params):#{{{
    cmpmethod = g_params['pairwise_comparison_method']
    if cmpmethod == 0:
        return PairwiseTopologyComparison_method0(topoRecordList, g_params)
    elif cmpmethod == 1:
        return PairwiseTopologyComparison_method1(topoRecordList, g_params)
    elif cmpmethod == 2: # no finished
        return PairwiseTopologyComparison_method2(topoRecordList, g_params)
    elif cmpmethod == 3:
        return PairwiseTopologyComparison_method3(topoRecordList, g_params)
    else:
        print >> sys.stderr, "Wrong method %d"%(cmpmethod)
        return 1
#}}}

def PairwiseTopologyComparison_method0(topoRecordList, g_params):#{{{
    """
    Do pairwise topology comparison.
    @params
        topoRecordList  A list of topology records of n-tuple
                        (seqID, anno, seq, seqIdentity, dgscore)
        processed topoRecored will be deleted from topoRecordList
    """

    fpout = g_params['fpout']
    fpLog = g_params['fpLog']
    dupPairSet = g_params['dupPairSet']
    signalpDict = g_params['signalpDict']

    DEBUG_TMMAPPING = g_params['DEBUG_TMMAPPING']

    numSeq = len(topoRecordList)
    idList = [r[0] for r in topoRecordList ]
    topoSeqList = [r[2] for r in topoRecordList ]
    seqIdentityList = [r[3] for r in topoRecordList ]
    dgScoreList = [r[4] for r in topoRecordList ]

    posTMList = [ myfunc.GetTMPosition(topo) for topo in topoSeqList]
    numTMList = [len(posTM) for posTM in posTMList]
    alignSeqLenList = [len(topo) for topo in topoSeqList ]

    # Get length of un-aligned sequences
    seqLenList=[ len(tp.replace(GAP,'')) for tp in topoSeqList]; 
    # Get topology at N terminal of each sequence 
    NtermStateList=[ lcmp.GetNtermState(topo) for topo in topoSeqList]; 
    numPair = numSeq / 2

    DGvalueTMList = []
    gapSeqList=[]
    for i in xrange(numSeq):
        topo = topoSeqList[i]
        dglist=dgScoreList[i]
        numTM = numTMList[i]
        if dglist and len(dglist) == numTM:
            DGvalueTMList.append(dglist)
        else:
            DGvalueTMList.append([INIT_DGVALUE]*numTM)
        gapSeqList.append ([(s==GAP) for s in topo]); 


    for i in xrange (numPair):
        idx1 = 2*i
        idx2 = 2*i+1
        if alignSeqLenList[idx1] != alignSeqLenList[idx2]:
            print >> sys.stderr, ("Length of aligned topology for %s (%d) " % 
                    (idList[idx1], alignSeqLenList[idx1]) 
                    + "and %s (%d) not equal.  Ignore." 
                    % (idList[idx2], alignSeqLenList[idx2]))
            continue
        if numTMList[idx1] < 1 or numTMList[idx2] < 1:
            print >> sys.stderr, ("numTM of aligned topology for %s (%d) or %s (%d) < 1. Ignore." % (idList[idx1], numTMList[idx1],
                        idList[idx2], numTMList[idx2]))
            continue
        if seqIdentityList[idx1] != seqIdentityList[idx2]:
            print >> sys.stderr, "%s - %s: sequence identity inconsistent. %.2f - %.2f. Ignore."%(idList[idx1], idList[idx2], seqIdentityList[idx1], seqIdentityList[idx2])
            continue
        commonMarray = [ (topoSeqList[idx1][j] == 'M' and 
            topoSeqList[idx2][j] == 'M') for j in
            range(alignSeqLenList[idx1])]; 

        topo1 = topoSeqList[idx1]
        topo2 = topoSeqList[idx2]
        seqLength1 = seqLenList[idx1]
        seqLength2 = seqLenList[idx2]
        (mapArray1, mapArray2) = MappingTM(posTMList[idx1], posTMList[idx2],
                commonMarray, idList[idx1], idList[idx2])
        class_global = ClassifyTopoComparison_pairwise(mapArray1, mapArray2,
                NtermStateList[idx1], NtermStateList[idx2],posTMList[idx1],
                posTMList[idx2], seqLength1, seqLength2, idList[idx1], idList[idx2],
                dupPairSet, signalpDict)

        if DEBUG_TMMAPPING and fpLog != None:
            PrintMappedArray(mapArray1, mapArray2, idList[idx1], idList[idx2],
                    fpLog)
        # for DIFF topology, do further analysis. Analyzing whether it has
        # N-terminal indels, C-terminal indels and the number of internal
        # indels
        ana1={}
        ana2={}
        if class_global=="DIFF":
            ana1=AnaDIFFTopology1(mapArray1)
            ana2=AnaDIFFTopology1(mapArray2); 
#check with DG, topcons_single and gap 
# 1. check if the compared topology are mostly gaps at the aligned region
            ana1 = CheckGapOfMSA(ana1, posTMList[idx1],gapSeqList[idx2])
            ana2 = CheckGapOfMSA(ana2, posTMList[idx2],gapSeqList[idx1])

# 2. check the DG value of the TM region to see if the prediction is reliable
# we should probably also check the DG value of the opposite topology, but it
# is not implemented at this stage
            ana1 = CheckDGOfMSA(ana1, DGvalueTMList[idx1])
            ana2 = CheckDGOfMSA(ana2, DGvalueTMList[idx2])
# 3. check with the predictions from topocons_single 
# not implemented yet

# Write overall information
        print >> fpout, "//Begin record", g_params['cntOutputPair']+1
        lcmp.WriteOverallInfo_pairwise(idList[idx1], idList[idx2],
                seqIdentityList[idx1],  class_global, numTMList[idx1],
                numTMList[idx2], seqLenList[idx1], seqLenList[idx2],  fpout,
                g_params['uniprot2pdbMap'], g_params['swissprotAcSet'])
        PrintMappedArray(mapArray1, mapArray2, idList[idx1], idList[idx2],
                fpout)
        # for DIFF topology, do further analysis. Analyzing whether it has

# Write detailed topology variation info
        if  ana1 != {}:
            lcmp.WriteAna(ana1,fpout, "1") 
        if  ana2 != {}:
            lcmp.WriteAna(ana2,fpout, "2") 

        print >> fpout, "//End record", g_params['cntOutputPair']+1
        g_params['cntOutputPair'] += 1

# remove records that are already compared
    for i in xrange (numPair):
        topoRecordList.pop(0)
        topoRecordList.pop(0)

    return 0
#}}}
def PairwiseTopologyComparison_method1(topoRecordList, g_params):#{{{
    """
    Do pairwise topology comparison.
    @params
        topoRecordList  A list of topology records of n-tuple
                        (seqID, anno, seq, seqIdentity, dgscore)
        processed topoRecored will be deleted from topoRecordList
    pairwise topology comparisons are classified according to 
    a) TM - TM
    b) TM - Seq  (SIGNALP should be included in this)
    c) TM - Gap  (duplications should be included in this)
    categories:
    1) Identical 
    2) Inverted
    3) Only b)
    4) Only c)
    5) b) + c)
    """

    fpout = g_params['fpout']
    fpLog = g_params['fpLog']
    dupPairSet = g_params['dupPairSet']
    signalpDict = g_params['signalpDict']
    localseqpairDict = g_params['localseqpairDict']

    DEBUG_TMMAPPING = g_params['DEBUG_TMMAPPING']

    numSeq = len(topoRecordList)
    idList = [r[0] for r in topoRecordList ]
    topoSeqList = [r[2] for r in topoRecordList ]
    seqIdentityList = [r[3] for r in topoRecordList ]
    dgScoreList = [r[4] for r in topoRecordList ]

    posTMList = [ myfunc.GetTMPosition(topo) for topo in topoSeqList]
    numTMList = [len(posTM) for posTM in posTMList]
    alignSeqLenList = [len(topo) for topo in topoSeqList ]

    # Get length of un-aligned sequences
    seqLenList=[ len(tp.replace(GAP,'')) for tp in topoSeqList]; 
    # Get topology at N terminal of each sequence 
    NtermStateList=[ lcmp.GetNtermState(topo) for topo in topoSeqList]; 
    numPair = numSeq / 2

    DGvalueTMList = []
    gapSeqList=[]
    for i in xrange(numSeq):
        topo = topoSeqList[i]
        dglist=dgScoreList[i]
        numTM = numTMList[i]
        if dglist and len(dglist) == numTM:
            DGvalueTMList.append(dglist)
        else:
            DGvalueTMList.append([INIT_DGVALUE]*numTM)
        gapSeqList.append ([(s==GAP) for s in topo]); 


    for i in xrange (numPair):
        idx1 = 2*i
        idx2 = 2*i+1
        if alignSeqLenList[idx1] != alignSeqLenList[idx2]:
            print >> sys.stderr, ("Length of aligned topology for %s (%d) " % 
                    (idList[idx1], alignSeqLenList[idx1]) 
                    + "and %s (%d) not equal.  Ignore." 
                    % (idList[idx2], alignSeqLenList[idx2]))
            continue
        if numTMList[idx1] < 1 or numTMList[idx2] < 1:
            print >> sys.stderr, ("numTM of aligned topology for %s (%d) "%(
                idList[idx1], numTMList[idx1]) +
                "or %s (%d) < 1. Ignore." % (idList[idx2], numTMList[idx2]))
            continue
        if seqIdentityList[idx1] != seqIdentityList[idx2]:
            print >> sys.stderr, ("%s - %s: sequence identity inconsistent. "%(
                idList[idx1], idList[idx2]) +
                "%.2f - %.2f. Ignore."%(seqIdentityList[idx1],
                    seqIdentityList[idx2]))
            continue
        commonMarray = [ (topoSeqList[idx1][j] == 'M' and 
            topoSeqList[idx2][j] == 'M') for j in
            range(alignSeqLenList[idx1])]; 

        topo1 = topoSeqList[idx1]
        topo2 = topoSeqList[idx2]
        seqLength1 = seqLenList[idx1]
        seqLength2 = seqLenList[idx2]
        (mapArray1, mapArray2) = MappingTM_method1(posTMList[idx1], posTMList[idx2],
                topo1, topo2, idList[idx1], idList[idx2])

        cmpclassSP = ""
        if g_params['isCompareSP'] == True:
            cmpclassSP = CompareSignalPeptide(idList[idx1], idList[idx2],
                    topo1, topo2, posTMList[idx1], posTMList[idx2],
                    mapArray1, mapArray2, localseqpairDict, signalpDict)


        class_global = ClassifyTopoComparison_pairwise_method1(mapArray1, mapArray2,
                NtermStateList[idx1], NtermStateList[idx2], idList[idx1],
                idList[idx2], cmpclassSP)

        # for DIFF topology, do further analysis. Analyzing whether it has
        # N-terminal indels, C-terminal indels and the number of internal
        # indels
        ana1={}
        ana2={}
        if class_global != "":
# Write overall information
            print >> fpout, "//Begin record", g_params['cntOutputPair']+1
            lcmp.WriteOverallInfo_pairwise(idList[idx1], idList[idx2],
                    seqIdentityList[idx1],  class_global, numTMList[idx1],
                    numTMList[idx2], seqLenList[idx1], seqLenList[idx2],  fpout,
                    g_params['uniprot2pdbMap'], g_params['swissprotAcSet'])
            PrintMappedArray_method1_1(mapArray1, mapArray2, idList[idx1], idList[idx2],
                    fpout)
            # for DIFF topology, do further analysis. Analyzing whether it has

# Write detailed topology variation info
            if  ana1 != {}:
                lcmp.WriteAna(ana1,fpout, "1") 
            if  ana2 != {}:
                lcmp.WriteAna(ana2,fpout, "2") 

            print >> fpout, "//End record", g_params['cntOutputPair']+1
            g_params['cntOutputPair'] += 1

# remove records that are already compared
    for i in xrange (numPair):
        topoRecordList.pop(0)
        topoRecordList.pop(0)

    return 0
#}}}
def PairwiseTopologyComparison_method2(topoRecordList, g_params):#{{{
    """
    Do pairwise topology comparison.
    not finished

    Distinguish TM helices from i->o and o->i, denote i->o as M and o->i as W

    @params
        topoRecordList  A list of topology records of n-tuple
                        (seqID, anno, seq, seqIdentity, dgscore)
        processed topoRecored will be deleted from topoRecordList
    pairwise topology comparisons are classified according to 
    a) TM - TM
    b) TM - Seq  (SIGNALP should be included in this)
    c) TM - Gap  (duplications should be included in this)
    categories:
    1) Identical 
    2) Inverted
    3) Only b)
    4) Only c)
    5) b) + c)
    """

    fpout = g_params['fpout']
    fpLog = g_params['fpLog']
    dupPairSet = g_params['dupPairSet']
    signalpDict = g_params['signalpDict']

    DEBUG_TMMAPPING = g_params['DEBUG_TMMAPPING']

    numSeq = len(topoRecordList)
    idList = [r[0] for r in topoRecordList ]
    topoSeqList = [r[2] for r in topoRecordList ]
    seqIdentityList = [r[3] for r in topoRecordList ]
    dgScoreList = [r[4] for r in topoRecordList ]

    posTMList = [ myfunc.GetTMPosition(topo) for topo in topoSeqList]
    numTMList = [len(posTM) for posTM in posTMList]
    alignSeqLenList = [len(topo) for topo in topoSeqList ]

    # Get length of un-aligned sequences
    seqLenList=[ len(tp.replace(GAP,'')) for tp in topoSeqList]; 
    # Get topology at N terminal of each sequence 
    NtermStateList=[ lcmp.GetNtermState(topo) for topo in topoSeqList]; 
    numPair = numSeq / 2

    DGvalueTMList = []
    gapSeqList=[]
    for i in xrange(numSeq):
        topo = topoSeqList[i]
        dglist=dgScoreList[i]
        numTM = numTMList[i]
        if dglist and len(dglist) == numTM:
            DGvalueTMList.append(dglist)
        else:
            DGvalueTMList.append([INIT_DGVALUE]*numTM)
        gapSeqList.append ([(s==GAP) for s in topo]); 


    for i in xrange (numPair):
        idx1 = 2*i
        idx2 = 2*i+1
        if alignSeqLenList[idx1] != alignSeqLenList[idx2]:
            print >> sys.stderr, ("Length of aligned topology for %s (%d) " % 
                    (idList[idx1], alignSeqLenList[idx1]) 
                    + "and %s (%d) not equal.  Ignore." 
                    % (idList[idx2], alignSeqLenList[idx2]))
            continue
        if numTMList[idx1] < 1 or numTMList[idx2] < 1:
            print >> sys.stderr, ("numTM of aligned topology for %s (%d) "%(
                idList[idx1], numTMList[idx1]) +
                "or %s (%d) < 1. Ignore." % (idList[idx2], numTMList[idx2]))
            continue
        if seqIdentityList[idx1] != seqIdentityList[idx2]:
            print >> sys.stderr, ("%s - %s: sequence identity inconsistent. "%(
                idList[idx1], idList[idx2]) +
                "%.2f - %.2f. Ignore."%(seqIdentityList[idx1],
                    seqIdentityList[idx2]))
            continue
        commonMarray = [ (topoSeqList[idx1][j] == 'M' and 
            topoSeqList[idx2][j] == 'M') for j in
            range(alignSeqLenList[idx1])]; 

        topo1 = topoSeqList[idx1]
        topo2 = topoSeqList[idx2]
        seqLength1 = seqLenList[idx1]
        seqLength2 = seqLenList[idx2]
        (mapArray1, mapArray2) = MappingTM_method1(posTMList[idx1], posTMList[idx2],
                topo1, topo2, idList[idx1], idList[idx2])
        class_global = ClassifyTopoComparison_pairwise_method1(mapArray1, mapArray2,
                NtermStateList[idx1], NtermStateList[idx2], idList[idx1],
                idList[idx2])

        # for DIFF topology, do further analysis. Analyzing whether it has
        # N-terminal indels, C-terminal indels and the number of internal
        # indels
        ana1={}
        ana2={}
        if class_global != "":
# Write overall information
            print >> fpout, "//Begin record", g_params['cntOutputPair']+1
            lcmp.WriteOverallInfo_pairwise(idList[idx1], idList[idx2],
                    seqIdentityList[idx1],  class_global, numTMList[idx1],
                    numTMList[idx2], seqLenList[idx1], seqLenList[idx2],  fpout,
                    g_params['uniprot2pdbMap'], g_params['swissprotAcSet'])
            PrintMappedArray_method1_1(mapArray1, mapArray2, idList[idx1], idList[idx2],
                    fpout)
            # for DIFF topology, do further analysis. Analyzing whether it has

# Write detailed topology variation info
            if  ana1 != {}:
                lcmp.WriteAna(ana1,fpout, "1") 
            if  ana2 != {}:
                lcmp.WriteAna(ana2,fpout, "2") 

            print >> fpout, "//End record", g_params['cntOutputPair']+1
            g_params['cntOutputPair'] += 1

# remove records that are already compared
    for i in xrange (numPair):
        topoRecordList.pop(0)
        topoRecordList.pop(0)

    return 0
#}}}
def PairwiseTopologyComparison_method3(topoRecordList, g_params):#{{{
    """
    Do pairwise topology comparison.
    @params
        topoRecordList  A list of topology records of n-tuple
                        (seqID, anno, seq, seqIdentity, dgscore)
        processed topoRecored will be deleted from topoRecordList
    pairwise topology comparisons are classified according to 
    a) TM - TM
    b) TM - GAP
    c) TM - SEQ
    d) TM - SP
    categories:
    1) IDT 
    2) INV
    3) TM2GAP
    4) TM2SEQ
    5) TM2SP       # the variation happened only by aligning a TM helix to signal peptide
    6) Mixed
    """

    fpout = g_params['fpout']
    fpLog = g_params['fpLog']
    dupPairSet = g_params['dupPairSet']
    signalpDict = g_params['signalpDict']
    pairalnStat = g_params['pairalnStat']

    DEBUG_TMMAPPING = g_params['DEBUG_TMMAPPING']

    numSeq = len(topoRecordList)
    idList = [r[0] for r in topoRecordList ]
    topoSeqList = [r[2] for r in topoRecordList ]
    seqIdentityList = [r[3] for r in topoRecordList ]
    dgScoreList = [r[4] for r in topoRecordList ]

    posTMList = [ myfunc.GetTMPosition(topo) for topo in topoSeqList]
    numTMList = [len(posTM) for posTM in posTMList]
    alignSeqLenList = [len(topo) for topo in topoSeqList ]

    # Get length of un-aligned sequences
    seqLenList=[ len(tp.replace(GAP,'')) for tp in topoSeqList]; 
    # Get topology at N terminal of each sequence 
    NtermStateList=[ lcmp.GetNtermState(topo) for topo in topoSeqList]; 
    numPair = numSeq / 2

    DGvalueTMList = []
    gapSeqList=[]
    for i in xrange(numSeq):
        topo = topoSeqList[i]
        dglist=dgScoreList[i]
        numTM = numTMList[i]
        if dglist and len(dglist) == numTM:
            DGvalueTMList.append(dglist)
        else:
            DGvalueTMList.append([INIT_DGVALUE]*numTM)
        gapSeqList.append ([(s==GAP) for s in topo]); 


    for i in xrange (numPair):
        idx1 = 2*i
        idx2 = 2*i+1
        if alignSeqLenList[idx1] != alignSeqLenList[idx2]:
            print >> sys.stderr, ("Length of aligned topology for %s (%d) " % 
                    (idList[idx1], alignSeqLenList[idx1]) 
                    + "and %s (%d) not equal.  Ignore." 
                    % (idList[idx2], alignSeqLenList[idx2]))
            continue
        if numTMList[idx1] < 1 or numTMList[idx2] < 1:
            print >> sys.stderr, ("numTM of aligned topology for %s (%d) "%(
                idList[idx1], numTMList[idx1]) +
                "or %s (%d) < 1. Ignore." % (idList[idx2], numTMList[idx2]))
            continue
        if seqIdentityList[idx1] != seqIdentityList[idx2]:
            print >> sys.stderr, ("%s - %s: sequence identity inconsistent. "%(
                idList[idx1], idList[idx2]) +
                "%.2f - %.2f. Ignore."%(seqIdentityList[idx1],
                    seqIdentityList[idx2]))
            continue
        commonMarray = [ (topoSeqList[idx1][j] == 'M' and 
            topoSeqList[idx2][j] == 'M') for j in
            range(alignSeqLenList[idx1])]; 

        topo1 = topoSeqList[idx1]
        topo2 = topoSeqList[idx2]
        seqLength1 = seqLenList[idx1]
        seqLength2 = seqLenList[idx2]
        seq2alignMap1 = lcmp.GetSeq2AlignMap(topo1.replace(GAP,""), topo1)
        seq2alignMap2 = lcmp.GetSeq2AlignMap(topo2.replace(GAP,""), topo2)

        # get sequence identity
        if pairalnStat != {}:
            keystr = idList[idx1] + '-' + idList[idx2]
            try:
                if g_params['seqidttype'] == 0:
                    seqidt = pairalnStat[keystr]['seqidt']
                elif g_params['seqidttype'] == 1:
                    seqidt = pairalnStat[keystr]['seqidt1']
                elif g_params['seqidttype'] == 2:
                    seqidt = pairalnStat[keystr]['seqidt2']
                else:
                    seqidt = pairalnStat[keystr]['seqidt']
            except KeyError:
                print >> sys.stderr, "%s not found in tableinfo"%(keystr)
                seqidt = INIT_SEQUENCE_IDENTITY
        else:
            seqidt = seqIdentityList[idx1]

        # signal peptide
        try:
            sp_pos1 = signalpDict[idList[idx1]]
        except KeyError:
            sp_pos1 = -1
            pass
        try:
            sp_pos2 = signalpDict[idList[idx2]]
        except KeyError:
            sp_pos2 = -1
            pass

        if sp_pos1 != -1:
            try:
                sp_pos1 = seq2alignMap1[sp_pos1]
            except KeyError:
                print >> sys.stderr, "idx %d not in seq2alignMap1 for %s"%(sp_pos1, idList[idx1])
                pass
        if sp_pos2 != -1:
            try:
                sp_pos2 = seq2alignMap2[sp_pos2]
            except KeyError:
                print >> sys.stderr, "idx %d not in seq2alignMap1 for %s"%(sp_pos2, idList[idx2])
                pass

        (mapArray1, mapArray2) = MappingTM_method3(posTMList[idx1], posTMList[idx2],
                topo1, topo2, sp_pos1, sp_pos2, idList[idx1], idList[idx2])
        class_global = ClassifyTopoComparison_pairwise_method3(mapArray1, mapArray2,
                NtermStateList[idx1], NtermStateList[idx2], idList[idx1],
                idList[idx2])

        # for DIFF topology, do further analysis. Analyzing whether it has
        # N-terminal indels, C-terminal indels and the number of internal
        # indels
        ana1={}
        ana2={}
        if class_global != "":
# Write overall information
            print >> fpout, "//Begin record", g_params['cntOutputPair']+1
            lcmp.WriteOverallInfo_pairwise(idList[idx1], idList[idx2], seqidt,
                    class_global, numTMList[idx1], numTMList[idx2],
                    seqLenList[idx1], seqLenList[idx2], fpout,
                    g_params['uniprot2pdbMap'], g_params['swissprotAcSet'])
            PrintMappedArray_method1_1(mapArray1, mapArray2, idList[idx1], idList[idx2],
                    fpout)
            # for DIFF topology, do further analysis. Analyzing whether it has

# Write detailed topology variation info
            if  ana1 != {}:
                lcmp.WriteAna(ana1,fpout, "1") 
            if  ana2 != {}:
                lcmp.WriteAna(ana2,fpout, "2") 

            print >> fpout, "//End record", g_params['cntOutputPair']+1
            g_params['cntOutputPair'] += 1

# remove records that are already compared
    for i in xrange (numPair):
        topoRecordList.pop(0)
        topoRecordList.pop(0)

    return 0
#}}}
def MultipleTopologyComparison(topoRecordList, g_params):#{{{
    print "method_topology_comparison=", g_params['method_topology_comparison']
    if g_params['method_topology_comparison'] == 0:
        return MultipleTopologyComparison_alignment(topoRecordList, g_params)
    elif g_params['method_topology_comparison'] == 9:
        return MultipleTopologyComparison_numTM(topoRecordList, g_params)
    else:
        print >> sys.stderr, "Wrong method_topology_comparison = %d"%(
                g_params['method_topology_comparison'])
        return 1
#}}}
def MultipleTopologyComparison_numTM(topoRecordList, g_params):#{{{
    """
    Compare multiply aligned topologies by number of TM helices
    """
    numSeq = len(topoRecordList)
    idList = [r[0] for r in topoRecordList ]
    topoSeqList = [r[2] for r in topoRecordList ]
    posTMList = [myfunc.GetTMPosition(topo) for topo in topoSeqList]
    numTMList = [myfunc.CountTM(topo) for topo in topoSeqList]

    Mcmp= [[0]*numSeq for x in xrange(numSeq)]
    for i in range(numSeq): 
        Mcmp[i][i] = 1
    # All-to-all comparison
    for i in range (numSeq):
        numTM1 = numTMList[i]
        for j in range(i+1, numSeq):
            numTM2 = numTMList[j]
            if numTM1 == numTM2:
                Mcmp[i][j] = 1
            else:
                Mcmp[i][j] = 0
    #copy the symmetric matrix
    for i in range (numSeq):
        for j in range(i+1, numSeq):
            Mcmp[j][i]=Mcmp[i][j]

    rootname = os.path.basename(os.path.splitext(g_params['inFile'])[0])
    if g_params['outpath'] != "":
        outpath = g_params['outpath']
    else:
        outpath = myfunc.my_dirname(g_params['inFile'])
    outSortedClusteredTopoMSAFile = "%s/%s.cluster.m%d.topomsa.fa"%(outpath,
            rootname, g_params['method_topology_comparison'])
    WriteSortedClusteredTopoMSA(outSortedClusteredTopoMSAFile, idList,
            topoSeqList, posTMList, Mcmp)

#}}}
def MultipleTopologyComparison_alignment(topoRecordList, g_params):#{{{
#   Check if it is a multiple sequence alignment by checking length
    """
    Compare the topology of multiple sequence alignemnt
    """
    fpLog = g_params['fpLog']
    fpout = g_params['fpout']
    fpoutGrouped = g_params['fpoutGrouped']; 
    numSeq = len(topoRecordList)
    topoSeqList = [r[2] for r in topoRecordList ]

    alignSeqLenList = [len(topo) for topo in topoSeqList ]
    if min(alignSeqLenList) != max(alignSeqLenList):
        print >> sys.stderr, "Length of aligned topologies are not equal."
        print >> sys.stderr, "You should probably run pairwise comparison."
        return -1

    idList = [r[0] for r in topoRecordList ]
    seqIdentityList = [r[3] for r in topoRecordList ]
    dgScoreList = [r[4] for r in topoRecordList ]
    posTMList=[myfunc.GetTMPosition(topo) for topo in topoSeqList]
    #topology at N terminal of each sequence 
    NtermStateList=[ lcmp.GetNtermState(topo) for topo in topoSeqList]; 

    if g_params['method_getIDTgroup'] == 0:
        if not IsTrimmedMSA(topoSeqList):
            print >> sys.stderr,("MSA are with gaps, when method_getIDTgroup" 
                    + "= 0, they should be trimmed. Exit.")
            return 1
        (consensusTopo, indexIDTTopo,  numTM_IDT, 
                Mcmp) = GetIDTTopoGroup(idList,topoSeqList)
    elif g_params['method_getIDTgroup'] == 1:
        (consensusTopo, indexIDTTopo,  gapSeqCons, numTM_IDT, 
                Mcmp) = GetIDTTopoGroup1(topoSeqList, NtermStateList,
                        posTMList)
    else:
        sys.stderr.write(("Sorry! method_getIDTgroup %d" %
            g_params['method_getIDTgroup']) + " not implemented.  Exit.\n" )
        return 1
    # number of sequences for the largest identity group
    numIDTTopo = len(indexIDTTopo); 
    numTMCons = myfunc.CountTM(consensusTopo)
    #print "numIDTTopo=", numIDTTopo

    if numTMCons <= 0:
    #if numTMCons != numTM_IDT:
        print >> sys.stderr, "Error! numTMCons (%d) <= 0. Exit." %numTMCons
        print >> sys.stderr, "\"%s\""%consensusTopo
        return 1

    method_comparison = g_params['method_comparison']
    if method_comparison == 0:
        if not IsTrimmedMSA(topoSeqList):
            print >> sys.stderr, ("method_comparison accept only trimmed "
                    + "topology MSA, but topology MSA in file %s " %
                    g_params['inFile'] 
                    + "are with gaps. Exit.")
            return 1
        (indexClass, infoDIFF) = CompareToConsensus(topoSeqList,
                indexIDTTopo,consensusTopo, idList, numTM_IDT,
                comparisonClassNameListAll[method_comparison])
    elif method_comparison == 1:
        (cmpToConsRecordList, indexClass, infoDIFF) = CompareToConsensus1(
                topoSeqList, NtermStateList, indexIDTTopo,consensusTopo,
                gapSeqCons,idList, posTMList, dgScoreList,
                comparisonClassNameListAll[method_comparison])
    elif method_comparison == 2:
        (cmpToConsRecordList, indexClass, infoDIFF) = CompareToConsensus2(
                topoSeqList, NtermStateList, indexIDTTopo,consensusTopo,
                gapSeqCons,idList, posTMList, dgScoreList,
                comparisonClassNameListAll[method_comparison])
    elif method_comparison == 11:
# this is a quick comparison, just compare if the other topology has the same
# Nterm topology. 2011-09-22
        (indexClass, infoDIFF) = CompareToConsensus11(topoSeqList,
                indexIDTTopo, consensusTopo, idList, numTM_IDT,
                comparisonClassNameListAll[method_comparison])
        lcmp.WriteOverallInfo_msa(comparisonClassNameListAll[method_comparison],
                g_params['rootname'], numSeq, numTM_IDT,
                myfunc.CountTM(consensusTopo), numIDTTopo, indexClass, fpout)
    else:
        sys.stderr.write("Sorry! Comparison method \"%d\"" %
                (g_params['method_comparison'])
                +" has not yet been developed. Exit!\n")
        return 1

    lengthAlignment=len(topoSeqList[0])

    if g_params['isPrintIndexIDTGroup']: 
        print ("There are %d sequences in the IDTgroup, indeces are:" %
                len(indexIDTTopo))
        for i in indexIDTTopo:
            print i

    if g_params['outSortedOrigTopoMSAFile'] != "":
        WriteSortedOrigTopoMSA(g_params['outSortedOrigTopoMSAFile'], idList,
                topoSeqList, posTMList, numIDTTopo, consensusTopo, numTM_IDT,
                indexIDTTopo, indexClass,
                comparisonClassNameListAll[method_comparison], infoDIFF)
    if g_params['outSortedOrigTopoMSAFileHTML'] != "":
        WriteSortedOrigTopoMSAHTML(g_params['outSortedTrimmedTopoMSAFileHTML'],
                topoSeqList, lengthAlignment,idList, topoSeqList,
                indexIDTTopo)
    if g_params['outSortedClusteredTopoMSAFile'] != "":
        WriteSortedClusteredTopoMSA(g_params['outSortedClusteredTopoMSAFile'],
                idList, topoSeqList, posTMList, Mcmp)

    if (g_params['outGroupedSortedOrigTopoMSAFile'] != "" or 
            g_params['outGroupedResultFile'] != ""):
        groupList = GetGroupedTopoMSA(numSeq, idList, topoSeqList, dgScoreList,
                cmpToConsRecordList, indexIDTTopo, consensusTopo, Mcmp,
                indexClass, comparisonClassNameListAll[method_comparison])
#         print groupList

        if g_params['outGroupedSortedOrigTopoMSAFile'] != "": # to be continues
            WriteGroupedSortedOrigTopoMSA(
                    g_params['outGroupedSortedOrigTopoMSAFile'], groupList,
                    idList, topoSeqList, numIDTTopo, consensusTopo,
                    indexIDTTopo,
                    comparisonClassNameListAll[method_comparison])
        if g_params['outGroupedResultFile'] != "":
            WriteGroupedTopoAnaResult(numSeq, cmpToConsRecordList, groupList,
                    numTM_IDT, indexClass, indexIDTTopo, dgScoreList,
                    comparisonClassNameListAll[method_comparison],
                    fpoutGrouped)
    if (g_params['outInvertedTopologyMSAFile'] != ""):
        WriteInvertedTopologyPairwise_MSA(g_params['outInvertedTopologyMSAFile'],
                idList, topoSeqList, posTMList, NtermStateList)

    if (g_params['outSortedTrimmedTopoMSAFile'] != "" or
            g_params['outSortedTrimmedTopoMSAFileHTML'] != ""):
        trimmedTopoSeqList=[ct.filterTopo(ct.trimTopo(topo)) for topo in
                topoSeqList]
        trimmedConsensusTopo = ct.filterTopo(ct.trimTopo(consensusTopo))
        if g_params['outSortedTrimmedTopoMSAFile'] != "":
            WriteSortedTrimmedTopoMSA(
                    g_params['outSortedTrimmedTopoMSAFile'],
                    idList, trimmedTopoSeqList, posTMList, numIDTTopo,
                    trimmedConsensusTopo, numTM_IDT, indexIDTTopo, indexClass,
                    comparisonClassNameListAll[method_comparison], infoDIFF)
        if g_params['outSortedTrimmedTopoMSAFileHTML'] != "":
            WriteSortedTrimmedTopoMSAHTML(
                    g_params['outSortedTrimmedTopoMSAFileHTML'],
                    lengthAlignment, idList, trimmedTopoSeqList, indexIDTTopo)
    return 0
#}}}

def main(g_params):#{{{
    argv = sys.argv
    numArgv=len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    tableinfoFile = ""
    pdbtospFile = ""
    sprotACListFile = ""

    # argument parsing#{{{
    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            isNonOptionArg=False
            i = i + 1
        elif argv[i] == "--":
            isNonOptionArg=True
            i = i + 1
        elif argv[i][0] == "-":
            if argv[i] ==  "-h" or  argv[i] == "--help":
                PrintHelp()
                return 1
            elif argv[i] in ["-mode", "--mode"]:
                g_params['mode_comparison'], i = myfunc.my_getopt_int(argv, i)
            elif argv[i] in [ "-i", "--infile"]:
                g_params['inFile'], i = myfunc.my_getopt_str(argv, i)
            elif (argv[i] in ["-o", "--o", "-outfile", "--outfile"]):
                g_params['outFile'], i = myfunc.my_getopt_str(argv, i)
            elif (argv[i] in ["-og", "--og"]):
                g_params['outGroupedResultFile'], i = myfunc.my_getopt_str(argv, i)
            elif (argv[i] in ["-outpath", "--outpath"]):
                g_params['outpath'], i = myfunc.my_getopt_str(argv, i)
            elif (argv[i] in ["-pdbtosp", "--pdbtosp"]):
                pdbtospFile, i = myfunc.my_getopt_str(argv, i)
            elif (argv[i] in ["-sprot", "--sprot"]):
                sprotACListFile, i = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-mcmp", "--mcmp"]:
                g_params['method_comparison'], i = myfunc.my_getopt_int(argv, i)
            elif argv[i] in ["-cmpsp", "--cmpsp"]:
                g_params['isCompareSP'] = True; i += 1
            elif argv[i] in ["-cmpdup", "--cmpdup"]:
                g_params['isCompareDup'] = True; i += 1
            elif argv[i] in ["-maxdgdiff", "--maxdgdiff"]:
                g_params['maxDGdifference'], i = myfunc.my_getopt_float(argv, i)
            elif argv[i] in ["-mcons", "--method-cons"]:
                g_params['method_consensus'], i = myfunc.my_getopt_int(argv, i)  
            elif argv[i] in ["-c-tm2tm", "--c-tm2tm"] :
                g_params['threshold_TM2TM'],  i = myfunc.my_getopt_float(argv, i)
            elif argv[i] in ["-c-tm2gap", "--c-tm2gap"] :
                g_params['threshold_TM2GAP'], i = myfunc.my_getopt_float(argv, i)
            elif argv[i] in [ "-midt", "--midt"] :
                g_params['method_getIDTgroup'], i = myfunc.my_getopt_int(argv, i)
            elif argv[i] in [ "-verbose", "--verbose" ]:
                g_params['verbose'], i = myfunc.my_getopt_int(argv, i)  
            elif argv[i] in [ "-log",  "--log"] :
                g_params['logFile'],  i = myfunc.my_getopt_str(argv, i)
            elif argv[i] in [ "-origmsa",  "--origmsa"] :
                g_params['origTopoMSAFile'],  i = myfunc.my_getopt_str(argv, i) 
            elif (argv[i] in ["-wo", "--wo", "-write-sorted-origmsa" ,
                "--write-sorted-origmsa"]) :
                g_params['outSortedOrigTopoMSAFile'], i = myfunc.my_getopt_str(argv, i) 
            elif (argv[i] in ["-obad", "-obad"]) :
                g_params['outBadMappedPairFile'], i = myfunc.my_getopt_str(argv, i) 
            elif (argv[i] in ["-woc","--woc"]):
                g_params['outSortedClusteredTopoMSAFile'], i=myfunc.my_getopt_str(argv, i)
            elif (argv[i] in ["-woinv","--woinv"]):
                (g_params['outInvertedTopologyMSAFile'], i) = myfunc.my_getopt_str(argv, i)
            elif (argv[i] in ["-wog", "--wog"] ) :
                g_params['outGroupedSortedOrigTopoMSAFile'], i = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-tableinfo", "--tableinfo"]:
                tableinfoFile, i = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-seqidttype", "--seqidttype"]:
                g_params['seqidttype'], i = myfunc.my_getopt_int(argv,i)
            elif (argv[i] in [ "-wt" , "-write-sorted-trimmedmsa" ,
                "--write-sorted-trimmedmsa"]) :
                g_params['outSortedTrimmedTopoMSAFile'],  i = myfunc.my_getopt_str(argv, i)
            elif (argv[i] in ["-dupfile", "-dupfile"]):
                g_params['dupfile'], i = myfunc.my_getopt_str(argv, i)
            elif (argv[i] in ["-mp", "--mp"]):
                g_params['pairwise_comparison_method'], i = myfunc.my_getopt_int(argv, i)
            elif (argv[i] in ["-mm", "--mm"]):
                g_params['method_topology_comparison'], i = myfunc.my_getopt_int(argv, i)
            elif (argv[i] in ["-signalp", "-signalp"]):
                (g_params['signalpfile'], i) = myfunc.my_getopt_str(argv, i)
            elif (argv[i] in ["-localali", "-localali"]):
                (g_params['localalnfile'], i) = myfunc.my_getopt_str(argv, i)
            elif (argv[i] in ["-rmsp", "-rmsp"]):
                g_params['isRemoveSignalP'] = True; i += 1
            elif argv[i] == "-wohtml":
                g_params['outSortedOrigTopoMSAFileHTML'], i = myfunc.my_getopt_str(argv, i)
            elif argv[i] == "-wthtml":
                g_params['outSortedTrimmedTopoMSAFileHTML'], i = myfunc.my_getopt_str(argv, i)
            elif (argv[i] in ["-pi", "-print-index-idt", "--print-index-idt"]):
                g_params['isPrintIndexIDTGroup']=True; i += 1
            elif (argv[i] in ['-pdgcons', '--pdgcons' ]):
                g_params['isPrintDGCons']=True; i += 1
            elif (argv[i] in [ "-mino" ,  "--mino", "-minoverlap", "--minoverlap"]):
                g_params['min_TM_overlap'], i = myfunc.my_getopt_int(argv, i)
            elif (argv[i] in ['-debug-grouping', '--debug-grouping' ]):
                g_params['DEBUG_GROUPING'] = True; i += 1
            elif (argv[i] in ['-debug-tmmapping', '--debug-tmmapping' ]):
                g_params['DEBUG_TMMAPPING'] = True; i += 1
            elif (argv[i] in ['-debug-consensus', '--debug-consensus' ]):
                g_params['DEBUG_CONSENSUS'] = True; i += 1
            elif (argv[i] in ['-debug', '--debug', "-debugall", "--debugall" ]):
                g_params['DEBUG_CONSENSUS'] = True
                g_params['DEBUG_TMMAPPING'] = True
                g_params['DEBUG_GROUPING'] = True
                i += 1
            else:
                print >> sys.stderr,"Error! Wrong argument:%s" % (argv[i])
                return 1
        else:
            g_params['inFile'] = argv[i]
            i += 1
#}}}
    g_params['fpout'] = myfunc.myopen(filename= g_params['outFile'],
            default_fp = sys.stdout, mode="w", isRaise=False); 
    g_params['fpout_badmap'] = myfunc.myopen(filename= g_params['outBadMappedPairFile'],
            default_fp = sys.stderr, mode="w", isRaise=False); 
    g_params['fpLog'] = myfunc.myopen(filename= g_params['logFile'],
            default_fp = None, mode="w", isRaise=False); 
    g_params['fpoutGrouped'] = myfunc.myopen(filename=
            g_params['outGroupedResultFile'], default_fp = None, mode="w",
            isRaise=False); 


    if g_params['inFile'] == "" or not os.path.exists(g_params['inFile']):
        print >> sys.stderr, ("input file not set or does not exist. Exit %s."
                % argv[0])
        return -1

    if g_params['isCompareSP'] == True and g_params['signalpfile'] == "":
        print >> sys.stderr, "signalpfile not set while -cmpsp is enabled. Exit."
        return 1
    if g_params['isCompareDup'] == True and g_params['dupfile'] == "":
        print >> sys.stderr, "dupfile not set while -cmpdup is enabled. Exit."
        return 1

    rootname = os.path.basename(os.path.splitext(g_params['inFile'])[0])
    rootname = rootname.split('.')[0]
    g_params['rootname'] = rootname

# Read in tableinfo file for pairwise sequence alignment info
    if tableinfoFile != "":
        g_params['pairalnStat'] = lcmp.ReadPairAlnTableInfo(tableinfoFile)

# Read in duplication file as idpair set
    if g_params['dupfile'] != "" :
        dupPairList = lcmp.ReadDupPairList(g_params['dupfile'])
        g_params['dupPairSet'] = set(dupPairList)

# Read in pdbtosp map
    if pdbtospFile != "":
        (pdb2uniprotMap, uniprot2pdbMap) =\
                myfunc.ReadPDBTOSP(pdbtospFile)
        g_params['uniprot2pdbMap'] = uniprot2pdbMap
    #debug
#     for pdbid in pdb2uniprotMap:
#         print pdbid, pdb2uniprotMap[pdbid]
#     for uniprotac in uniprot2pdbMap:
#         print uniprotac, uniprot2pdbMap[uniprotac]
#     return 0

# Read in swissprot ac list 
    if sprotACListFile != "":
        g_params['swissprotAcSet'] = set(myfunc.ReadIDList(sprotACListFile))

# Read in signalpeptide definition file as dictionary
    if g_params['signalpfile'] != "":
        g_params['signalpDict'] = lcmp.ReadSignalPDict(g_params['signalpfile'])

    if (g_params['pairwise_comparison_method'] == 3 
            and g_params['signalpDict'] == {}):
        print >> sys.stderr, "signalp is not set while mp is 3. Exit."
        return 1


    if g_params['isRemoveSignalP'] and g_params['signalpDict'] == {}:
        msg = "RMSP is enabled, but the signalp definition is empty. Exit."
        print >> sys.stderr, msg
        return 1

    if g_params['outpath'] != "":
        if not os.path.exists(g_params['outpath']):
            cmd = ["mkdir", "-p", g_params['outpath']]
            try:
                subprocess.check_call(cmd)
            except subprocess.CalledProcessError, e:
                print e
                g_params['outpath'] = ""

# read in the topology msa file (can be with dg scores)
    if g_params['mode_comparison'] == 0: # pairwise comparison
        # For pairwise comparison, read in a number of topoRecords at a time
        # to prevent the memory overflow problem. Changed 2011-11-03 
        if g_params['localalnfile'] != "":
            (idListLocal, annoListLocal, 
                    seqListLocal) = myfunc.ReadFasta(g_params['localalnfile'])
            nump = len(idListLocal)/2
            if nump > 0:
                for ii in xrange(nump):
                    id1 = idListLocal[2*ii]
                    id2 = idListLocal[2*ii+1]
                    unaligned_str = GetUnAlignedString(seqListLocal[2*ii],
                            seqListLocal[2*ii+1])
                    if unaligned_str != "":
                        g_params['localseqpairDict'][(id1, id2)] =\
                                [seqListLocal[2*ii], seqListLocal[2*ii+1], unaligned_str]

        try:
            fpin = open (g_params['inFile'], "rb")
        except IOError:
            msg =  "Failed to open input file %s. Exit."
            print >> sys.stderr,  msg%(g_params['inFile'])
            return -1
        unprocessedBuffer=""
        topoRecordList=[]
        isEOFreached = False
        while 1:
            buff = fpin.read(BLOCK_SIZE)
            if len(buff) < BLOCK_SIZE:
                isEOFreached=True
            buff = unprocessedBuffer + buff
            unprocessedBuffer = ReadTopoWithDGScoreFromBuffer(buff,
                    topoRecordList, isEOFreached)
            if len(topoRecordList) >= 2: 
                if g_params['isRemoveSignalP']:
                    topoRecordList = RemoveSignalPeptide(topoRecordList, g_params['signalpDict'])
                PairwiseTopologyComparison(topoRecordList, g_params)
            if isEOFreached == True:
                break
        fpin.close()
    elif g_params['mode_comparison'] == 1: # multiple sequence alignment
        # for multiple aligment, Read in topoRecordList at once, but in that
        # case, the size of the file is limited to 400 Mb. There should be a
        # way to deal with even larger files that are not limited by the
        # memory size. To do so, I can first index the fasta file and then
        # read in sequence data from the hard disk at needs.
        topoFileSize = os.path.getsize(g_params['inFile'])
        if topoFileSize > g_params['MAX_ALIGN_FILE_SIZE']:
            msg =  "Size of MSA is over the limit (%d). Exit."
            print >> sys.stderr, msg%(g_params['MAX_ALIGN_FILE_SIZE'])
            return -1

        topoRecordList = ReadTopoWithDGScore(g_params['inFile'])
#         #debug
#         print "Before removal of signal peptide"
#         for rd in topoRecordList:
#             print "%s %d\t%s"%(rd[0],len(rd[2]), rd[2]), rd[4]
#         #debug

        if g_params['isRemoveSignalP']:
            topoRecordList = RemoveSignalPeptide(topoRecordList, g_params['signalpDict'])

        numSeq = len(topoRecordList)
        if numSeq < 2:
            msg = "Too few sequences (%d). No comparison can be made. Exit."
            print >> sys.stderr, msg%(numSeq)
            return -1
        return MultipleTopologyComparison(topoRecordList, g_params) #debug
    else:
        msg = "mode %d has not been implemented yet. Exit %s."
        print >> sys.stderr, msg%(g_params['mode_comparison'], argv[0])
        return -1

    myfunc.myclose(g_params['fpout'])
    myfunc.myclose(g_params['fpout_badmap'])
    myfunc.myclose(g_params['fpLog'])
    myfunc.myclose(g_params['fpoutGrouped'])
#}}}


def InitGlobalParameter():#{{{
    g_params = {}
    g_params['rootname'] = ""
    g_params['outpath'] = ""
    g_params['isPrintIndexIDTGroup'] = False
    g_params['isPrintDGCons'] = False
    g_params['verbose'] =  1
    g_params['cntOutputPair'] = 0
    # method_comparison == 2 is corresponding to pairwise_comparison_method == 1
    g_params['method_comparison']  = 1 # just for multiple comparison method
    g_params['method_consensus']   = 2 # method to compare topology to the consensus
    g_params['method_getIDTgroup'] = 1 # method to get the largest identical group
    g_params['isCompareSP'] = False    # whether do signal peptide comparison
    g_params['mode_comparison'] = 1; # 0: pairwise, 1: multiple    
    g_params['pairwise_comparison_method'] = 0; #0: 
    g_params['MAX_ALIGN_FILE_SIZE']=400 * 1024 * 1024; # 400Mb
    g_params['outSortedOrigTopoMSAFile'] = ""
    g_params['outSortedClusteredTopoMSAFile'] = ""
    g_params['outSortedTrimmedTopoMSAFile'] = ""
    g_params['outSortedOrigTopoMSAFileHTML'] = ""
    g_params['outSortedTrimmedTopoMSAFileHTML']= ""
    g_params['outGroupedResultFile']= ""
    g_params['outGroupedSortedOrigTopoMSAFile']= ""
    g_params['outInvertedTopologyMSAFile']= ""
    g_params['fpout']= sys.stdout
    g_params['fpout_badmap'] = sys.stderr
    g_params['fpLog']= None
    g_params['fpoutGrouped']= None
    g_params['inFile'] = ""
    g_params['outFile'] = ""
    g_params['outBadMappedPairFile'] = ""
    g_params['logFile'] = ""
    g_params['origTopoMSAFile'] = ""
    g_params['min_TM_overlap'] = 5
    # maximum allowed DG difference for similar DGs
    g_params['maxDGdifference'] = 0.5; 

    g_params['dupfile'] = ""
    g_params['isCompareDup'] = False
    g_params['dupPairSet'] = set([])
    g_params['signalpfile'] = ""
    g_params['signalpDict'] = {}
    g_params['isRemoveSignalP'] = False
    g_params['localalnfile'] = ""
    g_params['localseqpairDict'] = {}

    g_params['seqidttype'] = 1
    g_params['pairalnStat'] = {} #info for pairwise alignment, e.g. seqidt
    g_params['uniprot2pdbMap'] = {} # dictionary from uniprotAC->pdbid
    g_params['swissprotAcSet'] = set([]) # set of accession numbers for swissprot

    g_params['threshold_TM2TM'] = 1/3.0
    g_params['threshold_TM2GAP'] = 1/2.0
    g_params['method_topology_comparison'] = 0 #0 or 9

    g_params['DEBUG_CONSENSUS'] = False
    g_params['DEBUG_GROUPING'] = False
    g_params['DEBUG_TMMAPPING'] = False
    #pTM = re.compile("([M-]*M[M-]*)")
    g_params['GAP'] = GAP
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))
#    cProfile.run("main()")

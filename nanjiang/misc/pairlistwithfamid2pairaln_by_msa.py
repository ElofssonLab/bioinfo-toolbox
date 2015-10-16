#!/usr/bin/env python
# ChangeLog 2013-05-15 
#   option -localali added
#   option -msaext added
# Description:
import os
import sys
import myfunc
import libtopologycmp as lcmp
import random
from math import sqrt
usage = """
usage:    pairlistwithfamid2pairaln_by_msa.py pairlistwithfamid-file
Description: Generate pairwise alignment file and tableinfo from multiple alignment

Options:
  -outname STR    rootname of the output file
  -all            When all is enabled, every pair will be output
  -msapath DIR    path for multiple alignment
  -msaext  STR    extension for msa file, (default: .msa.fa)
  -localali       MSA is local alignment created by dialign
  -maxsel  INT    maximum number of pairs selected for each sequence identity
                  group
  -seqidttype INT Set sequence identity type, (default: 1)
                  0: seqIDT = numIDTRes /alnLength
                  1: seqIDT = numIDTRes / min(len1, len2)
                  2: seqIDT = numIDTRes / (alnLength - NumGAP)
  -q              Quiet mode
  -verbose INT    verbose level, 0 or 1
  -h, --help      Print this help message and exit

Created 2012-08-20, updated 2013-11-28, Nanjiang Shu 

Examples:
     pairlistwithfamid2pairaln_by_msa.py -msapath msa -maxsel 30000 -seqidttype 1 pairlistwithclanid.txt
"""
SEQIDT_GROUP = [
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
        ];

def PrintHelp():
    print usage

def GetSeqIDTGroupIndex(seqidt, seqIDTGroupList):#{{{
    numGroup = len(seqIDTGroupList)/2;
    for i in xrange(numGroup):
        if seqidt >= seqIDTGroupList[i*2] and seqidt < seqIDTGroupList[i*2+1]:
            return i;
    return numGroup;
#}}}

def AddPairwiseAlignmentFactor(pairlistDict, msapath, msaext, #{{{
        isLocalAlignment):
    cntfamid = 0
    verbose = g_params['verbose']
    for famid in pairlistDict:
        cntfamid += 1
        if verbose >= 2:
            print "Add pairwise alignment factor for %d: %s"%(cntfamid, famid)
        msafile = msapath + os.sep + famid + msaext
        if not os.path.exists(msafile):
            print >> sys.stderr, "msafile %s does not exist. Ignore" % msafile
            continue
        (idList, annoList, seqList) = myfunc.ReadFasta(msafile)
        msaDict = {}
        for i in xrange(len(idList)):
            msaDict[idList[i]] = seqList[i]
        pairlist = pairlistDict[famid]
        #print "pairlist=", pairlist
        for i in xrange(len(pairlist)):
            pair = pairlist[i]
            #print "pair = ", pair
            seq1 = ""
            seq2 = ""
            id1 = pair[0]
            id2 = pair[1]
            if id1 in msaDict and id2 in msaDict:
                seq1 = msaDict[id1] 
                seq2 = msaDict[id2]
                [seq1, seq2] = lcmp.RemoveUnnecessaryGap([seq1, seq2])
                if len(seq1) != len(seq2):
                    print >> sys.stderr, "Bad alignment for %s and %s" %(id1,id2)
                else:
                    alignFactor = lcmp.GetAlignmentFactorFromPairAlignment(
                            seq1,seq2, isLocalAlignment)
                    pair.append(alignFactor)
            else:
                if id1 not in msaDict:
                    print >> sys.stderr, "%s not in msafile %s"%(id1, msafile)
                if id2 not in msaDict:
                    print >> sys.stderr, "%s not in msafile %s"%(id2, msafile)
    return 0
#}}}
def SelectPairListOfEachSeqIDTBin(pairlistDict):#{{{
# input is pairlist of one sequence identity group
    # families Na, Nb, Nc pairs
    # if Na + Nb + Nc > MAX
# introduce factors f1, f2, f3
# let 
#   f1 * Na + f2 * Nb + f3 * Nc = MAX
# where f1 = f / sqrt(Na)
#       f2 = f / sqrt(Nb)
#       f3 = f/ sqrt (Nc)
# therefore:
#       f * (sqrt(Na) + sqrt(Nb) + sqrt(Nc) = MAX
# ==>   f = MAX / (sqrt(Na) + sqrt(Nb) + sqrt(Nc)
# 
    maxsel = g_params['maxsel']
# first get numpairs in each family
    numpairList = []
    for famid in pairlistDict:
        numpairList.append(len(pairlistDict[famid]))
    if sum(numpairList) <= maxsel:
        return pairlistDict
    else:
        sum1 = 0.0
        for x in numpairList:
            sum1 += sqrt(float(x))
        f = float(maxsel) / sum1
        newPairlistDict = {}
        for famid in pairlistDict:
            pairlist = pairlistDict[famid]
            numPair = len(pairlist)
            newNumPair = int(f * sqrt(float(numPair)) + 0.5)
            if newNumPair >= 1:
                if newNumPair > numPair:
                    newNumPair = numPair
                newPairlistDict[famid] = random.sample(pairlist, newNumPair)
        return newPairlistDict
#}}}

def RandSelectPairList(pairlistDict):#{{{
    seqidttype = g_params['seqidttype']
    seqIDTGroupAll = SEQIDT_GROUP
    numSeqIDTGroupAll = len(seqIDTGroupAll)/2;
    pairlistDict_list = []
    for i in xrange(numSeqIDTGroupAll):
        pairlistDict_list.append({})
    for famid in pairlistDict:
        pairlist = pairlistDict[famid]
        for pair in pairlist:
            if len(pair) >= 3:
                id1 = pair[0]
                id2 = pair[1]
                seqidt = -1.0
                record = pair[2]
                seqidt = lcmp.GetSeqIDT(record, seqidttype)
                idxGroup = GetSeqIDTGroupIndex(seqidt, seqIDTGroupAll)
                if idxGroup < numSeqIDTGroupAll:
                    if not famid in pairlistDict_list[idxGroup]:
                        pairlistDict_list[idxGroup][famid] = []
                    pairlistDict_list[idxGroup][famid].append(pair)

    for i in xrange(numSeqIDTGroupAll):
        pairlistDict_list[i] = SelectPairListOfEachSeqIDTBin(pairlistDict_list[i])

    selectedPairListDict = {}
    for i in xrange(numSeqIDTGroupAll):
        for famid in pairlistDict_list[i]:
            if not famid in selectedPairListDict:
                selectedPairListDict[famid] = []
            selectedPairListDict[famid] += pairlistDict_list[i][famid]
    return selectedPairListDict
#}}}

def WritePairAln(pairlistDict, msapath, msaext, outname):#{{{
    verbose = g_params['verbose']
    outAlnFile = outname + ".pairaln"
    outTableFile = outname + ".tableinfo"
    outSelPairList = outname + ".pairlistwithpfamid"
    try:
        fpout_aln = open(outAlnFile, "w")
    except IOError:
        print >> sys.stderr, "Failed to write to file", outAlnFile
        return 1
    try:
        fpout_table = open(outTableFile, "w")
    except IOError:
        print >> sys.stderr, "Failed to write to file", outTableFile
        return 1

    try:
        fpout_list = open(outSelPairList, "w")
    except IOError:
        print >> sys.stderr, "Failed to write to file", outSelPairList
        return 1

    fpout_table.write("#%-15s %-15s %6s %6s %9s %6s %6s %9s %6s %6s %6s %6s %6s\n" % (
        "Seq1","Seq2", "IDT0", "SIM0", "AlnLength", "Len1","Len2",
        "Score","N_IDT", "N_SIM", "N_GAP", "IDT1", "IDT2"))

    for famid in pairlistDict:
        if verbose >= 2:
            print "Write pairwise alignment for %s"%(famid)
        msafile = msapath + os.sep + famid + msaext
        if not os.path.exists(msafile):
            print >> sys.stderr, "msafile %s does not exist. Ignore" % msafile
            continue
        (idList, annoList, seqList) = myfunc.ReadFasta(msafile)
        msaDict = {}
        annoDict = {}
        for i in xrange(len(idList)):
            msaDict[idList[i]] = seqList[i]
            annoDict[idList[i]] = annoList[i]
        pairlist = pairlistDict[famid]
        #print "pairlist2=", pairlist
        for pair in pairlist:
            #print "pair2 = ", pair
            seq1 = ""
            seq2 = ""
            id1 = pair[0]
            id2 = pair[1]
            if id1 in msaDict and id2 in msaDict:
                seq1 = msaDict[id1] 
                seq2 = msaDict[id2]
                [seq1, seq2] = lcmp.RemoveUnnecessaryGap([seq1, seq2])
                if len(seq1) != len(seq2):
                    print >> sys.stderr, "Bad alignment for %s and %s" %(id1,id2)
                else:
                    rd = pair[2]
                    fpout_aln.write(">%s aligned_to=%s seqIDT=%.1f seqIDT1=%.1f\n"%(
                        annoDict[id1], id2, rd['seqidt0'], rd['seqidt1']))
                    fpout_aln.write("%s\n"%seq1)
                    fpout_aln.write(">%s aligned_to=%s seqIDT=%.1f seqIDT1=%.1f\n"%(
                        annoDict[id2], id1, rd['seqidt0'], rd['seqidt1']))
                    fpout_aln.write("%s\n"%seq2)
                    fpout_table.write("%-16s %-15s %6.1f %6.1f %9d %6d %6d %9.1f %6d %6d %6d %6.1f %6.1f\n"% (
                        id1, id2, rd['seqidt0'], -1.0,
                        rd['alnLength'],
                        rd['seqLength1'], rd['seqLength2'],
                        -1.0,
                        rd['numIDT'], -1, rd['numGap'],
                        rd['seqidt1'], rd['seqidt2']))
                    fpout_list.write("%s %s %s\n"%(id1, id2, famid))
    fpout_aln.close()
    fpout_table.close()
    fpout_list.close()
    print "Result output to "
    print "\t%s"%outAlnFile
    print "\t%s"%outTableFile

    return 0
#}}}
def ReadPairListWithFamID(infile):#{{{
    """
    Input:
        file of the content  "id1 id2 pfamid" for each line
    Output:
        pairlistDict        {pfamid:[(id1,id2),(id1,id2),...]
    """
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
def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    msapath = ""
    outname = ""
    pairlistwithfamid_file = ""
    isLocalAlignment = False

    idListFile = None
    idList = []
    msaext = ".msa.fa"

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            pairlistwithfamid_file = argv[i]
            isNonOptionArg = False
            i += 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-outname", "--outname"]:
                (outname, i)  = myfunc.my_getopt_str(argv,i)
            elif argv[i] in ["-msapath", "--msapath"]:
                (msapath,i)  = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-msaext", "--msaext"]:
                (msaext,i)  = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-localali", "--localali"]:
                isLocalAlignment = True; i += 1
            elif argv[i] in ["-all", "--all"]:
                g_params['isOutputAll'] = True; i += 1
            elif argv[i] in ["-seqidttype", "--seqidttype"]:
                (g_params['seqidttype'],i) = myfunc.my_getopt_int(argv, i)
            elif argv[i] in ["-maxsel", "--maxsel"]:
                (g_params['maxsel'],i) = myfunc.my_getopt_int(argv,i)
            elif argv[i] in ["-verbose", "--verbose"]:
                (g_params['verbose'],i) =  myfunc.my_getopt_int(argv, i)
            elif argv[i] in ["-q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            pairlistwithfamid_file = argv[i]
            i += 1

    if pairlistwithfamid_file == "":
        print >> sys.stderr, "pairlistwithfamid_file not set"
        return 1

    if outname == "":
        addname = ""
        if isLocalAlignment is True:
            addname = ".local"
        outname = os.path.splitext(pairlistwithfamid_file)[0] + addname + ".selected"

    verbose = g_params['verbose']

    if verbose >= 1:
        print "Reading file", pairlistwithfamid_file
    # pairlistDict {pfamid:[(id1,id2),(id1,id2),...]    }
    pairlistDict = ReadPairListWithFamID(pairlistwithfamid_file)
    #print pairlistDict
    # obtain pairwise alignment factor
    if verbose >= 1:
        print "Obtaining pairwise alignment factor"
    AddPairwiseAlignmentFactor(pairlistDict, msapath, msaext, isLocalAlignment)
    if not g_params['isOutputAll']:
        if verbose >=1:
            print "Select pair list to limit to ", g_params['maxsel']
        pairlistDict = RandSelectPairList(pairlistDict)
    #print "selected pairlistDict = ", pairlistDict
    #print "pairlistDict_selected=",pairlistDict
    if verbose >= 1:
        print "Write out the alignment"
    WritePairAln(pairlistDict, msapath, msaext, outname)


#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['seqidttype'] = 1
    g_params['maxsel'] = 30000
    g_params['isOutputAll'] = False # when this is True, maxsel is not functioning
    g_params['verbose'] = 1
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

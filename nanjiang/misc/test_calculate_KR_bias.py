#!/usr/bin/env python
# Description:
# ChangeLog 
# ChangeLog 2013-01-28 
#     K, R residues that are >= MAX_DIST (e.g. 100) away from the TM helix are
#     not counted in statistics
# ChangeLog 2013-07-09 
#     Re-define the inter-distance and outer distance.
#     MAX_DIST should be much shorter, for SCAMPI, it uses 12 on both sides,
#     and the soft bonder within the TM boundary can be set as 5 instead of 10
import os
import sys
import myfunc
import libtopologycmp as lcmp
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))
usage = """
usage:  %s  -seqfile FILE -topofile
Description:
    Calculate KRbias given sequence
    output file name is $outpath/$rootname{seqfile}.krbias.txt

Options:
  -outpath DIR    Set ouput path
  -l      FILE    Set the idListFile
  -q              Quiet mode
  -maxdist INT    Maximal distant from the TM helix to be counted.
                  (default: 12)
  -flankwin INT   Flanking window of TM helix, (default: 5)
  -h, --help      Print this help message and exit
  -debug          Show debug information


Created 2009-06-08, updated 2014-11-05, Nanjiang Shu
"""%(progname)

def PrintHelp():
    print usage
def GetKRPos(seq, b, e):
    pos = []
    for j in range(b,e):
        if seq[j] in ["K","R"]:
            pos.append(j)
    return pos

def CalKRBias(seq, topo, flank_win, max_dist):
    kr_bias = None
    KR_pos_list = [] # [[1,5], [1,3]]
    
    posTM = myfunc.GetTMPosition(topo)
    NtermState = lcmp.GetNtermState(topo)
    numTM = len(posTM)
    length = len(seq)
    if numTM >= 1:
        sum_KR_odd = 0
        sum_KR_even = 0
        for i in xrange(numTM+1):
            if i == 0 or i == numTM:
                if i == 0:
                    b = max(0, posTM[i][0] - max_dist)
                    e = posTM[i][0] + flank_win
                else:
                    b = posTM[i-1][1] - flank_win
                    e = min(length, posTM[i-1][1] + max_dist)
                KRpos = GetKRPos(seq, b, e)
            else:
                if posTM[i][0] - posTM[i-1][1] > 2*max_dist:
                    b1 = posTM[i-1][1] - flank_win
                    e1 = posTM[i-1][1] + max_dist
                    b2 = posTM[i][0] - max_dist
                    e2 = posTM[i][0] + flank_win
                    KRpos = GetKRPos(seq, b1, e1)
                    KRpos += GetKRPos(seq, b2,e2)
                else:
                    b = posTM[i-1][1] - flank_win
                    e = posTM[i][0] + flank_win
#                     print (b,e)
#                     print len(seq)
#                     print "flank_win=",flank_win
#                     print "i=",i
#                     print "numTM=", len(posTM)
                    KRpos = GetKRPos(seq, b, e)
            KR_pos_list.append(KRpos)
            if i%2 == 0:
                sum_KR_odd += len(KRpos)
            else:
                sum_KR_even += len(KRpos)

        kr_bias = sum_KR_odd - sum_KR_even
#        print KR_pos_list
    return (kr_bias, KR_pos_list, numTM)

def WriteResult(idd, cmpclass, seq, numTM, kr_bias, KR_pos_list, fpout):
    fpout.write("%-5s %20s (%+d): "%(idd, cmpclass, kr_bias))
    for i in xrange(numTM+1):
        if i > 0:
            fpout.write(" %s%d%s "%("MMMM", i, "MMMM"))
        for j in KR_pos_list[i]:
            fpout.write(" %s%d "%(seq[j], j+1))
    fpout.write("\n")

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = "./"
    idListFile = None
    idList = []
    seqfile = ""
    topofile = ""
    max_dist = 12 # maximal distance to the TM helix so that K, R residues are counted
    flank_win = 5 # flanking window of the TM helix, residues at position
                   #TMbeg-flank_win and TMend+flank_win are also counted

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            idList.append(argv[i])
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
                (outpath,i) = myfunc.my_getopt_str(argv,i)
            elif argv[i] in ["-maxdist", "--maxdist"]:
                (max_dist,i) = myfunc.my_getopt_int(argv,i)
            elif argv[i] in ["-flankwin", "--flankwin"]:
                (flank_win,i) = myfunc.my_getopt_int(argv,i)
            elif argv[i] in ["-seqfile", "--seqfile"]:
                (seqfile,i) = myfunc.my_getopt_str(argv,i)
            elif argv[i] in ["-topofile", "--topofile"]:
                (topofile,i) = myfunc.my_getopt_str(argv,i)
            elif argv[i] in ["-l", "--l"] :
                (idListFile,i) = myfunc.my_getopt_str(argv,i)
            elif argv[i] in ["-q"]:
                g_params['isQuiet'] = True
                i += 1
            elif argv[i] in ["-debug"]:
                g_params['isDEBUG'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            idList.append(argv[i])
            i += 1
    (idListSeq, annoListSeq, seqList) = myfunc.ReadFasta(seqfile)
    (idListTopo, annoListTopo, topoList) = myfunc.ReadFasta(topofile)

    numSeq = len(idListSeq)
    numTopo = len(idListTopo)
    if numSeq < 1 or numTopo < 1: 
        print >> sys.stderr, "No seq set"
        return 1
    seqDict = {}
    for i in xrange(numSeq):
        seqDict[idListSeq[i]] = seqList[i]
    topoDict = {}
    for i in xrange(numTopo):
        topoDict[idListTopo[i]] = topoList[i]

    cmpclassDict = {}
    for anno in annoListTopo:
        anno = anno.lstrip(">")
        strs = anno.split()
        cmpclassDict[strs[0]] = strs[1]

    outpath = os.path.dirname(seqfile)
    if outpath == "":
        outpath = "."
    rootname = os.path.basename(os.path.splitext(seqfile)[0])
    outfile_kr_list = outpath + os.sep + rootname + ".krlist.txt"
    outfile_krbias = outpath + os.sep + rootname + ".krbias.txt"

    fpout_krlist = open(outfile_kr_list, "w")
    fpout_krbias = open(outfile_krbias, "w")

    for idd in idListSeq:
        if g_params['isDEBUG']:
            print "seqid: %s"%(idd)
        try:
            topo = topoDict[idd]
        except KeyError:
            print >> sys.stderr, "no topo for %s"%idd
            continue
        try:
            seq = seqDict[idd]
        except KeyError:
            print >> sys.stderr, "no seq for %s"%idd
            continue
        try:
            cmpclass = cmpclassDict[idd]
        except KeyError:
            cmpclass = "INV"
        (kr_bias, KR_pos_list, numTM) = CalKRBias(seq, topo, flank_win,
                max_dist)
        WriteResult(idd, cmpclass, seq, numTM, kr_bias, KR_pos_list, fpout_krlist)
        if cmpclass in ["IDT", "INV"]:
            fpout_krbias.write("%d\n"%kr_bias)
    fpout_krlist.close()
    fpout_krbias.close()

    #cmd = "drawhistogram -outstyle png -xlabel KR-bias -ylabel Occurrences %s"%(outfile_krbias)
    #os.system(cmd)

#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['isDEBUG'] = False
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

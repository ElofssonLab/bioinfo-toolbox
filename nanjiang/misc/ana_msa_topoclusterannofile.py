#!/usr/bin/python
# Description:

# ChangeLog
# ChangeLog 2014-03-10
#   option "-overwrite" added
import os
import sys
import myfunc
import re
import libtopologycmp as lcmp
import subprocess
#import bioinformatics
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))
DATADIR3 = os.environ['DATADIR3']
binpath = os.path.dirname(sys.argv[0])
if binpath  == "":
    binpath = "."

usage_short="""
Usage: %s clusterAnnoFile [FILE...]
"""%(progname)

usage_ext="""
Description:
    Analyze cluster file output by compareMSATopo.py

OPTIONS:
  -o FILE          Output the result to file
  -l LISTFILE      Set the listfile
  -pfamdef FILE    Pfam definition file,
                   (default:$DATADIR3/data/pfam/pfam26.0/Pfam-A.clans.tsv)
  -thncls2 INT     Threshold for the minimal number of sequences for the
                   second group, (default: 2)
  -thfrac2 FLOAT   Threshold for the minimal fraction of the
                   second group, (default: 0.05)
  -tableinfo FILE  tableinfo file, for getting sequence identity
  -seqidttype INT  Type of sequence identity, (default: 1)
  -pdbtosp FILE    PDB to swissprot id (uniprot ac) maplist 
  -sprot   FILE    Supply swissprot aclist
  -aapath   DIR    Supply path for amino acid sequences for each family
  -topoaln  FILE   Pairwise topology alignment
  -overwrite       Force overwrite the existing file, (default: no)
  -q               Quiet mode
  -h, --help       Print this help message and exit

Created 2013-12-16, updated 2014-03-10, Nanjiang Shu 
"""
usage_exp="""
Examples:
"""

def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}
def GetNumSeqInClusterFromAnnotation(line):#{{{
    if line:
        m=re.search('numSeqInCluster *=[^, ]*',line)
        if m: 
            rlty = m.group(0).split('=')[1]
            try:
                return int(rlty)
            except (ValueError, TypeError):
                return None
        else: 
            return None
    return None
#}}}
def GetNumTMFromAnnotation(line):#{{{
    if line:
        m=re.search('nTM *=[^, ]*',line)
        if m: 
            rlty = m.group(0).split('=')[1]
            try:
                return int(rlty)
            except (ValueError, TypeError):
                return None
        else: 
            return None
    return None
#}}}
def GetClusterNoFromAnnotation(line):#{{{
    if line:
        m=re.search('ClusterNo *=[^, ]*',line)
        if m: 
            rlty = m.group(0).split('=')[1]
            try:
                return int(rlty)
            except (ValueError, TypeError):
                return None
        else: 
            return None
    return None
#}}}
def GetSeqIDT(seqid1, seqid2, infoDict, seqidttype):#{{{
    """
    Get sequence identity according to the infoDict and seqidttype
    """
    if seqidttype == 0:
        try:
            seqidt = infoDict['seqidt']
        except KeyError:
            print >> sys.stderr, "seqidt not in pairinfo for %s %s"%(seqid1, seqid2)
            seqidt = -1
    elif seqidttype == 1:
        try:
            seqidt = infoDict['seqidt1']
        except KeyError:
            print >> sys.stderr, "seqidt1 not in pairinfo for %s %s"%(seqid1, seqid2)
            seqidt = -1
    elif seqidttype == 2:
        try:
            seqidt = infoDict['seqidt2']
        except KeyError:
            print >> sys.stderr, "seqidt2 not in pairinfo for %s %s"%(seqid1, seqid2)
            seqidt = -1
    else:
        seqidt = -1

    return seqidt
#}}}

def AnaClusterFile(infile, anaList):#{{{
    """
    Analyze
    Input:
        a file with a number of clustered topologies with the protein family,
        the topology is clustered by the number of TM helices of each topology
        The input file is in FASTA format, while 
        e.g.
        >Q81PI9, nTM=8 ClusterNo=1 numSeqInCluster=15
        iiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMoooooooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMooooooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMooooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiiiiiii
    Output:
        anaList   {'pfamid': pfamid; 'numseq': numseq; 'cluster': [[numTM,
        numseq, [seqid1, seqid2...]], ...]}

    """
    try:
        fpin = open(infile, "r")
        pfamid = os.path.basename(infile).split(".")[0]
        lines = fpin.read().split("\n")
        fpin.close()
        cntseq = 0
        addedClusterSet = set([])

        ana = {}
        ana['famid'] = pfamid
        ana['cluster'] = []
        idxCls = -1

        for line in lines:
            if not line or line[0] != ">":
                continue
            cntseq += 1
            numCls = GetClusterNoFromAnnotation(line)
            if not numCls in addedClusterSet: # for a new cluster
                numSeqCls = GetNumSeqInClusterFromAnnotation(line)
                numTM = GetNumTMFromAnnotation(line)
                seqid = myfunc.GetFirstWord(line).lstrip('>').rstrip(",")
                ana['cluster'].append([numTM, numSeqCls, [seqid]])
                addedClusterSet.add(numCls)
                idxCls += 1
            else:
                seqid = myfunc.GetFirstWord(line).lstrip('>').rstrip(",")
                ana['cluster'][idxCls][2].append(seqid)

        ana['numseq'] = cntseq

        if len(ana['cluster']) > 0:
            anaList.append(ana)
        else:
            print >> sys.stderr, "No cluster in file %s"%(infile)
    except IOError:
        print >> sys.stderr, "Failed to read file %s"%(infile)
#}}}
def CheckSwissProt(seqid1, seqid2, swissprotAcSet):#{{{
    """
    Check whether it has swissprot
    -1     Empty swissprotAcSet
    0      NO_SPROT
    1      just first seq is sprot
    2      just second seq is sprot
    3      both sprot
    """
    if swissprotAcSet != set([]):
        isSwissprot1 = False
        isSwissprot2 = False
        if seqid1 in swissprotAcSet:
            isSwissprot1 = True
        if seqid2 in swissprotAcSet:
            isSwissprot2 = True
        if isSwissprot1 and isSwissprot2:
            return 3
        elif isSwissprot1:
            return 1
        elif isSwissprot2:
            return 2
        else:
            return 0
    return -1

#}}}
def CheckPDB(seqid1, seqid2, uniprot2pdbMap):#{{{
    """
    Check whether it has swissprot
    -1     Empty swissprotAcSet
    0      NO_SPROT
    1      just first seq is sprot
    2      just second seq is sprot
    3      both sprot
    """
    if uniprot2pdbMap != {}:
        try:
            dt1 = uniprot2pdbMap[seqid1]
            numPDBID1 = len(dt1)
        except KeyError:
            numPDBID1 = 0

        try:
            dt2 = uniprot2pdbMap[seqid2]
            numPDBID2 = len(dt2)
        except KeyError:
            numPDBID2 = 0
        if numPDBID1 ==  0 and numPDBID2 == 0:
            return 0
        else:
            if numPDBID1 > 0 and numPDBID2 > 0:
                return 3
            elif numPDBID1 > 0:
                return 1
            elif numPDBID2 > 0:
                return 2
    else:
        return -1

#}}}
def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile = ""
    fileListFile = ""
    fileList = []
    pfamDefFile = "%s/data/pfam/pfam26.0/Pfam-A.clans.tsv"%(DATADIR3)
    threshold_Fraction_Group_2 = 0.05
    threshold_NumSeq_Group_2 = 2
    tableinfoFile = ""
    pdbtospFile = ""
    sprotACListFile = ""

    threshold_g12_seqidt = 20.0

    topoalnFile = ""
    aapath = ""


    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            fileList.append(argv[i])
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
            elif argv[i] in ["-l", "--l"] :
                (fileListFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-seqidttype", "--seqidttype"] :
                (g_params['seqidttype'], i) = myfunc.my_getopt_int(argv, i)
            elif argv[i] in ["-tableinfo", "--tableinfo"] :
                (tableinfoFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-topoaln", "--topoaln"] :
                (topoalnFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-aapath", "--aapath"] :
                (aapath, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-thncls2", "--thncls2"] :
                (threshold_NumSeq_Group_2, i) = myfunc.my_getopt_int(argv, i)
            elif argv[i] in ["-thfrac2", "--thfrac2"] :
                (threshold_Fraction_Group_2, i) = myfunc.my_getopt_float(argv, i)
            elif argv[i] in ["-pfamdef", "--pfamdef"] :
                (pfamDefFile, i) = myfunc.my_getopt_str(argv, i)
            elif (argv[i] in ["-pdbtosp", "--pdbtosp"]):
                pdbtospFile, i = myfunc.my_getopt_str(argv, i)
            elif (argv[i] in ["-sprot", "--sprot"]):
                sprotACListFile, i = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True; i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            fileList.append(argv[i])
            i += 1

    if fileListFile != "":
        fileList += myfunc.ReadIDList(fileListFile)
    if len(fileList) < 1:
        print >> sys.stderr, "No input set. exit"
        return 1

    if myfunc.checkfile(topoalnFile, "topoalnFile") != 0:
        return 1
    if myfunc.checkfile(aapath, "aapath") != 0:
        return 1
    if outfile == "":
        print >> sys.stderr, "outfile not set. Exit"
        return 1

    outpath = myfunc.my_dirname(outfile)
    if not os.path.exists(outpath):
        cmd = ["mkdir", "-p", outpath]
        try:
            subprocess.check_output(cmd)
        except subprocess.CalledProcessError, e:
            print e
            return 1

    (pfamidDefDict, clanidDefDict) = lcmp.ReadPfamDefFile(pfamDefFile)
# Read in pdbtosp map
    if pdbtospFile != "":
        (pdb2uniprotMap, uniprot2pdbMap) = myfunc.ReadPDBTOSP(pdbtospFile)
    else:
        (pdb2uniprotMap, uniprot2pdbMap) = ({}, {})
# Read in swissprot ac list 
    if sprotACListFile != "":
        swissprotAcSet = set(myfunc.ReadIDList(sprotACListFile))
    else:
        swissprotAcSet = set([])

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)

    if tableinfoFile != "":
        pairalnStat = lcmp.ReadPairAlnTableInfo(tableinfoFile)
    if pairalnStat != {}:
        outfile_pair_g12 = outfile + ".pair_group1_2.txt"
        fpout_pair = myfunc.myopen(outfile_pair_g12, sys.stdout, "w", False)
    else:
        fpout_pair = None

    selectedPairList = [] # select pairs for draw pairwise topology alignment
    # for each family, select the pair with highest sequence identity between
    # group 1 and group2 
    # and select the highest sequence identity pair from swissprot

    anaList = []
    for infile in fileList:
        AnaClusterFile(infile, anaList)

    anaList = sorted(anaList, key=lambda x:x['numseq'], reverse=True)
    # the structure of the new list is the same as the old list
    print >> fpout, "#%-6s %15s %10s %10s %6s %6s %10s"%("Pfamid", "PfamDef", "NumSeqCls2",  "NumSeqCls1",
            "NumSeq", "NumCls", "ClusterList")
    header_pair_file = """\
# isSP   Status of whether the pair of proteins are swissprot entries
#       -1. the swissprot definition file not supplied
#        0. none of the two proteins are from swissprot
#        1. the first protein is from swissprot
#        2. the second protein is from swissprot
#        3. both of the proteins are from swissprot
#
# isPDB  Status of whether the pair of proteins are PDB entries
#       -1. the PDB definition file not supplied
#        0. none of the two proteins are from PDB
#        1. the first protein is from PDB
#        2. the second protein is from PDB
#        3. both of the proteins are from PDB\
"""
    try:
        print >> fpout_pair, header_pair_file
    except IOError:
        pass
    cntFam = 0
    for ana in anaList:
        cntFam += 1
        numSeqCls1 = ana['cluster'][0][1]
        if len(ana['cluster']) >= 2:
            numSeqCls2 = ana['cluster'][1][1]
        else:
            numSeqCls2 = 0
        try:
            pfamdef = pfamidDefDict[ana['famid']]
        except KeyError:
            pfamdef = ""

        # output pair list
        pairlist = []
        if numSeqCls2 >= threshold_NumSeq_Group_2:
            for seqid1 in ana['cluster'][0][2]:
                for seqid2 in ana['cluster'][1][2]:
                    pairid1 = "%s-%s"%(seqid1, seqid2)
                    pairid2 = "%s-%s"%(seqid2, seqid1)
                    if pairid1 in pairalnStat:
                        dt = pairalnStat[pairid1]
                    elif pairid2 in pairalnStat:
                        dt = pairalnStat[pairid2]
                    else:
                        print >> sys.stderr, "No seqidt info for pair %s %s"%(seqid1, seqid2)
                        dt = {}
                    seqidt = GetSeqIDT(seqid1, seqid2, dt, g_params['seqidttype'])
                    if seqidt >= threshold_g12_seqidt:
                        pairlist.append([seqidt, seqid1, seqid2])
            pairlist = sorted(pairlist, key=lambda x:x[0], reverse=True)

        if len(pairlist) > 0:
            try:
                nTM1 =  ana['cluster'][0][0]
                nTM2 =  ana['cluster'][1][0]
                print >> fpout_pair, "#%-12s %8s %15s %6s %6s %6s %6s"%(
                        "FamilyNo", "PfamID", "PfamDef", "nSeq1", "nSeq2", "nTM1", "nTM2")
                print >> fpout_pair, "//Family %4d %8s %15s %6d %6d %6d %6d"%(
                        cntFam, ana['famid'], pfamdef, numSeqCls1, numSeqCls2,
                        nTM1, nTM2)
                print >> fpout_pair, "#%-6s %6s %6s %6s %6s"%(
                        "ID1", "ID2", "SeqIDT", "isSP", "isPDB")
                for j in xrange(len(pairlist)):
                    li = pairlist[j]
                    status_swissprot = CheckSwissProt(li[1], li[2], swissprotAcSet)
                    status_pdb = CheckPDB(li[1], li[2], uniprot2pdbMap)
                    print >> fpout_pair, "%-6s %6s %6.1f %6d %6d"%(
                            li[1],li[2],li[0],status_swissprot, status_pdb)
                    if j == 0 or status_swissprot == 3:
    # seqid1, seqid2, seqidt, famid, pfamdef, numSeqCls1, numSeqCls2, nTM1, nTM2, isSP, isPDB
                        selectedPairList.append([li[1], li[2], li[0],
                            ana['famid'], pfamdef, numSeqCls1, numSeqCls2, 
                            ana['numseq'], nTM1, nTM2, status_swissprot, status_pdb]) 
            except IOError:
                pass
            closest_pair = pairlist[0]
        else:
            closest_pair = [0, "-", "-"]

        clst = []
        for li in ana['cluster']:
            clst.append((li[0],li[1]))
        print  >> fpout, "%7s %15s %10d %10d %6d %6d %6s %6s %6.1f" %(ana['famid'], pfamdef,
                numSeqCls2, numSeqCls1, ana['numseq'], len(ana['cluster']),
                closest_pair[1], closest_pair[2], closest_pair[0],
                ), clst

    myfunc.myclose(fpout)
    myfunc.myclose(fpout_pair)

    cmpclassList = ["IDT"]

# get topN diff
    freqTopNList = [] # [ [topN, min, frac_DIFF]]
    sumTopNWithDiffTopoList = [0]*len(cmpclassList)
    for i in xrange(len(anaList)):
        ana = anaList[i]
        isHaveDiffTopo = False
        numSeqCls1 = ana['cluster'][0][1]
        if len(ana['cluster']) >= 2:
            numSeqCls2 = ana['cluster'][1][1]
        else:
            numSeqCls2 = 0

        fracCls2 = myfunc.FloatDivision(numSeqCls2, ana['numseq'])
        if (numSeqCls2 >= threshold_NumSeq_Group_2 and
                fracCls2 >= threshold_Fraction_Group_2):
            isHaveDiffTopo = True
        if isHaveDiffTopo:
            sumTopNWithDiffTopoList[0] += 1
        fracList = []
        for j in xrange(len(sumTopNWithDiffTopoList)):
            fracList.append(myfunc.FloatDivision(sumTopNWithDiffTopoList[j], i+1))

        freqTopNList.append( [i+1, ana['numseq']]+ fracList)


# output topN statistics
    outfileList = []
    item = "numseq_TMpro"
    outfile1 = outfile + ".topNdifffam.thN%d.thF%g.sortby_%s.txt"%(
            threshold_NumSeq_Group_2, threshold_Fraction_Group_2, item)
    outfileList.append(outfile1)
    fpout = myfunc.myopen(outfile1, sys.stdout, "w", False)

    ss_sort_item = "min_%s"%(item)
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


# output selected pairlist and draw pairwise topology alignment for selected
# pairs
    outpath = myfunc.my_dirname(outfile)
    ext_topoaln = ".topoaln.fa"
    outfile_selected_pair = outfile+ ".selected.pairlistwithfamid"
    fpout = myfunc.myopen(outfile_selected_pair, None, "w", True)
    print >> fpout, "#%-6s %6s %6s %7s %15s %10s %10s %6s %4s %4s %5s %5s"%(
            "seqid1", "seqid2", "seqidt", "famid", "pfamdef", "numSeqCls1",
            "numSeqCls2", "numSeq", "nTM1", "nTM2", "isSP", "isPDB")
    for li in selectedPairList:
        print >> fpout, "%-7s %6s %6.1f %7s %15s %10d %10d %6d %4d %4d %5d %5d"%(
                li[0], li[1], li[2], li[3], li[4], li[5], li[6], li[7], li[8],
                li[9], li[10], li[11])
    myfunc.myclose(fpout)

    cmd = ["%s/selectPairaln.py"%(binpath), "-pairaln", topoalnFile, "-l",
            outfile_selected_pair, "-split", "-outpath", outpath, "-ext",
            ext_topoaln]
    print '\n', " ".join(cmd), '\n'
    try:
        subprocess.check_output(cmd)
    except subprocess.CalledProcessError, e:
        print e
        return 1
# draw pairwise topology alignment
    for li in selectedPairList:
        single_topoaln_file = "%s%s%s%s"%(outpath, os.sep, "%s_%s"%(li[0],
            li[1]), ext_topoaln)
        aafile = "%s%s%s%s"%(aapath, os.sep, li[3], ".fa")

        if myfunc.checkfile(single_topoaln_file, "single_topoaln_file") != 0:
            continue
        if myfunc.checkfile(aafile, "aafile") != 0:
            continue
        method_shrink = "2"
        method_plot = "mat"
        shrinkrateTM = "3"
        maxHoldLoop = "4"
        cmd = ["%s/drawMSATopo.py"%(binpath), single_topoaln_file,  "-aaseq",
                aafile, "-pdg", "yes", "-shrink", "yes", "-m-shrink",
                method_shrink, "-method", method_plot, "-text", "yes",
                "-krbias", "-shrinkrateTM", shrinkrateTM, "-max-hold-loop",
                maxHoldLoop]
        try:
            subprocess.check_output(cmd)
        except subprocess.CalledProcessError, e:
            print "[Command failed]", " ".join(cmd)
            print e


#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['seqidttype'] = 1
    g_params['isForceOverWrite'] = False
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

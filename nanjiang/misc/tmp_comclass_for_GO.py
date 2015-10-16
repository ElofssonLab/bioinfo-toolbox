#!/usr/bin/python
# filename: tmp_comclass_for_GO.py 
# Description:  from ana5.cmpclass.pairinfo.txt calculate cmpclass for each GO
# term, the two proteins in the pair should have common GO term
import os
import sys
import myfunc
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

usage_short="""
Usage: %s ana5Files 
"""%(progname)

usage_ext="""
Description:

OPTIONS:
  -gomap FILE   GO mapping file, from uniprot AC to GO 
  -goterm FILE  GO term description file, from GO ID to GO term
  -q            Quiet mode
  -h, --help    Print this help message and exit

Created 2013-09-23, updated 2013-09-23, Nanjiang Shu 
"""
usage_exp="""
Examples:
"""

GOLevelOneSet = set([
"GO:0016209",
"GO:0005488",
"GO:0003824",
"GO:0016247",
"GO:0042056",
"GO:0045499",
"GO:0036370",
"GO:0009055",
"GO:0030234",
"GO:0016530",
"GO:0060089",
"GO:0016015",
"GO:0001071",
"GO:0045735",
"GO:0000988",
"GO:0031386",
"GO:0004872",
"GO:0030545",
"GO:0005198",
"GO:0045182",
"GO:0005215"
])

SEL_GOID_SET = set([
"GO:0060089",
"GO:0009055",
"GO:0005488",
"GO:0005215",
"GO:0005198",
"GO:0003824"
])
SEL_GOID_LIST = list(SEL_GOID_SET)

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
# seqIDTGroup= [
#         0 ,20,
#         20 ,30,
#         30 ,40,
#         40 ,100
#         ]
seqIDTGroup = seqIDTGroupAll
numSeqIDTGroup = len(seqIDTGroup)/2

cmpclassList = ["IDT", "INV", "TM2GAP","TM2SEQ","SP2TM", "Mixed"]
numCmpClass = len(cmpclassList)
def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}

def GetClassIndex(cls, classList):#{{{
    try:
        return classList.index(cls)
    except ValueError:
        return len(classList)
#}}}
def GetSeqIDTGroupIndex(seqidt, seqIDTGroupList):#{{{
    numGroup = len(seqIDTGroupList)/2
    for i in xrange(numGroup):
        if seqidt >= seqIDTGroupList[i*2] and seqidt < seqIDTGroupList[i*2+1]:
            return i
    return numGroup
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
def ReadGOTerm(infile): #{{{
    hdl = myfunc.ReadLineByBlock(infile)
    dt = {}
    if hdl.failure:
        return 1
    lines = hdl.readlines()
    while lines != None:
        for line in lines:
            if line != "" and line[0] != "#":
                strs = line.split("\t")
                if len(strs) >= 2:
                    goid = strs[0].strip()
                    dt[goid] = strs[1].strip()
        lines = hdl.readlines()
    hdl.close()
    return dt
#}}}
def ReadGOMap(infile): #{{{
    try:
        goMapDict = {}
        fpin = open(infile, "r")
        cntLine = 0
        while 1:
            gomap = {}
            line = fpin.readline()
            if line == "":
                break
            cntLine += 1
            strs = line.split()
            seqid = strs[0]
            gomap['seqid'] = seqid
            gomap['function'] = []
            gomap['process'] = []
            gomap['location'] = []
            nF = int(strs[2])
            nP = int(strs[4])
            nC = int(strs[6])
            for item in ["function", "process", "location"]:
                if item == "function":
                    num = nF
                elif item == "process":
                    num = nP
                else:
                    num = nC
                for i in range(num):
                    rd = {}
                    line = fpin.readline()
                    cntLine += 1
                    strs1 = line.split("\t")
                    if len(strs1) < 2:
                        errmsg =  "Bad line in gomapfile at line number %d:"
                        print >> sys.stderr, errmsg%(cntLine)
                        continue
                    thisGOterm = strs1[0]
                    ancGOterm = strs1[1]
                    strs_this = thisGOterm.split("; ")
                    strs_anc = ancGOterm.split("; ")
                    rd['GOID'] = strs_this[0].strip()
                    ancList = []
                    for tss in strs_anc:
                        ancList.append(tss.split())
                    ancList.reverse()
                    rd['ancestor'] = ancList
                    gomap[item].append(rd)
            goMapDict[seqid] = gomap
        fpin.close()
        return goMapDict
    except IOError:
        print >> sys.stderr, "Failed to read file %s"%(infile)
        return {}
#}}}
def ReadPairInfo(infile):#{{{
    hdl = myfunc.ReadLineByBlock(infile)
    lst = []
    if hdl.failure:
        return []
    lines = hdl.readlines()
    while lines != None:
        for line in lines:
            if line != "" and line[0] != "#"
                strs = line.split()
                seqid1 = strs[0]
                seqid2 = strs[1]
                NtermState1 = strs[2]
                NtermState2 = strs[3]
                numTM1 = int(strs[4])
                numTM2 = int(strs[5])
                seqLen1 = int(strs[6])
                seqLen2 = int(strs[7])
                seqidt = float(strs[8])
                lst.append([seqid1, seqid2, NtermState1, NtermState2, numTM1,
                    numTM2, seqLen1, seqLen2, seqidt])
        lines = hdl.readlines()
    hdl.close()
    return lst
#}}}
def GetAncenstorGOList_LevelOne(seqid, gomapDict, gotype):#{{{
    try:
        gomap = goMapDict[seqid]
    except KeyError:
        msg = "No GO info for seqid %s"
        print >> sys.stderr, msg%(seqid)
        return []
    ancGOList = []
    for di in gomap[gotype]:
        ancinfo = di['ancestor']
        if (len(ancinfo) > 2 and "all" in ancinfo[0] and "GO:0003674" in
                ancinfo[1]):
            ancGOList += di['ancestor'][2]
    ancGOList = list(SEL_GOID_SET & set(ancGOList))
    return ancGOList
#}}}
def WritePairInfoWithGO(pairinfo, outfile):#{{{
    fpout = open(outfile, fpout)
    for li in pairinfo:
        fpout.write("%s %s %s %s %2d %2d %4d %4d %6.1f"%(li[0], li[1], li[2], li[3], li[4], li[5], li|6], li[7], li[8]))
        if len(li) > 9:
            for goid in li[9:]:
                fpout.write(" %s"%goid)
        fpout.write("\n")
    fpout.close()
#}}}

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = "./"
    infileList = []
    gomapfile = "/data3/wk/MPTopo/pfamAna_refpro/pfammap_from_uniprot/Pfam-A-full.seqfrompfamfasta.percentTMpro_scampi.perTM75_nseq20.Family.nr100.filter.fragmented.uniq.pfam.goinfowithancestor.txt"
    gotermfile = "/data3/wk/MPTopo/pfamAna_refpro/GO_analysis/GO_term.txt"
    anclevel = 2
    gotype = "function"

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            infileList.append( argv[i])
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
            elif argv[i] in ["-gomap", "--gomap"]:
                (gomapfile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-goterm", "--goterm"]:
                (gotermfile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            infileList.append(argv[i])
            i += 1

#     print len(gomapfile), gomapfile
#     lines = open(gomapfile, "r").readlines()
#     print lines

    if myfunc.checkfile(gomapfile, "GO map file") != 0:
        return 1
    if myfunc.checkfile(gotermfile, "GO Term file") != 0:
        return 1


    numInfile = len(infileList)
    if numInfile < len(cmpclassList):
        print >> sys.stderr, "input file less than len(cmpclassList)=%d"%(len(cmpclassList))

    goMapDict = ReadGOMap(gomapfile)
    goTermDict = ReadGOTerm(gotermfile)

    pairinfoDict = {}
    for infile in infileList:
        tag = ""
        for cls in cmpclassList:
            if infile.find(".%s."%(cls)) != -1:
                tag = cls
                break
        if tag == "":
            print >> sys.stderr, "bad infile %s"%(infile)
            return 1
        pairinfoDict[tag] = ReadPairInfo(infile)

    for tag in pairinfoDict:
        pairinfo = pairinfoDict[tag]
        for j in xrange(len(pairinfo)):
            tup = pairinfo[j]
            ancList1 = GetAncenstorGOList_LevelOne(tup[0], gomapDict,
                    gotype)
            ancList2 = GetAncenstorGOList_LevelOne(tup[1], gomapDict,
                    gotype)
            common_ancList = list(set(ancList1) & set(ancList2))
            pairinfo[j].append(common_ancList)

# output pairinfo with GO common term
    stemname = os.path.splitext(infileList[0].replace(".pairoinfo.txt", ""))[0]
    for tag in pairinfoDict:
        outfile = outpath + os.sep + stemname + ".%s"%(tag) + ".pairinfowithGOterm.txt"
        WritePairInfoWithGO(pairinfoDict[tag], outfile)


    tableCmpclassDict = {}
    for goid in SEL_GOID_SET:
        tableCmpclassDict[goid] = {}
        InitTableCmpClass(tableCmpclassDict[goid], numSeqIDTGroup, numCmpClass)

    for tag in pairinfoDict.keys():
        cmpclass = tag
        idxClass = GetClassIndex(cmpclass, cmpclassList)
        pairinfo = pairinfoDict[tag]
        for li in pairinfo:
            if len(li) > 9:
                seqidt = li[8]
                for goid in li[9:]:
                    idxGOID = GetClassIndex(goid, SEL_GOID_LIST)
                    idxGroup = GetSeqIDTGroupIndex(seqidt, seqIDTGroup)
                    tableCmpclassDict[goid]['freq'][idxGroup][idxClass] += 1


#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

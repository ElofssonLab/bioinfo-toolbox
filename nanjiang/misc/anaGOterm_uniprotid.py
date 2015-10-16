#!/usr/bin/python
# filename: anaGOterm_uniprotid.py 
# Description:  Given uniprot ID, analyze GO term frequency
import os
import sys
import myfunc
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

usage_short="""
Usage: %s uniprot-aclist [-o OUTFILE]
"""%(progname)

usage_ext="""
Description:
    Analyze GO term frequency given a list of uniprot id

OPTIONS:
  -o OUTFILE    Output the result to OUTFILE
  -gomap FILE   GO mapping file, from uniprot AC to GO 
  -goterm FILE  GO term description file, from GO ID to GO term
  -q            Quiet mode
  -h, --help    Print this help message and exit

Created 2009-06-08, updated 2013-03-14, Nanjiang Shu
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
        

def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}
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

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = "./"
    outfile = ""
    infile = ""
    gomapfile = "/data3/wk/MPTopo/pfamAna_refpro/pfammap_from_uniprot/Pfam-A-full.seqfrompfamfasta.percentTMpro_scampi.perTM75_nseq20.Family.nr100.filter.fragmented.uniq.pfam.goinfowithancestor.txt"
    gotermfile = "/data3/wk/MPTopo/pfamAna_refpro/GO_analysis/GO_term.txt"
    anclevel = 2
    gotype = "function"
    
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
            infile = argv[i]
            i += 1

#     print len(gomapfile), gomapfile
#     lines = open(gomapfile, "r").readlines()
#     print lines

    if myfunc.checkfile(infile) != 0:
        return 1
    if myfunc.checkfile(gomapfile, "GO map file") != 0:
        return 1
    if myfunc.checkfile(gotermfile, "GO Term file") != 0:
        return 1


    goMapDict = ReadGOMap(gomapfile)
    goTermDict = ReadGOTerm(gotermfile)
    idList = myfunc.ReadIDList(infile)


    freqDict = {}
    for goid in GOLevelOneSet:
        freqDict[goid] = 0

    for seqid in idList:
        try:
            gomap = goMapDict[seqid]
        except KeyError:
            msg = "No GO info for seqid %s"
            print >> sys.stderr, msg%(seqid)
            continue
        for di in gomap[gotype]:
            ancinfo = di['ancestor']
            if (len(ancinfo) > 2 and "all" in ancinfo[0] and "GO:0003674" in
                    ancinfo[1]):
                ancGOList = di['ancestor'][anclevel]
                for idd in ancGOList:
                    #debuging
                    if idd not in GOLevelOneSet:
                        print >> sys.stderr, "seqid=", seqid, "ancID=", idd, gomap[gotype]
                        continue
                    #debuging
                    if not idd in freqDict:
                        freqDict[idd] = 0
                    freqDict[idd] += 1

    #freqList = sorted(freqDict.items(), key=lambda x:x[1], reverse=True)
    freqList = sorted(freqDict.items(), key=lambda x:x[0], reverse=True)

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)

    printTupList = []
    for tup in freqList:
        try:
            term = goTermDict[tup[0]]
        except KeyError:
            term = ""
        printTupList.append((tup[0], term, tup[1]))
    maxSizeTerm = max([len(x[1]) for x in printTupList])
    maxSizeGOID = max([len(x[0]) for x in printTupList])  
    total = sum([x[2] for x in printTupList])

    for tup in printTupList:
        fpout.write("%-*s\t%-*s\t%4d\t%6.6f\n"%(maxSizeGOID, tup[0],
            maxSizeTerm, tup[1], tup[2], myfunc.FloatDivision(tup[2], total)))

    myfunc.myclose(fpout)


#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

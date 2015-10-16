#!/usr/bin/env python
# File name: get_GOancenstor.py
# Description:
import os
import sys
import myfunc
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

DATADIR3 = os.environ['DATADIR3']

usage_short="""
Usage: %s goinfoFile -anc ancestorTable
             ID [ID ...]
"""%(progname)

usage_ext="""
Description:
    Given the goinfo file, append the ancenstor information of GO
    Format of the output file is:
    UniprotAC F #F  P #P  C #C 
    GO:xxxxx; F:xxxxxx; xxx \t GOancestor list
    GO:xxxxx; F:xxxxxx; xxx \t GOancestor list
    ...

OPTIONS:
  -o OUTFILE        Output the result to OUTFILE
  -wgoidlist  FILE  Write unique GO idlist
  -anc        FILE  GO ancenstor file
  -q                Quiet mode
  -h, --help        Print this help message and exit

Created 2013-09-10, updated 2013-09-10, Nanjiang Shu 

Examples: 
    
"""
usage_exp="""
Examples:
    %s test.goinfo.txt -anc GO_ancestor.table.txt -o test.goinfo.withancestor.txt
"""

def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}

def ReadGOAnc(infile):
    hdl = myfunc.ReadLineByBlock(infile)
    goAncDict = {}
    if hdl.failure:
        return 1
    lines = hdl.readlines()
    while lines != None:
        for line in lines:
            if line != "" and line[0] != "#":
                strs = line.split(";")
                goid = strs[0].strip()
                goAncDict[goid] = line
        lines = hdl.readlines()
    hdl.close()

    return goAncDict

def ScanfGOInfo(line):#{{{
    goinfo = {}
    goinfo['process'] = []
    goinfo['function'] = []
    goinfo['location'] = []
    strs = line.split("\t")
    if len(strs) < 2:
        return {}
    seqid = strs[0].strip()
    goinfo['seqid'] = seqid
    for i in range(1, len(strs)):
        slist = strs[i].split("; ")
        msg = strs[i]
        if len(slist) >= 2:
            goid = slist[0]
            goterm = slist[1].split(":")[1]
            if slist[1][0:2] == "P:": #process
                goinfo['process'].append([goid, goterm, msg])
            elif slist[1][0:2] == "C:": #location
                goinfo['location'].append([goid, goterm, msg])
            elif slist[1][0:2] == "F:": #location
                goinfo['function'].append([goid, goterm, msg])
            else:
                print >> sys.stderr, "bad GO:%s"%(msg)
    return goinfo
#}}}
def ReadGOInfo(infile):#{{{
    hdl = myfunc.ReadLineByBlock(infile)
    GOinfoList = []
    if hdl.failure:
        return 1
    lines = hdl.readlines()
    while lines != None:
        for line in lines:
            if line != "" and line[0] != "#":
                goinfo = ScanfGOInfo(line)
                if goinfo != {}:
                    GOinfoList.append(goinfo)
        lines = hdl.readlines()
    hdl.close()

    return GOinfoList
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
    goancfile = "%s%s%s"%(DATADIR3, os.sep, "wk/MPTopo/pfamAna_refpro/GO_analysis/GO_ancenstor.MF.txt")
    uniqGOIDListFile = ""

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
            elif argv[i] in ["-wgoidlist", "--wgoidlist"]:
                (uniqGOIDListFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-outpath", "--outpath"]:
                (outpath, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-anc", "--anc"]:
                (goancfile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-l", "--l"] :
                (idListFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-uniprotdb", "--uniprotdb"] :
                (uniprotDBname, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            infile = argv[i]
            i += 1

    if myfunc.checkfile(infile) != 0:
        return 1
    if myfunc.checkfile(goancfile, "ancestor file") != 0:
        return 1


    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)

    GOinfoList = ReadGOInfo(infile)
    GOAncDict = ReadGOAnc(goancfile)

    if uniqGOIDListFile != "":#{{{
        fp = open(uniqGOIDListFile, "w")
        for goinfo in GOinfoList:
            for item in ["function", "process", "location"]:
            #for item in ["function"]:
                if len(goinfo[item]) > 0:
                    for li in goinfo[item]:
                        goidSet.add(li[0])
        for goid in goidSet:
            print >> fp, goid
        fp.close()
#}}}
    for goinfo in GOinfoList:
        nF = len(goinfo["function"])
        nP = len(goinfo["process"])
        nC = len(goinfo["location"])
        fpout.write("%s  F %1d  P %1d   C %1d\n"%(goinfo['seqid'], nF, nP, nC))
        for item in ["function", "process", "location"]:
            for li in goinfo[item]:
                try:
                    goancinfo = GOAncDict[li[0]]
                except KeyError:
                    goancinfo = ""
                fpout.write("%s\t%s\n"%(li[2],goancinfo))

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

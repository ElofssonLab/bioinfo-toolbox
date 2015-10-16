#!/usr/bin/env python
# Description:
# run hhsearch given the pairwise tableinfo file
import os
import sys
import myfunc
import itertools
import subprocess
progname =  os.path.basename(sys.argv[0])
usage = """
Usage: %s id1 id2 [-l pairlistfile] -hhprofile DIR
Description:
    Align two sequence by hhalign 
OPTIONS:
  -outpath    DIR   Set ouput path, default = ./
  -hhprofile  DIR   Multiple path can be set by evoking multiple times, ordered
  -hhalign    DIR   When this is supplied, hhr files will first be searched in
                    the path before running
  -hhopt STR        Set hhalign option, enclosed by quote
  -outpng           Output png file of the alignment, (default: not)
  -outhfa           Output hfa file of the alignment, (default: not)
  -forcewrite       force write the already exist file (default: not)
  -h, --help        Print this help message and exit

Created 2013-05-16, updated 2013-12-09, Nanjiang Shu 
"""%(progname)

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
def GetProfileFileName(hhprofilepathList, hhprofilepathMapDictList, seqid, #{{{
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

def RunHHAlignPairwise(pairlist,   #{{{ 
        hhprofilepathList, hhprofilepathMapDictList, 
        hhalignpathList, hhalignpathMapDictList, 
        outpath):
    HHBIN = os.environ['HHBIN']
    isOutputHFA = g_params['isOutputHFA']
    isOutputPNG = g_params['isOutputPNG']
    for tup in pairlist:

        (seqid1, seqid2) = tup

        hhrfile = "%s%s%s_%s.hhr"%(outpath, os.sep, seqid1, seqid2)

        if g_params['isUsePreBuildHHalignResult']:
            keystr = "%s_%s"%(seqid1, seqid2)
            tmp_hhrfile = GetProfileFileName(hhsearchpathList,
                    hhsearchpathMapDictList, keystr, ".hhr")
            if os.path.exists(tmp_hhrfile):
                hhrfile = tmp_hhrfile
            else:
                print >> sys.stderr, "hhrfile %s does not exist in"\
                            " the prebuilt path"%(hhrfile)

        outfile_fasta = "%s%s%s_%s.hfa"%(outpath, os.sep, seqid1, seqid2)
        pngfile = "%s%s%s_%s.png"%(outpath, os.sep, seqid1, seqid2)

        a3mfile = GetProfileFileName(hhprofilepathList,
                hhprofilepathMapDictList, seqid1, ".a3m")
        hhmfile = GetProfileFileName(hhprofilepathList,
                hhprofilepathMapDictList, seqid2, ".hhm")
        if a3mfile == "" or not os.path.exists(a3mfile):
            print >> sys.stderr, "a3mfile not found for %s. Ignore." %(seqid1)
        elif hhmfile == "" or not os.path.exists(hhmfile):
            print >> sys.stderr, "hhmfile not found for %s. Ignore." %(seqid2)
        else:
            cmd = ["%s/hhalign"%(HHBIN), "-i", a3mfile, "-t",  hhmfile,  "-o", hhrfile]
            if isOutputHFA:
                cmd += ["-ofas" , outfile_fasta]
            if isOutputPNG:
                cmd += ["-png", pngfile]

            cmd += g_params['hhalignopt'].split()
# Note, when using subprocess, a string will be quoted, therefore,
# hhalignopt should be splitted
            if not os.path.exists(hhrfile) or g_params['isForceOverWrite']:
                subprocess.call(cmd)
                if os.path.exists(hhrfile):
                    print hhrfile, "output"
    return 0
#}}}
def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = "./"
    listfile = ""
    hhprofilepathList = []
    hhalignpathList  = []

    nonPosArgList = []

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            nonPosArgList.append(argv[i])
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
            elif argv[i] in ["-hhalign", "--hhalign"] :
                (ss, i) = myfunc.my_getopt_str(argv, i)
                hhalignpathList.append(ss)
            elif argv[i] in ["-l", "--l"]:
                (listfile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-hhopt", "--hhopt"]:
                g_params['hhalignopt'] = argv[i+1]
                i += 2
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True
                i += 1
            elif argv[i] in ["-outpng", "--outpng"]:
                g_params['isOutputPNG'] = True; i += 1
            elif argv[i] in ["-outhfa", "--outhfa"]:
                g_params['isOutputHFA'] = True; i += 1
            elif argv[i] in ["-overwrite", "-forcewrite", "--forcewrite"]:
                g_params['isForceOverWrite'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            nonPosArgList.append(argv[i])
            i += 1

    if len(hhprofilepathList) < 1:
        print >> sys.stderr, "hhprofilepath not set. exit"
        return 1

    pairlist = []
    if len(nonPosArgList) == 2:
        pairlist.append((nonPosArgList[0], nonPosArgList[1]))
    elif len(nonPosArgList) > 0:
        msg = "Wrong number of non positional argument (%d)"
        print >> sys.stderr,msg%(len(nonPosArgList))
        return 1


    if listfile != "":
        pairlist += myfunc.ReadPairList(listfile)

    numpair = len(pairlist)
    if numpair <= 0:
        print >> sys.stderr, "input pair is 0, exit"
        return 1

    if not os.path.exists(outpath):
        try:
            os.makedirs(outpath)
        except OSError:
            pass
        if not os.path.exists(outpath):
            print >> sys.stderr, "failed to created outpath %s"%(outpath)
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


# read in index dictionary for hhalign result file
    hhalignpathMapDictList = []
    if len(hhalignpathList) > 0:
        g_params['isUsePreBuildHHalignResult'] = True
        for hhalignpath in hhalignpathList:
            hhalignmapfile = hhalignpath + os.sep + "id2pathmap.txt"
            if not os.path.exists(hhalignmapfile):
                print >> sys.stderr, "hhalignmapfile not exist. exit"
                hhalignpathMapDictList.append({})
            else:
                hhalignpathMapDictList.append(ReadSeqPathMapDict(hhalignmapfile))

    RunHHAlignPairwise(pairlist, 
            hhprofilepathList, hhprofilepathMapDictList, 
            hhalignpathList, hhalignpathMapDictList, 
            outpath)

    return 0
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['hhprofilepath'] = ""
    g_params['hhprofilepathMapDict'] = {}
    g_params['isForceOverWrite'] = False
    g_params['isOutputPNG'] = False
    g_params['isOutputHFA'] = False
    g_params['hhalignopt'] = ""
    g_params['isUsePreBuildHHalignResult'] = False
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

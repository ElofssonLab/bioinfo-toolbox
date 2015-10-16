#!/usr/bin/python
# Description:
import os
import sys
import myfunc
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))
DATADIR3 = os.environ['DATADIR3']

usage_short="""
Usage: %s FILE [FILE ...] [-o OUTFILE]
"""%(progname)

usage_ext="""
Description:
    Count the unique pair given file *.inverted.info.txt

OPTIONS:
  -o OUTFILE    Output the result to OUTFILE
  -q            Quiet mode
  -h, --help    Print this help message and exit

Created 2013-07-11, updated 2013-07-11, Nanjiang Shu
"""
usage_exp="""
Examples:
"""

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
def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}
def CountUniquePairInvertedInfo(infile, pfamidDefDict, fpout):
    idset1 = set([])
    idset2 = set([])
    numTMSet = set([])
    numInvPair = 0
    numAllPair = 0
    ratio = 0.0
    hdl = myfunc.ReadLineByBlock(infile)
    if hdl.failure:
        return 1
    lines = hdl.readlines()
    while lines != None:
        for line in lines:
            if line.find("General") == 0:
                strs = line.split()
                numInvPair = int(strs[1])
                numAllPair = int(strs[2])
                ratio = float(strs[3])
            if line.find("Pair") == 0:
                strs = line.split()
                id1 = strs[1]
                id2 = strs[2]
                NtermState1 = strs[3]
                NtermState2 = strs[4]
                numTM = int(strs[5])
                numTMSet.add(numTM)
                if NtermState1 == 'i':
                    idset1.add(id1)
                    idset2.add(id2)
                else:
                    idset1.add(id2)
                    idset2.add(id1)
        lines = hdl.readlines()
    hdl.close()
    pfamid = os.path.basename(infile).split(".")[0]
    try:
        pfamdef = pfamidDefDict[pfamid]
    except KeyError:
        pfamdef = "N/A"
    if len(idset1) > 0 or len(idset2) > 0:
        fpout.write("%-8s %20s %4d %4d %2d %8s   %5d %5d %6.3f\n"%(
            pfamid,
            pfamdef,
            len(idset1),
            len(idset2),
            len(numTMSet),
            str(list(numTMSet)),
            numInvPair,
            numAllPair,
            ratio
            ))


def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = "./"
    outfile = ""
    fileListFile = ""
    fileList = []
    pfamDefFile = "%s/data/pfam/pfam26.0/Pfam-A.clans.tsv"%(DATADIR3)

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
            elif argv[i] in ["-outpath", "--outpath"]:
                (outpath, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-l", "--l"] :
                (fileListFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            fileList.append(argv[i])
            i += 1

    if fileListFile != "":
        fileList += myfunc.ReadIDList(fileListFile)

    (pfamidDefDict, clanidDefDict) = ReadPfamDefFile(pfamDefFile)

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    for i in xrange(len(fileList)):
        CountUniquePairInvertedInfo(fileList[i], pfamidDefDict, fpout)

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

#!/usr/bin/python
# Description:
import os
import sys
import myfunc
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

usage_short="""
Usage: %s -map clanid2pfamid_mapfile -i pfam_percentTM_datafile [-o OUTFILE]
"""%(progname)

usage_ext="""
Description:
    Calculate clan percent TM from pfam percent TM

OPTIONS:
  -o OUTFILE    Output the result to OUTFILE
  -q            Quiet mode
  -h, --help    Print this help message and exit

Created 2013-05-02, updated 2013-05-02, Nanjiang Shu 
"""
usage_exp="""
Examples:
"""

def ReadPercentTM(infile):#{{{
    hdl = myfunc.ReadLineByBlock(infile)
    if hdl.failure:
        return {}
    percentTMDict = {}
    lines = hdl.readlines()
    cntline = 0
    while lines != None:
        for line in lines:
            cntline += 1
            if not line or line[0] == "#":
                continue
            strs = line.split()
            if len(strs) == 6:
                try:
                    pfamid = strs[0]
                    numTM = int(strs[1])
                    numSeq = int(strs[3])
                    percentTMDict[pfamid] = [numTM, numSeq]
                except (IndexError, ValueError):
                    msg = "Error in mapfile %s at line %d: \"%s\""
                    print >> sys.stderr, msg%(infile, cntline, line)
                    pass
        lines = hdl.readlines()
    hdl.close()
    return percentTMDict
#}}}
def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}

def GetPercentTMOfClan(pfamPercentTMDict, clanid2pfamidDict, fpout):
    for clanid in clanid2pfamidDict:
        sumTMSeq = 0
        sumSeq = 0
        pfamidlist = clanid2pfamidDict[clanid]
        for pfamid in pfamidlist:
            try:
                li = pfamPercentTMDict[pfamid]
                sumTMSeq += li[0]
                sumSeq += li[1]
            except (KeyError, IndexError):
                print >> sys.stderr, "pfamid %s not found in percentTMDict"%(pfamid)
                pass
        fpout.write("%s %d / %d = %.1f\n"%(clanid, sumTMSeq, sumSeq,
            float(sumTMSeq)/sumSeq*100))

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile = ""
    infile = ""
    mapfile = ""

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
            elif argv[i] in ["-map", "--map", "-mapfile"]:
                (mapfile, i) = myfunc.my_getopt_str(argv, i)
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
    if myfunc.checkfile(mapfile) != 0:
        return 1

    clanid2pfamidDict = myfunc.ReadFam2SeqidMap(mapfile)
    pfamPercentTMDict = ReadPercentTM(infile)

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    GetPercentTMOfClan(pfamPercentTMDict, clanid2pfamidDict, fpout)
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

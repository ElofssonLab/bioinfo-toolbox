#!/usr/bin/env python
# Get the Position of TM helices by given topology file
import os
import sys
import myfunc
usage = """
usage:  GetTMPosition.py [-i] topofile 
Description:
Get the Position of TM helices by given topology file

Options:
  -o      FILE    Set the output file
  -q              Quiet mode
  -h, --help      Print this help message and exit

Created 2012-11-30, updated 2012-11-30, Nanjiang Shu  
"""

def PrintHelp():
    print usage

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    topofile  = ""
    outfile = ""
    isGapLess = False

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            topofile = argv[i]
            isNonOptionArg = False
            i += 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-o", "--o"] :
                outfile = argv[i+1]
                i += 2
            elif argv[i] in ["-i", "--i"] :
                topofile = argv[i+1]
                i += 2
            elif argv[i] in ["-gapless", "--gapless"] :
                isGapLess = True
                i += 1
            elif argv[i] in ["-q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            topofile = argv[i]
            i += 1
    if topofile == "":
        print >> sys.stderr, "topofile not set. exit"
        return 1
    try: 
        (idList, annoList, seqList) =  myfunc.ReadFasta(topofile)
        fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
        for i in xrange(len(idList)):
            topo = seqList[i]
            seqid = idList[i]
            if isGapLess:
                topo = topo.replace("-","").replace(".","")
            posTMList = myfunc.GetTMPosition(topo)
            print >> fpout, seqid, posTMList
        myfunc.myclose(fpout)
    except (IOError, IndexError):
        pass
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

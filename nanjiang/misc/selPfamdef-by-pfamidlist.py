#!/usr/bin/env python
# Description:
import os
import sys
usage = """
usage:   selPfamdef-by-pfamidlist.py [-l FILE]  [-q]
                      ID [ID ...]
Description:

Options:
  -l      FILE    Set the idListFile
  -q              Quiet mode
  -h, --help      Print this help message and exit

Created 2009-06-08, updated 2011-11-07, Nanjiang Shu
"""

def PrintHelp():
    print usage

def ReadPfamDEList(infile):#{{{
    pfamDefDict={};
    try:
        fpin=open(infile);
        lines=fpin.readlines();
        fpin.close()
        for line in lines:
            strs=line.split(':');
            pfamid=strs[0].strip();
            pfamDefDict[pfamid]=strs[1].strip();
        return pfamDefDict;
    except IOError: 
        print >> sys.stderr, "%s: Error reading %s"%(sys.argv[0], infile);
        return 1;
#}}}

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = "./"
    outfile = ""
    idListFile = None
    idList = []
    pfamdeffile="/data3/wk/MPTopo/pfamAna/pfamA.seed.ac-delist"

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
                outpath = argv[i+1]
                i += 2
            elif argv[i] in ["-l", "--l"] :
                idListFile = argv[i+1]
                i += 2
            elif argv[i] in ["-q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            idList.append(argv[i])
            i += 1

    if idListFile != None:
        try:
            fp = open(idListFile,"r")
            idList += fp.read().split()
            fp.close()
        except IOError:        
            print >> sys.stderr, "file %s does not exist." %idListFile
    pfamDefDict = ReadPfamDEList(pfamdeffile)
    for pfamid in idList:
        if pfamid in pfamDefDict:
            print "%s : %s"%(pfamid, pfamDefDict[pfamid])

#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

#!/usr/bin/env python
# Description:
# get pairlist from cd-hit cluster file
import os
import sys
import myfunc
usage = """
usage:   cluster2pair_cdhit.py cdhit-clstr-file
                               
Description: get pairlist from cd-hit cluster file
Options:
  -o      FILE   Output the result to file
  -q             quiet mode
  -h, --help     Print this help message and exit

Created 2012-03-22, updated 2012-03-22, Nanjiang Shu  
"""

def PrintHelp():
    print usage

def ClstrToPairList(infile, fpout):
    try:
        fpin = open(infile, "r")
        lines = fpin.readlines()
        fpin.close()
    except IOError:
        print >> sys.stderr, "Failed to open file %s"%infile
        raise
    pairlist = []
    numLine = len(lines)
    i = 0
    while i < numLine:
        if lines[i][0] == ">":
            mbrlist = []
            j = 1
            while i+j < numLine and lines[i+j][0] != ">":
                mbrlist.append(lines[i+j].split(">")[1].split()[0].rstrip("."))
                j += 1
            if len(mbrlist) > 1:
                n_mbr = len(mbrlist)
#                fpout.write("%s\n"%lines[i])
                for m in xrange(n_mbr):
                    for n in xrange(m+1,n_mbr):
                        fpout.write("%s %s\n"%(mbrlist[m],mbrlist[n]))
            i += j
        else:
            i += 1

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile = ""
    clstrfile = ""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            clstrfile = argv[i]
            isNonOptionArg = False
            i += 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-o", "--o"]:
                outfile = argv[i+1]
                i += 2
            elif argv[i] in ["-q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            clstrfile = argv[i]
            i += 1

    if clstrfile == "" or not os.path.exists(clstrfile):
        print >> sys.stderr, "clstrfile file not set or not exist"
        return 1

    fpout = sys.stdout
    if outfile != "":
        try:
            fpout = open(outfile, "w")
        except IOError:
            print >> sys.stderr, "Failed to write to file %s"%outfile
            fpout = sys.stdout
    ClstrToPairList(clstrfile, fpout)
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

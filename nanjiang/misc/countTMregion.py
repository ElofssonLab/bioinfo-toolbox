#!/usr/bin/env python
# 
import sys, os
import myfunc

usage="""
Usage:  countTMregion.py topoFile
Options:
  -o FILE     Set the output file, default is stdout
  -q          Quiet mode
  -ni         Do not print the id name
  -h|--help   Print this help message and exit
Created 2011-03-31, updated 2012-06-10, Nanjiang Shu
"""

def PrintHelp():
    print usage

def main(g_params):#{{{
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    isQuiet = False
    isPrintIDName = True
    outfile = ""
    topofile = ""
    i = 1
    isNonOptionArg = False
    while i < numArgv:
        if isNonOptionArg == True:
            isNonOptionArg=False
            topofile=sys.argv[i]
            i = i + 1
        elif sys.argv[i] == "--":
            isNonOptionArg=True
            i = i + 1
        elif sys.argv[i][0] == "-":
            if sys.argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif sys.argv[i] in [ "-o" , "--o" , "-out"]:
                outfile=sys.argv[i+1]
                i = i + 2
            elif sys.argv[i] in [ "-q" , "--q"]:
                isQuiet=True
                i = i + 1
            elif sys.argv[i] in [ "-ni", "--ni", "-noid"] :
                isPrintIDName=False
                i = i + 1
            else:
                print "Error! Wrong argument:", sys.argv[i]
                return 1
        else:
            topofile = sys.argv[i]
            i = i + 1

    if topofile == "":
        print >> sys.stderr , "topofile not set. Exit."
        return 1
    elif not os.path.exists(topofile):
        print >> sys.stderr , "topofile %s doe not exist. Exit." %topofile
        return 1

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)

    fpin = open (topofile, "rb")
    if not fpin:
        print >> sys.stderr, "Failed to open input file %s"%(topofile)
        return 1
    unprocessedBuffer=""
    isEOFreached = False
    while 1:
        buff = fpin.read(BLOCK_SIZE)
        if len(buff) < BLOCK_SIZE:
            isEOFreached=True
        buff = unprocessedBuffer + buff
        recordList = []
        unprocessedBuffer = myfunc.ReadFastaFromBuffer(buff,recordList, isEOFreached)
        if len(recordList) > 0: 
            for rd in recordList:
                if isPrintIDName:
                    fpout.write("%s\t"%rd[0]);
                fpout.write("%d\n"% myfunc.CountTM(rd[2]))
        if isEOFreached == True:
            break
    fpin.close()

    myfunc.myclose(fpout)

    return 0
#}}}

if __name__ == '__main__' :
    g_params = {}
    BLOCK_SIZE=100000
    sys.exit(main(g_params))

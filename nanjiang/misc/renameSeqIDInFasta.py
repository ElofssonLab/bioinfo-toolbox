#!/usr/bin/env python
import sys
import re
import os
import myfunc 

usage="""
usage: renameSeqIDInFasta.py [-i] infile [-o OUTFILE]

Description: Rename annotation line in fasta file, e.g. using
             GetSeqIDFromAnnotation
             if map file supplied, replace it by mapped id

  -i        FILE  Input file
  -map      Map file, format, id TAB mappedid
  -o        FILE  Outputfile
  -h|--help       Print this help message and exit

Created 2012-06-08, updated 2013-02-11, Nanjiang Shu 

Examples:
    renameSeqIDInFasta.py input.fa -o renamedid.fa
"""

BLOCK_SIZE=100000

def PrintHelp():
    print usage

def ReadMapFile(infile):#{{{
    mapDict = {}
    hdl = myfunc.ReadLineByBlock(infile)
    if hdl.failure:
        return mapDict

    lines = hdl.readlines()
    while lines != None:
        for line in lines:
            if not line or line[0] == "#":
                continue
            strs = line.split("\t")
            try:
                mapDict[strs[0]] = strs[1]
            except IndexError:
                pass
        lines = hdl.readlines()
    hdl.close()
    return mapDict
#}}}
def main(g_params):#{{{
    numArgv = len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile  = ""
    infile  = ""
    mapfile = ""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            infile = sys.argv[i+1]
            isNonOptionArg=False
            i = i + 1
        elif sys.argv[i] == "--":
            isNonOptionArg=True
            i = i + 1
        elif sys.argv[i][0] == "-":
            if sys.argv[i] ==  "-h" or  sys.argv[i] == "--help":
                PrintHelp()
                return 1
            elif sys.argv[i] in [ "-i", "--i", "--infile"]:
                infile = sys.argv[i+1]
                i = i + 2
            elif sys.argv[i] in [ "-o", "--o", "-outfile", "--outfile"]:
                outfile =sys.argv[i+1]
                i = i + 2
            elif sys.argv[i] in [ "-map", "--map"]:
                mapfile =sys.argv[i+1]
                i = i + 2
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i])
                return 1
        else:
            infile = sys.argv[i]
            i+=1

    if infile == "":
        print >> sys.stderr,"Error! Topology file not set."
        return 1


    isMapSupplied = False
    mapDict = {}
    if mapfile != "" and os.path.exists(mapfile):
        mapDict = ReadMapFile(mapfile)
        isMapSupplied = True

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)

    try: 
        fpin = open (infile, "rb")
        unprocessedBuffer=""
        isEOFreached = False
        processedTopoIDSet = set([])
        while 1:
            buff = fpin.read(BLOCK_SIZE)
            if len(buff) < BLOCK_SIZE:
                isEOFreached=True
            buff = unprocessedBuffer + buff
            recordList = []
            unprocessedBuffer = myfunc.ReadFastaFromBuffer(buff,recordList, isEOFreached)
            if len(recordList) > 0: 
                for rd in recordList:
# if not isMapSupplied use seqid as the first word
# else, using mapped id, and if keyerror, use the original annotation line
                    if not isMapSupplied:
                        renamedID = rd[0]
                        fpout.write(">%s %s\n"%(renamedID, rd[1]))
                    else:
                        try:
                            renamedID = mapDict[rd[0]]
                            fpout.write(">%s %s\n"%(renamedID, rd[1]))
                        except KeyError:
                            msg = "ID %s not found in mapDict"
                            print >> sys.stderr, msg%(rd[0])
                            fpout.write(">%s\n"%(rd[1]))
                    fpout.write("%s\n"%(rd[2]))
            if isEOFreached == True:
                break
        fpin.close()
    except IOError:
        print >> sys.stderr, "Failed to read input file %s"%(infile)
        return 1

    myfunc.myclose(fpout)

#}}}
def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    # Check argv
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

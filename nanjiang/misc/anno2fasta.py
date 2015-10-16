#!/usr/bin/env python
# Description:
import os
import sys
import myfunc
usage = """
usage:  anno2fasta.py annotationFile -seqdb FILE [-o OUTFILE] 

Description: given the annotation file, append the fasta sequence (aligned or
             non-aligned). Everything in the annotation is kept but the
             sequence is retrieved from seqdb (recognized by seqID)

Options:
  -o FILE    Output the result to FILE, (default: stdout)
  -h, --help      Print this help message and exit

Created 2012-08-30, updated 2012-08-30, Nanjiang Shu 
"""

def PrintHelp():
    print usage

def GetSeqDict(fastafile):#{{{
    seqDict = {}
    (idList, seqList) = myfunc.ReadFasta_without_annotation(fastafile)
    for i in xrange(len(idList)):
        seqDict[idList[i]] = seqList[i]
    return seqDict
#}}}

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile = ""
    seqdbfile = ""
    infile = ""

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
            elif argv[i] in ["-outfile", "--outfile"]:
                outfile = argv[i+1]
                i += 2
            elif argv[i] in ["-seqdb", "--seqdb"]:
                seqdbfile = argv[i+1]
                i += 2
            elif argv[i] in ["-q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            infile = argv[i]
            i += 1
    if infile == "":
        print >> sys.stderr, "annotation file not set"
        return 1
    elif not os.path.exists(infile):
        print >> sys.stderr, "annotation file %s does not exist"%(infile)
        return 1
    if seqdbfile == "":
        print >> sys.stderr, "seqdbfile file not set"
        return 1
    elif not os.path.exists(seqdbfile):
        print >> sys.stderr, "seqdbfile file %s does not exist"%(seqdbfile)
        return 1
    seqDict = GetSeqDict(seqdbfile)
    if seqDict == {}:
        print >> sys.stderr, "Failed to read seqdbfile %s"%(seqdbfile)
        return 1
    (idList, annoList, contentList) = myfunc.ReadFasta(infile)
    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    for i in xrange(len(idList)):
        seqid = idList[i]
        try:
            seq = seqDict[seqid]
            fpout.write(">%s\n"%(annoList[i]))
            fpout.write("%s\n"%(seq))
            if contentList[i] != "":
                fpout.write("%s\n"%(contentList[i]))
        except KeyError:
            print >> sys.stderr, "seqid %s not found in seqdb"%(seqid)

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

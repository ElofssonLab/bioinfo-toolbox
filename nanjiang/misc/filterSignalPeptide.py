#!/usr/bin/env python

# Filename: filterSignalPeptide.py 
# filter signal peptide given the predicted TM topology

# ChangeLog 2013-07-11 
#   add option -deleteseq
# ChangeLog 2013-08-22
#   if newtopo has no TM helices, delete the entry

import os
import sys
import libtopologycmp as lcmp
import myfunc

BLOCK_SIZE = myfunc.BLOCK_SIZE
progname = os.path.basename(sys.argv[0])
usage="""
usage: %s -topo topoFile -sig sigpepFile [-o OUTFILE]

Description: Filter signal peptide given the topoFile and signal peptide
             prediction file
Options:
  -o OUTFILE   Output the result to outfile
  -deleteseq   Delete protein with signal peptide predicted, by default remove
               only the signal peptide
  -q           Quiet mode
  -h, --help   Print this help message and exit

Selection control options:

Created 2011-12-16, updated 2013-08-22, Nanjiang Shu 
"""%(progname)

def PrintHelp():
    print usage

def RemoveSigPep(topo, posTM):
    newtopo = ""
    if len(posTM) > 1:
        state = topo[posTM[1][0]-1]
        sizeToReplace = posTM[1][0]
        li = [state]*sizeToReplace
        newtopo += "".join(li)
        newtopo += topo[posTM[1][0]:]
    else:
        state = topo[posTM[0][1]]
        li = [state] * len(topo)
        newtopo = "".join(li)
    return newtopo

def FilterSignalPeptide(topofile, sigpepDict, outfile,
        isDeleteSeqWithSignalPeptide):
    hdl = myfunc.ReadFastaByBlock(topofile)
    if hdl.failure:
        return 1

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    recordList = hdl.readseq()
    while recordList != None:
        for rd in recordList:
            try:
                sp_pos = sigpepDict[rd.seqid]
            except KeyError:
                sp_pos = -1
            if sp_pos != -1:
                if isDeleteSeqWithSignalPeptide:
                    newtopo = ""
                else:
                    newtopo = lcmp.FilterSignalPeptideInTopology(rd.seq, sp_pos)
            else:
                newtopo = rd.seq
            if newtopo != "" and myfunc.CountTM(newtopo) > 0:
                fpout.write(">%s\n"%(rd.description))
                fpout.write("%s\n"%(newtopo))
        recordList = hdl.readseq()
    hdl.close()
    myfunc.myclose(fpout)
    return 0

def main():#{{{
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    argv = sys.argv
    topofile = ""
    sigpepfile = ""
    outfile = ""
    isQuiet = False
    isDeleteSeqWithSignalPeptide = False
    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            isNonOptionArg = False
            i += 1
        elif sys.argv[i] == "--":
            isNonOptionArg=True
            i += 1
        elif sys.argv[i][0] == "-":
            if sys.argv[i] in [ "-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in [ '-o',  '--o', "-outfile", "--outfile"]:
                outfile, i = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ['-topo', '--topo']:
                topofile, i = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ['-sig', '--sig']:
                sigpepfile, i = myfunc.my_getopt_str(argv, i)
            elif argv[i] == "-q":
                isQuiet=True; i += 1
            elif argv[i] in ["-deleteseq", "--deleteseq"]:
                isDeleteSeqWithSignalPeptide=True; i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            print >> sys.stderr, "Error! Wrong argument:", argv[i]
            return 1

    if myfunc.checkfile(topofile, 'topofile') != 0:
        return 1
    if myfunc.checkfile(sigpepfile, 'sigpepfile') != 0:
        return 1

    sigpepDict = lcmp.ReadSignalPDict(sigpepfile)
    FilterSignalPeptide(topofile, sigpepDict, outfile, isDeleteSeqWithSignalPeptide)

    return 0

#}}}
if __name__ == '__main__' :
    sys.exit(main())

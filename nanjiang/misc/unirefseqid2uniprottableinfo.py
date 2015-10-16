#!/usr/bin/env python
# Description:  given the list of uniref seqid, match uniprot table info 

# 
import os
import sys
import myfunc
usage = """
usage:  unirefseqid2uniprottableinfo.py -l seqidlistfile [ID [ID...]]
Description: Given the list of uniref seqid, match uniprot tableinfo
             # output format
             # uniref_seqid uniprot_AC Length GN OS OC
             tab delimited
Options:
  -datafile FILE  Set uniprot table info file (default: 
                  /data3/data/uniprot/uniprot_trembl.tableinfo)
  -q              Quiet mode
  -h, --help      Print this help message and exit

Created 2012-06-01, updated 2012-06-01, Nanjiang Shu
"""

def PrintHelp():
    print usage

def ReadUniprotInfo(infile):
    try:
        uniprotInfoDict = {}
        fpin = open(infile, "r")
        line = fpin.readline()
        while line:
            line = line.rstrip("\n")
            if line[0] != "#":
                strs = line.split("\t")
                accession = strs[0]
                aclist = accession.split(";")
                for ac in aclist:
                    if ac != "":
                        uniprotInfoDict[ac] = line
            line = fpin.readline()
        fpin.close()
        return uniprotInfoDict
    except IOError:
        print >> sys.stderr, "Failed to read file ",infile
        return {}

def UnirefSeqid2UniprotTableInfo(idList, uniprotInfoDict, fpout):#{{{
    for seqid in idList:
        ac = ""
        if seqid.find("_") != -1:
            ac = seqid.split("_")[1]
            if ac.find("-") != -1:
                ac = ac.split("-")[0]
        else:
            ac = seqid
        if ac in uniprotInfoDict:
            fpout.write("%s\t%s\n" % (seqid, uniprotInfoDict[ac]))
    return 0
#}}}

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile = ""
    datafile = "/data3/data/uniprot/uniprot_trembl.tableinfo"
    idList = []
    idListFile = ""

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
            elif argv[i] in ["-o", "--o", "-outfile", "--outfile"]:
                outfile = argv[i+1]
                i += 2
            elif argv[i] in ["-datafile", "--datafile"] :
                datafile = argv[i+1]
                i += 2
            elif argv[i] in ["-l", "--l", "-listfile", "--listfile"] :
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

    if os.path.exists(idListFile):
        idList += myfunc.ReadIDList(idListFile)

    if (len(idList)) < 1:
        print >> sys.stderr, "id not set. Exit"
        return 1

    if not os.path.exists(datafile):
        print >> sys.stderr, "datafile %s not set or not exists. Exit" %(datafile)
        return 1

    uniprotInfoDict = ReadUniprotInfo(datafile)
    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)

    UnirefSeqid2UniprotTableInfo(idList, uniprotInfoDict, fpout)

    myfunc.myclose(fpout)
    return 0
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    BLOCK_SIZE=100000
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

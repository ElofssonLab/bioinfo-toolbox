#!/usr/bin/env python
# Description:
import os
import sys
import tempfile
import myfunc
usage = """
usage:  buildPfamid_seqid_list_from_uniprot.py uniprotTableinfofile -outpath DIR
    Description:
        Create two lists
    list 1:
        Pfamid1 seqid1 seqid2 seqid3 ...
        pfamid2 seqid1 seqid2 seqid3 ...
    list 2:
        seqid1 pfamid1 pfamid2 ...
        seqid2 pfamid1 pfamid2 ...     
Options:
  -pfamclandef FILE Set pfamclan definition file, (default:
                    /data3/data/pfam/pfam26.0/Pfam-A.clans.tsv)
  -q                Quiet mode
  -h, --help        Print this help message and exit

Created 2012-06-07, updated 2012-06-07, Nanjiang Shu 

Examples:
"""

rundir = os.path.dirname(sys.argv[0])
binpath = rundir

def PrintHelp():
    print usage

def ReadPfamClanDBFile(dbfile):#{{{
    try:
        pfamid2clanidDict = {}
        fpin = open(dbfile,"r")
        lines = fpin.readlines()
        fpin.close()
        for line in lines:
            if line and line[0] != "#":
                strs = line.split("\t")
                if len(strs) >= 2:
                    pfamid = strs[0]
                    clanid = strs[1]
                    if clanid == "\N":
                        clanid = pfamid
                    pfamid2clanidDict[pfamid] = clanid
        return pfamid2clanidDict
    except IOError:
        print >> sys.stderr, "Failed to read dbfile ", dbfile
        return {}
#}}}
def ReadUniprotTableInfo(infile):#{{{
    try:
        seqid2pfamidDict = {}
        fpin = open(infile, "r")
        line = fpin.readline()
        while line:
            if line != "" and line[0] != "#":
                strs = line.split("\t")
                if len(strs) >= 6:
                    seqidList = strs[0].split(";")
                    seqidList = filter(None, seqidList)
                    seqid = seqidList[0]
                    if seqid == "":
                        print >> sys.stderr, "Bad record, line =", line
                        continue
                    pfamidList = strs[5].split(";")
                    pfamidList = filter(None, pfamidList)
                    if len(pfamidList) > 0:
                        seqid2pfamidDict[seqid] = pfamidList
            line = fpin.readline()
        fpin.close()
        return seqid2pfamidDict
    except IOError:
        print >> sys.stderr, "Failed to read file %s" %infile
        return {}
#}}}
def WriteIDMapFile(outfile, idmapdict):#{{{
    try:
        fpout = open(outfile, "w")
        for idd in idmapdict:
            fpout.write("%s %d"%(idd, len(idmapdict[idd])))
            for idd2 in idmapdict[idd]:
                fpout.write(" %s"%(idd2))
            fpout.write("\n")
        fpout.close()
        return 0
    except IOError:
        print >> sys.stderr, "Failed to write to file %s" %outfile
        return 1
#}}}

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = "./"
    uniprotTableInfoFile = ""
    pfamclandefFile = "/data3/data/pfam/pfam26.0/Pfam-A.clans.tsv"

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            isNonOptionArg = False
            uniprotTableInfoFile = argv[i]
            i += 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-outpath", "--outpath"] :
                outpath = argv[i+1]
                i += 2
            elif argv[i] in ["-pfamclandef", "--pfamclandef"] :
                pfamclandefFile = argv[i+1]
                i += 2
            elif argv[i] in ["-q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            uniprotTableInfoFile = argv[i]
            i += 1
            
    if uniprotTableInfoFile == "":
        print >> sys.stderr, "Error! uniprotTableInfoFile not set. Exit."
        return 1
    elif not os.path.exists(uniprotTableInfoFile):
        print >> sys.stderr, "Error! tablefile %s does not exist. Exit." % (
                uniprotTableInfoFile)
        return 1

    if not os.path.exists(outpath):
        os.system("mkdir -p %s"%outpath)

    rootname=os.path.basename(os.path.splitext(uniprotTableInfoFile)[0]);

    seqid2pfamidDict = ReadUniprotTableInfo(uniprotTableInfoFile)
    pfamid2clanidDict = ReadPfamClanDBFile(pfamclandefFile)

    seqid2clanidDict = {}
    pfamid2seqidDict = {}
    clanid2seqidDict = {}

    for seqid in seqid2pfamidDict:
        clanidSet = set([])
        pfamidList = seqid2pfamidDict[seqid]
        if len(pfamidList) <= 0:
            continue
        for pfamid in pfamidList:
            if pfamid in pfamid2clanidDict:
                clanidSet.add(pfamid2clanidDict[pfamid])

            if pfamid in pfamid2seqidDict:
                pfamid2seqidDict[pfamid].append(seqid)
            else:
                pfamid2seqidDict[pfamid] = []
                pfamid2seqidDict[pfamid].append(seqid)
        if len(clanidSet) > 0:
            seqid2clanidDict[seqid] = list(clanidSet)
        for clanid in clanidSet:
            if clanid in clanid2seqidDict:
                clanid2seqidDict[clanid].append(seqid)
            else:
                clanid2seqidDict[clanid] = []
                clanid2seqidDict[clanid].append(seqid)


    outListFile1 = outpath + os.sep + rootname + ".seqid2pfamid"
    outListFile2 = outpath + os.sep + rootname + ".pfamid2seqid"
    outListFile3 = outpath + os.sep + rootname + ".seqid2clanid"
    outListFile4 = outpath + os.sep + rootname + ".clanid2seqid"

    WriteIDMapFile(outListFile1, seqid2pfamidDict)
    WriteIDMapFile(outListFile2, pfamid2seqidDict)
    WriteIDMapFile(outListFile3, seqid2clanidDict)
    WriteIDMapFile(outListFile4, clanid2seqidDict)

    print "result output to"
    print "\t", outListFile1
    print "\t", outListFile2
    print "\t", outListFile3
    print "\t", outListFile4

#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :

    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

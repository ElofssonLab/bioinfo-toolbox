#!/usr/bin/env python
# Description:
import os
import sys
import myfunc
usage = """
usage:   pfamid2clanid.py [-l FILE] 
                      pfamid [pfamid ...]
Description: Mapping pfamid to clanid, if no_clan, set as pfamid

Options:
  -o      FILE    Output the result to file
  -omap   FILE    Output the clanid to pfamid map file
  -l      FILE    Set the idListFile
  -db     FILE    Set dbfile, (default:
                  /data3/data/pfam/pfam26.0/Pfam-A.clans.tsv)
  -q              Quiet mode
  -h, --help      Print this help message and exit

Created 2012-05-29, updated 2012-05-29, Nanjiang Shu 
"""

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
def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    dbfile = "/data3/data/pfam/pfam26.0/Pfam-A.clans.tsv"
    outfile = ""
    outfile_map = ""
    idListFile = None
    idList = []

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
            elif argv[i] in ["-omap", "--omap"]:
                outfile_map = argv[i+1]
                i += 2
            elif argv[i] in ["-db", "--db"]:
                dbfile = argv[i+1]
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

    pfamid2clanidDict = ReadPfamClanDBFile(dbfile)

    clanidSet = set([])
    clanid2pfamidDict = {}
    for pfamid in idList:
        if pfamid in pfamid2clanidDict:
            clanid =  pfamid2clanidDict[pfamid]
            clanidSet.add(clanid)
            if clanid in clanid2pfamidDict:
                clanid2pfamidDict[clanid].append(pfamid)
            else:
                clanid2pfamidDict[clanid] = []
                clanid2pfamidDict[clanid].append(pfamid)
        else:
            print >> sys.stderr, "pfamid %s not in dict"%(pfamid)

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False) 
    for clanid in clanidSet:
        print >> fpout, clanid
    myfunc.myclose(fpout)

    if outfile_map != "":
        fpout_map = myfunc.myopen(outfile_map, None, "w", False)
        if fpout_map:
            for clanid in clanid2pfamidDict:
                fpout_map.write("%s %d"%(clanid, len(clanid2pfamidDict[clanid])))
                for pfamid in clanid2pfamidDict[clanid]:
                    fpout_map.write(" %s"%(pfamid))
                fpout_map.write("\n")
            myfunc.myclose(fpout_map)
    return 0
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

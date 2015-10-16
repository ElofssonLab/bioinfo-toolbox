#!/usr/bin/env python
# Description:
import os
import sys
import tempfile
import myfunc
usage = """
usage:  findPfam-for-seqpair.py seqpairListFile

Description: Find Pfam ID given a pair of sequence id
The pair should be from the same family
        
Options:
  -q              Quiet mode
  -h, --help      Print this help message and exit

Created 2012-05-14, updated 2012-05-14, Nanjiang Shu 

Examples:
"""

rundir = os.path.dirname(sys.argv[0])
binpath = rundir

def ReadIDMap2(infile):#{{{
    try:
        idmap2Dict = {}
        fpin = open(infile,"r")
        lines = fpin.readlines()
        for line in lines:
            if line:
                strs = line.split()
                if len(strs) > 2:
                    idmap2Dict[strs[0]] = strs[2:]
                else:
                    print >> sys.stderr, "broken item in file %s: line \"%s\"" \
                            % (infile, line)
        fpin.close()
        return idmap2Dict
    except IOError:
        print "Failed to read listfile ", infile
        return {}
#}}}

def PrintHelp():
    print usage

def ReadIDList(infile):
    try:
        fpin = open(infile,"r")
        li = fpin.read().split()
        fpin.close()
        return li
    except IOError:
        print "Failed to read listfile ", infile
        return []

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    idmap2file = "/data3/wk/MPTopo/pfamAna/pfamfullseq.selTM_uniq.idmap2"
    pairlistfile =  ""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            pairlistfile = argv[i]
            isNonOptionArg = False
            i += 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            pairlistfile = argv[i]
            i += 1
    if pairlistfile == "":
        print >> sys.stderr, "Error! pairlistfile not set. Exit."
        return 1

    idmap2Dict = ReadIDMap2(idmap2file)
    seqidPairList = myfunc.ReadPairList(pairlistfile)
    for pair in seqidPairList:
        id1 = pair[0]
        id2 = pair[1]
        pfamidlist1 = idmap2Dict[id1]
        pfamidlist2 = idmap2Dict[id2]
        interlist = list(set(pfamidlist1) & set(pfamidlist2))
        unionlist = list(set(pfamidlist1) | set(pfamidlist2))
        if len(unionlist) == 1:
            print "%-10s %-10s %s" %(id1, id2, interlist[0])
        else:
            print >> sys.stderr, "%-10s %-10s" % (id1, id2), \
                    "inter", interlist, "union", unionlist

#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :

    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

#!/usr/bin/env python
# Description:
import os
import sys
import tempfile
usage = """
usage:   build-singlefam-idmap.py

Description: Build idmap of which one seqid only belongs to one family
        
Options:
  -q              Quiet mode
  -h, --help      Print this help message and exit

Created 2012-05-14, updated 2012-05-14, Nanjiang Shu 

Examples:
    ./build-singlefam-idmap.py  > result.txt
"""

rundir = os.path.dirname(sys.argv[0])
binpath = rundir

def ReadIDMap2(infile):
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


def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)

    idmap2file = "/data3/wk/MPTopo/pfamAna/pfamfullseq.selTM_uniq.idmap2"

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

    idmap2Dict = ReadIDMap2(idmap2file)

    singlefamIDMapDict = {}
    for seqid in idmap2Dict:
        pfamidlist = idmap2Dict[seqid]
        if len(pfamidlist) == 1:
            pfamid = pfamidlist[0]
            if pfamid in singlefamIDMapDict:
                singlefamIDMapDict[pfamid].append(seqid)
            else:
                singlefamIDMapDict[pfamid] = []
                singlefamIDMapDict[pfamid].append(seqid)
    for pfamid in singlefamIDMapDict:
        sys.stdout.write("%s %d"%(pfamid, len(singlefamIDMapDict[pfamid])))
        for seqid in singlefamIDMapDict[pfamid]:
            sys.stdout.write(" %s"%(seqid))
        sys.stdout.write("\n")

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

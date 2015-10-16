#!/usr/bin/env python
# Description:
import os
import sys
import tempfile
usage = """
usage: buildPfamid_seqid_list.py -dbname seqdbname
    Description:
        Create two lists
    list 1:
        Pfamid1 seqid1 seqid2 seqid3 ...
        pfamid2 seqid1 seqid2 seqid3 ...
    list 2:
        seqid1 pfamid1 pfamid2 ...
        seqid2 pfamid1 pfamid2 ...     
Options:
  -q              Quiet mode
  -l  FILE        Set pfamID list file
  -h, --help      Print this help message and exit

Created 2012-05-14, updated 2012-05-14, Nanjiang Shu 

Examples:
    buildPfamid_seqid_list.py -dbname pfamfullseq.selTM_uniq -l  pfamAfull.selTM.pfamidlist
"""

rundir = os.path.dirname(sys.argv[0])
binpath = rundir

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

    outpath = ""
    pfamidListFile = ""
    dbname =  ""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            isNonOptionArg = False
            i += 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-dbname", "--dbname"]:
                dbname = argv[i+1]
                i += 2
            elif argv[i] in ["-l", "--l"]:
                pfamidListFile = argv[i+1]
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
            print >> sys.stderr, "Error! Wrong argument:", argv[i]
            return 1
    if dbname == "":
        print >> sys.stderr, "Error! dbname not set. Exit."
        return 1
    if pfamidListFile == "":
        print >> sys.stderr, "Error! pfamidListFile not set. Exit."
        return 1

    outpath = os.path.dirname(dbname)
    if outpath == "":
        outpath = "."
    basename = os.path.basename(dbname)

    outListFile1 = outpath + os.sep + basename + ".idmap1"
    outListFile2 = outpath + os.sep + basename + ".idmap2"
    fpout1 = open(outListFile1, 'w')
    fpout2 = open(outListFile2, 'w')

    seqidToPfamIDDict = {}

    pfamidList = ReadIDList(pfamidListFile)
    tmpdir = tempfile.mkdtemp()
    for pfamid in pfamidList:
        tmpseqfile = tmpdir + os.sep + "seq.fa"
        tmpidlistfile = tmpdir + os.sep + "seq.idlist"
        cmd = "python %s/my_extractdb.py %s -dbname %s -o %s"%(binpath, pfamid,
                dbname, tmpseqfile)
        os.system(cmd)
        if os.path.getsize(tmpseqfile) <= 0:
            print >> sys.stderr, pfamid + " ignored because of emtpy seqfile"
            continue
        cmd = "python %s/getfastaid.py %s -o %s" % (binpath, tmpseqfile,
                tmpidlistfile)
        os.system(cmd)
        if os.path.getsize(tmpidlistfile) <= 0:
            print >> sys.stderr, pfamid + " ignored because of emtpy seqidlistfile"
            continue

        seqIDList = ReadIDList(tmpidlistfile)
        if len(seqIDList) > 0:
            fpout1.write("%s %d"%(pfamid, len(seqIDList)))
            for seqid in seqIDList:
                fpout1.write(" %s"%(seqid))
            fpout1.write("\n")

            # add to seqidToPfamIDDict
            for seqid in seqIDList:
                if seqid in seqidToPfamIDDict:
                    seqidToPfamIDDict[seqid].append(pfamid)
                else:
                    seqidToPfamIDDict[seqid] = []
                    seqidToPfamIDDict[seqid].append(pfamid)

        os.remove(tmpseqfile)
        os.remove(tmpidlistfile)

    for seqid in seqidToPfamIDDict:
        fpout2.write("%s %d"%(seqid, len(seqidToPfamIDDict[seqid])))
        for pfamid in seqidToPfamIDDict[seqid]:
            fpout2.write(" %s"%(pfamid))
        fpout2.write("\n")

    close(fpout1)
    close(fpout2)

    print "result output to"
    print "\t", outListFile1
    print "\t", outListFile2

            

#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :

    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

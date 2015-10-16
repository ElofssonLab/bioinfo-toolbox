#!/usr/bin/python
# Filename: get_TMproList_per_family.py
# Description: Get the fasta file of TM proteins for each protein family
# 
import os
import sys
import subprocess
import myfunc

progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

usage_short="""
Usage: %s -pfamid2seqid FILE -topodb DBNAME -seqdb DBNAME -outpath DIR
"""%(progname)

usage_ext="""
Description:
    Get fasta file of TM proteins (defined in the topology database) for each 
    protein family. The protein family and its member proteins are listed in 
    the pfamid2seqid file.

OPTIONS:
  -outpath DIR          Output the fasta file to outpath
  -pfamid2seqid  FILE   A mapping file defining pfam families and their member proteins
  -topodb DBNAME        Database name for membrane protein topology
  -seqdb  DBNAME        Database name for amino acid sequences
  -q                    Quiet mode
  -h, --help            Print this help message and exit

Created 2013-11-27, updated 2013-11-27, Nanjiang Shu
"""
usage_exp="""
Examples:
"""

def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}

def GetTMProList_per_family(pfamid2seqidDict, idSet_topo, hdl_seq,#{{{ 
        hdl_topo, outpath):
    """
    Get fasta file of TM proteins (defined in the topology database) for each 
    protein family. The protein family and its member proteins are listed in 
    the pfamid2seqid file.

    Input:
        pfamid2seqidDict    {pfamid: [seqid1, seqid2, ...], pfamid:[]
        idSet_topo          A set of seqids of TM proteins with topologies
        outpath             path to output the result
    Output:
        outpath/pfamid.fa
    """
    for pfamid in pfamid2seqidDict:
        idSet_pfam = set(pfamid2seqidDict[pfamid]) & idSet_topo
        outfile_seq = "%s%s%s.fa"%(outpath, os.sep, pfamid)
        outfile_topo = "%s%s%s.topo"%(outpath, os.sep, pfamid)
        try:
            fpout = open(outfile_seq, "w")
            for seqid in idSet_pfam:
                record = hdl_seq.GetRecord(seqid)
                if record != None:
                    fpout.write("%s"%(record))
            fpout.close()
            if not g_params['isQuiet']:
                print "Sequence file output to %s"%(outfile_seq)
        except IOError:
            print >> sys.stderr, "Failed to write to %s"%(outfile_seq)

        try:
            fpout = open(outfile_topo, "w")
            for seqid in idSet_pfam:
                record = hdl_topo.GetRecord(seqid)
                if record != None:
                    fpout.write("%s"%(record))
            fpout.close()
            if not g_params['isQuiet']:
                print "Topology file output to %s"%(outfile_topo)
        except IOError:
            print >> sys.stderr, "Failed to write to %s"%(outfile_topo)
    #}}}

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = ""
    pfamid2seqidFile = ""
    topodb = ""
    seqdb = ""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            print >> sys.stderr, "Error! Wrong argument:", argv[i]
            return 1
        elif argv[i] == "--":
            isNonOptionArg = True; i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-outpath", "--outpath"]:
                (outpath, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-pfamid2seqid", "--pfamid2seqid"] :
                (pfamid2seqidFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-seqdb", "--seqdb"] :
                (seqdb, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-topodb", "--topodb"] :
                (topodb, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True; i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            print >> sys.stderr, "Error! Wrong argument:", argv[i]
            return 1

    if outpath == "":
        print >> sys.stderr, "Error! outpath not set. Exit"
        print usage_short
        return 1
    elif not os.path.exists(outpath):
        try:
            subprocess.check_output(["mkdir", "-p", outpath])
        except subprocess.CalledProcessError, e:
            print e
            return 1

    if myfunc.checkfile(pfamid2seqidFile, "pfamid2seqidFile") != 0:
        return 1

    if myfunc.checkfile("%s0.db"%topodb, "topodb") != 0:
        return 1

    if myfunc.checkfile("%s0.db"%seqdb, "seqdb") != 0:
        return 1

    pfamid2seqidDict = myfunc.ReadFam2SeqidMap(pfamid2seqidFile)

    hdl_topo = myfunc.MyDB(topodb)
    if not hdl_topo.failure:
        idSet_topo = set(hdl_topo.indexedIDList)
    else:
        idSet_topo = set([])
        print >> sys.stderr, "Failed to open topology database %s"%(topodb)
        return 1

    hdl_seq = myfunc.MyDB(seqdb)
    if hdl_seq.failure:
        print >> sys.stderr, "Failed to open sequence database %s"%(seqdb)
        return 1

    GetTMProList_per_family(pfamid2seqidDict, idSet_topo, hdl_seq, hdl_topo,
            outpath)

    if not hdl_topo.failure:
        hdl_topo.close()
    if not hdl_seq.failure:
        hdl_seq.close()
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

#!/usr/bin/env python
# Description: This script just build seqid2pfamid table,
#              thus it will not be limited by the RAM size
# ChangeLog 2014-12-12 
#   for seqid2pfamid table, the pfamid is unique by uniquelist() now

import os
import sys
import myfunc
import subprocess
from collections import deque
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

usage_short="""
Usage: %s FILE [FILE ...] -outpath DIR
"""%(progname)

usage_ext="""
Description:
    Generate seqid2pfamid table file from the hmmscan result file
    two files will be output to out folder
    1. $outpath/$rootname.seqid2pfamid
    2. $outpath/$rootname.domainlistperseq

OPTIONS:
  -l LISTFILE   Set the listfile
  -outpath DIR  Output the result to DIR, if not set, the outpath will be the
                same as the input file 
  -e   FLOAT    Set e-value threshold, (default: 1e-3)
  -overwrite    Overwite the existing result file, (default: not)
  -h, --help    Print this help message and exit

Created 2014-12-11, updated 2014-12-12, Nanjiang Shu 
"""
usage_exp="""
Examples:
    %s uniref100_1.domtbl 
"""%(progname)

def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}

def ScanfHmmscanRecord(line):#{{{
    """
    get record from line
    """
    strs = line.split()
    if len(strs) < 23:
        return None
    try:
        pfamid = strs[0]
        target_length = int(strs[2])
        seqid = strs[3]
        query_length = int(strs[5])
        evalue = float(strs[6])
        env_begin = int (strs[19])
        env_end = int(strs[20])
        rd = (seqid, pfamid, evalue, query_length, target_length, env_begin, env_end)
        return rd
    except ValueError:
        return None
#}}}
def Build_seqid2pfamid(infile, g_outpath):#{{{
    outpath = ""
    dirname_infile = myfunc.my_dirname(infile)
    if g_outpath != "":
        outpath = g_outpath
    else:
        outpath = dirname_infile

    rootname = os.path.basename(os.path.splitext(infile)[0]);

    domainfile = "%s/%s.domainlistperseq"%(outpath, rootname)
    seqid2pfamidfile = "%s/%s.seqid2pfamid"%(outpath, rootname)
    if os.path.exists(domainfile) and os.path.exists(seqid2pfamidfile) and not g_params['isOverwrite']:
        print >> sys.stderr, "result file %s and %s exist. Ignore"%(domainfile, seqid2pfamidfile)
        return 1

    fpout_domain = myfunc.myopen(domainfile, None, "w", True)
    fpout_table = myfunc.myopen(seqid2pfamidfile, None, "w", True)

    evalue_threshold = g_params['evalue_threshold']
    hdl = myfunc.ReadLineByBlock(infile)
    queue = deque([])
    if hdl.failure:
        return 1
    lines = hdl.readlines()
    while lines != None:
        for line in lines:
            if not line or line[0] == "#":
                continue
            rd = ScanfHmmscanRecord(line)
#             print rd
            if rd == None:
                print >> sys.stderr, "%s: bad record. line=\"%s\""%(infile, line)
            else:
                evalue = rd[2]
                if evalue <= evalue_threshold:
                    queue.append(rd)
        # scan queue and output result
        idlist = [x[0] for x in queue]
#         print idlist
        # the top records is complete for one query if there are more than one unique seqids
        idlist_unique = myfunc.uniquelist(idlist)
        idlist_complete = idlist_unique[:-1] # remove the last item
        idset_complete = set(idlist_complete)
        if len(idset_complete) > 0:
            recordDict = {}
            cnt_used_rd = 0
            for i in xrange(len(queue)):
                seqid = queue[i][0]
                if not seqid in idset_complete:
                    continue
                if not seqid in recordDict:
                    recordDict[seqid] = []
                recordDict[seqid].append(queue[i])
                cnt_used_rd += 1
            #output
            for seqid in idlist_complete:
                try:
                    li = recordDict[seqid]
                except KeyError:
                    print "seqid=%s"%(seqid), "idlist_complete=", idlist_complete
                    raise
                famidlist = [x[1] for x in li]
                famidlist = myfunc.uniquelist(famidlist)
                fpout_table.write("%s %d"%(seqid, len(famidlist)))
                for pfamid in famidlist:
                    fpout_table.write(" %s"%(pfamid))
                fpout_table.write("\n")


                fpout_domain.write("%s %d"%(seqid, len(li)))
                for rd in li:
                    fpout_domain.write(" %s,%d,%d"%(rd[1], rd[5],rd[6]))
                fpout_domain.write("\n")
            # pop up queue
            for i in xrange(cnt_used_rd):
                queue.popleft()
        lines = hdl.readlines()

    if len(queue) > 0: # output the last item
        seqid = queue[0][0]
        li = queue
        famidlist = [x[1] for x in li]
        famidlist = myfunc.uniquelist(famidlist)
        fpout_table.write("%s %d"%(seqid, len(famidlist)))
        for pfamid in famidlist:
            fpout_table.write(" %s"%(pfamid))
        fpout_table.write("\n")

        fpout_domain.write("%s %d"%(seqid, len(li)))
        for rd in li:
            fpout_domain.write(" %s,%d,%d"%(rd[1], rd[5],rd[6]))
        fpout_domain.write("\n")

    hdl.close()
    myfunc.myclose(fpout_domain)
    myfunc.myclose(fpout_table)

    return 0
#}}}
def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = ""
    fileListFile = ""
    fileList = []

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            fileList.append(argv[i])
            isNonOptionArg = False
            i += 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-outpath", "--outpath"]:
                (outpath, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-e", "--e"]:
                (g_params['evalue_threshold'], i) = myfunc.my_getopt_float(argv, i)
            elif argv[i] in ["-l", "--l"] :
                (fileListFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True
                i += 1
            elif argv[i] in ["-overwrite", "--overwrite", "-force", "--force"]:
                g_params['isOverwrite'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            fileList.append(argv[i])
            i += 1


    if fileListFile != "":
        fileList += myfunc.ReadIDList(fileListFile)

    if len(fileList) < 1:
        print >> sys.stderr, "%s: no input file is set. exit"%(sys.argv[0])
        return 1

    if outpath != "" and not os.path.exists():
        cmd =  ["mkdir", "-p", outpath]
        subprocess.check_output(cmd)

    for i in xrange(len(fileList)):
        Build_seqid2pfamid(fileList[i], outpath)
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['evalue_threshold'] = 1e-3
    g_params['isOverwrite'] = False
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

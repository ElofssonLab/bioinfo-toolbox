#!/usr/bin/env python
# Description:
# Build pfamid2seqid table from a number of seqid2pfamid files
import os
import sys
import myfunc
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

usage_short="""
Usage: %s FILE [FILE ...] [-o OUTFILE]
"""%(progname)

usage_ext="""
Description:
    build pfamid2seqid table from a number of seqid2pfamid files

OPTIONS:
  -o OUTFILE    Output the result to OUTFILE
  -l LISTFILE   Set the listfile
  -nmax    INT  Set the maximum number of pfamids to handel at a time
                May read multiple times to solve memory overflow
                (default: 999999999, unlimited)
  -h, --help    Print this help message and exit

Created 2014-12-12, updated 2014-12-12, Nanjiang Shu
"""
usage_exp="""\
Examples:
    %s -l seqid2pfamid.list -o rst.pfamid2seqid
""" %(progname)

def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}

def IsUniqueSeq(infile, method):#{{{
    hdl = myfunc.ReadFastaByBlock(infile)
    if hdl.failure:
        return -1

    isunique = 1 #init value
    myset = set([])

    recordList = hdl.readseq()
    while recordList != None:
        for rd in recordList:
            if method == "id":
                if rd.seqid in myset:
                    isunique = 0
                    break
                myset.add(rd.seqid)
            elif method == "seq":
                if rd.seq in myset:
                    isunique = 0
                    break
                myset.add(rd.seq)
        recordList = hdl.readseq()

    hdl.close()
    return isunique
#}}}
def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile = ""
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
            elif argv[i] in ["-o", "--o", "-outfile"]:
                (outfile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-nmax", "--nmax"]:
                (g_params['nmax'], i) = myfunc.my_getopt_int(argv, i)
            elif argv[i] in ["-l", "--l"] :
                (fileListFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True
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

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    pfamidset_all = set([])
    pfamidset_output = set([])
    nmax = g_params['nmax']

    cnt_round = 0
    while 1:
        cnt_round += 1
        famid2seqidDict = {}
        for i in xrange(len(fileList)):
            hdl = myfunc.ReadLineByBlock(fileList[i])
            if hdl.failure:
                continue
            lines = hdl.readlines()
            while lines != None:
                for line in lines:
                    line = line.strip()
                    if not line or line[0] == "#":
                        continue
                    strs = line.split()
                    if len(strs) > 2:
                        seqid = strs[0]
                        pfamidlist = strs[2:]
                        for pfamid in pfamidlist:
                            if cnt_round == 1:
                                pfamidset_all.add(pfamid)
                            if pfamid in pfamidset_output:
                                continue
                            if not pfamid in famid2seqidDict:
                                if len(famid2seqidDict) < nmax:
                                    famid2seqidDict[pfamid] = []
                            if pfamid in famid2seqidDict:
                                famid2seqidDict[pfamid].append(seqid)
                    else:
                        msg="broken item in file %s: line \"%s\""
                        print >> sys.stderr, msg%(fileList[i], line)
                lines = hdl.readlines()
            hdl.close()

        for pfamid in famid2seqidDict:
            pfamidset_output.add(pfamid)
            seqidlist = famid2seqidDict[pfamid]
            seqidlist = myfunc.uniquelist(seqidlist)
            fpout.write("%s %d"%(pfamid, len(seqidlist)))
            for seqid in seqidlist:
                fpout.write(" %s"%(seqid))
            fpout.write("\n")

        if len(pfamidset_output) == len(pfamidset_all):
            break
        else:
            print " %d / %d "%(len(pfamidset_output), len(pfamidset_all))

    myfunc.myclose(fpout)
    if outfile != "":
        print "result output to %s"%(outfile)

    return 0
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['nmax'] = 999999999
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

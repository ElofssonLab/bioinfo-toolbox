#!/usr/bin/python
# Description:
import os
import sys
import myfunc
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))
GAP = "-"

usage_short="""
Usage: %s -msa seqMSAFile -topo topoFile [-o OUTFILE]
"""%(progname)

usage_ext="""
Description: 
    Output the topology alignment given sequence alignment.  All files should
    be in Fasta format.  Sequences are matched by seqID identified from the
    annotation line

OPTIONS:
  -o OUTFILE     Output the result to OUTFILE
  -topodb DBNAME Input topology is a formatted db
  -localaln      Treat the alignment as local alignment, in this case,
                 Lower case letters will be considered as non-aligned and 
                 in the matched topology replace it with X
  -ig, -ignore-badseq y|n
                 whether ignore bad sequence matching (default:yes)
                 if not enabled, the badly matched topology will be output as 
                 \"BADSEQ\"
  -m  0|1        method match seq to topo,(default: 1)
                 0: iterate sequence one by one
                 1: using segment matching
  -q             Quiet mode
  -h, --help     Print this help message and exit

Created 2013-05-15, updated 2013-05-15, Nanjiang Shu
"""
usage_exp="""
Examples:
    %s -msa test.mfa -topo test.topo -o test.topomsa.fa
"""%(progname)

def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}
def GetTopoDict(topofile):#{{{
    (idList, topoList) = myfunc.ReadFasta_without_annotation(topofile)
    topoDict = {}
    for i in xrange(len(idList)):
        topoDict[idList[i]] = topoList[i]
    return topoDict
#}}}
def MatchSeqToTopo(alignedseq, topo, method_match):
    gaplessseq = alignedseq.replace(GAP, "")
    if len(topo) != len(gaplessseq):
        return "BADSEQ"

    lengthAlignedSeq = len(alignedseq)
    lengthTopo = len(topo)

    if method_match == 0:
        j = 0 
        topoAliList = []
        for i in xrange(lengthAlignedSeq):
            if alignedseq[i] == GAP:
                topoAliList.append(GAP)
            else:
                topoAliList.append(topo[j])
                j += 1
        return "".join(topoAliList)

    elif method_match == 1:
        gapPosList = myfunc.GetSegPos(alignedseq, GAP)
#         print gapPosList
        num = len(gapPosList)
        if num < 1 :
            return topo
        else:
            topoAliList = []
            jTopo = 0
            if gapPosList[0][0] > 0:
                topoAliList.append(topo[jTopo:(jTopo+gapPosList[0][0])])
                jTopo += gapPosList[0][0]
            for i in xrange(num-1):
                topoAliList.append(GAP*(gapPosList[i][1] - gapPosList[i][0]))
                topoAliList.append(topo[jTopo:(jTopo + gapPosList[i+1][0] - gapPosList[i][1])])
                jTopo += (gapPosList[i+1][0] - gapPosList[i][1])
            topoAliList.append(GAP*(gapPosList[num-1][1] - gapPosList[num-1][0]))
            if gapPosList[num-1][1] < lengthAlignedSeq:
                topoAliList.append(topo[jTopo:])
            return "".join(topoAliList)

def MatchMSATopo_using_topofile(msafile, topofile, isIgnoreBadseq, #{{{
        method_match, outfile):
    topoDict = GetTopoDict(topofile)
    hdl = myfunc.ReadFastaByBlock(msafile)
    if hdl.failure:
        return 1

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)

    recordList = hdl.readseq()
    while recordList != None:
        for rd in recordList:
            try:
                topo = topoDict[rd.seqid]
            except KeyError:
                print >> sys.stderr, "topo not found for ID %s"%(rd.seqid)
                topo = ""
            matchedtopo = MatchSeqToTopo(rd.seq, topo, method_match)
            if not (matchedtopo == "BADSEQ" and isIgnoreBadseq):
                print >> fpout, ">%s"%(rd.description)
                print >> fpout, "%s"%(matchedtopo)
        recordList = hdl.readseq()

    myfunc.myclose(fpout)
    hdl.close()

    return 0
#}}}
def MatchMSATopo_using_topodb(msafile, topodb, isIgnoreBadseq, #{{{
        method_match, outfile):
    hdl_topo = myfunc.MyDB(topodb)
    if hdl_topo.failure:
        return 1

    hdl = myfunc.ReadFastaByBlock(msafile)
    if hdl.failure:
        return 1

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)

    recordList = hdl.readseq()
    while recordList != None:
        for rd in recordList:
            topowithanno = hdl_topo.GetRecord(rd.seqid)
            if topowithanno != None:
                (topoid, topoanno, topo) = myfunc.ExtractFromSeqWithAnno(topowithanno)
            else:
                print >> sys.stderr, "topo not found for ID %s"%(rd.seqid)
                topo = ""
            matchedtopo = MatchSeqToTopo(rd.seq, topo, method_match)
            if not (matchedtopo == "BADSEQ" and isIgnoreBadseq):
                print >> fpout, ">%s"%(rd.description)
                print >> fpout, "%s"%(matchedtopo)
        recordList = hdl.readseq()

    myfunc.myclose(fpout)
    hdl.close()
    hdl_topo.close()

    return 0
#}}}
def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile = ""
    msafile = ""
    topofile = ""
    topodb = ""
    isIgnoreBadseq = True
    method_match = 1


    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            print >> sys.stderr, "Error! Wrong argument:", argv[i]
            return 1
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
            elif argv[i] in ["-msa", "--msa"]:
                (msafile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-topo", "--topo"]:
                (topofile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-m", "--m"]:
                (method_match, i) = myfunc.my_getopt_int(argv, i)
                if method_match not in [0,1]:
                    print >> sys.stderr, "method_match %d not in [0,1]"%method_match
                    return 1
            elif argv[i] in ["-topodb", "--topodb"]:
                (topodb, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-ig", "--ig", "-ignore-badseq","--ignore-badseq"]:
                (tmpss, i) = myfunc.my_getopt_str(argv, i)
                if tmpss[0].lower() == "y":
                    isIgnoreBadseq = True
                else:
                    isIgnoreBadseq = False
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            print >> sys.stderr, "Error! Wrong argument:", argv[i]
            return 1

    if myfunc.checkfile(msafile) != 0:
        return 1

    if topodb != "":
        if not os.path.exists(topodb+'0.db'):
            print >> sys.stderr, "topodb %s does not exist"%(topodb)
            return 1
        else:
            return MatchMSATopo_using_topodb(msafile, topodb, isIgnoreBadseq,
                    method_match, outfile)
    elif topofile != "":
        if not os.path.exists(topofile):
            print >> sys.stderr, "topofile %s does not exist"%(topofile)
            return 1
        else:
            return MatchMSATopo_using_topofile(msafile, topofile, isIgnoreBadseq, method_match, outfile)
    else:
        print >> sys.stderr, "neither topofile nor topodb is set. exit"
        return 1

    return 0

#}}}


def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['method_match'] = 1
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

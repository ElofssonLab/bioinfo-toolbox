#!/usr/bin/python
# Description:
import os
import sys
sys.path.append("%s/wk/MPTopo/src"%(os.environ['DATADIR3']))
import myfunc
import libtopologycmp as lcmp

#ChangeLog 2014-10-14 
#   1. count also topologies that have been incorrectly predicted with inverted
#   orientation
#   2. Also separate single spanning and multi-spanning TM proteins
#ChangeLog 2015-08-26 
#   1. output wrongly predicted topology

progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

usage_short="""
Usage: %s -realtopo FILE -path_predtopo -basename STR
                [-mode tp|tps] [-o OUTFILE] [-rmsp]
"""%(progname)

usage_ext="""
Description:
    Benchmark the topology prediction by TOPCONS or TOPCONS-single

OPTIONS:
  -path_predtopo DIR   Set directory with prediction results
  -basename STR        Set the basename of prediction files
  -rmsp                Using predictions with signal peptide filtered
  -debug               Debug mode
  -restrictidlist FILE SeqIDs used in statistics
  -o OUTFILE           Output the result to OUTFILE
  -seqfile  FILE       Supply AA sequence file
  -owrong FILE         Output the wrongly predicted topology
  
  -h, --help           Print this help message and exit

Created 2014-10-10, updated 2015-08-26, Nanjiang Shu
"""
usage_exp="""
Examples:
    %s -realtopo opm.topo -path_predtopo pred_topcons -basename opm -rmsp
"""%(progname)

def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}

def CountIdenticalTopology(pred_topodict, real_topodict, agreement, TM_type,
        fpout_wrong, seqDict, pred_method_item):#{{{
    """
    return (cntIDT, cntINV)
    """
    numPredTopo = len(pred_topodict)

    cntIDT = 0
    cntINV = 0
    cntDIFF = 0
    for seqid in pred_topodict:
        predtopo = pred_topodict[seqid]
        try: 
            realtopo = real_topodict[seqid]
        except KeyError:
            print >> sys.stderr, "%s a nonTM protein predicted as TM protein"%(seqid)
            realtopo = "i"*len(predtopo)
            pass

        pred_NtermState = lcmp.GetNtermState(predtopo)
        real_NtermState = lcmp.GetNtermState(realtopo.replace('.', '-'))
        pred_posTM = myfunc.GetTMPosition(predtopo)
        real_posTM = myfunc.GetTMPosition(realtopo)
        pred_numTM = len(pred_posTM)
        real_numTM = len(real_posTM)

#         if g_params['isDEBUG'] and seqid == "3cx5I":
#             print "pred_NtermState = <%s>"%pred_NtermState
#             print "real_NtermState = <%s>"% real_NtermState
#             print "pred_posTM = ", pred_posTM
#             print "real_posTM = ", real_posTM



        if lcmp.IsIdenticalTopology(pred_NtermState, real_NtermState, pred_numTM,
                real_numTM, pred_posTM, real_posTM, predtopo, realtopo,
                g_params['min_TM_overlap']):
            cntIDT += 1
        else:
            if fpout_wrong != None:
                # output the wrongly predict topology
                fpout_wrong.write(">%s Number %d mtd_%s\n"%(seqid, cntDIFF+1, pred_method_item))
                try:
                    seq = seqDict[seqid]
                    fpout_wrong.write("%-10s %s\n"%("AASeq",seq))
                except KeyError:
                    seq = ""
                fpout_wrong.write("%-10s %s\n"%("RealTop",realtopo))
                fpout_wrong.write("%-10s %s\n"%("PredTop",predtopo))
                fpout_wrong.write("\n")

            if lcmp.IsInvertedTopology(pred_NtermState, real_NtermState, pred_numTM,
                    real_numTM, pred_posTM, real_posTM, predtopo, realtopo,
                    g_params['min_TM_overlap']):
                cntINV += 1
            if g_params['isDEBUG']:
                print >> sys.stderr, "%-7s(real %3s) nTM=%2d %s"%(seqid, agreement, real_numTM, realtopo)
                print >> sys.stderr, "%-7s(pred %3s) nTM=%2d %s"%(seqid, agreement, pred_numTM, predtopo)
                print >> sys.stderr
            cntDIFF += 1
    return (cntIDT, cntINV)
#}}}

def Benchmark(real_topodict, idSet_single, idSet_multi, TM_type, fpout, fpout_wrong, seqDict):#{{{
    if g_params['mode'] == "tps":
        itemlist = ["40","41","42","43","44", "All"]
    elif g_params['mode'] == "tp":
        itemlist = ["50","51","52","53","54","55", "All"]

    isRestrictIDList = g_params['isRestrictIDList']
    addname = ""
    if g_params['isRMSP']:
        addname = ".RMSP"

    numRealTopo = len(real_topodict)

    if isRestrictIDList:
        numRealTopo = len(g_params['restrictIDset'] & set(real_topodict.keys()))

    pred_topofile_list = []
    pred_topodict_list = []
# Step 1, read in predicted topology
    for item in itemlist:
        pred_topofile = ""
        if item.upper() ==  "ALL":
            if g_params['mode'] == "tps":
                pred_topofile = "%s/%s.topcons-single_topcons_single%s.topo"%(g_params['path_predtopo'], g_params['basename'], addname)
            elif g_params['mode'] == "tp":
                pred_topofile = "%s/%s.topcons.result_TOPCONS%s.topo"%(g_params['path_predtopo'], g_params['basename'], addname)

        else:
            if g_params['mode'] == "tps":
                pred_topofile = "%s/%s_topcons_single.m1.agree-%s%s.topo"%(g_params['path_predtopo'], g_params['basename'], item, addname)
            elif g_params['mode'] == "tp":
                pred_topofile = "%s/%s.topcons.result_TOPCONS.m1.agree-%s%s.topo"%(g_params['path_predtopo'], g_params['basename'], item, addname)

        (pred_idlist, pred_annolist, pred_topolist) = myfunc.ReadFasta(pred_topofile)
        if len(pred_idlist) <= 0:
            print >> sys.stderr, "Failed to read pred_topofile %s"%(pred_topofile)
        pred_topodict = {}
        for i in xrange(len(pred_idlist)):
            if ((not isRestrictIDList)  or pred_idlist[i] in g_params['restrictIDset']) :
                #if (TM_type == "All_Alpha" or (TM_type == "Single" and pred_idlist[i] in idSet_single) or (TM_type == "Multi" and pred_idlist[i] in idSet_multi)):
                pred_topodict[pred_idlist[i]] = pred_topolist[i]
        pred_topodict_list.append(pred_topodict)

# Step 2, calculate precision of the prediction
    #header line
    fpout.write("#%s\n"%(TM_type))
    fpout.write("#%2s %7s %8s %8s %8s %8s %8s %8s %8s\n"%("No", "Group", "nIDT", "nINV", "nPred", "PPV(%)", "NPV_INV", "NPV_Other", "nAllReal"))
    for i in xrange(len(itemlist)):
        item = itemlist[i]
        pred_topodict = pred_topodict_list[i]
        numPredTopo = len(pred_topodict)

        (numIDTtopo, numINVtopo) = CountIdenticalTopology(pred_topodict,
                real_topodict, item, TM_type, fpout_wrong, seqDict, item)

        ss = "%-3d %7s %8d %8d %8d %8.1f %8.1f %8.1f %8d"%(i, item, numIDTtopo, numINVtopo, numPredTopo,
                myfunc.FloatDivision(numIDTtopo, numPredTopo)*100.0,
                myfunc.FloatDivision(numINVtopo, numPredTopo)*100.0,
                myfunc.FloatDivision(numPredTopo-numIDTtopo-numINVtopo, numPredTopo)*100.0,
                numRealTopo)
        fpout.write("%s\n"%(ss))
    fpout.write("\n")
#}}}


def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = "./"
    outfile = ""
    real_topofile = ""
    seqfile = ""
    restrictIDListFile = ""
    outfile_wrong_predtopo = ""

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
            elif argv[i] in ["-owrong", "--owrong"]:
                (outfile_wrong_predtopo, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-realtopo", "--realtopo"]:
                (real_topofile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-seqfile", "--seqfile"]:
                (seqfile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-mode", "--mode"]:
                (g_params['mode'], i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-path_predtopo", "--path_predtopo"]:
                (g_params['path_predtopo'], i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-basename", "--basename"]:
                (g_params['basename'], i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-restrictidlist", "--restrictidlist"]:
                (restrictIDListFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True
                i += 1
            elif argv[i] in ["-rmsp", "--rmsp"]:
                g_params['isRMSP'] = True
                i += 1
            elif argv[i] in ["-debug", "--debug"]:
                g_params['isDEBUG'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            print >> sys.stderr, "Error! Wrong argument:", argv[i]
            return 1
            i += 1

    if myfunc.checkfile(g_params['path_predtopo'], "path_predtopo") != 0:
        return 1
    if g_params['basename'] == "":
        print >> sys.stderr, "%s: basename not set. exit"%(argv[0])
        return 1
    if myfunc.checkfile(real_topofile, "real_topofile") != 0:
        return 1

    if restrictIDListFile != "":
        g_params['restrictIDset'] = set(myfunc.ReadIDList(restrictIDListFile))
        g_params['isRestrictIDList'] = True


    if g_params['mode'] == "":
        if g_params['path_predtopo'].find("topcons_single") >= 0:
            g_params['mode'] = "tps"
        elif  g_params['path_predtopo'].find("topcons") >= 0: 
            g_params['mode'] = "tp"
        else:
            print >> sys.stderr, "mode not set, and can not be recognized from path_predtopo=%s"%(path_predtopo)
            return 1

    if not g_params['mode'] in ["tp", "tps"]:
        print >> sys.stderr, "Unrecognized mode = %s"%(g_params['mode'])
        return 1

    (real_idlist, real_annolist, real_topolist) = myfunc.ReadFasta(real_topofile)
    seqDict = {}
    if seqfile != "" and os.path.exists(seqfile):
        (seq_idlist, seq_annolist, seqlist) = myfunc.ReadFasta(seqfile)
        for i in xrange(len(seq_idlist)):
            seqDict[seq_idlist[i]] = seqlist[i]

    if len(real_idlist) <= 0:
        print >> sys.stderr, "Failed to read real_topofile %s"%(real_topofile)
        return 1

    real_topodict = {}
    for i in xrange(len(real_idlist)):
        real_topodict[real_idlist[i]] = real_topolist[i]


    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    fpout_wrong = myfunc.myopen(outfile_wrong_predtopo, None, "w", False)


    idSet_single = set([])
    idSet_multi = set([])
    for seqid in real_topodict:
        topo = real_topodict[seqid]
        numTM = myfunc.CountTM(topo)
        if numTM == 1:
            idSet_single.add(seqid)
        elif numTM > 1:
            idSet_multi.add(seqid)

#     print "len(real_topodict)", len(real_topodict)
#     print "len(idSet_single)", len(idSet_single)
#     print "len(idSet_multi)", len(idSet_multi)

    #for TM_type in ["All_Alpha", "Single", "Multi"]:
    for TM_type in ["All_Alpha"]:
        if TM_type == "All_Alpha":
            sub_real_topodict = real_topodict
        else:
            sub_real_topodict = {}
            for seqid in real_topodict:
                topo = real_topodict[seqid]
                numTM = myfunc.CountTM(topo)
                if TM_type == "Single" and numTM == 1:
                    sub_real_topodict[seqid] = topo
                elif TM_type == "Multi" and numTM > 1:
                    sub_real_topodict[seqid] = topo
        Benchmark(sub_real_topodict, idSet_single, idSet_multi, TM_type, fpout,
                fpout_wrong, seqDict)


    myfunc.myclose(fpout)


#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['path_predtopo'] = ""
    g_params['basename'] = ""
    g_params['mode'] = ""
    g_params['isRMSP'] = False
    g_params['isDEBUG'] = False
    g_params['restrictIDset'] = set([])
    g_params['isRestrictIDList'] = False
    g_params['min_TM_overlap'] = 5
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))


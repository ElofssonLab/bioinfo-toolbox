#!/usr/bin/python
# Description:
import os
import sys
import myfunc
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

usage_short="""
Usage: %s topoalnfile -localali FILE [-o OUTFILE]
"""%(progname)

usage_ext="""
Description:

OPTIONS:
  -o OUTFILE    Output the result to OUTFILE
  -q            Quiet mode
  -h, --help    Print this help message and exit

Created 2013-06-17, updated 2013-06-17, Nanjiang Shu 
"""
usage_exp="""
Examples:
"""

def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}
def GetUnAlignedString(alnseq1, alnseq2):#{{{
    # for local alignment
    sslist = []
    if len(alnseq1) != len(alnseq2):
        print >> sys.stderr, "alnseq1 and alnseq2 length not matched"
        return ""
    else:
        for j in xrange(len(alnseq1)):
            if alnseq1[j].islower() or alnseq2[j].islower():
                sslist.append('0')
            else:
                sslist.append('1')
        return "".join(sslist)
#}}}
def GetUnglianedTermStatus(posTM_term1, posTM_term2):#{{{
    if len(posTM_term1) >0 and len(posTM_term2) > 0:
        s_term = 2
    elif len(posTM_term1) >0 or len(posTM_term2) > 0:
        s_term = 1
    else:
        s_term = 0
    return s_term

#}}}
def GetUnglianedTermStatus1(posTM_term1, posTM_term2):#{{{
    if len(posTM_term1) == 0:
        if len(posTM_term2) == 0:
            s1_term = 0
        else:
            s1_term = 1
    elif len(posTM_term1) == 1:
        if len(posTM_term2) == 0:
            s1_term = 1
        else:
            if posTM_term1[0][1]-posTM_term1[0][0] < 5:
                s1_term = 1.5
            else:
                if len(posTM_term2) == 1 and posTM_term2[0][1] - posTM_term2[0][0] < 5:
                    s1_term = 1.5
                else:
                    s1_term = 2
    else:
        if len(posTM_term2) == 0:
            s1_term = 1
        elif len(posTM_term2) == 1 and posTM_term2[0][1] - posTM_term2[0][0] < 5:
            s1_term = 1.5
        else:
            s1_term = 2
    return s1_term
#}}}
def AnaLocalTopoAln(idList, topoList, localseqpairDict, fpout, fpout1):
    # fpout: write result for those one end with no TM region
# fpout1: write result for those one end with one TM region but less then 5
    # residues of TM region
    numseq = len(idList)
    numpair = numseq/2
    for i in xrange(numpair):
        id1 = idList[2*i]
        id2 = idList[2*i+1]
        topo1 = topoList[2*i]
        topo2 = topoList[2*i+1]
        lengthAln = len(topo1)
        try:
            unaligned_str = localseqpairDict[(id1,id2)][2]
        except (KeyError, IndexError):
            print >> sys.stderr, "no local alignment found for %s %s"%(id1,id2)
            continue
        alignedPosList = myfunc.GetSegPos(unaligned_str, "1")
        if len(alignedPosList) != 1:
            print >> sys.stderr, "aligned region not equal 1 for %s %s"%(id1,id2)
            continue
        else:
            alignedPos = alignedPosList[0]
        if alignedPos[0] == 0 and alignedPos[1] == lengthAln:
            print "%s %s: Full aligned"%(id1, id2)
        else:
            alignedPos = alignedPosList[0]
            topo_Nterm1 = topo1[:alignedPos[0]]
            topo_Cterm1 = topo1[alignedPos[1]:]
            topo_Nterm2 = topo2[:alignedPos[0]]
            topo_Cterm2 = topo2[alignedPos[1]:]
            posTM_Nterm1 = myfunc.GetTMPosition(topo_Nterm1)
            posTM_Cterm1 = myfunc.GetTMPosition(topo_Cterm1)
            posTM_Nterm2 = myfunc.GetTMPosition(topo_Nterm2)
            posTM_Cterm2 = myfunc.GetTMPosition(topo_Cterm2)

            s_Nterm = GetUnglianedTermStatus(posTM_Nterm1, posTM_Nterm2)
            s_Cterm = GetUnglianedTermStatus(posTM_Cterm1, posTM_Cterm2)
            # s1_Nterm, s2_Nterm  is used to record the status of those
            # unaligned terminals with one has a splitted TM helices and
            # another has >= 1 TM helix
            s1_Nterm = GetUnglianedTermStatus1(posTM_Nterm1, posTM_Nterm2)
            s1_Cterm = GetUnglianedTermStatus1(posTM_Cterm1, posTM_Cterm2)


            if s_Nterm < 2 and s_Cterm < 2 and (s_Nterm + s_Cterm) > 0:
                if len(posTM_Nterm1) > 0:
                    num_res_unaligned_Nterm = len(topo_Nterm2.replace("-",""))
                    numTM_unaligned_Nterm = len(posTM_Nterm1)
                    num_res_to_TM_Nterm = len(topo_Nterm1) - posTM_Nterm1[len(posTM_Nterm1)-1][1]
                elif len(posTM_Nterm2) > 0:
                    num_res_unaligned_Nterm = len(topo_Nterm1.replace("-",""))
                    numTM_unaligned_Nterm = len(posTM_Nterm2)
                    num_res_to_TM_Nterm = len(topo_Nterm2) - posTM_Nterm2[len(posTM_Nterm2)-1][1]
                else:
                    num_res_unaligned_Nterm = 0
                    numTM_unaligned_Nterm = 0
                    num_res_to_TM_Nterm = 0

                if len(posTM_Cterm1) > 0:
                    num_res_unaligned_Cterm = len(topo_Cterm2.replace("-",""))
                    numTM_unaligned_Cterm = len(posTM_Cterm1)
                    num_res_to_TM_Cterm =  posTM_Cterm1[0][0]
                elif len(posTM_Cterm2) > 0:
                    num_res_unaligned_Cterm = len(topo_Cterm1.replace("-",""))
                    numTM_unaligned_Cterm = len(posTM_Cterm2)
                    num_res_to_TM_Cterm = posTM_Cterm2[0][0]
                else:
                    num_res_unaligned_Cterm = 0
                    numTM_unaligned_Cterm = 0
                    num_res_to_TM_Cterm = 0


                ss = "%s %s %4d %4d %4d           %4d %4d %4d"
                print >> fpout, ss%(id1, id2, 
                        num_res_unaligned_Nterm,
                        num_res_to_TM_Nterm,
                        numTM_unaligned_Nterm,
                        num_res_unaligned_Cterm,
                        num_res_to_TM_Cterm,
                        numTM_unaligned_Cterm
                        )
            if ((s1_Nterm == 1.5 or s1_Cterm == 1.5) 
                    and s1_Nterm < 2 
                    and s1_Cterm < 2):

                num_res_unaligned_Nterm = -1
                numRes_PartHelix_Nterm =  -1
                numTM_unaligned_Nterm = -1
                num_res_to_TM_Nterm = -1
                num_res_unaligned_Cterm = -1
                numRes_PartHelix_Cterm =  -1
                numTM_unaligned_Cterm = -1
                num_res_to_TM_Cterm = -1
                if s1_Nterm == 1.5:
                    if len(posTM_Nterm1) == 1 and posTM_Nterm1[0][1] - posTM_Nterm1[0][0] < 5:
                        num_res_unaligned_Nterm = len(topo_Nterm1.replace("-",""))
                        numRes_PartHelix_Nterm =  posTM_Nterm1[0][1] - posTM_Nterm1[0][0]
                        numTM_unaligned_Nterm = len(posTM_Nterm2)
                        num_res_to_TM_Nterm = len(topo_Nterm2) - posTM_Nterm2[len(posTM_Nterm2)-1][1]
                    elif len(posTM_Nterm2) == 1 and posTM_Nterm2[0][1] - posTM_Nterm2[0][0] < 5:
                        num_res_unaligned_Nterm = len(topo_Nterm2.replace("-",""))
                        numRes_PartHelix_Nterm =  posTM_Nterm2[0][1] - posTM_Nterm2[0][0]
                        numTM_unaligned_Nterm = len(posTM_Nterm1)
                        num_res_to_TM_Nterm = len(topo_Nterm1) - posTM_Nterm1[len(posTM_Nterm1)-1][1]

                if s1_Cterm == 1.5:
                    if len(posTM_Cterm1) == 1 and posTM_Cterm1[0][1] - posTM_Cterm1[0][0] < 5:
                        num_res_unaligned_Cterm = len(topo_Cterm1.replace("-",""))
                        numRes_PartHelix_Cterm =  posTM_Cterm1[0][1] - posTM_Cterm1[0][0]
                        numTM_unaligned_Cterm = len(posTM_Cterm2)
                        num_res_to_TM_Cterm = len(topo_Cterm2) - posTM_Cterm2[len(posTM_Cterm2)-1][1]
                    elif len(posTM_Cterm2) == 1 and posTM_Cterm2[0][1] - posTM_Cterm2[0][0] < 5:
                        num_res_unaligned_Cterm = len(topo_Cterm2.replace("-",""))
                        numRes_PartHelix_Cterm =  posTM_Cterm2[0][1] - posTM_Cterm2[0][0]
                        numTM_unaligned_Cterm = len(posTM_Cterm1)
                        num_res_to_TM_Cterm = len(topo_Cterm1) - posTM_Cterm1[len(posTM_Cterm1)-1][1]


                ss = "%s %s %4d %4d %4d %4d          %4d %4d %4d %4d"
                print >> fpout1, ss%(id1, id2, 
                        num_res_unaligned_Nterm,
                        numRes_PartHelix_Nterm,
                        num_res_to_TM_Nterm,
                        numTM_unaligned_Nterm,
                        num_res_unaligned_Cterm,
                        numRes_PartHelix_Cterm,
                        num_res_to_TM_Cterm,
                        numTM_unaligned_Cterm
                        )

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile = ""
    topoalnfile = ""
    localalifile = ""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            topoalnfile = argv[i]
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
            elif argv[i] in ["-localali", "--localali"]:
                (localalifile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            topoalnfile = argv[i]
            i += 1

    if myfunc.checkfile(topoalnfile) != 0:
        return  1
    if myfunc.checkfile(localalifile) != 0:
        return  1

    (idList, topoList) = myfunc.ReadFasta_without_annotation(topoalnfile)
    (idListLocal, seqListLocal) = myfunc.ReadFasta_without_annotation(localalifile)
    numseqLocal = len(idListLocal)
    numpairLocal = numseqLocal/2
    localseqpairDict = {}
    for i in xrange(numpairLocal):
        id1 = idListLocal[2*i]
        id2 = idListLocal[2*i+1]
        unaligned_str = GetUnAlignedString(seqListLocal[2*i], seqListLocal[2*i+1])
        if unaligned_str != "":
            localseqpairDict[(id1, id2)] = [seqListLocal[2*i],
                    seqListLocal[2*i+1], unaligned_str]
    del idListLocal,seqListLocal

    if outfile != "":
        outfile1 = outfile + ".1"
    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    fpout1 = myfunc.myopen(outfile1, sys.stdout, "w", False)

    AnaLocalTopoAln(idList, topoList, localseqpairDict, fpout, fpout1)

    myfunc.myclose(fpout)
    myfunc.myclose(fpout1)


#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

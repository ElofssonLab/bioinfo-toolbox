#!/usr/bin/env python

# Draw pairwise topology comparison by given a number of pairwise sequence
# alignment, dg values profiles are added
import os
import sys
rundir = os.path.dirname(sys.argv[0])
binpath = rundir
dgscanprog = rundir + '/../../../program/dgpred_standalone/myscanDG.pl'
topofile = rundir + '/../pfamAna/pairwise/all/pred_topcons_single_method4/pfamfullseq.selTM_uniq.topcons-single.fa'
sys.path.append(binpath)
import myfunc

usage="""
usage:    drawpairwisetopo.py seq_pairaln_file  [-outpath outpath]

Description: 
    Draw pairwise topology comparison by given a number of pairwise
    sequence alignment, dg values profiles are added, 
Options:
  -outpath DIR    Output the result to DIR, (default: ./)   
  -topofile FILE  Supply topology file in FASTA format. (default: 
                  pred_topcons_single_method4/pfamfullseq.selTM_uniq.topcons-single.fa)
  -cmpclass STR   Draw only selected cmpclass (OK, IDT, INV, DIFF)
  -q              Quiet mode
  -h, --help      Print this help message and exit

Created 2012-03-22, updated 2013-04-19, Nanjiang Shu  

Examples:
    drawpairwisetopo.py allpfam.pairaln.fa 
    drawpairwisetopo.py allpfam.pairaln.fa  -topofile topo.fa -outpath outdir
"""

def DrawPairwiseTopo(pairtopoAlnFile, aaSeqDict, pairCmpclassDict,
        outpath):
    (idList, annoList, seqList) = myfunc.ReadFasta(pairtopoAlnFile);
    numSeq = len(idList);
    numPair = numSeq/2;

    print "numSeq = ", numSeq
    print "numPair = ", numPair

    for i in range(numPair):
        id1 = idList[2*i]
        id2 = idList[2*i+1]
        if len(seqList[2*i]) != len(seqList[2*i+1]):
            print "Error for %s - %s "%(idList[2*i], idList[2*i+1])
            continue;
        basename = "%s-%s"%(id1, id2)

        isSatisfied = True
#         if basename in pairCmpclassDict:
#             if g_params['cmpclassList'] != []:
#                 if (not pairCmpclassDict[basename] in
#                         g_params['cmpclassList']):
#                     isSatisfied = False
#             elif pairCmpclassDict[basename] == 'OK':
#                 isSatisfied = False
        if isSatisfied:
            outPairAlnFile = outpath + os.sep + "%s.topoaln.fa"%(basename)
            fpout = open (outPairAlnFile, 'w')
            print >> fpout, ">%s" %annoList[2*i];
            print >> fpout, "%s" %seqList[2*i];
            print >> fpout, ">%s" %annoList[2*i+1];
            print >> fpout, "%s" %seqList[2*i+1];
            fpout.close()
            outAASeqFile = outpath + os.sep + "%s.fa" % (basename)
            fpout = open (outAASeqFile, "w")
            if id1 in aaSeqDict:
                print >> fpout, ">%s"%id1
                print >> fpout, "%s"%aaSeqDict[id1]
            if id2 in aaSeqDict:
                print >> fpout, ">%s"%id2
                print >> fpout, "%s"%aaSeqDict[id2]
            fpout.close()
  
# Output dgscan file
            dgpfile = outpath + os.sep + basename + '.dgscan'
            cmd = "%s %s -lmin 21 -lmax 21 -o %s" %(dgscanprog, outAASeqFile, dgpfile)
            os.system(cmd)

            outpngfile = outpath + os.sep + "%s.topoaln.png" %basename
            outShrinkedFile = (outpath + os.sep +
                    "%s.topoaln.shrinked.png"%basename)
            thumb_outShrinkedFile = (outpath + os.sep + 'thumb.' +
                    "%s.topoaln.shrinked.png"%basename)
            outNonShrinkedFile = (outpath + os.sep +
                    "%s.topoaln.nonshrinked.png"%basename)
            thumb_outNonShrinkedFile = (outpath + os.sep + 'thumb.' +
                    "%s.topoaln.nonshrinked.png"%basename)
            os.system("python %s/drawMSATopo.py %s -pfm no -shrink yes -method mat" %
                    (binpath, outPairAlnFile))
            os.system("mv %s %s"%(outpngfile, outShrinkedFile))
            os.system("python %s/drawMSATopo.py %s -pfm no -shrink no -pdg yes -method yes -dgpfile %s" % (binpath, outPairAlnFile, dgpfile))
            os.system("mv %s %s"%(outpngfile, outNonShrinkedFile))
            os.system("convert -thumbnail 200 %s %s" %(outShrinkedFile, thumb_outShrinkedFile))
            os.system("convert -thumbnail 200 %s %s" %(outNonShrinkedFile, thumb_outNonShrinkedFile))
            os.system("rm -f %s %s"%(outAASeqFile, dgpfile))

def ReadPaircmpCmpclass(infile):
    try:
        fpin = open(infile, "r")
        lines = fpin.readlines()
        fpin.close()
        pairCmpclassDict = {}
        for line in lines:
            strs = line.split()
            if strs[0] == "PairwiseComparison:":
                id1 = strs[1]
                id2 = strs[2]
                cmpclass = strs[4]
                pairid = id1 + '-' + id2
                pairCmpclassDict[pairid] = cmpclass
        return pairCmpclassDict
    except IOError:
        print >> sys.stderr, "Failed to read file %s"%infile
        return {}

def PrintHelp():
    print usage

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1;

    outpath="./";
    isQuiet=False;
    diffseqidtgroup = "0"

    pairseqAlnFile = ''
    cmpclassList = []
    topofile = ""

    i = 1;
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            pairseqAlnFile = argv[i];
            isNonOptionArg=False;
            i += 1;
        elif argv[i] == "--":
            isNonOptionArg=True;
            i += 1;
        elif argv[i][0] == "-":
            if argv[i] ==  "-h" or  sys.argv[i] == "--help":
                PrintHelp();
                return 1
            elif (argv[i] in ["-outpath", "--outpath"]):
                outpath = argv[i+1];
                i += 2;
            elif (argv[i] in ["-topofile", "--topofile"]):
                topofile = argv[i+1];
                i += 2;
            elif (argv[i] in ["-cmpclass", "--cmpclass"]):
                cmpclassList.append(argv[i+1]);
                i += 2;
            elif sys.argv[i] == "-q":
                isQuiet=True;
                i += 1;
            else:
                print >> sys.stderr, "Error! Wrong argument:", sys.argv[i];
                return -1;
        else:
            pairseqAlnFile = argv[i];
            i += 1
    g_params['outpath'] = outpath
    g_params['cmpclassList'] = cmpclassList

    if pairseqAlnFile == "":
        print >> sys.stderr, "pairseqAlnFile not set. Exit.";
        return 1;
    if not os.path.exists(pairseqAlnFile):
        print >> sys.stderr, "pairseqAlnFile %s does not exists. Exit."%pairseqAlnFile;
        return 1

    rootname = os.path.basename(os.path.splitext(pairseqAlnFile)[0]);
# Read in aaSeqDict
    print "Read in aaSeqDict"
    os.system("mkdir -p %s"%outpath);
    (idList, seqList) = myfunc.ReadFasta_without_annotation(pairseqAlnFile)

    # create seqdbfile
    seqdbfile = outpath + os.sep + rootname + ".seqdb.fa"
    fo = open(seqdbfile,"w")
    for i in range(len(idList)):
        print >> fo, ">%s"%(idList[i])
        print >> fo, "%s"%(seqList[i].replace('-', ''))
    fo.close()
    cmd = "%s/indexfasta.py %s"
    os.system(cmd%(binpath, seqdbfile))
    seqdbname = outpath + os.sep + rootname + ".seqdb"

    aaSeqDict = {}
    for i in xrange(len(idList)):
        aaSeqDict[idList[i]] = seqList[i].replace('-','')
# # Output uniqid included in pairtopoAlnFile
#     uniqid_set = set(idList)
# # output uniqid seqfile
#     uniqidAASeqFile = outpath + os.sep + rootname + '.uniqid.aaseq.fa'
#     print "Output uniqid seqfile to %s"%uniqidAASeqFile
#     fpout = open(uniqidAASeqFile, "w")
#     for idd in uniqid_set:
#         fpout.write(">%s\n"%idd)
#         fpout.write("%s\n"%aaSeqDict[idd])
#     fpout.close()
# # Output dgscan file
#     dgpfile = outpath + os.sep + rootname + '.uniqid.dgscan'
#     print "Output dgscan file to %s"%dgpfile
#     cmd = "%s %s -lmin 21 -lmax 21 -o %s" %(dgscanprog, uniqidAASeqFile,
#             dgpfile)
#     os.system(cmd)
# Output pairwise topology comparison
    if not os.path.exists(topofile):
        print >> sys.stderr, "topofile %s not exist. exit." % topofile
        return 1
    cmd = "%s/seqpairaln_to_topopaircmp.sh %s -outpath %s -topofile %s -seqdb %s" % (
            binpath, pairseqAlnFile, outpath, topofile, seqdbname)
    print "Output paircmp file to %s"%outpath
    os.system(cmd)
    paircmpFile = outpath + os.sep + rootname + '.paircmp'
    pairtopoAlnFile = outpath + os.sep + rootname + '.topoaln.fa'
# Read in paircmp file 
    pairCmpclassDict = ReadPaircmpCmpclass(paircmpFile)

    print "Draw pairwise topology comparison ..."
    DrawPairwiseTopo(pairtopoAlnFile, aaSeqDict, pairCmpclassDict, outpath)

    return 0;

#}}}
def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['outpath'] = ""
    return g_params
#}}}

if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))


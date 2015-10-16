#!/usr/bin/env python
# Description:
# select pairwise alignment by given the id of pairs
# ChangeLog 2014-03-10 
#       the name of the outfile is changed so that the output file name id1_id2
#       will be in the same order as in the input pairlist
#           outfile = "%s%s%s%s"%(g_params['outpath'], os.sep,
#               "%s_%s"%(seqid1, seqid2), g_params['ext']) 
#       to
#           outfile = "%s%s%s%s"%(g_params['outpath'], os.sep,
#               "%s_%s"%(p_rd1[0], p_rd2[0]), g_params['ext']) 
#   
import os
import sys
import myfunc
import subprocess
usage = """
usage:  selectPairaln.py -pairaln pairalnFile id1 id2 [-o OUTFILE]
Description:
# select pairwise alignment by given the id of pairs
Options:
  -l LISTFILE     Set list file, one record per line
  -split          Whether split the pairwise alignment to individual files
                  named as $outpath/id1_id2$ext
  -outpath DIR    Output path, when split is enabled, (default: ./)
  -ext     STR    File extension, (default: .pairaln.fa)
  -h, --help      Print this help message and exit

Created 2012-11-30, updated 2014-03-10, Nanjiang Shu  
"""
BLOCK_SIZE=100000;

def PrintHelp():
    print usage
def SelectRecord(recordList, pairlistSet, fpout_global):#{{{
    numRecord = len(recordList)
    numPair = numRecord/2
    cntSelectedPair = 0
    for i in xrange(numPair):
        rd1 = recordList[2*i]
        rd2 = recordList[2*i+1]
        seqid1 = rd1[0]
        seqid2 = rd2[0]
        isFound = False
        if (seqid1, seqid2) in pairlistSet:
            p_rd1 = rd1
            p_rd2 = rd2
            isFound = True
        elif (seqid2, seqid1) in pairlistSet:
            p_rd1 = rd2
            p_rd2 = rd1
            isFound = True
        if isFound:
            if g_params['isSplit']:
                outfile = "%s%s%s%s"%(g_params['outpath'], os.sep,
                        "%s_%s"%(p_rd1[0], p_rd2[0]), g_params['ext']) #changed 2014-03-10, so that the output file name id1_id2 will be in the same order as in the input pairlist
                try:
                    fpout = open(outfile, "w")
                except IOError:
                    print >> sys.stderr, "Failed to write to %s"%(outfile)
                    fpout = None
            else:
                fpout = fpout_global

            if fpout != None:
                fpout.write(">%s\n" %(p_rd1[1]))
                fpout.write("%s\n"  %(p_rd1[2]))
                fpout.write(">%s\n" %(p_rd2[1]))
                fpout.write("%s\n"  %(p_rd2[2]))
            if fpout != None and g_params['isSplit']:
                print "%s output"%(outfile)
                fpout.close()
            cntSelectedPair += 1
    return cntSelectedPair
#}}}
def SelectPairaln(pairalnFile, pairlistSet, fpout):#{{{
    try:
        numPairToSelect = len(pairlistSet)
        cntSelectedPair = 0
        fpin = open(pairalnFile, "r")
        unprocessedBuffer="";
        isEOFreached = False;
        processedTopoIDSet = set([]);
        remainedRd = None
        while 1:
            buff = fpin.read(BLOCK_SIZE);
            if len(buff) < BLOCK_SIZE:
                isEOFreached=True;
            buff = unprocessedBuffer + buff;
            recordList = [];
            if remainedRd != None:
                recordList.append(remainedRd)
            unprocessedBuffer = myfunc.ReadFastaFromBuffer(
                    buff, recordList, isEOFreached);
            numRecord = len(recordList)
            if numRecord > 0: 
                cntSelectedPair += SelectRecord(recordList, pairlistSet,
                        fpout)
                if cntSelectedPair >= numPairToSelect:
                    break
                if numRecord%2 == 1:
                    remainedRd = recordList[numRecord-1]
                else:
                    remainedRd = None

            if isEOFreached == True:
                break;
        fpin.close();
    except IOError:
        print >> sys.stderr, "Failed to open file %s for read"%(pairalnFile);
        return 1
#}}}
def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile = ""
    pairalnFile = ""
    pair = []
    listfile = ""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            pair.append(argv[i])
            isNonOptionArg = False
            i += 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-o", "--o"]:
                outfile, i = myfunc.my_getopt_str(argv,i)
            elif argv[i] in ["-pairaln", "--pairaln"]:
                pairalnFile, i = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-outpath", "--outpath"]:
                g_params['outpath'], i = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-ext", "--ext"]:
                g_params['ext'], i = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-l", "--l"]:
                listfile, i = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q"]:
                g_params['isQuiet'] = True; i += 1
            elif argv[i] in ["-split", "--split"]:
                g_params['isSplit'] = True; i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            pair.append(argv[i])
            i += 1

    if myfunc.checkfile(pairalnFile, "pairalnFile") != 0:
        return 1

    if g_params['isSplit']:
        if g_params['outpath'] == "":
            print >> sys.stderr, "Error! outpath string is empty when 'split'"\
                    " is enabled. exit"
            return 1
        elif not os.path.exists(g_params['outpath']):
            cmd = ["mkdir" , "-p", g_params['outpath']]
            subprocess.check_call(cmd)

    pairlist = []
    if len(pair) >= 2:
        pairlist.append((pair[0],pair[1]))
    if listfile != "":
        pairlist += myfunc.ReadPairList(listfile)

    pairlistSet = set([])
    for pair in pairlist:
        pairlistSet.add(pair)
    del pairlist

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    SelectPairaln(pairalnFile, pairlistSet, fpout)
    myfunc.myclose(fpout)
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['keywordList'] = []
    g_params['isCaseSensitive'] = False
    g_params['isSplit'] = False
    g_params['outpath'] = "."
    g_params['ext'] = ".pairaln.fa"
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

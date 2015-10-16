#!/usr/bin/env python
# Description:
import os
import sys
import myfunc
import tempfile
usage = """
usage: run_calEvoDist.py pairalnFile [-o OUTFILE]
Description:

Options:
  -q              Quiet mode
  -h, --help      Print this help message and exit

Created 2012-08-24, updated 2012-08-24, Nanjiang Shu
"""

def PrintHelp():
    print usage



def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    pairalnFile = ""
    outfile = ""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            pairalnFile = argv[i]
            isNonOptionArg = False
            i += 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-outfile", "--outfile"]:
                outfile = argv[i+1]
                i += 2
            elif argv[i] in ["-q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            pairalnFile = argv[i]
            i += 1
    if pairalnFile == "":
        print >> sys.stderr, "pairalnFile not set"
        return 1
    elif not os.path.exists(pairalnFile):
        print >> sys.stderr, "pairalnFile %s does not exist" %pairalnFile
        return 1
    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    (idList, seqList) = myfunc.ReadFasta_without_annotation(pairalnFile)
    numSeq = len(idList)
    numPair = numSeq / 2
    for i in xrange(numPair):
        id1 = idList[2*i]
        id2 = idList[2*i+1]
        seq1 = seqList[2*i]
        seq2 = seqList[2*i+1]
        if len(seq1) != len(seq2):
            print >> sys.stderr, "Bad alignment, seq length conflicts, %d (%s) = %d (%s)" %(
                    len(seq1), id1, len(seq2), id2)
            continue
        tmpfile = tempfile.mktemp()
        fpout = open(tmpfile, "w")
        fpout.write(">%s\n"%(id1))
        fpout.write("%s\n"%(seq1))
        fpout.write(">%s\n"%(id2))
        fpout.write("%s\n"%(seq2))
        fpout.close()
        cmd = "%s/calEvoDist.sh -f 1 %s" %(binpath, tmpfile)
        os.system(cmd)
        os.remove(tmpfile)


    myfunc.myclose(fpout)

#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    rundir = os.path.dirname(sys.argv[0])
    binpath = rundir
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

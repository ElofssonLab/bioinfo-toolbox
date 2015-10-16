#!/usr/bin/python
# Description:
import os
import sys
sys.path.append("%s/wk/MPTopo/src"%(os.environ['DATADIR3']))
import myfunc
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

usage_short="""
Usage: %s seqid_to_taxid_file -gram GRAM-P-N-FILE [-o OUTFILE]
"""%(progname)

usage_ext="""
Description:
    format of GRAM-P-N-FILE 
        TAXID gram+
        TAXID euk
        TAXID gram-
    format of seqid_to_taxid_file
        seqid taxid

OPTIONS:
  -o OUTFILE    Output the result to OUTFILE
  -q            Quiet mode
  -h, --help    Print this help message and exit

Created 2013-09-02, updated 2013-09-02, Nanjiang Shu 
"""
usage_exp="""
Examples:
    
"""

def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile = ""
    infile = ""
    gramfile = ""
    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            infile = argv[i]
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
            elif argv[i] in ["-gram", "--gram"]:
                (gramfile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            infile = argv[i]
            i += 1

    if myfunc.checkfile(infile) != 0:
        return 1
    if myfunc.checkfile(gramfile) != 0:
        return 1

    grampairlist = myfunc.ReadPairList(gramfile)
    gramMapDict = {}
    for tup in grampairlist:
        gramMapDict[tup[0]] = tup[1]

    gi2taxidpairlist = myfunc.ReadPairList(infile)

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)

    for tup in gi2taxidpairlist:
        try:
            fpout.write("%s\t%s\t%s\n"%(tup[0], tup[1], gramMapDict[tup[1]]))
        except KeyError:
            fpout.write("%s\t%s\t%s\n"%(tup[0], tup[1], "NA"))

    myfunc.myclose(fpout)


#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

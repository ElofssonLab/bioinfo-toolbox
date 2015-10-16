#!/usr/bin/python
# Description:
import os
import sys
sys.path.append("%s/wk/MPTopo/src"%(os.environ['DATADIR3']))
import myfunc
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

usage_short="""
Usage: %s IDLISTFILE -euk EUK_IDLISTFILE -gram+ GRAM+_IDLISTFILE -gram- GRAM-_IDLISTIFILE [-o OUTFILE]
"""%(progname)

usage_ext="""
Description:
    create a file with the format
    SEQID\tNCBI_TaxID\tclass
    class can be euk, gram+, gram- or NA
    NCBI_TaxID can be empty

OPTIONS:
  -o OUTFILE    Output the result to OUTFILE
  -q            Quiet mode
  -h, --help    Print this help message and exit

Created 2014-09-10, updated 2014-09-10, Nanjiang Shu 
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
    idListFile = ""
    euk = ""
    gram_pos =""
    gram_neg = ""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            idListFile  = argv[i]
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
            elif argv[i] in ["-euk", "--euk"]:
                (euk, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-gram+", "--gram+"] :
                (gram_pos, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-gram-", "--gram-"] :
                (gram_neg, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            idListFile = argv[i]
            i += 1

    if myfunc.checkfile(idListFile, "idListFile") != 0:
        return 1
    if myfunc.checkfile(euk, "euk") != 0:
        return 1
    if myfunc.checkfile(gram_pos, "gram_pos") != 0:
        return 1
    if myfunc.checkfile(gram_neg, "gram_neg") != 0:
        return 1



    idList = myfunc.ReadIDList(idListFile)
    set_euk_idlist = set(myfunc.ReadIDList(euk))
    set_gram_pos_idlist = set(myfunc.ReadIDList(gram_pos))
    set_gram_neg_idlist = set(myfunc.ReadIDList(gram_neg))

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)

    NCBI_TaxID = ""
    for i in xrange(len(idList)):
        seqid = idList[i]
        cls = ""
        if seqid in set_euk_idlist:
            cls = "euk"
        elif seqid in set_gram_pos_idlist:
            cls = "gram+"
        elif seqid in set_gram_neg_idlist:
            cls = "gram-"
        else:
            cls = "NA"
        print >> fpout, "%s\t%s\t%s"%(seqid, NCBI_TaxID, cls)
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

#!/usr/bin/python
# Description:
import os
import sys
sys.path.append("%s/wk/MPTopo/src"%(os.environ['DATADIR3']))
import myfunc
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

usage_short="""
Usage: %s taxidwithtaxoFile -gram+ FILE -gram- FILE -euk FILE
"""%(progname)

usage_ext="""
Description:

OPTIONS:
  -o OUTFILE    Output the result to OUTFILE
  -q            Quiet mode
  -h, --help    Print this help message and exit
  -debug

Created 2013-09-02, updated 2014-10-09, Nanjiang Shu 
"""
usage_exp="""
Examples:
"""

def ReadSignalPFile(infile):#{{{
    hdl = myfunc.ReadLineByBlock(infile)
    dt = {}
    if hdl.failure:
        return 1
    lines = hdl.readlines()
    while lines != None:
        for line in lines:
            if line == "" or line[0] == "#":
                continue
            #seqid = myfunc.GetFirstWord(line)
            seqid = myfunc.GetSeqIDFromAnnotation(line)
            dt[seqid] = line
        lines = hdl.readlines()
    hdl.close()
    return dt
#}}}
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
    gramPositiveFile = ""
    gramNegativeFile = ""
    eukFile = ""

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
            elif argv[i] in ["-gram+", "--gram+"]:
                (gramPositiveFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-gram-", "--gram-"]:
                (gramNegativeFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-euk", "--euk"]:
                (eukFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True
                i += 1
            elif argv[i] in ["-debug", "--debug"]:
                g_params['isDEBUG'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            infile = argv[i]
            i += 1

    if myfunc.checkfile(infile, "taxidwithtaxoFile") != 0:
        return 1
    if myfunc.checkfile(gramPositiveFile, "gramPositiveFile") != 0:
        return 1
    if myfunc.checkfile(gramNegativeFile, "gramNegativeFile") != 0:
        return 1
    if myfunc.checkfile(eukFile, "eukFile") != 0:
        return 1

    gramPositiveDict = ReadSignalPFile(gramPositiveFile)
    gramNegativeDict = ReadSignalPFile(gramNegativeFile)
    eukDict = ReadSignalPFile(eukFile)


    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)

    hdl = myfunc.ReadLineByBlock(infile)
    if hdl.failure:
        return 1
    lines = hdl.readlines()
    while lines != None:
        for line in lines:
            strs = line.split("\t")
            if len(strs) == 3:
                seqid = strs[0].strip()
                taxo = strs[2].strip()
                info = ""
                try:
                    if taxo == "Gram+" or taxo == "gram+":
                        info = gramPositiveDict[seqid]
                    elif taxo == "Gram-" or taxo == "gram-":
                        info = gramNegativeDict[seqid]
                    elif taxo == "Euk" or taxo == "euk":
                        info = eukDict[seqid]

                    if g_params['isDEBUG']:
                        print >> sys.stderr, "%s: %s"%(seqid, taxo)
                except KeyError:
                    info = ""
                if info != "":
                    fpout.write("%s\n"%info)

        lines = hdl.readlines()
    hdl.close()
    myfunc.myclose(fpout)

#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['isDEBUG'] = False
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

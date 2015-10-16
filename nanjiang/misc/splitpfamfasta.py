#!/usr/bin/python
# Description:
import os
import sys
#sys.path.append("%s/wk/MPTopo/src"%(os.environ['DATADIR3']))
import myfunc
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

usage_short="""
Usage: %s PfamFastaFile -outpath DIR
"""%(progname)

usage_ext="""
Description:
  Give Pfam fasta file, split sequences named by each Pfam ID

OPTIONS:
  -outpath DIR  Set output path, compulsory
  -q            Quiet mode
  -h, --help    Print this help message and exit

Created 2013-04-23, updated 2013-04-23, Nanjiang Shu 
"""
usage_exp="""
Examples:
    %s Pfam-A.fasta -outpath pfamA-full-aa
"""%(progname)

def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}

def ExtractPfamIDFromDescription(description):#{{{
    try:
        return description[description.find(" PF")+1:].split(";")[0].split(".")[0]
    except IndexError:
        print >> sys.stderr, "Bad description: %s"%(description)
        raise
        return ""
#}}}
def ExportFastaForPfamID(pfamid, li, outpath, cnt):#{{{
    outfile = "%s%s%s.fa"%(outpath, os.sep, pfamid)
    try:
        fpout = open(outfile, "w")
        for rd in li:
            print >> fpout, ">%s"%(rd.description)
            print >> fpout, "%s"%(rd.seq)
        fpout.close()
        print "%5d %s output"%(cnt, outfile)
        return 0
    except IOError:
        print >> sys.stderr, "%5d Failed to write to file %s"%(cnt, outfile)
        return 1
#}}}
def SplitPfamFasta(infile, outpath):
    hdl = myfunc.ReadFastaByBlock(infile, 0, 0)
    if hdl.failure:
        return 1
    pfamIDThis = ""
    liThis = []
    recordList = hdl.readseq()
    cnt = 0
    while recordList != None:
        for rd in recordList:
            pfamid = ExtractPfamIDFromDescription(rd.description)
            if pfamIDThis == "":
                pfamIDThis = pfamid
            if pfamid == pfamIDThis:
                liThis.append(rd)
            else: # starting a new record, export this record
                ExportFastaForPfamID(pfamIDThis, liThis, outpath, cnt)
                cnt += 1
                pfamIDThis = pfamid
                liThis = []
                liThis.append(rd)
        recordList = hdl.readseq()
    if len(liThis) > 0:
        ExportFastaForPfamID(pfamIDThis, liThis, outpath, cnt)
        cnt += 1
    hdl.close()
    return 0


def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = ""
    infile = ""

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
            elif argv[i] in ["-outpath", "--outpath"]:
                (outpath, i) = myfunc.my_getopt_str(argv, i)
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
    if outpath == "":
        print >> sys.stderr, "outpath not set"
    elif not os.path.exists(outpath):
        os.system("mkdir -p %s"%(outpath))

    SplitPfamFasta(infile, outpath)

#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

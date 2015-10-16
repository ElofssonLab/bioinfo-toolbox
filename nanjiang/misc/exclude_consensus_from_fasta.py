#!/usr/bin/python
# Description:
# exclude the consensus sequence from the fasta file
import os
import sys
#sys.path.append("%s/wk/MPTopo/src"%(os.environ['DATADIR3']))
import myfunc
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

usage_short="""
Usage: %s FILE [FILE ...] [-l LISTFILE] [-outpath DIR]
"""%(progname)

usage_ext="""
Description:
    Exclude the censensus sequence from the fatsa file 

OPTIONS:
  -outpath DIR  Output the result to outpath, (default: the same folder as the
                input file), outname = $rootname.nocons.fasta)
  -l LISTFILE   Set the listfile
  -q            Quiet mode
  -h, --help    Print this help message and exit

Created 2014-11-25, updated 2014-11-25, Nanjiang Shu
"""
usage_exp="""
Examples:
"""

def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}

def ExcludeConsensus(infile, g_outpath):
    if g_outpath == "":
        outpath = os.path.dirname(infile)
        if outpath == "":
            outpath = "."
    else:
        outpath = g_outpath
    rootname = os.path.basename(os.path.splitext(infile)[0])
    outfile = "%s%s%s.nocons.fasta"%(outpath, os.sep, rootname)
    try:
        fpout = open(outfile, "w")

        hdl = myfunc.ReadFastaByBlock(infile)
        if hdl.failure:
            return 1
        recordList = hdl.readseq()
        while recordList != None:
            for record in recordList:
                if record.seqid.lower().find("consensus") == -1:
                    fpout.write(">%s\n%s\n"%(record.description, record.seq))
            recordList = hdl.readseq()
        fpout.close()

        hdl.close()
    except IOError:
        print sys.stderr,"Failed to write to file %s"%(outfile)
        return 1


def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = ""
    outfile = ""
    fileListFile = ""
    fileList = []

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            fileList.append(argv[i])
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
            elif argv[i] in ["-outpath", "--outpath"]:
                (outpath, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-l", "--l"] :
                (fileListFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            fileList.append(argv[i])
            i += 1


    if fileListFile != "":
        fileList += myfunc.ReadIDList(fileListFile)

    if outpath != "" and not os.path.exists(outpath):
        os.system("mkdir -p %s"%(outpath))

    for i in xrange(len(fileList)):
        ExcludeConsensus(fileList[i], outpath)
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

#!/usr/bin/env python
# rewite fasta file, sequence in one line
import os,sys
import myfunc
progname =  os.path.basename(sys.argv[0])
usage="""
Usage: %s fastafile [-o OUTFILE]

Description: rewrite fasta file, so that sequence is written in one line

Created 2013-04-16, updated 2013-04-18, Nanjiang Shu 
"""%(progname)

def PrintHelp():
    print usage

def ReWriteFasta(infile, outfile):#{{{
    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    hdl = myfunc.ReadFastaByBlock(infile, 0, 1)
    if hdl.failure:
        return 1
    recordList = hdl.readseq()
    while recordList != None:
        for rd in recordList:
            fpout.write(">%s\n"%rd.description)
            fpout.write("%s\n"%rd.seq)
        recordList = hdl.readseq()
    hdl.close()
    myfunc.myclose(fpout)
    return 0
#}}}
def main():#{{{
    numArgv = len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    argv = sys.argv

    outfile = ""
    infile = ""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            infile = sys.argv[i]
            isNonOptionArg=False
            i = i + 1
        elif sys.argv[i] == "--":
            isNonOptionArg=True
            i = i + 1
        elif sys.argv[i][0] == "-":
            if sys.argv[i] in [ "-h", "--help"]:
                PrintHelp()
                return 1
            elif sys.argv[i] in [ "-o", "--o"]:
                outfile, i = myfunc.my_getopt_str(argv, i)
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % argv[i])
                return 1
        else:
            infile = argv[i]
            i += 1

    if myfunc.checkfile(infile) != 0:
        return 1

    return ReWriteFasta(infile, outfile)
#}}}
if __name__ == '__main__' :
    sys.exit(main())

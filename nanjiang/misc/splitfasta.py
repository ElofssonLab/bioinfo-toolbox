#!/usr/bin/env python
# split the given fasta sequence file in pieces
# ChangeLog 2012-05-23
#   options added
#       -maxseq INT  specified maximum number of sequences in each splitted
#       file, default 1
#   fasta file can also be splitted in pieces with specified numbers
#   or pieces with specified number of sequences
# ChangeLog 2013-02-21
#   Set verbose mode, by default, splitted file name do not output
import sys
import re
import os
import myfunc
from math import ceil

BLOCK_SIZE = 100000
progname =  os.path.basename(sys.argv[0])

usage="""
Usage:  %s [OPTIONS] [-i] fastafile

Description: Split the given sequence file in pieces.
             By default, it splits the fasta file to individual
             sequence files, one sequence per file.
OPTIONS:
  -i       FILE   Input file
  -outpath  DIR   Output the result to dir, (default: ./)
  -ext      STR   Set the file extension, default = aa
  -nseq     INT   Split fasta file into pieces with (at most) N sequences,
                  (default: 1). When nseq > 1, namemode is forced to 1
  -nfile    INT   Split fasta file into N pieces. When this is set, namemode
                  is forced to 1
  -namemode INT   Name mode 0 or 1, (default: 0)
                  0: file name from annotation line
                  1: file named by $rootname_i, i= 0,1,2...
  -q              Quiet mode, do not write any messages
  -v              Verbose mode, print each output file
  -h, --help      Print this help message and exit

Outdated options:
  -nameseq        name the splitted files by $rootname_i, i = 0,1,2...
                  this is useful when it is the filename extracted from the 
                  annotation line having a duplicated name

Created 2010-10-22, updated 2012-07-05, Nanjiang

Examples:
# split the fasta file in to one sequence per file, name obtained from
# annotation line
    %s example.fasta
    %s example.fasta -outpath outdir

# split the fasta file into (at most) 10 pieces
    %s example.fasta -outpath outdir -nfile 10

# split the fasta file into pieces, each piece has (at most) 15 sequences
    %s example.fasta -outpath outdir -nseq 15
"""%(progname, progname, progname, progname, progname)

def PrintHelp():
    print usage


def OutputSplittedSeq(seqWithAnno, rootname, cntsplit, #{{{
        cntseq_of_split, fpout):
# return (cntsplit, cntseq_of_split, fpout)
    begseq = seqWithAnno.find("\n")
    seq = seqWithAnno[begseq:]
    seqID = myfunc.GetSeqIDFromAnnotation(seqWithAnno[0:begseq])
    outfile = ""

    if cntseq_of_split < g_params['numseq_per_split']:
        if fpout == None:
            if g_params['isNameFileSequentially']:
                outfile = (g_params['outpath']+os.sep+ rootname + "_%d" %
                        cntsplit  + "." + g_params['file_ext'])
            else:
                if seqID == "":
                    seqID = rootname + "_%d" % cntsplit
                outfile = (g_params['outpath']+os.sep+ seqID + "." +
                        g_params['file_ext'])

            try:
                fpout = open(outfile,"w")
            except IOError:
                print >> sys.stderr ,"Failed to write to file %s"%outfile
                return 1

        fpout.write("%s"%seqWithAnno[0:begseq])
        fpout.write("%s\n" % seq)

        cntseq_of_split += 1
        if cntseq_of_split >= g_params['numseq_per_split']: 
            fpout.close()
            fpout = None
            cntseq_of_split = 0
            cntsplit += 1
        if g_params['verbose'] >= 2:
            if outfile != "":
                print >> sys.stdout, "split %d\t%s output"%(cntsplit, outfile)
    else:
        msg = "Error! cntseq_of_split (%d) >= numseq_per_split (%d)"
        print >> sys.stderr,  msg%(cntseq_of_split,
                g_params['numseq_per_split'])

    return (cntsplit, cntseq_of_split, fpout)
#}}}
def SplitFasta(inFile):#{{{
# The faster version
    if 'numsplit' in g_params or g_params['numseq_per_split'] > 1:
        g_params['isNameFileSequentially'] = True
    if 'numsplit' in g_params and g_params['numsplit'] > 1:
        numTotalSeq = myfunc.CountFastaSeq(inFile)
        g_params['numseq_per_split'] = int(ceil(
            numTotalSeq/float(g_params['numsplit'])))
        if g_params['verbose'] >= 1:
            msg = "file %s (with %d sequences) is going to"\
                    "be splitted into %d files"
            print  msg%(inFile, numTotalSeq, g_params['numsplit'])

    rootname = os.path.basename(os.path.splitext(inFile)[0])

    cntTotalSeq = 0
    fpout = None
    cntsplit = 0
    cntseq_of_split = 0

    fpin = open(inFile, "r")
    buff = fpin.read(BLOCK_SIZE)
    brokenSeqWithAnnoLine = ""; ##for the annotation line broken by BLOCK read
    while buff:
        beg=0
        end=0
        while 1:
            if brokenSeqWithAnnoLine:
                if brokenSeqWithAnnoLine[len(brokenSeqWithAnnoLine)-1] == "\n":
                    end=buff.find(">")
                else:
                    end=buff.find("\n>")
                if end >= 0:
                    seqWithAnno = brokenSeqWithAnnoLine + buff[0:end]
                    (cntsplit, cntseq_of_split, fpout) = OutputSplittedSeq(
                            seqWithAnno, rootname, cntsplit, cntseq_of_split,
                            fpout)
                    brokenSeqWithAnnoLine = ""
                    cntTotalSeq += 1
                    beg=end
                else:
                    brokenSeqWithAnnoLine += buff
                    break

            beg=buff.find(">",beg)
            end=buff.find("\n>",beg+1)
            if beg >= 0:
                if end >=0:
                    seqWithAnno=buff[beg:end]
                    (cntsplit, cntseq_of_split, fpout) = OutputSplittedSeq(
                            seqWithAnno, rootname, cntsplit, cntseq_of_split,
                            fpout)
                    cntTotalSeq += 1
                    beg=end
                else:
                    brokenSeqWithAnnoLine=buff[beg:]
                    break
            else:
                break

        buff = fpin.read(BLOCK_SIZE)
    fpin.close()
    if brokenSeqWithAnnoLine:
        seqWithAnno=brokenSeqWithAnnoLine
        (cntsplit, cntseq_of_split, fpout) = OutputSplittedSeq(
                seqWithAnno, rootname, cntsplit, cntseq_of_split, fpout)
        cntTotalSeq += 1

    return cntTotalSeq

#}}}
def main(g_params): #{{{
    argv = sys.argv
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    inFile = ""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            isNonOptionArg=False
            i = i + 1
        elif sys.argv[i] == "--":
            isNonOptionArg=True
            i = i + 1
        elif sys.argv[i][0] == "-":
            if sys.argv[i] in  [ "-h",  "--help"]:
                PrintHelp()
                return 1
            elif sys.argv[i] in ["-i", "--infile"]:
                inFile=sys.argv[i+1]
                i = i + 2
            elif sys.argv[i] == "-ext" or sys.argv[i] == "--ext":
                g_params['file_ext'] = sys.argv[i+1]
                i = i + 2
            elif sys.argv[i] in ["-outpath", "--outpath"]:
                g_params['outpath'] = sys.argv[i+1]
                i = i + 2
            elif sys.argv[i] in ["-q" , "--q" , "--quiet" ]:
                g_params['isQuiet'] = True
                i = i + 1
            elif sys.argv[i] in ["-v" , "--v" ]:
                g_params['verbose'] = 2
                i = i + 1
            elif sys.argv[i] in ["-nameseq", "--nameseq" ]:
                g_params['namemode'] = 1
                i = i + 1
            elif sys.argv[i] in ["-namemode", "--namemode" ]:
                g_params['namemode'] = int(argv[i+1])
                i = i + 2
            elif sys.argv[i] in  ["-nfile", "--nfile" ]:
                g_params['numsplit'] = int(argv[i+1])
                i = i + 2
            elif sys.argv[i] in  ["-nseq", "--nseq" ]:
                g_params['numseq_per_split'] = int(argv[i+1])
                i = i + 2
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i])
                return 1
        else:
            inFile = sys.argv[i]
            i+=1

    if g_params['isQuiet']: 
        g_params['verbose'] = 0

    if inFile == "":
        print sys.stderr, "infile not set. Exit."
        return 1
    elif not os.path.exists(inFile):
        print >> sys.stderr,"infile %s does not exist. Exit."%(inFile)
        return 1

    if not os.path.exists(g_params['outpath']):
        os.system("mkdir -p %s"%g_params['outpath'])

    if g_params['namemode'] == 0:
        g_params['isNameFileSequentially'] = False
    elif g_params['namemode'] == 1:
        g_params['isNameFileSequentially'] = True
    else:
        msg = "Wrong namemode (%d). Exit."
        print >> sys.stderr, msg%(g_params['namemode'])
        return 1

#     print "namemode", g_params['namemode']
#     print g_params['isNameFileSequentially']
    numSeq = SplitFasta(inFile)
    if g_params['verbose'] >= 1:
        print "%d sequences are splitted and output to %s"%(
            numSeq, g_params['outpath'])

    return 0
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = False
    # verbose level,
    # 0: do not print anything except error message
    # 1: print minimum information
    # 2. print more information
    g_params['verbose'] = 1
    g_params['outpath'] = "./"
    g_params['namemode'] = 0
    g_params['numseq_per_split'] = 1
    g_params['file_ext'] = "aa"
    g_params['isNameFileSequentially'] = False
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))
    # Check argv

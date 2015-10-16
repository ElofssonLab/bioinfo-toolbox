#!/usr/bin/env python
# Description:
import os
import sys
import myfunc
import tempfile
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))
usage = """
Usage:  %s FILE [-o OUTFILE]
Description: Replace non standard amino acid in sequence by 'X'
Options:
  -o OUTFILE    Output the result to OUTFILE
  -l LISTFILE   Set the listfile
  -overwrite    Overwrite the original file, in that case, 
                outfile is not needed
  -q            Quiet mode
  -h, --help    Print this help message and exit

Created 2013-03-11, updated 2013-03-11, Nanjiang Shu
"""%(progname)


STD1CharAA_alphabet = "ARNDCQEGHILKMFPSTWYVX"

def PrintHelp():
    print usage
def IsHasNonStandardAminoAcid(infile):#{{{
    isHasNonStdAA = 0
    hdl = myfunc.ReadFastaByBlock(infile, 0, 0)
    if hdl.failure:
        return -1
    recordList = hdl.readseq()
    while recordList != None:
        for rd in recordList:
            for i in xrange(len(rd.seq)):
                if STD1CharAA_alphabet.find(rd.seq[i]) == -1:
                    isHasNonStdAA = 1
                    break
        recordList = hdl.readseq()
    hdl.close()
    return isHasNonStdAA
#}}}
def ReplaceNonStandardAminoAcid(seq):#{{{
    seq = seq.upper()
    li = []
    apd = li.append
    for i in xrange(len(seq)):
        if STD1CharAA_alphabet.find(seq[i]) != -1:
            apd(seq[i])
        else:
            apd('X')
    return ''.join(li)
#}}}
def CleanAASeq(infile, isOverWrite, fpout):#{{{
    isHasNonStandardAminoAcid = IsHasNonStandardAminoAcid(infile)
    if isHasNonStandardAminoAcid == 0:
        if isOverWrite:
            msg = "seqfile %s is already cleaned. Ignore"
            print >> sys.stderr, msg%(infile)
        else:
            fpin = open(infile, "r")
            BLOCK_S = 10000
            buff = fpin.read(BLOCK_S)
            while buff:
                fpout.write(buff)
                buff = fpin.read(BLOCK_S)
            fpin.close()
        return 0
    elif isHasNonStandardAminoAcid == -1:
        return -1

    fpout_local = None
    if not isOverWrite:
        fpout_local = fpout
    else:
        try:
            fpout_local = tempfile.NamedTemporaryFile(delete=False)
        except IOError:
            msg = "Failed to write to temporary file for running seqfile %s"
            print >> sys.stderr, msg%(infile)
            return -1
    hdl = myfunc.ReadFastaByBlock(infile,0,0)
    if hdl.failure:
        return -1
    recordList = hdl.readseq()
    while recordList != None:
        for rd in recordList:
            seq = ReplaceNonStandardAminoAcid(rd.seq)
            fpout_local.write(">%s\n"%(rd.description))
            fpout_local.write("%s\n"%(seq))
        recordList = hdl.readseq()
    hdl.close()
    if isOverWrite:
        fpout_local.close()
        #print "tmpfile=",fpout_local.name
        os.system("/bin/mv -f %s %s"%(fpout_local.name, infile))
        os.system("chmod 644 %s"%(infile))
        if not g_params['isQuiet']:
            msg = "seqfile %s cleaned"
            print msg%(infile)
#}}}

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile=""
    fileListFile = ""
    fileList = []
    isOverWrite=0

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
            elif argv[i] in ["-outpath", "--outpath"]:
                outpath = argv[i+1]
                i += 2
            elif argv[i] in ["-l", "--l"] :
                fileListFile = argv[i+1]
                i += 2
            elif argv[i] in ["-overwrite", "--overwrite"] :
                isOverWrite = 1
                i += 1
            elif argv[i] in ["-q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            fileList.append(argv[i])
            i += 1

    if fileListFile != "":
        try:
            fp = open(fileListFile,"r")
            fileList += fp.read().split()
            fp.close()
        except IOError:
            msg = "Failed to read idlistfile {}."
            print >> sys.stderr, msg.format(fileListFile)
    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    for i in xrange(len(fileList)):
        CleanAASeq(fileList[i], isOverWrite, fpout)
    myfunc.myclose(fpout)
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = False
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

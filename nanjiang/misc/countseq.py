#!/usr/bin/env python
# get the number of sequences in the fastafile
import sys
import os
import myfunc

BLOCK_SIZE = 100000

progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))
usage="""
Usage:  %s [-l LISTFILE] [-o OUTFILE]
        %s[FILE [FILE ...]]
Options:
  -l  LISTFILE   Set list of filenames
  -o   OUTFILE   Output the result to OUTFILE
  -bs      INT   Size for blocks when reading file, (default: 50000)
  -nf            Do not print file name
  -h, --help     Print this help message and exit

Created 2010-11-03, updated 2013-02-16, Nanjiang Shu

Examples:
    %s test.fa
"""%(progname, wspace, progname)

def PrintHelp():
    print usage

def CountSeq(inFile, fpout):#{{{
# The faster version
    isPrintFileName = g_params['isPrintFileName']
    try:
        isFirstSeq=True
        isPreviousBuffEndWithNewLine=False
        fpin = open(inFile, "r")
        buff = fpin.read(BLOCK_SIZE)
        cntSeq = 0
        while buff:
            if isFirstSeq and buff[0] == '>':
                cntSeq +=1
                isFirstSeq = False
            if isPreviousBuffEndWithNewLine and buff[0] == '>':
                cntSeq += 1
                isPreviousBuffEndWithNewLine = False
            cntSeq += buff.count("\n>")
            if buff[len(buff)-1] == '\n':
                isPreviousBuffEndWithNewLine = True
            buff = fpin.read(BLOCK_SIZE)
        fpin.close()

        if isPrintFileName:
            fpout.write("%s\t"%inFile)
        fpout.write("%d\n"%(cntSeq))
        return cntSeq
    except IOError:
        print >> sys.stderr, "Failed to read file %s"%(inFile)
        return -1
#}}}

def main(g_params):#{{{
    # Check argv
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outFile=""
    fileList=[]
    listFile=""

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
            if sys.argv[i] ==  "-h" or  sys.argv[i] == "--help":
                PrintHelp()
                return 1
            elif sys.argv[i] in ["-o", "--outfile"]:
                try:
                    outFile = sys.argv[i+1]
                except IndexError:
                    msg = "option %s should be followed by a string"
                    print >> sys.stderr, msg%(sys.argv[i])
                    return 1
                i += 2
            elif sys.argv[i] in [ "-nf", "--nf"]:
                g_params['isPrintFileName'] = False
                i += 1
            elif sys.argv[i] in ["-l", "--list"]:
                listFile = sys.argv[i+1]
                i +=  2
            elif sys.argv[i] in [ "-bs" ,  "--block-size" , "-block-size"]:
                BLOCK_SIZE = int(sys.argv[i+1])
                if BLOCK_SIZE < 0:
                    print >> sys.stderr,"Error! BLOCK_SIZE should >0"
                    return 1
                i = i + 2
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i])
                return 1
        else:
            fileList.append(sys.argv[i])
            i += 1

    if listFile and os.path.exists(listFile):
        fileList += myfunc.ReadListFile(listFile, "\n")

    if len(fileList) <= 0:
        print >> sys.stderr,"Error! No input set."
        return 1


    fpout = myfunc.myopen(outFile, sys.stdout, "w", False)

    for f in fileList:
        CountSeq(f, fpout)

    myfunc.myclose(fpout)

#}}}
def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isPrintFileName'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))


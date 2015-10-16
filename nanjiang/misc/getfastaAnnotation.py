#!/usr/bin/env python
# get the annotation line of the fasta file, the leading '>' is kept, one record per line
import sys,re,os
import myfunc

BLOCK_SIZE = myfunc.BLOCK_SIZE

usage="""
Usage:  getfastaAnnotation.py [Options] [-i] fastafile
Options:
  -i         <file>    : input file
  -o         <file>    : outputfile
  -m|--method 1|2      : method for reading, default = 2
  -bs|--block-size int : size for blocks when reading file, default = 50000
  -h|--help            : print this help message and exit
Created 2010-08-20, updated 2010-08-24, Nanjiang
"""

def PrintHelp():
    print usage

def PrintFastaAnnotation(inFile, fpout):#{{{
    fpin = open(inFile, "r")
    line = fpin.readline()
    while line:
        line = line.rstrip('\n').strip()
        if line and line[0] == ">":
            print >>fpout,"%s"%line
        line = fpin.readline()
    fpin.close()
    return 0
#}}}
def PrintFastaAnnotation2(inFile, fpout):#{{{
# The faster version
    fpin = open(inFile, "r")
    buff = fpin.read(BLOCK_SIZE)
    brokenAnnoLine=""; ##for the annotation line broken by BLOCK read
    while buff:
        beg=0
        end=0
        while 1:
            if brokenAnnoLine:
                end=buff.find("\n")
                if end >= 0:
                    line = brokenAnnoLine + buff[0:end]
                    print >> fpout, line
                    brokenAnnoLine = ""
                    beg=end
                else:
                    brokenAnnoLine += buff
                    break

            beg=buff.find(">",beg)
            end=buff.find("\n",beg+1)
            if beg >= 0:
                if end >=0:
                    line=buff[beg:end]
                    print >> fpout,line
                    beg=end
                else:
                    brokenAnnoLine=buff[beg:]
                    break
            else:
                break

        buff = fpin.read(BLOCK_SIZE)
    fpin.close()
    return 0
#}}}


if __name__ == '__main__' :
    # Check argv
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        sys.exit(1)

    outFile=""
    inFile=""
    method=2


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
                sys.exit(0)
            elif sys.argv[i] == "-i" or sys.argv[i] == "--infile":
                inFile=sys.argv[i+1]
                i = i + 2
            elif sys.argv[i] == "-o" or sys.argv[i] == "--outfile":
                outFile=sys.argv[i+1]
                i = i + 2
            elif sys.argv[i] == "-m" or sys.argv[i] == "--method" or sys.argv[i] == "-method":
                method=int(sys.argv[i+1])
                if method < 1 or method > 2:
                    print >> sys.stderr,"Error! method should be 1 or 2"
                    sys.exit(1)
                i = i + 2
            elif sys.argv[i] == "-bs" or sys.argv[i] == "--block-size" or sys.argv[i] == "-block-size":
                BLOCK_SIZE=int(sys.argv[i+1])
                if BLOCK_SIZE < 0:
                    print >> sys.stderr,"Error! BLOCK_SIZE should >0"
                    sys.exit(1)
                i = i + 2
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i])
                sys.exit(1)
        else:
            inFile=sys.argv[i]
            i+=1


    if inFile == "":
        print >> sys.stderr,"Error! Input file not set."

    fpout = sys.stdout
    if outFile != "":
        fpout = open(outFile,"w")

    try :
        if method ==1:
            PrintFastaAnnotation(inFile,fpout)
        else:
            PrintFastaAnnotation2(inFile,fpout)
        if fpout != sys.stdout:
            fpout.close()

    except :
        print >>sys.stderr, "except for the input file: %s" % inFile
        raise 

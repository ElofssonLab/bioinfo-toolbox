#!/usr/bin/env python
# get the length of the sequences in the fasta file
# Note: this python code is about 7 times faster than the perl script getseqlen.pl
import sys
import os
import myfunc


progname = os.path.basename(sys.argv[0])
usage="""
Usage: %s fastafile [-o OUTFILE]
Description:
    Get the length sequences given a fasta seq file
OPTIONS:
  -o OUTFILE         Output the result to file
  -just-print-sum    Just print the total number of amino acids
  -i, -printid       Print sequence IDs, (default: no)
  -bs INT            Set block size when reading file, (default: 100000)
  -h, --help         print this help message and exit

Created 2010-09-02, updated 2013-05-21, Nanjiang

Examples:
    %s test.fa           #print the length of sequences, one record by line
    %s -printid test.fa  #print seqid and length of sequences
    %s -just-print-sum test.fa # print the total number of amino acids of the file

"""%(progname, progname, progname, progname)

def PrintHelp():
    print usage

def GetFromRawSeq(seqWithAnno, isPrintID, isJustPrintSum, fpout):#{{{
#     begseq=seqWithAnno.find("\n")
#     seq=seqWithAnno[begseq:]
#     seq=seq.replace('\n','').replace(' ','')
    length=len(seqWithAnno[seqWithAnno.find("\n"):].replace('\n','').replace(' ',''))
    if not isJustPrintSum:
        if isPrintID:
            seqID=myfunc.GetSeqIDFromAnnotation(seqWithAnno)
            fpout.write("%s\t"%seqID)
        fpout.write("%d\n" % length)
    return length
#}}}
def Getseqlen(infile, isPrintID, isJustPrintSum, BLOCK_SIZE, fpout):#{{{
# The faster version
    try:
        isFirstSeq=True
        totalLength=0
        fpin = open(infile, "r")
        buff = fpin.read(BLOCK_SIZE)
        brokenSeqWithAnnoLine=""; ##for the annotation line broken by BLOCK read
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
                        length = GetFromRawSeq(seqWithAnno, isPrintID,
                                isJustPrintSum, fpout)
                        totalLength +=length
                        brokenSeqWithAnnoLine = ""
                        beg=end
                    else:
                        brokenSeqWithAnnoLine += buff
                        break

                beg=buff.find(">",beg)
                end=buff.find("\n>",beg+1)
                if beg >= 0:
                    if end >=0:
                        seqWithAnno=buff[beg:end]
                        length=GetFromRawSeq(seqWithAnno, isPrintID, isJustPrintSum, fpout)
                        totalLength +=length
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
            length=GetFromRawSeq(seqWithAnno, isPrintID, isJustPrintSum, fpout)
            totalLength +=length

        if isJustPrintSum:
            fpout.write("%d\n"%totalLength)

        return 0
    except IOError:
        print >> sys.stderr, "Failed to read infile %s"%(infile)
        return 1

#}}}
def main(g_params):#{{{#{{{
    # Check argv
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile = ""
    infile = ""
    BLOCK_SIZE = 100000
    isPrintID = False
    isJustPrintSum = False

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            isNonOptionArg=False
            infile = sys.argv[i]
            i = i + 1
        elif sys.argv[i] == "--":
            isNonOptionArg=True
            i = i + 1
        elif sys.argv[i][0] == "-":
            if sys.argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 0
            elif sys.argv[i] in ["-i", "--i", "-printid", "--printid"]:
                isPrintID = True
                i += 1
            elif sys.argv[i] in ["-just-print-sum", "--just-print-sum"]:
                isJustPrintSum = True
                i += 1
            elif sys.argv[i] in [ "-o", "--o", "-outfile", "--outfile"]:
                outfile, i = myfunc.my_getopt_str(sys.argv, i)
            elif sys.argv[i] in [ "-bs", "--bs", "-block-size", "--block-size"]:
                BLOCK_SIZE, i = myfunc.my_getopt_int(sys.argv, i)
                if BLOCK_SIZE < 0:
                    print >> sys.stderr,"Error! BLOCK_SIZE should >0"
                    return 1
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i])
                return 1
        else:
            infile=sys.argv[i]
            i += 1

    if myfunc.checkfile(infile) != 0:
        return 1

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    status = Getseqlen(infile, isPrintID, isJustPrintSum, BLOCK_SIZE, fpout)
    myfunc.myclose(fpout)

    return status
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    return g_params
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

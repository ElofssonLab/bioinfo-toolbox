#!/usr/bin/env python
# convert msa to mpa format
# mpa format description:
# Length: 34503
# >description
# 1 A  5 C 9 E
# 
# that is list only positions of residues, gaps are not shown
# 
# 
# 
import sys ,re,os
import myfunc

progname =  os.path.basename(sys.argv[0])

usage="""
Usage: %s MSA-FASTA-FILE [-o OUTFILE]

Description: Convert multiple sequence alignment to mpa format

  -o  OUTFILE   Output the result to OUTFILE
  -if STR       input format, (default: mfa)
  -h, --help    Print this help message and exit

Created 2013-04-16, updated 2013-04-16, Nanjiang Shu
"""%(progname)

GAP='-'
BLOCK_SIZE=100000

def PrintHelp():
    print usage

def MFA2MPA_obsolete(infile, fpout):#{{{
    hdl = myfunc.ReadFastaByBlock(infile, 0, 1)
    if hdl.failure:
        return 1
    lengthList = []
    recordList = hdl.readseq()
    while recordList != None:
        for rd in recordList:
            print >> fpout, ">%s"%(rd.description)
            print >> fpout, "Length: %d"%(len(rd.seq))
            seq = rd.seq
            firstrespos = -1
            for i in xrange(len(seq)):
                if seq[i] != GAP:
                    if firstrespos == -1:
                        firstrespos = i
                        fpout.write("%d %s "%(i, seq[i]))
                    else:
                        fpout.write("%d %s "%(i-firstrespos, seq[i]))
            fpout.write("\n")
            lengthList.append(len(seq))
        recordList = hdl.readseq()
    hdl.close()
    if len(set(lengthList)) > 1:
        msg = "Warning! Length of the MSA file %s are not equal!"
        print >> sys.stderr, msg%(infile)

    return 0
#}}}
def MFA2MPA_old2(infile, fpout):#{{{
    hdl = myfunc.ReadFastaByBlock(infile, 0, 1)
    if hdl.failure:
        return 1
    lengthList = []
    recordList = hdl.readseq()
    while recordList != None:
        for rd in recordList:
            print >> fpout, ">%s"%(rd.description)
            seq = rd.seq
            gapPosList = myfunc.GetSegPos(seq, GAP)
            num = len(gapPosList)
            length = len(seq)
            if num < 1 :
                fpout.write("%s\n"%(seq))
            else:
                if gapPosList[0][0] > 0:
                    fpout.write("%s "%(seq[0:gapPosList[0][0]]))
                for i in xrange(num-1):
                    fpout.write("%d-%d "%(gapPosList[i][0], gapPosList[i][1])) 
                    fpout.write("%s "%(seq[gapPosList[i][1]:gapPosList[i+1][0]]))
                fpout.write("%d-%d "%(gapPosList[num-1][0], gapPosList[num-1][1])) 
                if gapPosList[num-1][1] < length:
                    fpout.write("%s"%(seq[gapPosList[num-1][1]:length]))

                fpout.write("\n")
            lengthList.append(length)
        recordList = hdl.readseq()
    hdl.close()
    if len(set(lengthList)) > 1:
        msg = "Warning! Length of the MSA file %s are not equal!"
        print >> sys.stderr, msg%(infile)
    return 0
#}}}
def MFA2MPA(infile, fpout):#{{{
    hdl = myfunc.ReadFastaByBlock(infile, 0, 1)
    if hdl.failure:
        return 1
    lengthList = []
    recordList = hdl.readseq()
    while recordList != None:
        for rd in recordList:
            print >> fpout, ">%s"%(rd.description)
            seq = rd.seq
            gapPosList = myfunc.GetSegPos(seq, GAP)
            num = len(gapPosList)
            length = len(seq)
            if num < 1 :
                fpout.write("%s\n"%(seq))
            else:
                if gapPosList[0][0] > 0:
                    fpout.write("%s "%(seq[0:gapPosList[0][0]]))
                for i in xrange(num-1):
                    fpout.write("%d "%(gapPosList[i][1] - gapPosList[i][0]))
                    fpout.write("%s "%(seq[gapPosList[i][1]:gapPosList[i+1][0]]))
                fpout.write("%d "%(gapPosList[num-1][1] - gapPosList[num-1][0]))
                if gapPosList[num-1][1] < length:
                    fpout.write("%s"%(seq[gapPosList[num-1][1]:length]))

                fpout.write("\n")
            lengthList.append(length)
        recordList = hdl.readseq()
    hdl.close()
    if len(set(lengthList)) > 1:
        msg = "Warning! Length of the MSA file %s are not equal!"
        print >> sys.stderr, msg%(infile)
    return 0
#}}}

def MSA2MPA(infile, input_format, fpout):
    if input_format in ["mfa", "fa", "fasta"]:
        MFA2MPA(infile, fpout)
        return 0
    else:
        return 1


def main():#{{{
    numArgv = len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    argv = sys.argv

    outfile = ""
    infile = ""
    input_format = "mfa"

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
            elif sys.argv[i] in [ "-if", "--if"]:
                input_format,i  = myfunc.my_getopt_str(argv, i)
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % argv[i])
                return 1
        else:
            infile = argv[i]
            i += 1

    if myfunc.checkfile(infile, "MSA file") != 0:
        return 1

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)

    MSA2MPA(infile, input_format, fpout)

    myfunc.myclose(fpout)


#}}}
if __name__ == '__main__' :
    sys.exit(main())

#!/usr/bin/env python
# convert mpa in the format of 
# mpa format description:
# Length: 34503
# >description
# 1 A  5 C 9 E
# 
# that is list only positions of residues, gaps are not shown
# 
# to msa
# 
# 
import sys ,re,os
import myfunc

progname =  os.path.basename(sys.argv[0])

usage="""
Usage: %s MPA_FILE [-o OUTFILE]

Description: Convert multiple sequence alignment to mpa format

  -o  OUTFILE   Output the result to OUTFILE
  -of STR       input format, (default: mfa)
  -h, --help    Print this help message and exit

Created 2013-04-16, updated 2013-04-16, Nanjiang Shu
"""%(progname)

GAP='-'
BLOCK_SIZE=100000

def PrintHelp():
    print usage

def MPA2MSA_old(infile, output_format, fpout):#{{{
    hdl = myfunc.ReadLineByBlock(infile)
    if hdl.failure:
        return 1
    lengthList = []
    remainLineList = []
    lines = hdl.readlines()
    while lines != None:
        lines = remainLineList + lines
        numLine = len(lines)
        numRD = numLine/2
        for i in xrange(numRD):
            fpout.write("%s\n"%lines[2*i])
            strs = lines[2*i+1].split()
            for ss in strs:
                if ss.find("-") != -1:
                    strs1 = ss.split("-")
                    b = int(strs1[0])
                    e = int(strs1[1])
                    li = ["-"]*(e-b)
                    fpout.write("%s"%(''.join(li)))
                else:
                    fpout.write("%s"%(ss))
            fpout.write("\n")
        if numRD*2 < numLine:
            remainLineList = [lines[numLine-1]]
        else:
            remainLineList = []
        lines = hdl.readlines()
    hdl.close()

    return 0
#}}}
def MPA2MSA(infile, output_format, fpout):#{{{
    hdl = myfunc.ReadLineByBlock(infile)
    if hdl.failure:
        return 1
    lengthList = []
    remainLineList = []
    lines = hdl.readlines()
    while lines != None:
        lines = remainLineList + lines
        numLine = len(lines)
        numRD = numLine/2
        for i in xrange(numRD):
            li = []
            fpout.write("%s\n"%lines[2*i])
            strs = lines[2*i+1].split()
            for ss in strs:
                if ss[0].isdigit():
                    lgap = int(ss)
                    li.append("-"*lgap)
                else:
                    li.append(ss)
            fpout.write("%s\n"%("".join(li)))
        if numRD*2 < numLine:
            remainLineList = [lines[numLine-1]]
        else:
            remainLineList = []
        lines = hdl.readlines()
    hdl.close()

    return 0
#}}}
# def MPA2MSA(infile, input_format, fpout):
#     if  in ["mfa", "fa", "fasta"]:
#         MFA2MPA(infile, fpout)
#         return 0
#     else:
#         return 1


def main():#{{{
    numArgv = len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    argv = sys.argv

    outfile = ""
    infile = ""
    output_format = "mfa"

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
            elif sys.argv[i] in [ "-of", "--of"]:
                output_format,i  = myfunc.my_getopt_str(argv, i)
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % argv[i])
                return 1
        else:
            infile = argv[i]
            i += 1

    if myfunc.checkfile(infile, "MSA file") != 0:
        return 1

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    # detect the format of mpa files, the old format

    MPA2MSA(infile, output_format, fpout)

    myfunc.myclose(fpout)


#}}}
if __name__ == '__main__' :
    sys.exit(main())

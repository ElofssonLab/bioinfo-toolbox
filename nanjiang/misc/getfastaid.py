#!/usr/bin/env python
# get the idlist from the annotation line in the fasta file
import sys,re,os
import myfunc 

BLOCK_SIZE=100000

usage="""
Usage:  getfastaid.py [Options] [-i] fastafile
Options:
  -i         <file>      input file
  -o         <file>      outputfile
  -printanno             Print also annotation line after seqid, tab delimited
  -m|--method 1|2        method for reading, default = 2
  -bs|--block-size int   size for blocks when reading file, default = 50000
  -h|--help              print this help message and exit
Created 2010-08-20, updated 2012-06-07, Nanjiang
"""

def PrintHelp():
    print usage

def GetFastaID(infile, fpout):#{{{
    fpin = open(infile, "r")
    line = fpin.readline()
    while line:
        line = line.rstrip('\n').strip()
        if line and line[0] == ">":
            idd = myfunc.GetSeqIDFromAnnotation(line)
            print >> fpout, idd
        line = fpin.readline()
    fpin.close()
    return 0
#}}}
def GetFastaID2(infile, fpout):#{{{
# The faster version
    isPrintAnnoLine = g_params['isPrintAnnoLine']
    fpin = open(infile, "r")
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
                    line = line.lstrip(">").rstrip("\n")
                    idd=myfunc.GetSeqIDFromAnnotation(line)
                    if not isPrintAnnoLine:
                        print >> fpout, idd
                    else:
                        fpout.write("%s\t%s\n"%(idd, line))
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
                    line = line.lstrip(">").rstrip("\n")
                    idd=myfunc.GetSeqIDFromAnnotation(line)
                    if not isPrintAnnoLine:
                        print >> fpout, idd
                    else:
                        fpout.write("%s\t%s\n"%(idd, line))
                    beg=end
                else:
                    brokenAnnoLine=buff[beg:]
                    break
            else:
                break

        buff = fpin.read(BLOCK_SIZE)
    fpin.close()
#}}}

def main(g_params):#{{{
    # Check argv
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile = ""
    infile=""
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
                return 1
            elif sys.argv[i] == "-i" or sys.argv[i] == "--infile":
                infile=sys.argv[i+1]
                i = i + 2
            elif sys.argv[i] in ["-printanno", "-pa", "--printanno", "--pa"]:
                g_params['isPrintAnnoLine'] = True
                i = i + 1
            elif sys.argv[i] in ["-o", "--o", "-outfile", "--outfile"]:
                outfile = sys.argv[i+1]
                i = i + 2
            elif sys.argv[i] == "-m" or sys.argv[i] == "--method" or sys.argv[i] == "-method":
                method=int(sys.argv[i+1])
                if method < 1 or method > 2:
                    print >> sys.stderr,"Error! method should be 1 or 2"
                    return 1
                i = i + 2
            elif sys.argv[i] == "-bs" or sys.argv[i] == "--block-size" or sys.argv[i] == "-block-size":
                BLOCK_SIZE=int(sys.argv[i+1])
                if BLOCK_SIZE < 0:
                    print >> sys.stderr,"Error! BLOCK_SIZE should >0"
                    return 1
                i = i + 2
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i])
                return 1
        else:
            infile = sys.argv[i]
            i+=1
           

    if infile == "":
        print >> sys.stderr,"Error! Input file not set."
        return 1
    elif not os.path.exists(infile):
        print >> sys.stderr,"Error! input file %s does not exist." %(infile)
        return 1


    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    if method == 1:
        GetFastaID(infile, fpout)
    else:
        GetFastaID2(infile, fpout)
    myfunc.myclose(fpout)

#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['isPrintAnnoLine'] = False
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

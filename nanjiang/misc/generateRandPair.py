#!/usr/bin/env python
# 
import os,sys
import random
import myfunc
progname =  os.path.basename(sys.argv[0])
usage="""
usage: %s [ID [ID...]] [-l IDLISTFILE]
  
Options:
  -l  LISTFILE   Set list file
  -o      FILE   Output to FILE
  -m 0|1         Method for generating pairs, (default: 0)
  -maxpair INT   Set maximum number of pairs to output, (default: 10Mb)
  -seed    INT   Set random seed, default is set by time
  -q             quiet mode
  -h, --help     print this help message and exit

Created 2011-10-20, updated 2013-04-23, Nanjiang Shu 
"""%(progname)

def PrintHelp():
    print usage


def main():#{{{
    argv = sys.argv
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        sys.exit()

    max_numpair = 10*1000*1000
    isQuiet = False
    rand_seed = None
    idList = []
    idListFile = ""
    outfile=""
    method = 0

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            idList.append(argv[i])
            isNonOptionArg=False
            i += 1
        elif sys.argv[i] == "--":
            isNonOptionArg=True
            i += 1
        elif sys.argv[i][0] == "-":
            if sys.argv[i] in [  "-h" , "--help"]:
                PrintHelp()
                sys.exit()
            elif sys.argv[i] in [ "-o" , "--o", "-outfile" , "--outfile"]:
                outfile, i = myfunc.my_getopt_str(argv,i)
            elif argv[i] in ["-m", "--m", "-method", "--method"]:
                method, i = myfunc.my_getopt_int(argv,i)
            elif sys.argv[i] in [ "-l" , "--l", "-listfile" , "--listfile"]:
                idListFile, i = myfunc.my_getopt_str(argv,i)
            elif sys.argv[i]  in [ "-maxpair" , "--maxpair"]:
                max_numpair, i = myfunc.my_getopt_int(argv,i)
            elif sys.argv[i]  in [ "-seed" , "--seed"]:
                rand_seed, i = myfunc.my_getopt_int(argv,i)
            elif sys.argv[i] == "-q":
                isQuiet=True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", sys.argv[i]
                return 1
        else:
            idList.append(argv[i])
            i += 1

    if idListFile != "":
        idList += myfunc.ReadIDList(idListFile)

    numseqid = len(idList)
    if numseqid <= 0:
        print >> sys.stderr, "List file is empty."
        return 1
    elif numseqid < 2:
        print >> sys.stderr, "Too few items. At least 2 are required."
        return 1

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)

    if method == 0:
        pairlist = myfunc.GenerateRandomPair(len(idList), max_numpair,
                rand_seed)
    elif method == 1:
        pairlist = myfunc.GenerateRandomPair_no_repeat_use(len(idList),
                max_numpair, rand_seed)

    for pair in pairlist:
        print >> fpout, "%s %s" %(idList[pair[0]], idList[pair[1]])

    myfunc.myclose(fpout)
    return 0
#}}}
if __name__ == '__main__' :
    exit(main())

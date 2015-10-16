#!/usr/bin/python
# Description:
import os
import sys
import myfunc
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

usage_short="""
Usage: %s FILE -class classFile [-outpath DIR]
"""%(progname)

usage_ext="""
Description:
    classify idlistfile

    format of class file: ID class

OPTIONS:
  -outpath DIR  Set output dir, default: the same as input file
  -q            Quiet mode
  -h, --help    Print this help message and exit

Created 2013-07-01, updated 2013-07-01, Nanjiang Shu 
"""
usage_exp="""
Examples:
"""

def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}

def ReadClassDict(infile):
    try:
        fpin = open(infile, "r")
        lines = fpin.readlines()
        fpin.close()
        id2ClassDict = {}
        classSet = set([])
        for line in lines:
            if not line or line[0] == "#":
                continue
            strs = line.split()
            if len(strs) == 2:
                id2ClassDict[strs[0]] = strs[1]
                classSet.add(strs[1])
        return (id2ClassDict, list(classSet))
    except IOError:
        print >> sys.stderr, "failed to read infile %s"%(infile)
        return ({}, [])

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = ""
    infile = ""
    classfile = ""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            infile = argv[i]
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
                (outpath, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-class", "--class"]:
                (classfile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            infile = argv[i]
            i += 1
    if myfunc.checkfile(infile) != 0:
        return 1
    if myfunc.checkfile(classfile, "Class File") != 0:
        return 1
    if outpath == "":
        outpath = os.path.dirname(infile)
        if outpath == "":
            outpath = "."
    (id2ClassDict, classList) = ReadClassDict(classfile)
    idList = myfunc.ReadIDList(infile)
    rootname = os.path.basename(os.path.splitext(infile)[0])
    ext = os.path.splitext(infile)[1]

    fpoutList = {}
    for i in range(len(classList)):
        outfile = outpath + os.sep + rootname + ".%s"%classList[i] + ext
        fpoutList[classList[i]] = open(outfile, "w")

    for idd in idList:
        try:
            cls = id2ClassDict[idd]
        except:
            print >> sys.stderr, "id %s not in classDict"%idd
            continue
        fpoutList[cls].write("%s\n"%idd)


    for i in range(len(classList)):
        fpoutList[classList[i]].close()

#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

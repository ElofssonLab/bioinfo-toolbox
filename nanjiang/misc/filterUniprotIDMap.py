#!/usr/bin/env python
# Description:
# This scirpt is not working well. 
import os
import sys
import myfunc
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))
usage = """
usage:  %s FILE [-o OUTFILE]
Description:

Options:
  -o      FILE    Output the result to file
  -q              Quiet mode
  -h, --help      Print this help message and exit

Created 2013-02-06, updated 2013-02-06, Nanjiang Shu 
"""%(progname)

def PrintHelp():
    print usage

def GetUniprotIDFromLongName(longname):
    if longname.find("_") != -1:
        strs = longname.split("_")
        if longname.find("UniRef") == 0:
            try:
                return strs[1]
            except IndexError:
                return ""
        elif strs[0] == "NX":
            try:
                return strs[1]
            except IndexError:
                return ""
        else:
            try:
                if len(strs) == 2 and len(strs[0]) == 6:
                    return strs[0]
            except IndexError:
                return ""
    elif longname.find(":") != -1:
        strs = longname.split(":")
        if len(strs) == 2 and strs[0].isdigit() and len(strs[1]) == 6:
            return strs[1]
        else:
            return ""

def FilterUniprotIDMap(infile, fpout):
    hdl = myfunc.ReadLineByBlock(infile)
    if not hdl:
        return 1
    lines = hdl.readlines()
    while lines != None:
        for line in lines:
            isIgnore = False
            strs = line.split("\t")
            try: 
                if strs[2].find(strs[0]) != -1:
                    uniprotid = GetUniprotIDFromLongName(strs[2])
                    if uniprotid != "":
                        if uniprotid == strs[0]:
                            isIgnore = True
                        else:
                            print >> sys.stderr, "Error\t", line
                    else:
                        print >> sys.stderr, "Null\t", line
                if not isIgnore:
                    print >> fpout, line
            except IndexError:
                print >> sys.stderr, "IndexError\t",line
        lines = hdl.readlines()

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    infile = ""
    outfile = ""

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
            elif argv[i] in ["-o", "--o", "-outfile", "--outfile"]:
                outfile = argv[i+1]
                i += 2
            elif argv[i] in ["-q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            infile = argv[i]
            i += 1

    if infile == "":
        print >> sys.stderr, "infile not set. exit"
        return 1
    elif not os.path.exists(infile):
        print >> sys.stderr, "infile %s does not exist. exit"%(infile)
        return 1
    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)

    FilterUniprotIDMap(infile, fpout)

    myfunc.myclose(fpout)

#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

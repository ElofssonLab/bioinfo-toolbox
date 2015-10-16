#!/usr/bin/python
# Description:
import os
import sys
import myfunc
#import bioinformatics
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

usage_short="""
Usage: %s [ID [ID...]] -datapath DIR -ext STR -l IDLISTFILE -o OUTFILE
"""%(progname)

usage_ext="""
Description:
   generate file list given idlist, datafiles may stored under 
   subfolders and in that case, id2pathmap.txt will be read

OPTIONS:
  -o  OUTFILE      output the pairaln file
  -l LISTFILE      Set the listfile
  -ext STR         file extention
  -q               Quiet mode
  -h, --help       Print this help message and exit

Created 2013-06-14, updated 2013-06-14, Nanjiang Shu 
"""
usage_exp="""
Examples:
    %s -datapath hhalign -ext .hhr -l idlistfile.txt -o hhrfilelist.txt
"""%(progname)

def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}
def ID2File(idd, datapath, id2pathMapDict, ext):#{{{
    try:
        subdir = id2pathMapDict[idd]
        filename = "%s%s%s%s%s%s"%(datapath, os.sep, subdir, os.sep, idd, ext)
    except KeyError:
        filename = "%s%s%s%s"%(datapath, os.sep, idd, ext)
    if not os.path.exists(filename):
        print >> sys.stderr, "File not found for ID %s"%(idd)
        return ""
    else:
        return filename
#}}}
def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile = ""
    datapath = ""
    idListFile = ""
    idList = []
    ext = ""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            idList.append(argv[i])
            isNonOptionArg = False
            i += 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-o", "--o", "-outfile"]:
                (outfile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-ext", "--ext"]:
                (ext, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-datapath", "--datapath"]:
                (datapath, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-l", "--l"] :
                (idListFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            idList.append(argv[i])
            i += 1

    if ext == "":
        print >> sys.stderr, "file extension not set. exit"
        return 1
    if datapath == "":
        print >> sys.stderr, "datapath not set. exit"
        return 1
    elif not os.path.exists(datapath):
        print >> sys.stderr, "datapath %s does not exist. exit"%(datapath)
        return 1

    if idListFile != "":
        idList += myfunc.ReadIDList(idListFile)
    if len(idList) < 1:
        print >> sys.stderr, "No input set. exit"
        return 1

    fpout = sys.stdout
    if outfile != "":
        fpout = myfunc.myopen(outfile, sys.stdout, "w", False)

    id2pathmapfile = "%s%s%s"%(datapath, os.sep, "id2pathmap.txt")
    id2pathMapDict = myfunc.ReadIDPathMapDict(id2pathmapfile)
    for idd in idList:
        filename = ID2File(idd, datapath, id2pathMapDict, ext)
        if filename != "":
            print >> fpout, filename

    if outfile != "":
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

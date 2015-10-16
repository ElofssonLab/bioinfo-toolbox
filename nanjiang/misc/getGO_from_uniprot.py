#!/usr/bin/python
# File name: getGO_from_uniprot.py  
# Description:
import os
import sys
import myfunc
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

usage_short="""
Usage: %s -l uniprotACList -uniprotdb DBNAME [-o OUTFILE]
             ID [ID ...]
"""%(progname)

usage_ext="""
Description:

OPTIONS:
  -o OUTFILE     Output the result to OUTFILE
  -l LISTFILE    Set the listfile
  -uniprotdb STR UniprotDBname
  -q             Quiet mode
  -h, --help     Print this help message and exit

Created 2013-09-02, updated 2013-09-02, Nanjiang Shu 
"""
usage_exp="""
Examples:
"""

def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}
def WriteGOInfo(seqid, goinfo, fpout):
    fpout.write("%s"%(seqid))
    for item in ["function", "process", "location"]:
        for msg in goinfo[item]:
            fpout.write("\t%s"%msg)
    fpout.write("\n")


def GetGOInfoFromUniprotData(data):#{{{
    lines = data.split("\n")
    goinfo = {}
    goinfo['process'] = []
    goinfo['function'] = []
    goinfo['location'] = []
    for line in lines:
        if line[0:7] == "DR   GO":
            msg = line[9:]
            strs = msg.split("; ")
            if len(strs) >= 2:
                if strs[1][0:2] == "P:": #process
                    goinfo['process'].append(msg)
                elif strs[1][0:2] == "C:": #location
                    goinfo['location'].append(msg)
                elif strs[1][0:2] == "F:": #location
                    goinfo['function'].append(msg)
                else:
                    print >> sys.stderr, "bad GO:%s"%(msg)
    return goinfo
#}}}

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = "./"
    outfile = ""
    idListFile = ""
    uniprotDBname = ""
    idList = []

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
            elif argv[i] in ["-outpath", "--outpath"]:
                (outpath, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-l", "--l"] :
                (idListFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-uniprotdb", "--uniprotdb"] :
                (uniprotDBname, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            idList.append(argv[i])
            i += 1

    if idListFile != "":
        idList += myfunc.ReadIDList(idListFile)

    if uniprotDBname == "":
        print >> sys.stderr, "uniprotdb not set"
        return 1
    uniprotdbfile = "%s0.db"%uniprotDBname
    if myfunc.checkfile(uniprotdbfile, "uniprotdbfile") != 0:
        return 1

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)

    hdl = myfunc.MyDB(uniprotDBname)
    if hdl.failure:
        return 1

    for seqid in idList:
        data = hdl.GetRecord(seqid)
        if data != None:
            goinfo = GetGOInfoFromUniprotData(data)
            WriteGOInfo(seqid, goinfo, fpout)
    hdl.close()
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

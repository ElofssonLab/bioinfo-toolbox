#!/usr/bin/env python
# Description:
import os
import sys
import myfunc

# Format of the input file
# SEQID1 2 FAM1 FAM2
# SEQID2 1 FAM3

progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))
usage = """
Usage:  %s FILE [-idlist1 FILE] [-idlist2 FILE] [-o OUTFILE]
Description: Filter seqid2fam map file by setting restriction idlist
Options:
  -idlist1 FILE  Set restriction idlist for keys
  -idlist2 FILE  Set restriction idlist for content ids
  -o OUTFILE     Output the result to OUTFILE
  -q             Quiet mode
  -h, --help     Print this help message and exit

Created 2013-03-12, updated 2013-03-12, Nanjiang Shu
"""%(progname)

def PrintHelp():
    print usage

def Filter_seqid2fam_map(infile, keyIDSet, contentIDSet, isKeyIDSet,
        isContentIDSet, fpout):
    hdl = myfunc.ReadLineByBlock(infile)
    if hdl.failure:
        return 1
    lines = hdl.readlines()
    while lines != None:
        for line in lines:
            if line:
                strs = line.split()
                try:
                    key = strs[0]
                    num = int(strs[1])
                    idlist = strs[2:]
                    tmp_idlist = []
                    if (not isKeyIDSet) or (key in keyIDSet):
                        for idd in idlist:
                            if (not isContentIDSet) or (idd in contentIDSet):
                                tmp_idlist.append(idd)
                        if len(tmp_idlist) > 0:
                            fpout.write("%s %d"%(key, len(tmp_idlist)))
                            for idd in tmp_idlist:
                                fpout.write(" %s"%(idd))
                            fpout.write("\n")
                except (IndexError):
                    msg = "Error in infile %s with line \"%s\""
                    print >> sys.stderr, msg%(infile, line)
                    return 1
        lines = hdl.readlines()
    hdl.close()
    return 0

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    infile = ""
    outfile = ""
    keyIDListFile = ""
    contentIDListFile = ""
    isKeyIDSet = False
    isContentIDSet = False

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
            elif argv[i] in ["-o", "--o"]:
                outfile = argv[i+1]
                i += 2
            elif argv[i] in ["-idlist1", "--idlist1"] :
                keyIDListFile = argv[i+1]
                isKeyIDSet = True
                i += 2
            elif argv[i] in ["-idlist2", "--idlist2"] :
                contentIDListFile = argv[i+1]
                isContentIDSet = True
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
        print >> sys.stderr, "infile not set. Exit"
        return 1
    elif not os.path.exists(infile):
        print >> sys.stderr, "infile %s does not exist. Exit"%(infile)
        return 1

    keyIDSet = {}
    contentIDSet = {}
    if keyIDListFile != "":
        keyIDSet = set(myfunc.ReadIDList(keyIDListFile))
    if contentIDListFile != "":
        contentIDSet = set(myfunc.ReadIDList(contentIDListFile))

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)

# if idlist1 and idlist2 is empty, output nothing
    if isKeyIDSet or isContentIDSet:
        Filter_seqid2fam_map(infile, keyIDSet, contentIDSet, isKeyIDSet,
                isContentIDSet, fpout)
    else:
        os.system("cat %s"%(infile))

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

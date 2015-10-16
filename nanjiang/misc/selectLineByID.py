#!/usr/bin/env python
# Description:
import os
import sys
import myfunc
progname = os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

# ChangeLog 2013-12-06 
#   enabled selection of lines by pairs of IDs, this is done by using multiple
#   specified selected field

usage_short = """
Usage: %s -i infile [-l LISTFILE] [-o OUTFILE] [-q]
       %s ID [ID ...]
"""%(progname, wspace)

usage_ext="""
Description: Select lines according to the given idlist
             White space is not allowed in the id

OPTIONS:
  -i FILE  Set input table file
  -l FILE  Set the idListFile, one line per record
           If multiple IDs are supplied in one line, all IDs will be matched
           On application is for pair match
  -m, -methodid 0|1|2|3
           Method to get id from the line record, (default: 0)
           0: id = firstword
           1: id = firstword delimited by [; ]
           2: id = GetSeqIDFromAnnotation(line)
           3: using specified field
  -se, -selfield INT:
           use field N
  -q       Quiet mode
  -h       Print this help message and exit

Created 2012-05-31, updated 2013-12-06, Nanjiang Shu
"""
usage_exp="""
Examples:
    # select lines in table.txt by IDs supplied in file idlist
    %s -i table.txt -l idlist -o selected.table.txt

    # select line in table.txt by IDs in the second column
    %s -i table.txt -l idlist -se 2 -o selected.table.txt

    # select lines in table.txt by pairs of ID supplied in file t1.pairlist
    %s -i table.txt -l t1.pairlist -se 1 -se 2 -o selected.table.txt
"""%(progname, progname, progname)

def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}

def SelectLineByID(infile, idListSet, fpout):#{{{
    hdl = myfunc.ReadLineByBlock(infile)
    if hdl.failure:
        return 1
    method_getid = g_params['method_getid']
    sel_field_list = g_params['sel_field_list']
    if method_getid == 3:
        if len(sel_field_list) == 0:
            sel_field = 0
        elif len(sel_field_list) == 1:
            sel_field = sel_field_list[0]


    lines = hdl.readlines()
    while lines != None:
        for line in lines:
            if not line or line[0] == "#":
                fpout.write("%s\n"%line)
            else:
                try:
                    if method_getid == 0:
                        idd = line.split(None, 1)[0]
                    elif method_getid == 1:
                        idd = (line.split(None, 1)[0]).partition(";")[0]
                    elif method_getid == 2:
                        idd = myfunc.GetSeqIDFromAnnotation(line)
                    elif method_getid == 3:
                        if len(sel_field_list) < 2:
                            idd = line.split()[sel_field-1]
                        else:
                            strs = line.split()
                            tmpli = []
                            for ff in sel_field_list:
                                tmpli.append(strs[ff-1])
                            idd = tuple(tmpli)
                    else:
                        print method_getid
                except (IndexError):
                    print >> sys.stderr, ("Bad line \"%s\"\n"%line)
                if idd in idListSet:
                    fpout.write("%s\n"%line)
        lines = hdl.readlines()
    hdl.close()
    return 0
#}}}
def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile = ""
    idListFile = ""
    idList = []
    infile = ""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            strs = argv[i].split()
            if len(strs) == 1:
                idList.append(argv[i])
            else:
                idList.append(tuple(strs))
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
            elif argv[i] in ["-l", "--l"] :
                (idListFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-m", "--m", "-methodid", "--methodid"] :
                (g_params['method_getid'], i) = myfunc.my_getopt_int(argv, i)
            elif argv[i] in ["-se,", "--se", "-selfield", "--selfield"] :
                (tmpint, i) = myfunc.my_getopt_int(argv, i)
                g_params['sel_field_list'].append(tmpint)
                g_params['method_getid'] = 3
            elif argv[i] in ["-i", "--i"] :
                (infile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            strs = argv[i].split()
            if len(strs) == 1:
                idList.append(argv[i])
            else:
                idList.append(tuple(strs))
            i += 1

    if infile == "":
        print >> sys.stderr, "infile not set"
        print >> sys.stderr, usage_short
        return 1
    elif not os.path.exists(infile):
        print >> sys.stderr, "infile %s does not exist"%(infile)
        return 1

    if idListFile != "":
        hdl = myfunc.ReadLineByBlock(idListFile)
        if not hdl.failure:
            lines = hdl.readlines()
            while lines != None:
                for line in lines:
                    if line:
                        strs = line.split()
                        if len(strs) == 1:
                            idList.append(strs[0])
                        else:
                            idList.append(tuple(strs))
                lines = hdl.readlines()
            hdl.close()


    if len(idList) < 1:
        print >> sys.stderr, "ID not set"
        print >> sys.stderr, usage_short
    idListSet = set(idList)
    fpout = myfunc.myopen(outfile, sys.stdout, "w", False) 
    status = SelectLineByID(infile, idListSet, fpout)
    myfunc.myclose(fpout)
    return status
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['method_getid'] = 0
    g_params['sel_field_list'] = []
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

#!/usr/bin/env python
# Description:
# Check whether the fasta sequence contains duplicated sequences
import os
import sys
import myfunc
import md5
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

usage_short="""
Usage: %s FILE [FILE ...] [-o OUTFILE]
"""%(progname)

usage_ext="""
Description:
    Check whether the fasta sequence contains duplicated sequences

OPTIONS:
  -o OUTFILE    Output the result to OUTFILE
  -m id|seq     Set method, by unique seqid or unique seq, (default: id)
  -md5 yes|no   Whether use md5 encoding, (default: yes)
                When enabled, it may use fewer memory if -m = seq
  -l LISTFILE   Set the listfile
  -h, --help    Print this help message and exit

Created 2014-12-12, updated 2014-12-12, Nanjiang Shu 
"""
usage_exp="""
Examples:
    %s test.fa -m id    #check whether the test.fa has duplicated seqIDs

""" %(progname)

def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}

def IsUniqueSeq(infile, method, isUseMD5):#{{{
    """
    return value 
    yes     :  1
    no      :  0
    failed  : -1
    """
    hdl = myfunc.ReadFastaByBlock(infile)
    if hdl.failure:
        return -1

    myset = set([])

    recordList = hdl.readseq()
    while recordList != None:
        for rd in recordList:
            if method == "id":
                key = rd.seqid
            elif method == "seq":
                if isUseMD5:
                    key = md5.new(rd.seq).digest()
                else:
                    key = rd.seq
            if key in myset: # duplicated
                return 0     # not unique
            myset.add(key)
        recordList = hdl.readseq()
    hdl.close()
    return 1  #unique
#}}}
def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile = ""
    fileListFile = ""
    fileList = []

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            fileList.append(argv[i])
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
            elif argv[i] in ["-m", "--m"]:
                (g_params['method'], i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-md5", "--md5"]:
                (tmpstr, i) = myfunc.my_getopt_str(argv, i)
                if tmpstr.lower() == "yes":
                    g_params['isUseMD5'] = True
                elif tmpstr.lower() == "no":
                    g_params['isUseMD5'] = False
                else:
                    print >> sys.stderr, "Bad syntax. option -md5 must be followed by yes or no"
                    return 1
            elif argv[i] in ["-l", "--l"] :
                (fileListFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            fileList.append(argv[i])
            i += 1


    if fileListFile != "":
        fileList += myfunc.ReadIDList(fileListFile)

    if len(fileList) < 1:
        print >> sys.stderr, "%s: no input file is set. exit"%(sys.argv[0])

    if not g_params['method'] in ["id","seq"]:
        print >> sys.stderr, "%s: bad method \"%s\""%(sys.argv[0], g_params['method'])


    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    for i in xrange(len(fileList)):
        status =  IsUniqueSeq(fileList[i], g_params['method'], g_params['isUseMD5'])
        if status >= 0:
            if status == 1:
                yes_or_no =  "yes"
            else:
                yes_or_no = "no"
            print >> fpout, "%s\t%s" %(fileList[i], yes_or_no)
        else:
            print >> sys.stderr, "%s: Failed to read file %s" %(sys.argv[0], fileList[i])

    myfunc.myclose(fpout)
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['method'] = "id"
    g_params['isUseMD5'] = False
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

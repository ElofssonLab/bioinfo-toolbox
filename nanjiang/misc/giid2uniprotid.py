#!/usr/bin/env python

import string
import os
import sys
sys.path.append("%s/wk/MPTopo/src"%(os.environ['DATADIR3']))
import myfunc
import libtopologycmp as lcmp
progname =  os.path.basename(sys.argv[0])


import urllib,urllib2

url = 'http://www.uniprot.org/mapping/'


def GIID2UniprotID(params, fpout):
    data = urllib.urlencode(params)
    request = urllib2.Request(url, data)
    contact = "" # Please set your email address here to help us debug in case of problems.
    request.add_header('User-Agent', 'Python %s' % contact)
    response = urllib2.urlopen(request)
    page = response.read(200000)
    print >> fpout, page

usage_short="""
Usage: %s ID [ID...] [-l IDLISTFILE] [-o OUTFILE]
"""%(progname)

usage_ext="""
Description:
    Get the species name from fasta file annotation line

OPTIONS:
  -h, --help    Print this help message and exit

Created 2013-10-29, updated 2013-10-29, Nanjiang Shu 
"""
usage_exp="""
Examples:
"""
def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}
def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile = ""
    idListFile = ""
    idList = []

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            idList.append(argv[i])
            i += 1
            isNonOptionArg = False
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-o", "--o", "-outfile"]:
                (outfile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-l", "--l", "-list"]:
                (idListFile, i) = myfunc.my_getopt_str(argv, i)
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            idList.append(argv[i])
            i += 1

    if idListFile != "":
        idList += myfunc.ReadIDList(idListFile)

    numID = len(idList)

    if numID < 1:
        print >> sys.stderr, "No ID set. exit"
        return 1

    params = {}
    params ['from'] = 'P_GI'
    params ['to'] = 'ID' # to uniprot id
    params ['format'] = 'tab'
    params['query'] = " ".join(idList)
    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)

    GIID2UniprotID(params, fpout)
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

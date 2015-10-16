#!/usr/bin/env python

import os
import sys
import urllib
import urllib2
import time
#import uniprot

# import socket
# nodename = socket.gethostname()

uniprot_url = "http://www.uniprot.org/uniprot/"
srcpath = "/data3/wk/MPTopo/src"
binpath = "/data3/bin"

sys.path.append(srcpath)
import myfunc

progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))
usage = """
Usage:  {} ID [ID...] [-o OUTFILE]
Description:
Options:
  -l LISTFILE   Set IDList File
  -o OUTFILE    Output the result to OUTFILE
  -sp           Show progress
  -q            Quiet mode
  -h, --help    Print this help message and exit

Created 2013-02-11, updated 2013-02-11, Nanjiang Shu 
"""

def PrintHelp(progname, numparallel):
    print usage.format(progname, numparallel)

def UniprotID2Seq_HTTP(uniprot_url, idList, addedIDset, fpout):#{{{
    cntRetrieved = 0
    cnt = 0
    isShowProgress = g_params['isShowProgress']
    for uniprotid in idList:
        cnt += 1
        if isShowProgress:
            if (cnt%100) == 1:
                print >> sys.stderr, "Running %d ..."%(cnt)
        if uniprotid in addedIDset:
            continue
        try:
            url = uniprot_url + uniprotid + ".fasta"
            filename = uniprotid + ".fasta"
            data = urllib2.urlopen(url).read()
            if data:
#                fpout.write(data)
                (tmpid, tmpanno, tmpseq
                        ) = myfunc.ExtractFromSeqWithAnno(data)
                fpout.write(">%s %s\n"%(uniprotid, tmpanno))
                fpout.write("%s\n"%(tmpseq))
                addedIDset.add(uniprotid)
                cntRetrieved += 1
        except urllib2.HTTPError, e:
            print >> sys.stderr, "HTTPError for ID %s"%(uniprotid)
#             print e.code
#             print e.message
#             print e.headers
    return cntRetrieved
#}}}
        

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)

    numparallel = 5

    outfile = ""
    listfile = ""
    idList = []
    i = 1

    if numArgv < 2:
        PrintHelp(progname, numparallel)
        return 1
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
                PrintHelp(progname, numparallel)
                return 1
            elif argv[i] in ["-o", "--o"] :
                outfile = argv[i+1]
                i += 2
            elif argv[i] in ["-l", "--l"] :
                listfile = argv[i+1]
                i += 2
            elif argv[i] in ["-p", "--p"] :
                numparallel = int(argv[i+1])
                i += 2
            elif argv[i] in ["-sp", "--sp"]:
                g_params['isShowProgress'] = True
                i += 1
            elif argv[i] in ["-q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            idList.append(argv[i])
            i += 1

    if listfile != "":
        try:
            fp = open(listfile,"r")
            idList += fp.read().split()
            fp.close()
        except IOError:        
            print >> sys.stderr, "file %s does not exist." %listfile
    numID = len(idList)

    if numID < 1:
        print >> sys.stderr, "No ID set. Exit"
        return 1
    
    outdir = os.path.dirname(outfile)
    if outdir != "" and not os.path.exists(outdir):
        os.system("mkdir -p %s"%(outdir))

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    idList = myfunc.uniquelist(idList)
    addedIDset = set([])

    isShowProgress = g_params['isShowProgress']
    if isShowProgress:
        print >> sys.stderr, "%d ID to be retrieved" %(len(idList))

    

    start = time.time()
    numRetrieved = UniprotID2Seq_HTTP(uniprot_url, idList, addedIDset, fpout)
    end = time.time()

    msg =  "Retrieving %d (out of %d) sequences costs %.3fs seconds"
    print  msg%(numRetrieved, len(idList), (end-start))
    
    numRemain = len(idList) - numRetrieved
    if numRemain > 0:
        msg = "%d records are failed to retrieve. They are:"
        print msg%(numRemain) 
        for seqid in idList:
            if not seqid in addedIDset:
                print seqid

    myfunc.myclose(fpout)


#}}}
def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['isShowProgress'] = False
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

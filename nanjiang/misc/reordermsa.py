#!/usr/bin/env python
import sys
import os
import myfunc

# ChangeLog 2012-08-28
#   Can output only annotation line or fasta format, annotation line is one
#   line per record

BLOCK_SIZE=100000

usage="""
Usage:  reordermsa.py [Options] msafile -orderlist orderfile
Options:
  -o FILE   Output the result to file, (default: stdout)
  -of STR   Set the output format, anno or fasta, (default: fasta)
            anno: output only the annotation line
            fasta: output the fasta file

Created 2012-03-19, updated 2012-08-28, Nanjiang Shu 

Examples:
    reordermsa.py t1.topomsa.fa -orderlist t1-listorder.txt > t1.reordered.topomsa.fa
    reordermsa.py t1.topomsa.fa -orderlist t1-listorder.txt -of anno > t1.reordered.topomsa.anno
"""

def PrintHelp():
    print usage

def ReadOrderList(infile):#{{{
    try:
        fpin = open(infile, "r")
        lines = fpin.readlines()
        fpin.close()
        orderlist = []
        for line in lines:
            line = line.strip()
            if line and line != "0.1":
                orderlist.append(line.strip())
        return orderlist
    except IOError:
        print >> sys.stderr, "Failed to read orderlist file ", infile
        return []
#}}}
def main(g_params):

    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outFile = ""
    orderlistfile = ""
    msafile = ""
    outformat = "fasta" # fasta or anno

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            msafile = sys.argv[i]
            isNonOptionArg=False
            i = i + 1
        elif sys.argv[i] == "--":
            isNonOptionArg=True
            i = i + 1
        elif sys.argv[i][0] == "-":
            if sys.argv[i] ==  "-h" or  sys.argv[i] == "--help":
                PrintHelp()
                return 1
            elif sys.argv[i] in [ "-o", "--o"] :
                outFile=sys.argv[i+1]
                i = i + 2
            elif sys.argv[i] == "-orderlist" or sys.argv[i] == "--orderlist":
                orderlistfile = sys.argv[i+1]
                i = i + 2
            elif sys.argv[i] == "-msafile" or sys.argv[i] == "--msafile":
                msafile = sys.argv[i+1]
                i = i + 2
            elif sys.argv[i] in ["-of", "--of", "-outformat", "--outformat"]:
                outformat = sys.argv[i+1].lower()
                i += 2
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i])
                return 1
        else:
            msafile = sys.argv[i]
            i+=1

    if not outformat in ["anno", "fasta"]:
        print >> sys.stderr, "Unrecognized outformat \"%s\","%(
                outformat) + " should be either \"anno\" or \"fasta\"."
        return 1

    if orderlistfile == "":
        print >> sys.stderr, "orderlist file not set. Exit"
        return 1
    if msafile == "":
        print >> sys.stderr, "msafile not set. Exit"
    orderList = ReadOrderList(orderlistfile)  
    (idList, annoList, seqList) = myfunc.ReadFasta(msafile)

    if len(orderList) > 0  and len(idList) > 0:
        fpout = sys.stdout
        fpout = myfunc.myopen(outFile, sys.stdout, "w", False)

        seqDict = {}
        annoDict = {}
        numSeq = len(idList)
        for i in xrange(numSeq):
            annoDict[idList[i]] = annoList[i]
        if outformat != "anno":
            for i in xrange(numSeq):
                seqDict[idList[i]] = seqList[i]
        for sid in orderList:
            if sid in annoDict:
                fpout.write(">%s\n"%annoDict[sid])
                if outformat != "anno":
                    fpout.write("%s\n"%seqDict[sid])
            else:
                print >> sys.stderr, "seqid %s not in msafile %s"%(
                        sid, msafile)
        myfunc.myclose(fpout)

    return 0

if __name__ == '__main__' :
    g_params = {}
    sys.exit(main(g_params))

#!/usr/bin/env python
# Filename: topo2TMfrag.py 
# Description:
#   From the topology file output by mySCAMPI_run.pl, get the peptides of TM
#   fragments.
#   The output result can be used by calc_dG.pl to calculate the DG (apparent
#   free energy) values. 
# Author:
#   Nanjiang Shu    nanjiang.shu@scilifelab.se
import sys,re,os
import myfunc 

#ChangeLog 2011-10-26 
#   1. support topology with gaps
#   2. option -printid must be appened with yes or no option
#      and by default sequence id is output
#ChangeLog 2011-10-30
#   redundant topologies (recognized by seqID) in the topologyFile are ignored
#      if seqID in  processedTopoIDSet:
#           continue
#ChangeLog 2011-11-04 
#   topologyFile is read in by BLOCK reading using ReadFastaFromBuffer
#   so that the program will not be limited by the size of topology file. 
#ChangeLog 2011-11-07 
#   Read Fasta block by block updated
#   

usage="""
usage:  topo2TMfrag.py  [-i] topologyFile -f fastaFile

Description:
    Output amino acid fragments of all TM regions in the topologyFile

OPTIONS:
  -i        FILE    Input file
  -o        FILE    Output file
  -printid  y|n     Print sequence ID of each fragment, (default: yes)
  -h, --help        Print this help message and exit

Created 2010-09-03, updated 2011-11-07,  Nanjiang Shu
"""

GAP='-'
MAX_FASTA_AA_FILE_SIZE = 2*1024*1024*1024;  # 2GB
BLOCK_SIZE=100000

def PrintHelp():
    print usage
def GetTMPosition_gapless(topo):#{{{
# Get the position of TM helices given the topology (without gaps)
# The return value is a list of 2-tuples
# [ (beg, end), (beg, end)...]
    posTM=[]
    m=re.finditer("(M+)",topo)
    for i in m:
        posTM.append((i.start(0), i.end(0)))
    return posTM
#}}}
def Topo2TMFrag(idList, topoList, aaSeqDict, processedTopoIDSet, fpout):#{{{
    for iSeq in xrange (len(topoList)):
        seqID= idList[iSeq]
        if seqID not in aaSeqDict:
            print >> sys.stderr, "Warning! not aaSeq found for ID %s"%seqID
            continue
        if seqID in  processedTopoIDSet:
            continue
        aaSeq=aaSeqDict[seqID]
        topo=topoList[iSeq]
        gapless_topo=topo.replace(GAP, '')
        #gapless_topo = ''.join (ch for ch in topo if ch not in GAP)
        #gapless_topo = ''.join (topo.split(GAP))
        if len(aaSeq) == len(gapless_topo):
            processedTopoIDSet.add(seqID)
            posTM = GetTMPosition_gapless(gapless_topo)
            for (b,e) in posTM:
                if isPrintSeqID:
                    print >> fpout, "%s \t %s" % (seqID,aaSeq[b:e])
                else:
                    print >> fpout, "%s" % (aaSeq[b:e])
        else:
            print >> sys.stderr, "length of aaSeq (%d) and topology (%d) does not match for ID %s. Ignore." %(len(aaSeq), len(gapless_topo), seqID)
    return 0
#}}}

def main():#{{{
    numArgv = len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    global isPrintSeqID
    outFile=""
    inFile=""
    fastaFile=""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            isNonOptionArg=False
            i = i + 1
        elif sys.argv[i] == "--":
            isNonOptionArg=True
            i = i + 1
        elif sys.argv[i][0] == "-":
            if sys.argv[i] ==  "-h" or  sys.argv[i] == "--help":
                PrintHelp()
                return 1
            elif sys.argv[i] == "-i" or sys.argv[i] == "--infile":
                inFile=sys.argv[i+1]
                i = i + 2
            elif sys.argv[i] == "-f" or sys.argv[i] == "--fasta":
                fastaFile=sys.argv[i+1]
                i = i + 2
            elif sys.argv[i] == "-printid" or sys.argv[i] == "--printid":
                if (sys.argv[i+1].lower())[0] == "y": 
                    isPrintSeqID=True
                else:
                    isPrintSeqID=False
                i = i + 2
            elif sys.argv[i] == "-o" or sys.argv[i] == "--outfile":
                outFile=sys.argv[i+1]
                i = i + 2
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i])
                return 1
        else:
            inFile=sys.argv[i]
            i+=1

    if inFile == "":
        print >> sys.stderr,"Error! Topology file not set."
        return 1
    if fastaFile == "":
        print >> sys.stderr,"Error!  amino acid fasta file not set."
        return 1

    fpout = sys.stdout
    if outFile != "":
        fpout = open(outFile,"w")
        if not fpout:
            print >> sys.stderr, "Failed to write to outfile %s. "%(outFile)
            print >> sys.stderr, "Reset output to stdout."
            fpout = sys.stdout
    sizeAASeqFile = os.path.getsize(fastaFile)

    if sizeAASeqFile > MAX_FASTA_AA_FILE_SIZE:
        print >> sys.stderr, ("size (%d)"%sizeAASeqFile 
                + " of fasta sequence file (%s)"%fastaFile
                + " is over the limit (%d). Exit."% MAX_FASTA_AA_FILE_SIZE)
        return 1

    (idListSeq, annotationListSeq, seqList) = myfunc.ReadFasta(fastaFile)
    if idListSeq == None:
        print >> sys.stderr, "%s exit with error."%sys.argv[0]
        return 1
    elif idListSeq < 1:
        print >> sys.stderr, ("Warning! zero aa sequences have" 
                + " been read in for file %s" %fastaFile)
    aaSeqDict={}
    for i in xrange (len(idListSeq)):
        aaSeqDict[idListSeq[i]] = seqList[i]


    fpin = open (inFile, "rb")
    if not fpin:
        print >> sys.stderr, "Failed to open input file %s"%(inFile)
        return -1
    unprocessedBuffer=""
    isEOFreached = False
    processedTopoIDSet = set([])
    while 1:
        buff = fpin.read(BLOCK_SIZE)
        if len(buff) < BLOCK_SIZE:
            isEOFreached=True
        buff = unprocessedBuffer + buff
        recordList = []
        unprocessedBuffer = myfunc.ReadFastaFromBuffer(buff,recordList, isEOFreached)
        if len(recordList) > 0: 
            idListTopo = [r[0] for r in recordList]
            topoList = [r[2] for r in recordList]
            Topo2TMFrag(idListTopo, topoList,aaSeqDict, processedTopoIDSet, fpout)
        if isEOFreached == True:
            break
    fpin.close()

    if fpout != None and fpout != sys.stdout:
        fpout.close()
#}}}
if __name__ == '__main__' :
    # Check argv
    isPrintSeqID=True
    sys.exit(main())

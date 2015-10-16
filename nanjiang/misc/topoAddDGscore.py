#!/usr/bin/env python
# better to be called by the master script  "getDGvalueTMOfTopo.sh"

# add DG scores to the topology file, 
# output two files:
#   one is the topology file adding with one line with DG scores under each TM
# segment
#   Another is a list 
#   ID  TMsegment(aa) DGscore
# dG file in the format
# ====
# ID1 FVSLGFLLIIIVPAFISCHAR       -0.530
# ID1 PYIYLGGAILAEVIGTTLMKF       1.464 
# ID1 VGTIICYCASFWLLAQTLAYI       0.739 
# ID1 AYAIWSGVGIVLISLLSWGFF       -1.146
# ====
#
#ChangeLog 2011-10-26 
#   1. support topology with gaps
#   2. input dg file should be with seqID for correct matching
#ChangeLog 2011-11-04 
#   topologyFile is read in by BLOCK reading using ReadFastaFromBuffer
#   so that the program will not be limited by the size of topology file. 
#ChangeLog 2011-11-07 
#   Read Fasta block by block updated

import sys,re,os
import myfunc 
import libtopologycmp as lcmp

BLOCK_SIZE=100000
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

usage="""
usage:  %s [-i] topologyFile -dg dgFile

Add deltaG values of TM fragments to the topology file in Fasta format
deltaG values are written below the topology sequence and enclosed by {dgscore }

  -i         FILE  Input file
  -outpath   DIR   Output path, (default: ./), output file is
                   $rootname.topowithdgscore
  -h, --help       Print this help message and exit

Created 2010-09-03, updated 2012-12-12,  Nanjiang Shu
"""%(progname)

def PrintHelp():
    print usage

def TopoAddDGscore(idListTopo, annotationListTopo, topoList, dgScoreDict, #{{{
        fpTopoDG):
    for iSeq in xrange (len(topoList)):
        seqID=idListTopo[iSeq]
        topo=topoList[iSeq]
        annoLine=annotationListTopo[iSeq]
        numTM = myfunc.CountTM(topo)
        #write topoDG
        fpTopoDG.write(">%s\n"%annoLine)
        fpTopoDG.write("%s\n"%topo)
        fpTopoDG.write("{dgscore ")
        if seqID in dgScoreDict:
            dglist = dgScoreDict[seqID]
            numDGscore = len(dglist)
            if numDGscore != numTM: 
                print >> sys.stderr, ("num DGscores for seqID %s (%d) "%(seqID,numDGscore)
                + "!= numTM (%d) for the topology. dglist = "%(numTM)), dglist
            else:
                for i in range(numTM):
                    fpTopoDG.write("%s "%dglist[i])
        fpTopoDG.write("}\n")
    return 0
#}}}

def main():#{{{
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1
    outpath="./"
    inFile=""
    dgFile=""
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
            elif sys.argv[i] == "-dg" or sys.argv[i] == "--dg":
                dgFile=sys.argv[i+1]
                i = i + 2
            elif sys.argv[i] == "-outpath" or sys.argv[i] == "--outpath":
                outpath=sys.argv[i+1]
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
    if dgFile == "":
        print >> sys.stderr,"Error!  dgFile not set."
        return 1

    os.system("mkdir -p %s" % outpath)
    dgScoreDict = lcmp.ReadDGScore(dgFile)
    if dgScoreDict == {}:
        print >> sys.stderr, "Read DG score failed. exit"
        return 1

    rootname = os.path.basename(os.path.splitext(inFile)[0])
    outFileTopoDG = outpath+ os.sep + rootname + ".topowithdgscore"

    fpTopoDG = myfunc.myopen(outFileTopoDG, None, "w", False)
    if fpTopoDG == None:
        return 1

    fpin = open (inFile, "rb")
    if not fpin:
        print >> sys.stderr, "Failed to open input file %s"%(inFile)
        return -1
    unprocessedBuffer=""
    isEOFreached = False
    while 1:
        buff = fpin.read(BLOCK_SIZE)
        if len(buff) < BLOCK_SIZE:
            isEOFreached=True
        buff = unprocessedBuffer + buff
        recordList = []
        unprocessedBuffer = myfunc.ReadFastaFromBuffer(buff,recordList, isEOFreached)
        if len(recordList) > 0: 
            idListTopo = [r[0] for r in recordList]
            annotationListTopo = [r[1] for r in recordList]
            topoList = [r[2] for r in recordList]
            TopoAddDGscore(idListTopo, annotationListTopo, topoList,
                    dgScoreDict, fpTopoDG)
        if isEOFreached == True:
            break
    myfunc.myclose(fpTopoDG)
    print  outFileTopoDG, "output"
#}}}

if __name__ == '__main__' :
    sys.exit(main())

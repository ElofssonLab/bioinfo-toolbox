#!/usr/bin/env python
# Select One TM segment by location and DG value
# Input 
#   dgFile and topologyFile
# output two files
#   1. each seq with the selected TM frag (amino acid seq) in Fasta format
#   2. segment
#   Another one is a list 
#   ID  TMsegment(aa) DGscore

# dG file in the format
# ====
# Sequence             Predicted_DG 
#                                   
# FVSLGFLLIIIVPAFISCHAR       -0.530
# PYIYLGGAILAEVIGTTLMKF       1.464 
# VGTIICYCASFWLLAQTLAYI       0.739 
# AYAIWSGVGIVLISLLSWGFF       -1.146
# ====
#
import sys,re,os
import myfunc 

usage="""
Usage:  selectOneTMsegment.py [Options] [-i] topologyFile -dg dgFile
Options:
  -i               <file>    : input file
  -outpath         <dir>     : output path, default=./ ,two files will be output
                             : 1. $rootname.oneTMfrag.fa
                             : 2. $rootname.dgscorelist
  -method|-m 0|1|2|3|4       : select method
  -h|--help            : print this help message and exit
Created 2011-06-16, updated 2011-06-16,  Nanjiang
"""

def PrintHelp():
    print usage

def Read_TMSeg_DGScore(dgFile):#{{{
    fpin = open(dgFile,"r")
    lines=fpin.readlines()
    fpin.close()
    TMSegList=[]
    dgScoreList=[]
    for line in lines:
        strs=line.split()
        if len(strs) >= 2 and strs[0] != "Sequence":
            TMSegList.append(strs[0])
            dgScoreList.append(strs[1])
    return (TMSegList,dgScoreList)
#}}}
def SelectTMFrag(seqID, tmseqList, dgScoreList):#{{{
    dgScoreList=[float(x) for x in dgScoreList]
    numTM=len(tmseqList)
    threshold=0.0

    #debug
    print 
    print "##############################################"
    print seqID
    print tmseqList
    print dgScoreList
    print "numTM=", numTM
#

    selTMseq=""
    selTMdgscore=9999.0
    if numTM == 0: 
        sys.stdout.write("%s, numTM=%d, neglected\n"%(seqID, numTM))
    elif numTM == 1:
        if dgScoreList[0] <= threshold:
            selTMseq=tmseqList[0]
            selTMdgscore=dgScoreList[0]
            print selTMseq, selTMdgscore
        else:
            sys.stdout.write("%s, numTM=%d, but the dgscore(%f) > threshold(%f), neglected\n"%(seqID, numTM, dgScoreList[0],threshold))
    else:
        for i in range(numTM):
            if dgScoreList[i] < selTMdgscore:
                selTMseq = tmseqList[i]
                selTMdgscore=dgScoreList[i]
                print i, tmseqList[i],dgScoreList[i], "sel", selTMseq, selTMdgscore; 
                
        if selTMdgscore > threshold:
            sys.stdout.write("%s, numTM=%d, but the dgscore (%f)> threshold(%f), neglected\n"%(seqID, numTM, selTMdgscore, threshold))
            selTMseq=""
            selTMdgscore=9999.0
        else:
            print selTMseq, selTMdgscore
    print "##############################################"
    print 
    return (selTMseq, selTMdgscore)
#}}}
def SelectTMFrag1(seqID, tmseqList, dgScoreList):#{{{
    dgScoreList=[float(x) for x in dgScoreList]
    numTM=len(tmseqList)
    threshold=0.0

    selTMseq=""
    selTMdgscore=9999.0
    if numTM != 1: 
        sys.stdout.write("%s, numTM=%d, neglected\n"%(seqID, numTM))
    else:
        if dgScoreList[0] <= threshold:
            selTMseq=tmseqList[0]
            selTMdgscore=dgScoreList[0]
            print selTMseq, selTMdgscore
        else:
            sys.stdout.write("%s, numTM=%d, but the dgscore(%f) > threshold(%f), neglected\n"%(seqID, numTM, dgScoreList[0],threshold))
    return (selTMseq, selTMdgscore)
#}}}

def SelectOneTMsegment(idListTopo, annotationListTopo, topoList, tmseqList , dgScoreList, outFileSeloneTM, outFileDGList):#{{{

    fpSelTMFrag=open(outFileSeloneTM,"w")
    fpDGList=open(outFileDGList , "w")

#     print "lengthdgscore = ", len(dgScoreList)

    cntDGScore=0
    for iSeq in range (len(topoList)):
        seqID=idListTopo[iSeq]
        topo=topoList[iSeq]
        annoLine=annotationListTopo[iSeq]
        newstr=re.sub("M+","M", topo)
        numTM=newstr.count("M")
        #write topoDG
        selTMseq = ""
        if method==0:
            (selTMseq, selTMdgscore) = SelectTMFrag(seqID, tmseqList[cntDGScore:cntDGScore+numTM], dgScoreList[cntDGScore:cntDGScore+numTM] )
        elif method == 1:
            (selTMseq, selTMdgscore) = SelectTMFrag1(seqID, tmseqList[cntDGScore:cntDGScore+numTM], dgScoreList[cntDGScore:cntDGScore+numTM] )
        elif method == 2:
            (selTMseq, selTMdgscore) = SelectTMFrag2(seqID, tmseqList[cntDGScore:cntDGScore+numTM], dgScoreList[cntDGScore:cntDGScore+numTM] )
        elif method == 3:
            (selTMseq, selTMdgscore) = SelectTMFrag3(seqID, tmseqList[cntDGScore:cntDGScore+numTM], dgScoreList[cntDGScore:cntDGScore+numTM] )
        else:
            (selTMseq, selTMdgscore) = SelectTMFrag(seqID, tmseqList[cntDGScore:cntDGScore+numTM], dgScoreList[cntDGScore:cntDGScore+numTM] )
        if selTMseq != "":
            fpSelTMFrag.write(">%s dgscore= %f\n"%(annoLine, selTMdgscore))
            fpSelTMFrag.write("%s\n"%selTMseq)

        for i in range(numTM):
            fpDGList.write("%s %s %s\n"%(seqID, tmseqList[cntDGScore], dgScoreList[cntDGScore]))
            cntDGScore+=1
    if cntDGScore != len(dgScoreList):
        sys.stderr.write("Error! cntDGscore and dgScoreList does not match. cntDGscore= %d, len(dgScoreList)=%d\n" %(cntDGscore, len(dgScoreList)))
    fpSelTMFrag.close()
    fpDGList.close()
    return 0
#}}}

if __name__ == '__main__' :
    # Check argv
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        sys.exit(1)

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
                sys.exit(0)
            elif sys.argv[i] == "-i" or sys.argv[i] == "--infile":
                inFile=sys.argv[i+1]
                i = i + 2
            elif sys.argv[i] == "-dg" or sys.argv[i] == "--dg":
                dgFile=sys.argv[i+1]
                i = i + 2
            elif sys.argv[i] == "-m" or sys.argv[i] == "-method":
                method=int(sys.argv[i+1])
                i = i + 2
            elif sys.argv[i] == "-outpath" or sys.argv[i] == "--outpath":
                outpath=sys.argv[i+1]
                i = i + 2
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i])
                sys.exit(1)
        else:
            inFile=sys.argv[i]
            i+=1
           

    if inFile == "":
        print >> sys.stderr,"Error! Topology file not set."
        sys.exit(1)
    if dgFile == "":
        print >> sys.stderr,"Error!  dgFile not set."
        sys.exit(1)

    os.system("mkdir -p %s" % outpath)

    try :
        (idListTopo,annotationListTopo, topoList) = myfunc.ReadFasta(inFile)
        (TMSegList, dgScoreList)= Read_TMSeg_DGScore(dgFile)
        rootname=os.path.basename(os.path.splitext(inFile)[0])
        outFileDGList=outpath+"/"+rootname+".dgscorelist"
        outFileSeloneTM=outpath+"/"+rootname+".oneTMfrag.fa"

        SelectOneTMsegment(idListTopo, annotationListTopo, topoList, TMSegList
                , dgScoreList, outFileSeloneTM, outFileDGList)
        print  outFileDGList, "output"
        print  outFileSeloneTM, "output"

    except :
        print >> sys.stderr, "except for the input file: %s" % inFile
        raise 

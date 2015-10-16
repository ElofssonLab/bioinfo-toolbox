#!/usr/bin/env python
# add the gap penalty information to the fasta sequence according to the
# predicted dg values by myscanDG.pl
import sys,re,os,math
import myfunc
import subprocess

DEBUG=False

weight_M      = 0.9
p_shift       = 0.1

progname =  os.path.basename(sys.argv[0])

usage="""
Usage: %s [Options] fasta-seq-files ...

Note: the output file will be named as $outpath/$rootname.gp.fa

Options:
  -l         <file>: set the list file for fasta files
  -outpath   <dir> : set the output path for fasta files with gap penalties added, 
                   : default=the same as the fasta-seq-file
  -dgpath  <dir>   : path for the predicted topology files
                   : it must be set, and the dg files is stored as $dgpath/$rootname.dgscan
  -q               : quiet mode
  -h|--help        : print this help message and exit
  
  Other parameters related to the model
  -weight-m <float>
  -p-shift <float>

Created 2010-11-09, updated 2013-06-05, Nanjiang Shu
"""%(progname)

def PrintHelp():
    print usage

def AddFileList(fileList, inFile):#{{{
# return fileList
    fpin = open(inFile, "r")
    lines = fpin.readlines()
    fpin.close()
    for line in lines:
        line=line.strip()
        if line:
            fileList.append(line)
    return (fileList)
#}}}
def ReadDGScan(inFile):#{{{
# return dgscanList indexID
# indexID is a hash table to quickly locate the index of
# dgscanList given the id
    dgscanList=[]
    indexID={}
    fpin = open(inFile, "r")
    lines = fpin.readlines()
    fpin.close()
    i = 0
    cntRecord = 0
    while i < len(lines):
        line = lines[i]
        if line.find("#SeqID") == 0:
            dgscanList.append({})
            dgscanList[cntRecord]['id'] = line.split(':')[1].strip()
            i+=1

            while i < len(lines):
                line = lines[i]
                if line.find("#SeqLength") ==0:
                    dgscanList[cntRecord]['length']=int(line.split(':')[1].strip()); 
                elif line.find("#L") ==0:
                    dgscanList[cntRecord]['L']=int(line.split(':')[1].strip()); 
                elif line.find("#Number of sliding windows") ==0:
                    dgscanList[cntRecord]['numdata']=int(line.split(':')[1].strip()); 
                elif line[0] != '#':
                    break
                i+=1

            dgscanList[cntRecord]['seqindex'] = []
            dgscanList[cntRecord]['dgscores'] = []
            for j in range(dgscanList[cntRecord]['numdata']):
                line = lines[i+j]
                strs=line.split()
                dgscanList[cntRecord]['seqindex'].append(int(strs[0]))
                dgscanList[cntRecord]['dgscores'].append(float(strs[1]))
            i+= dgscanList[cntRecord]['numdata']
            
            indexID[dgscanList[cntRecord]['id']]=(cntRecord)
            cntRecord+=1
        i+=1
    return (dgscanList, indexID)
#}}}

def GetGapPenaltyFromDG(dgscan, seqLength):#{{{
# return (gpoArray, gpeArray, tgpeArray)
# increase gapopens for TM regions 
    gpoArray=[1]*seqLength
    gpeArray=[1]*seqLength
    tgpeArray=[1]*seqLength
    if DEBUG:
        print >> sys.stdout, "dgscan\n", dgscan

    seqIndex=dgscan['seqindex']
    dgscores=dgscan['dgscores']

    for j in range(len(seqIndex)):
        dg = dgscores[j]
        pinsertion = 1-1/(1+math.exp(-(dg-0.5)))-p_shift
        incPercent=weight_M*pinsertion
        idx=seqIndex[j]
        gpoArray[idx]  = gpoArray[idx]*(1.0+incPercent) 
        gpeArray[idx]  = gpeArray[idx]*(1.0+incPercent) 
        tgpeArray[idx] = tgpeArray[idx]*(1.0+incPercent)

    return (gpoArray, gpeArray, tgpeArray)
#}}}

def AddGapPenaltyByPredDG(fastaFile, ssPath, outPath):#{{{
# add gap penalties
    rootname = os.path.basename(os.path.splitext(fastaFile)[0])
    inFilePath=os.path.dirname(fastaFile)
    if inFilePath == "":
        inFilePath='./'
    (annotationList, seqList) = myfunc.ReadFasta_without_id(fastaFile)

    if outPath == "":
        localOutPath=inFilePath
    else:
        localOutPath=outPath

    outFile="%s/%s.gp.fa"%(localOutPath,rootname)
    fpout = open(outFile, "w")

    dgscanFile="%s/%s.dgscan"%(dgPath, rootname)
    (dgscanList, indexID)= ReadDGScan(dgscanFile)
    if len(seqList) != len(dgscanList):
        print >> sys.stderr, "Error! numseq not match, numSeq=%d, numDGList=%d, for file %s"(len(seqList), len(dgscanList), fastaFile)
    for i in range(len(seqList)) :
        seqLength = len(seqList[i])
        dgscan=dgscanList[i]

        if DEBUG: 
            print >> sys.stdout, "%d:%s" %(i, annotationList[i])
        (gpoArray, gpeArray, tgpeArray) = GetGapPenaltyFromDG(dgscan, seqLength)
# write out the result
        fpout.write(">%s\n"% annotationList[i])
        fpout.write("%s\n"% seqList[i])
        fpout.write("{gpo: ")
        for j in range(seqLength):
            fpout.write("%f "% gpoArray[j])
        fpout.write("}\n")
        fpout.write("{gpe: ")
        for j in range(seqLength):
            fpout.write("%f "% gpeArray[j])
        fpout.write("}\n")
        fpout.write("{tgpe: ")
        for j in range(seqLength):
            fpout.write("%f "% tgpeArray[j])
        fpout.write("}\n")
            
    fpout.close()
    return len(seqList)
#}}}

if __name__ == '__main__' :
    # Check argv
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        sys.exit(1)

    outPath=""
    listFile="";  
    dgPath=""
    fastaFileList=[]
    isQuiet=False

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
            elif sys.argv[i] ==  "-outpath" or  sys.argv[i] == "--outpath":
                outPath=sys.argv[i+1]
                i+=2
            elif sys.argv[i] ==  "-dgpath" or  sys.argv[i] == "--dgpath":
                dgPath=sys.argv[i+1]
                i+=2
            elif sys.argv[i] ==  "-l" or  sys.argv[i] == "--list":
                listFile=sys.argv[i+1]
                i+=2
            elif sys.argv[i] ==  "-q" or  sys.argv[i] == "--q" or sys.argv[i] == "--quiet":
                isQuiet=True
                i+=1
            elif sys.argv[i] ==  "-weight-m" or  sys.argv[i] == "--weight-m":
                weight_M=float(sys.argv[i+1])
                i+=2
            elif sys.argv[i] ==  "-p-shift" or  sys.argv[i] == "--p-shift":
                p_shift=float(sys.argv[i+1])
                i+=2
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i])
                sys.exit(1)
        else:
            fastaFileList.append(sys.argv[i])
            i+=1

    if dgPath == "":
        print >> sys.stderr,"Error! dgPath not set."
        sys.exit(1)
    if len(fastaFileList) == 0 and listFile == "":
        print >> sys.stderr,"Error! no input is set"
        sys.exit(1)

    try :
        if outPath != "" :
            try:
                subprocess.check_call(["mkdir", "-p", outPath])
            except subprocess.CalledProcessError, e:
                sys.exit(1)
        if listFile != "":
            fastaFileList = AddFileList(fastaFileList, listFile)

        cntFile= 0
        for fastaFile in fastaFileList:
            numseq = AddGapPenaltyByPredDG(fastaFile, dgPath, outPath)
            if not isQuiet:
                print >> sys.stdout, "%d\t%d sequences in the file %s have been added gap penalties" %(cntFile, numseq, fastaFile)
            cntFile+=1

    except :
        print >>sys.stderr, "except for the function:%s"%sys._getframe().f_code.co_name 
        raise 

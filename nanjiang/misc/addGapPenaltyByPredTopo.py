#!/usr/bin/env python
# add the gap penalty information to the fasta sequence according to the
# predicted topology by scampi
import sys,re,os,math;
import myfunc;

DEBUG=False;

weight_M      = 0.9;

usage="""
Usage: addGapPenaltyByPredTopo.py [Options] fasta-seq-files ...

Note: the output file will be named as $outpath/$rootname.gp.fa

Options:
  -l         <file>: set the list file for fasta files
  -outpath   <dir> : set the output path for fasta files with gap penalties added, 
                   : default=the same as the fasta-seq-file
  -topopath  <dir> : path for the predicted topology files
                   : it must be set, and the secondary structure files is stored as $topopath/$rootname.topo
  -q               : quiet mode
  -h|--help        : print this help message and exit
  
  Other parameters related to the model
  -weight-m <float>

Created 2010-11-09, updated 2010-11-09, Nanjiang
"""

def PrintHelp():
    print usage;

def AddFileList(fileList, inFile):#{{{
# return fileList
    fpin = open(inFile, "r");
    lines = fpin.readlines();
    fpin.close();
    for line in lines:
        line=line.strip();
        if line:
            fileList.append(line);
    return (fileList);
#}}}


def GetGapPenaltyFromTopo(topo, seqLength):#{{{
# return (gpoArray, gpeArray, tgpeArray)
# increase gapopens for TM regions 
    gpoArray=[1]*seqLength;
    gpeArray=[1]*seqLength;
    tgpeArray=[1]*seqLength;
    if DEBUG:
        print >> sys.stdout, "%s" % topo;
    i = 0;
    while i < seqLength:
        if topo[i] == 'M':
            state = topo[i];
            j = 0;
            while  i+j < seqLength and topo[i+j] == state:
                j += 1;
            if j>=1:
                mid=(j+1)/2;
                for jj in range (1,j-1):
                    posShift=1.0-abs(jj+1-mid)/float(mid); # posShift (0.0-1.0)
                    if jj == 0 or jj == j-1:
                        weightPos = 0.05;
                    else:
                        weightPos=posShift**1;
                    incPercent=weight_M*weightPos;
                    gpoArray[i+jj] = gpoArray[i+jj]*(1.0+incPercent) ; 
                    gpeArray[i+jj] = gpeArray[i+jj]*(1.0+incPercent) ; 
                    tgpeArray[i+jj] = tgpeArray[i+jj]*(1.0+incPercent) ; 
            i+= j;
        else:
            i+= 1;
    return (gpoArray, gpeArray, tgpeArray);
#}}}

def AddGapPenaltyByPredTopo(fastaFile, ssPath, outPath):#{{{
# add gap penalties
    rootname = os.path.basename(os.path.splitext(fastaFile)[0]);
    inFilePath=os.path.dirname(fastaFile);
    if inFilePath == "":
        inFilePath='./';
    (annotationList, seqList) = myfunc.ReadFasta_without_id(fastaFile);

    if outPath == "":
        localOutPath=inFilePath;
    else:
        localOutPath=outPath;

    outFile="%s/%s.gp.fa"%(localOutPath,rootname)
    fpout = open(outFile, "w");

    topoFile="%s/%s.topo"%(topoPath, rootname);
    (tmp, topoList)=ReadFasta(topoFile);
    if len(seqList) != len(topoList):
        print >> sys.stderr, "Error! numseq not match, numSeq=%d, numTopo=%d, for file %s"(len(seqList), len(topoList), fastaFile);
    for i in range(len(seqList)) :
        seqLength = len(seqList[i]);
        topo=topoList[i];
        if seqLength != len(topo):
            print >> sys.stderr, "Error! aaSeqLength = %d, but the topoSeqLength = %d, for file %s" %(seqLength, len(topo), topoFile);
            sys.exit(1);

        if DEBUG: 
            print >> sys.stdout, "%d:%s" %(i, annotationList[i])
        (gpoArray, gpeArray, tgpeArray) = GetGapPenaltyFromTopo(topo, seqLength);
# write out the result
        fpout.write(">%s\n"% annotationList[i]);
        fpout.write("%s\n"% seqList[i]);
        fpout.write("{gpo: ");
        for j in range(seqLength):
            fpout.write("%f "% gpoArray[j]);
        fpout.write("}\n");
        fpout.write("{gpe: ");
        for j in range(seqLength):
            fpout.write("%f "% gpeArray[j]);
        fpout.write("}\n");
        fpout.write("{tgpe: ");
        for j in range(seqLength):
            fpout.write("%f "% tgpeArray[j]);
        fpout.write("}\n");
            
    fpout.close();
    return len(seqList);
#}}}

if __name__ == '__main__' :
    # Check argv
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp();
        sys.exit(1);

    outPath="";
    listFile="";  
    topoPath="";
    fastaFileList=[];
    isQuiet=False;

    i = 1;
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            isNonOptionArg=False;
            i = i + 1;
        elif sys.argv[i] == "--":
            isNonOptionArg=True;
            i = i + 1;
        elif sys.argv[i][0] == "-":
            if sys.argv[i] ==  "-h" or  sys.argv[i] == "--help":
                PrintHelp();
                sys.exit(0);
            elif sys.argv[i] ==  "-outpath" or  sys.argv[i] == "--outpath":
                outPath=sys.argv[i+1];
                i+=2;
            elif sys.argv[i] ==  "-topopath" or  sys.argv[i] == "--topopath":
                topoPath=sys.argv[i+1];
                i+=2;
            elif sys.argv[i] ==  "-l" or  sys.argv[i] == "--list":
                listFile=sys.argv[i+1];
                i+=2;
            elif sys.argv[i] ==  "-q" or  sys.argv[i] == "--q" or sys.argv[i] == "--quiet":
                isQuiet=True;
                i+=1;
            elif sys.argv[i] ==  "-weight-m" or  sys.argv[i] == "--weight-m":
                weight_M=float(sys.argv[i+1]);
                i+=2;
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i]);
                sys.exit(1);
        else:
            fastaFileList.append(sys.argv[i]);
            i+=1;
           

    if topoPath == "":
        print >> sys.stderr,"Error! topoPath not set.";
        sys.exit(1);
    if len(fastaFileList) == 0 and listFile == "":
        print >> sys.stderr,"Error! no input is set";
        sys.exit(1);

        

    try :
        if outPath != "" :
            os.system("mkdir -p %s" % outPath);
        if listFile != "":
            fastaFileList = AddFileList(fastaFileList, listFile);

        cntFile= 0;
        for fastaFile in fastaFileList:
            numseq = AddGapPenaltyByPredTopo(fastaFile, topoPath, outPath);
            if not isQuiet:
                print >> sys.stdout, "%d\t%d sequences in the file %s have been added gap penalties" %(cntFile, numseq, fastaFile);
            cntFile+=1;

    except :
        print >>sys.stderr, "except for the function:%s"%sys._getframe().f_code.co_name ;
        raise ;

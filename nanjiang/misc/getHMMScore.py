#!/usr/bin/env python
# get HMM scores from the mySCAMPI output file 
import sys,re,os;
import myfunc ;


usage="""
Usage:   getHMMScore.py [Options] [-i] scampioutputfile
Note: scampi output file has the extension xml.res
Options:
  -o         <file>    : outputfile
  -h|--help            : print this help message and exit
Created 2010-08-29, updated 2010-08-29, Nanjiang
"""

def PrintHelp():
    print usage;

def GetHMMScore(inFile):#{{{
# return (idList, lengthList, normLogList, logoddsList, reversiList, isTMProList);
    idList=[];
    lengthList=[];
    normLogList=[];
    logoddsList=[];
    reversiList=[];
    isTMProList=[];
    fpin = open(inFile, "r");
    line = fpin.readline();
    while line:
        if line.find("SeqID:")>=0:
            seqID=".";
            normLog=".";
            logodds=".";
            reversi=".";
            isTMPro=".";
            length=0;
            strs=line.split(":");
            if strs[1].strip():
                seqID=strs[1].strip();
            while 1: 
                line = fpin.readline();
                if not line.strip():
                    break;
                strs=line.split(":");
                if strs[0] == "SeqLength":
                    length=int(strs[1]);
                elif strs[0] == "NormalizedLogLikelihood":
                    normLog=strs[1].strip();
                elif strs[0] == "Logodds":
                    logodds=strs[1].strip();
                elif strs[0] == "Reversi":
                    reversi=strs[1].strip();
                elif strs[0] == "IsTMProtein":
                    isTMPro=strs[1].strip();
                elif strs[0][0] == ">":
                    if seqID=='.':
                        seqID=myfunc.GetSeqIDFromAnnotation(strs[0]);
            idList.append(seqID);
            lengthList.append(length);
            normLogList.append(normLog);
            logoddsList.append(logodds);
            reversiList.append(reversi);
            isTMProList.append(isTMPro);
        line = fpin.readline();

    fpin.close();
    return (idList, lengthList, normLogList, logoddsList, reversiList, isTMProList);
#}}}

if __name__ == '__main__' :
    # Check argv
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp();
        sys.exit(1);

    outFile="";
    inFile="";
    begin=0;
    end=999999999;

    method=2;

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
            elif sys.argv[i] == "-i" or sys.argv[i] == "--infile":
                inFile=sys.argv[i+1];
                i = i + 2;
            elif sys.argv[i] == "-o" or sys.argv[i] == "--outfile":
                outFile=sys.argv[i+1];
                i = i + 2;
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i]);
                sys.exit(1);
        else:
            inFile=sys.argv[i];
            i+=1;
           

    if inFile == "":
        print >> sys.stderr,"Error! Input file not set.";

    fpout = sys.stdout;
    if outFile != "":
        fpout = open(outFile,"w");

    try :
        (idList, lengthList, normLogList, logoddsList, reversiList, isTMProList) = GetHMMScore(inFile);
        # print out the result
        print >> fpout, "%-16s %8s %8s %8s %8s %6s" % ("#seqID", "NormLog", "Logodds", "reversi","isTMPro","length");
        for i in range (len(idList)):
            print >> fpout, "%-16s %8s %8s %8s %8s %6d" % (idList[i], normLogList[i], logoddsList[i], reversiList[i], isTMProList[i], lengthList[i]);

        if fpout != sys.stdout:
            fpout.close();
    except :
        print >>sys.stderr, "except for the input file: %s" % inFile;
        raise ;

#!/usr/bin/env python

# Modified by Nanjiang 2010-08-11
# Three methods for pairwise topology comparison are used.
# 1. Gapless
# 2. Local 
# 3. Global


import sys,re,os;
# sys.path.append('../../topology_comparison/bin');
import comptopo as ct;

#class_gapless: one of the four (OK, SHIFT, INV, DIFF) relationships between two topologies
#num1_gapless : number of trans membrane segments of protein1, filtered by the gapless algorithm
#num2_gapless : number of trans membrane segments of protein2, filtered by the gapless algorithm 

usage="""
Usage:  compare_topos.py argument-list
Options:
  -i         <file> : input file
  -f|--format <int> : set the file format, default, auto detect
                    : 0. multiple topology alignment in Fasta format, the first one is the query sequence
                    : 1. multiple topology alignment in alignment format, one record per line
                    :    in the format
                    :    text:aligned-seq
                    : 2  pairwise alignment in the format
                    :    #Number of alignments: 10
                    :    #Topology alignment 1
                    :    >seq1
                    :    topology in one line
                    :    >seq2
                    :    topology in one line
                    :
  --comp      <int> : comparison method, default = 0 
                    : 0. query (first sequence) to others
                    : 1. all-to-all
  --log      <file> : output the logfile 
  -o         <file> : output the result to file
  -h|--help         : print this help message and exit
Created 2010-08-11, updated 2010-09-07, Nanjiang
"""
def PrintHelp():
    print usage

def ReadPairwiseTopoAlignment(inFile):#{{{
# return (idList1, idList2, topoSeqList1,topoSeqList2, pidList, evalueList)
    idList1=[];
    idList2=[];
    topoSeqList1=[];
    topoSeqList2=[];
    pidList=[];
    evalueList=[];
    fpin = open(inFile, "r");
    lines = fpin.readlines();
    fpin.close();
    i = 0;
    numAlign=0;
    while i < len(lines):
        line = lines[i];
        if line.find("#Number of alignments") >=0 :
            numAlign=int(line.split(":")[1]);
        elif line.find("#Topology alignment") >=0:
            id1=ct.GetSeqIDFromAnnotation(lines[i+1]);
            topo1=lines[i+2].strip();
            id2=ct.GetSeqIDFromAnnotation(lines[i+3]);
            pid=ct.GetPIDFromAnnotation(lines[i+3]);
            evalue=ct.GetEvalueFromAnnotation(lines[i+3]);
            topo2=lines[i+4].strip();
            idList1.append(id1);
            idList2.append(id2);
            topoSeqList1.append(topo1);
            topoSeqList2.append(topo2);
            pidList.append(pid);
            evalueList.append(evalue);
            i=i+5
        i+=1;
    #check
    if numAlign != len(idList1):
        print >> sys.stderr, "The number of alignment read in (%d) does not match the annotation (%d)" % (len(idList1), numAlign);
    return (idList1, idList2, topoSeqList1,topoSeqList2, pidList, evalueList);
#}}}

if __name__ == '__main__' :
    # Check argv
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp();
        sys.exit(1);

    logFile = "";
    isOutPutLogFile = False;
    isFileFormatForcedSet = False;
    outFile="";
    isQuiet=False;
    fileFormat=0;
    compareMethod=0;
    inFile="";

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
            elif sys.argv[i] == "-f" or sys.argv[i] == "--format" :
                fileFormat=int(sys.argv[i+1]);
                isFileFormatForcedSet = True;
                i = i + 2;
            elif sys.argv[i] == "-comp" or sys.argv[i] == "--comp" :
                compareMethod=int(sys.argv[i+1]);
                i = i + 2;
            elif sys.argv[i] == "-log" or sys.argv[i] == "--log" :
                logFile=sys.argv[i+1];
                i = i + 2;
            elif sys.argv[i] == "-q":
                isQuiet=True;
                i = i + 1;
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i]);
                sys.exit(1);
        else:
            print >> sys.stderr,("Error! Wrong argument:%s"% sys.argv[i]);
            sys.exit(1);

    if inFile == "":
        print >> sys.stderr,"Error! Input file not set.";

    fpout = sys.stdout;
    if outFile != "":
        fpout = open(outFile,"w");

    fpLog = 0;
    if logFile != "":
        fpLog = open(logFile,"w");
        isOutPutLogFile = True;

    try :
        if not isFileFormatForcedSet:
            fileFormat = ct.DetectFileFormat(inFile);

        if fileFormat != 2:
            (idList, topoSeqList, pidList, evalueList) = ct.ReadAliTopoNew(inFile, fileFormat);
            numSeq = len(idList);
            # print header line
            print >> fpout, "#%-9s %-15s %8s %5s %5s %8s %5s %5s %8s %5s %5s %5s %6s" % ("protein1", "protein2", "gapless", "NrTM1", "NrTM2", "local", "NrTM1", "NrTM2", "global", "NrTM1","NrTM2", "PID", "Evalue");
            for i in range (1, numSeq):
                strTop1 = topoSeqList[0];
                strTop2 = topoSeqList[i];
                strProtein1 = idList[0];
                strProtein2 = idList[i];

                if fpLog != 0:
                    print >> fpLog, "========   original, removeUngaps, gapless, local, global  ==================";
                    print >> fpLog, "%-20s:%s"%(strProtein1, strTop1);
                    print >> fpLog, "%20s:%s"%(strProtein2, strTop2);
                    print >> fpLog;

                class_gapless, num1_gapless, num2_gapless = ct.CompareToposGaplesslyNew(strTop1, strTop2, strProtein1, strProtein2, fpLog);
                class_local, num1_local, num2_local = ct.CompareToposLocallyNew(strTop1, strTop2, strProtein1, strProtein2, fpLog);
                class_global, num1_global, num2_global = ct.CompareToposGloballyNew(strTop1, strTop2, strProtein1, strProtein2, fpLog);

                print >> fpout, "%-9s %-16s %8s %5s %5s %8s %5s %5s %8s %5s %5s %5s %6s"  % (strProtein1, strProtein2, class_gapless, num1_gapless, num2_gapless, class_local, num1_local, num2_local, class_global, num1_global, num2_global, pidList[i], evalueList[i]);
        else:
            (idList1, idList2, topoSeqList1,topoSeqList2, pidList, evalueList) = ReadPairwiseTopoAlignment(inFile);
            numSeq = len(idList1);
            # print header line
            print >> fpout, "#%-9s %-15s %8s %5s %5s %8s %5s %5s %8s %5s %5s %5s %6s" % ("protein1", "protein2", "gapless", "NrTM1", "NrTM2", "local", "NrTM1", "NrTM2", "global", "NrTM1","NrTM2", "PID", "Evalue");
            for i in range (numSeq):
                strTop1 = topoSeqList1[i];
                strTop2 = topoSeqList2[i];
                strProtein1 = idList1[i];
                strProtein2 = idList2[i];

                if fpLog != 0:
                    print >> fpLog, "========   original, removeUngaps, gapless, local, global  ==================";
                    print >> fpLog, "%-20s:%s"%(strProtein1, strTop1);
                    print >> fpLog, "%-20s:%s"%(strProtein2, strTop2);
                    print >> fpLog;

                class_gapless, num1_gapless, num2_gapless = ct.CompareToposGaplesslyNew(strTop1, strTop2, strProtein1, strProtein2, fpLog);
                class_local, num1_local, num2_local = ct.CompareToposLocallyNew(strTop1, strTop2, strProtein1, strProtein2, fpLog);
                class_global, num1_global, num2_global = ct.CompareToposGloballyNew(strTop1, strTop2, strProtein1, strProtein2, fpLog);

                print >> fpout, "%-9s %-16s %8s %5s %5s %8s %5s %5s %8s %5s %5s %5s %6s"  % (strProtein1, strProtein2, class_gapless, num1_gapless, num2_gapless, class_local, num1_local, num2_local, class_global, num1_global, num2_global, pidList[i], evalueList[i]);

        
        if fpout != sys.stdout:
            fpout.close();
        if fpLog != 0:
            fpLog.close();

    except :
        print >>sys.stderr, "except for the input file: %s" % inFile;
        raise ;
    

#!/usr/bin/env python
# get the histogram of topology comparison vs pid or evalue

import sys,re,os;

usage="""
Usage:  histoTopocmp.py [Options] -i topoCmpFile
Options:
  -i         <file> : input file
  -bins      <file> : set the bins file
  -pc 0|1           : whether print the count, defualt = 0
  -pf 0|1           : whether print normalized frequencies, default =1
  -o         <file> : output the result to file
  -not-ignore-100p  : do not ignore the one with 100% sequence identity
  -h|--help         : print this help message and exit
Created 2010-08-16, updated 2010-08-18, Nanjiang
"""

def PrintHelp():
    print usage;

pidBins =[#{{{
        5,
        10,
        15,
        20,
        25,
        30,
        35,
        40,
        45,
        50,
        55,
        60,
        65,
        70,
        75,
        80,
        85,
        90,
        95
        ];#}}}

def GetBinIndex(pid, pidBins):
    numBins = len(pidBins);
    index = numBins;
    for i in range (0, numBins):
        if pid < pidBins[i]:
            index = i;
            break;
    return index;

def ReadPIDBins(inFile):
    pidBins=[];
    fpin = open(inFile, "r");
    line = fpin.readline();
    while line:
        line = line.rstrip('\n');
        if line and line[0] != "#":
            pidBins.append(float(line));
        line = fpin.readline();
    fpin.close();
    return pidBins;

def ReadTopoCmpFile(inFile):#{{{
#return (gaplessCmpList, localCmpList, globalCmpList, pidList, evalueList) 
    gaplessCmpList=[];
    localCmpList=[];
    globalCmpList=[];
    pidList=[];
    evalueList=[];

    fpin = open(inFile, "r");
    line = fpin.readline();
    while line:
        line = line.rstrip('\n');
        if line and line[0] != "#":
            strs=line.split();
            gaplessCmpList.append(strs[2]);
            localCmpList.append(strs[5]);
            globalCmpList.append(strs[8]);
            pidList.append(float(strs[11]));
            evalueList.append(float(strs[12]));
        line = fpin.readline();
    fpin.close();
    return (gaplessCmpList, localCmpList, globalCmpList, pidList, evalueList);
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
    binFile="";
    isPrintFreq = True;
    isPrintCount = False;
    isIgnore100p = True; 

    states=['DIFF','SHIFT','INV','OK'];

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
            elif sys.argv[i] == "-bins" or sys.argv[i] == "--bins":
                binFile=sys.argv[i+1];
                i = i + 2;
            elif sys.argv[i] == "-pf" or sys.argv[i] == "--pf":
                isPrintFreq=bool(int(sys.argv[i+1]));
                i = i + 2;
            elif sys.argv[i] == "-pc" or sys.argv[i] == "--pc":
                isPrintCount=bool(int(sys.argv[i+1]));
                i = i + 2;
            elif sys.argv[i] == "-not-ignore-100p":
                isIgnore100p=False;
                i = i + 1;
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
        gaplessHist = [{}];
        localHist =[{}];
        globalHist =[{}];
        gaplessHistNorm = [{}];
        localHistNorm =[{}] ;
        globalHistNorm =[{}];

        if binFile != "":
            pidBins = ReadPIDBins(binFile);

        # initialization
        numBins = len(pidBins);
        for i in range(0,numBins+1):
            for s in states:
                gaplessHist[i][s]=0;
                localHist[i][s]=0;
                globalHist[i][s]=0;
                gaplessHistNorm[i][s]=0.0;
                localHistNorm[i][s]=0.0;
                globalHistNorm[i][s]=0.0;
            gaplessHist.append({});
            localHist.append({});
            globalHist.append({});
            gaplessHistNorm.append({});
            localHistNorm.append({});
            globalHistNorm.append({}) ;

#         print gaplessHist ;
#         print localHist ;
#         print globalHist;
#         print gaplessHistNorm ;
#         print localHistNorm ;
#         print globalHistNorm ;

        (gaplessCmpList, localCmpList, globalCmpList, pidList, evalueList) = ReadTopoCmpFile(inFile);
        #print pidList;
        #print evalueList;
#         print gaplessCmpList;
#         print len(gaplessCmpList);
        
        numRecord = len(gaplessCmpList);

        for i in range (numRecord):
            if isIgnore100p and pidList[i] == 100.0:
                continue;
            binIndex = GetBinIndex(pidList[i],pidBins);
            gaplessHist[binIndex][gaplessCmpList[i]] += 1;
            localHist[binIndex][localCmpList[i]] += 1;
            globalHist[binIndex][globalCmpList[i]] += 1;
        
#        print gaplessHist;
        # normalized histogram
        for i in range(0,numBins+1):
            total = float(sum(gaplessHist[i][key] for key in gaplessHist[i].keys()));
            if total == 0: total = 1e-6;
            for key in gaplessHist[i].keys():
                gaplessHistNorm[i][key] = gaplessHist[i][key]/float(total);

            total = float(sum(localHist[i][key] for key in localHist[i].keys()));
            if total == 0: total = 1e-6;
            for key in localHist[i].keys():
                localHistNorm[i][key] = localHist[i][key]/float(total);

            total = float(sum(globalHist[i][key] for key in globalHist[i].keys()));
            if total == 0: total = 1e-6;
            for key in globalHist[i].keys():
                globalHistNorm[i][key] = globalHist[i][key]/float(total);

        # print out the result
        # print header line
        print >> fpout, "#topology comparison vs pid (sequence identity)"
        print >> fpout, "#";

        if isPrintFreq:
            print >> fpout, "#normalized histogram"
            print >> fpout, "%-10s %20s %8s %20s %8s %20s" % ("#", "gapless", "", "local", "", "global");
            print >> fpout, "%-10s %5s %5s %5s %5s %4s %8s %5s %5s %5s %5s %4s %8s %5s %5s %5s %5s %4s" %("pidBins","DIFF","SHIFT","INV","OK","sum","","DIFF","SHIFT","INV","OK","sum","","DIFF","SHIFT","INV","OK", "sum");

            for i in range(0,numBins+1):
                if i == 0:
                    fpout.write("%-10s"%("<%d" % pidBins[i]));
                elif i == numBins:
                    fpout.write("%-10s"%(">%d" %pidBins[i-1]));
                else:
                    fpout.write("%-10s"%("%d-%d" %(pidBins[i-1], pidBins[i])));

                for key in states:
                    fpout.write(" %5.3f"% (gaplessHistNorm[i][key]));
                fpout.write(" %4s"% ("%d"% (sum(gaplessHist[i][key] for key in gaplessHist[i].keys()))));
                fpout.write(" %8s"% "");

                for key in states:
                    fpout.write(" %5.3f"%localHistNorm[i][key]);
                fpout.write(" %4s"% ("%d"% (sum(localHist[i][key] for key in localHist[i].keys()))));
                fpout.write(" %8s"% "");

                for key in states:
                    fpout.write(" %5.3f"%globalHistNorm[i][key]);
                fpout.write(" %4s"% ("%d"% (sum(globalHist[i][key] for key in globalHist[i].keys()))));
                fpout.write("\n");

        if isPrintCount:
            print >> fpout, "#";
            print >> fpout, "#frequency histogram"
            print >> fpout, "%-10s %20s %8s %20s %8s %20s" % ("#", "gapless", "", "local", "", "global");
            print >> fpout, "%-10s %5s %5s %5s %5s %4s %8s %5s %5s %5s %5s %4s %8s %5s %5s %5s %5s %4s" %("pidBins","DIFF","SHIFT","INV","OK","sum","","DIFF","SHIFT","INV","OK","sum","","DIFF","SHIFT","INV","OK", "sum");
            for i in range(0,numBins+1):
                if i == 0:
                    fpout.write("%-10s"%("<%d" % pidBins[i]));
                elif i == numBins:
                    fpout.write("%-10s"%(">%d" %pidBins[i-1]));
                else:
                    fpout.write("%-10s"%("%d-%d" %(pidBins[i-1], pidBins[i])));

                for key in states:
                    fpout.write(" %5d"%gaplessHist[i][key]);
                fpout.write(" %4s"% ("%d"% (sum(gaplessHist[i][key] for key in gaplessHist[i].keys()))));
                fpout.write(" %8s"% "");

                for key in states:
                    fpout.write(" %5d"%localHist[i][key]);
                fpout.write(" %4s"% ("%d"% (sum(localHist[i][key] for key in localHist[i].keys()))));
                fpout.write(" %8s"% "");

                for key in states:
                    fpout.write(" %5d"%globalHist[i][key]);
                fpout.write(" %4s"% ("%d"% (sum(globalHist[i][key] for key in globalHist[i].keys()))));
                fpout.write("\n");



        if fpout != sys.stdout:
            fpout.close();
        if fpLog != 0:
            fpLog.close();

    except :
        print >>sys.stderr, "except for the input file: %s" % inFile;
        raise ;


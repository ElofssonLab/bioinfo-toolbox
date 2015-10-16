#!/usr/bin/env python
# clean single span and non TM proteins
# Input: topoFile
import sys,re,os;
import myfunc ;

usage="""
Usage:  topoAddDGscore.py [Options] [-i] topologyFile 
Options:
  -i   <file>    : input file
  -o   <file>    : output file
  -mintm <int>   : set the minimal number of TM helices, default = 2
  -h|--help      : print this help message and exit
Created 2010-09-03, updated 2011-09-02,  Nanjiang
"""

def PrintHelp():
    print usage;

def CleanSingleSpanTMPro(idListTopo, annotationListTopo, topoList):#{{{
    fpout=sys.stdout;
    if outFile != "":
        fpout = open(outFile,"w");
    
    for i in range(len(idListTopo)):
        seqID=idListTopo[i];
        topo=topoList[i];
        annoLine=annotationListTopo[i];
        newstr=re.sub("-+","", topo);
        newstr=re.sub("M+","M", newstr);
        numTM=newstr.count("M");
        if numTM >= MIN_NUMTM :
            fpout.write(">%s\n"%annoLine);
            fpout.write("%s\n"%topo);
        else:
            sys.stderr.write("Seq %s, numTM (%d) <%d, removed. \n" %(seqID, numTM, MIN_NUMTM));
        
    if fpout != sys.stdout:
        fpout.close();

#}}}

if __name__ == '__main__' :
    # Check argv
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp();
        sys.exit(1);

    outFile="";
    inFile="";
    MIN_NUMTM=2; #minimal allowed number of TM regions

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
            elif sys.argv[i] == "-mintm" or sys.argv[i] == "--mintm":
                MIN_NUMTM=int(sys.argv[i+1]);
                i = i + 2;
            elif sys.argv[i] == "-o" or sys.argv[i] == "--out":
                outFile=sys.argv[i+1];
                i = i + 2;
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i]);
                sys.exit(1);
        else:
            inFile=sys.argv[i];
            i+=1;
           

    if inFile == "":
        print >> sys.stderr,"Error! Topology file not set.";
        sys.exit(1);


    try :
        (idListTopo,annotationListTopo, topoList) = myfunc.ReadFasta(inFile);
        CleanSingleSpanTMPro(idListTopo, annotationListTopo, topoList);
    except :
        print >>sys.stderr, "except for the input file: %s" % inFile;
        raise ;

#!/usr/bin/env python
# Given topology multiple alignment, fill in gaps and so on
# Input: topoFile
import sys,re,os;
import myfunc, comptopo ;

usage="""
Usage:  trimTopoMSA.py [Options] [-i] topologyFile 
Options:
  -i   <file>    : input file
  -o   <file>    : output file
  -h|--help      : print this help message and exit
Created 2011-08-15, updated 2011-08-15,  Nanjiang
"""

def PrintHelp():
    print usage;

def TrimTopoMSA(idListTopo, annotationListTopo, topoList):#{{{
    fpout=sys.stdout;
    if outFile != "":
        fpout = open(outFile,"w");
    
    for i in range(len(idListTopo)):
        seqID=idListTopo[i];
        topo=topoList[i];
        annoLine=annotationListTopo[i];
        topo=comptopo.trimTopo(topo);
        #topo=comptopo.trimTopo(topo);
        topo=comptopo.filterTopo(topo);
        fpout.write(">%s\n" %(annoLine));
        fpout.write("%s\n"%(topo));
        
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
        TrimTopoMSA(idListTopo, annotationListTopo, topoList);
    except :
        print >>sys.stderr, "except for the input file: %s" % inFile;
        raise ;

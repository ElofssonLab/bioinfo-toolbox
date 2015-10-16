#!/usr/bin/env python
# match the pairwise topology alignment by the pairwise sequence alignment
import sys,re,os;
import myfunc ;

usage="""
Usage:   matchTopoPairAln.py [Options] -aln alignFile -qtopo queryTopoFile -ttopo targetsTopologyFile
Options:
  -o         <file>    : outputfile
  -h|--help            : print this help message and exit
Created 2010-08-29, updated 2010-08-29, Nanjiang

output format
#Number of alignments: 2
#Topology alignment 1
>seq1
topology in one line
>seq2
topology in one line

#Topology alignment 2

"""

def PrintHelp():
    print usage;

def ReadNeedleAlignment(inFile):#{{{
    alns=[];
    try:
        cntAln=0;
        fpin = open(inFile, "r");
        line = fpin.readline();
        while line:
            if line.find("Aligned_sequences") >=0:
                alns.append({});
                line = fpin.readline();
                id1=line.split(":")[1].strip();
                line = fpin.readline();
                id2=line.split(":")[1].strip();
                while 1:
                    line = fpin.readline();
                    if line.find("#=====")>=0:
                        break;
                line = fpin.readline(); # neglect one blank line
                # now start to read the alignment 
                alnseq1="";
                alnseq2="";
                alnrel="";
                while 1:
                    line = fpin.readline();
                    if not line.strip():
                        break;
                    alnseqend=71;
                    indexws = line.find(" ", 21); # find the index in case the alignseq does not reach 50
                    if (indexws < 71):
                        alnseqend=indexws;
                    alnseq1+=line[21:alnseqend];
                    line = fpin.readline();
                    alnrel+=line[21:alnseqend];
                    line = fpin.readline();
                    alnseq2+=line[21:alnseqend];
                    line = fpin.readline();
                alns[cntAln]['alnseq1']=alnseq1;
                alns[cntAln]['seqid1']=id1;
                alns[cntAln]['alnrel']=alnrel;
                alns[cntAln]['alnseq2']=alnseq2;
                alns[cntAln]['seqid2']=id2;
                cntAln+=1;
            line = fpin.readline();
        fpin.close();
    except IOError:
        print >>sys.stderr, "except for the function:%s"%sys._getframe().f_code.co_name ;
    return alns;
#}}}
def MatchTopoPairAln(queryTopoFile,alignFile, targetsTopologyFile, fpout):#{{{
#     fptmp=open(queryTopoFile);
#     print fptmp.readlines();
#     fptmp.close();
    try:
        (queryID, queryAnnotation, queryTopology) = myfunc.ReadSingleFasta(queryTopoFile);
        # read in alignment
        alns = ReadNeedleAlignment(alignFile);

        # read in topologys
        (targetIDList, targetAnnotationList, targetTopoList) = myfunc.ReadFasta(targetsTopologyFile);

        # match and print the result
        print >> fpout, "#Number of alignments: %d" % len(targetIDList);

        for i in range (len(targetIDList)):
            seqID=targetIDList[i];
            alnseq1=alns[i]['alnseq1'];
            alnseq2=alns[i]['alnseq2'];
            topoaln1="";
            topoaln2="";

            if seqID != alns[i]['seqid2']:
                print >> sys.stderr, "seqID does not match, record %d" %i;

            cnt1=0;
            cnt2=0;
            for j in range(len(alnseq1)):
                if alnseq1[j] != '-':
                    if alnseq2[j] != '-':
                        topoaln1+=queryTopology[cnt1];
                        topoaln2+=targetTopoList[i][cnt2];
                    else:
                        topoaln1+=queryTopology[cnt1];
                        topoaln2+='-';
                else:
                    if alnseq2[j] != '-':
                        topoaln1+='-';
                        topoaln2+=targetTopoList[i][cnt2];
                    else:
                        topoaln1+='-';
                        topoaln2+='-';
                if alnseq1[j] != '-':
                    cnt1 +=1;
                if alnseq2[j] != '-':
                    cnt2 += 1;
            #print the result
            print >> fpout, "#Topology alignment %d" %( i+1);
            print >> fpout, ">%s" % queryAnnotation;
            print >> fpout, "%s" % topoaln1;
            print >> fpout, ">%s" % targetAnnotationList[i];
            print >> fpout, "%s" % topoaln2;
            print >> fpout;
    except: 
        print >>sys.stderr, "except for the function:%s"%sys._getframe().f_code.co_name ;
        raise ;
    return 0;
#}}}
if __name__ == '__main__' :
    # Check argv
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp();
        sys.exit(1);

    outFile="";
    targetsTopoFile="";
    queryTopoFile="";
    alignFile="";   #aln file in needle format


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
            elif sys.argv[i] == "-qtopo" or sys.argv[i] == "--qtopo":
                queryTopoFile=sys.argv[i+1];
                i = i + 2;
            elif sys.argv[i] == "-ttopo" or sys.argv[i] == "--ttopo":
                targetsTopoFile=sys.argv[i+1];
                i = i + 2;
            elif sys.argv[i] == "-aln" or sys.argv[i] == "--aln":
                alignFile=sys.argv[i+1];
                i = i + 2;
            elif sys.argv[i] == "-o" or sys.argv[i] == "--outfile":
                outFile=sys.argv[i+1];
                i = i + 2;
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i]);
                sys.exit(1);
        else:
            print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i]);
            sys.exit(1);
           

    if targetsTopoFile == "":
        print >> sys.stderr,"Error! targetsTopoFile not set.";
    if alignFile == "":
        print >> sys.stderr,"Error! alignFile not set.";
    if queryTopoFile == "":
        print >> sys.stderr,"Error! queryTopoFile not set.";

    fpout = sys.stdout;
    if outFile != "":
        fpout = open(outFile,"w");

    try :
        MatchTopoPairAln(queryTopoFile,alignFile, targetsTopoFile, fpout);
        if fpout != sys.stdout:
            fpout.close();

    except :
        print >>sys.stderr, "except for the function:%s"%sys._getframe().f_code.co_name ;
        raise ;

#!/usr/bin/env python
# label the fasta sequence according to the predicted topology and pairwise
# sequence alignment
import sys,re,os;
import myfunc ;

usage="""
Usage:  labeltopologyfastaseq.py [Options] -topo queryTopoFile -aln alignFile -fasta fastaFile_to_be_Labelled
Options:
  -o         <file>    : outputfile
  -h|--help            : print this help message and exit
Created 2010-08-29, updated 2010-08-29, Nanjiang
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
    except:
        print >>sys.stderr, "except for the function:%s"%sys._getframe().f_code.co_name ;
    return alns;
#}}}

def GetTopologyLabels(queryTopology, alns):
    topologyLabels={};
    try:
        for i in range(len(alns)):
            jTopoSeq=0;
            seqID=alns[i]['seqid2'];
            labelseq="";
            for j in range(len(alns[i]['alnseq1'])):
                if alns[i]['alnseq2'][j] != '-':
                    if alns[i]['alnseq1'][j] != '-':
                        labelseq+=queryTopology[jTopoSeq];
                    else:
                        labelseq+='.';
                if alns[i]['alnseq1'][j] != '-':   # if the query seq is not a gap, increament the topology seq index
                    jTopoSeq+=1;
            labelseq=re.sub( '[^M]', '.', labelseq);
            topologyLabels[seqID]=labelseq;
    except:
        print >>sys.stderr, "except for the function:%s"%sys._getframe().f_code.co_name ;
    return topologyLabels;

def Labeltopologyfastaseq(queryTopoFile,alignFile, fastaFile, fpout):#{{{
#     fptmp=open(queryTopoFile);
#     print fptmp.readlines();
#     fptmp.close();
    try:
        (queryID, queryAnnotation, queryTopology) = myfunc.ReadSingleFasta(queryTopoFile);
        # read in alignment
        alns = ReadNeedleAlignment(alignFile);
#         print alns;
        topologyLabels = GetTopologyLabels(queryTopology, alns);

        fpin = open(fastaFile, "r");
        lines = fpin.readlines();
        fpin.close();

        i = 0;
        while i < len(lines):
            line = lines[i];
            if line[0]== '>':
                seqID=myfunc.GetSeqIDFromAnnotation(line);
                aaSeq="";
                fpout.write("%s"%line);
                i = i + 1;
                while i < len(lines) and lines[i][0] != '>':
                    fpout.write("%s"%lines[i]);
                    aaSeq+=lines[i].strip();
                    i=i+1;
                fpout.write("/%s/\n"%topologyLabels[seqID]);
                if len(aaSeq) != len(topologyLabels[seqID] ):
                    print >> sys.stderr,"%s: length not match" % seqID;
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
    fastaFile="";
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
            elif sys.argv[i] == "-topo" or sys.argv[i] == "--topo":
                queryTopoFile=sys.argv[i+1];
                i = i + 2;
            elif sys.argv[i] == "-fasta" or sys.argv[i] == "--fasta":
                fastaFile=sys.argv[i+1];
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
           

    if fastaFile == "":
        print >> sys.stderr,"Error! fastaFile not set.";
    if alignFile == "":
        print >> sys.stderr,"Error! alignFile not set.";
    if queryTopoFile == "":
        print >> sys.stderr,"Error! queryTopoFile not set.";

    fpout = sys.stdout;
    if outFile != "":
        fpout = open(outFile,"w");

    try :
        Labeltopologyfastaseq(queryTopoFile,alignFile, fastaFile, fpout);
        if fpout != sys.stdout:
            fpout.close();

    except :
        print >>sys.stderr, "except for the function:%s"%sys._getframe().f_code.co_name ;
        raise ;

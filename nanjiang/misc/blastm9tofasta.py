#!/usr/bin/env python
# get the fasta sequence from blast output file m9
import sys,re,os;
import myfunc ;
import tempfile;

EvalueThreshold=1e-3;
usage="""
Usage:   blastm9tofasta.py [Options] [-i] blastfile
Note: blast file should be output by the option -m 9
Options:
  -i         <file>    : input file
  -o         <file>    : outputfile
  -blastdb   <str>     : blastdb for fastacmd to obtain the amino acid sequences
  -evalue    <float>   : set evalue threshold, default = 1e-3
  -round     <int>     : to use the hits of which iteration, default is using the last iteration
  -h|--help            : print this help message and exit
Created 2010-08-29, updated 2011-11-09, Nanjiang
"""

def PrintHelp():
    print usage;

def BlastM9toFasta(inFile, blastdb, iteration, fpout):#{{{
#first get line number of the Iteration record
    linenumber_iteration=[];
    fpin = open(inFile, "r");
    lines = fpin.readlines();
    fpin.close();
    cntLine=0;
    for line in lines:
        if line[0] == "#" and line.find("Iteration") >=0:
            linenumber_iteration.append(cntLine);
        cntLine +=1;

    if iteration > len(linenumber_iteration):
        iteration = len(linenumber_iteration);

    cntLine=0;
    records={};
    for line in lines:
        if cntLine > linenumber_iteration[iteration-1]:
            if line.find("BLASTP") >=0 :
                break;
            elif line[0] != "#":
                strs=line.split() ;
                hitID=strs[1];
                evalue=float(strs[10]);
                if not hitID in records and evalue <= EvalueThreshold:
                    records[hitID]={};
                    pid=float(strs[2]);
                    alnLength=int(strs[3]);
                    misMatch=int(strs[4]);
                    gapOpens=int(strs[5]);
                    qstart=int(strs[6]);
                    qend=int(strs[7]);
                    sstart=int(strs[8]);
                    send=int(strs[9]);
                    evalue=float(strs[10]);
                    bitscore=float(strs[11]);
                    records[hitID]['pid']=pid;
                    records[hitID]['alnLength']=alnLength;
                    records[hitID]['misMatch']=misMatch;
                    records[hitID]['gapOpens']=gapOpens;
                    records[hitID]['qstart']=qstart;
                    records[hitID]['qend']=qend;
                    records[hitID]['sstart']=sstart;
                    records[hitID]['send']=send;
                    records[hitID]['evalue']=evalue;
                    records[hitID]['bitscore']=bitscore;
        cntLine+=1;

    #create a new dict to sort the idlist by evalue in  ascending order
    tmpdict={};
    for key in records.keys():
        tmpdict[key]=records[key]['evalue'];
    keys=sorted(tmpdict);
    keys.sort(key=tmpdict.get, reverse=False);

    tmpFile=tempfile.mktemp()
    fptmp=open(tmpFile, "w");
    for key in keys:
        print >> fptmp, "%s"%key;
    fptmp.close();
    
    tmpFastaFile=tempfile.mktemp()
    os.system ("fastacmd -i %s -d %s > %s" %(tmpFile, blastdb, tmpFastaFile));

#replace annotations in the fasta file
    (annotationList,seqList)= myfunc.ReadFasta_without_id(tmpFastaFile);
    numSeq = len(seqList);
    for i in xrange(numSeq):
        anno = annotationList[i];
        seq = seqList[i];
        firstword=anno.split()[0];
        if firstword.find('|') >= 0:
            strs=firstword.split('|');
#modified according to _get_seq_identifier in bioperl
            if strs[0] == "lcl" : 
                seqID=strs[1];
            elif strs[0] == 'gnl':
                seqID=strs[2];
            else:
                seqID=firstword;
        else:
            seqID=firstword;
        if seqID in records:
            annoLine="%s evalue=%g pid=%g alnLength=%d bitscore=%g" %(seqID, records[seqID]['evalue'], records[seqID]['pid'], records[seqID]['alnLength'],records[seqID]['bitscore']);
            fpout.write(">%s\n"%annoLine);
            fpout.write("%s\n"%seq)
        else:
            print >> sys.stderr, "seqID", seqID, "does not found in records. Ignore.";

    os.remove(tmpFile);
    os.remove(tmpFastaFile);
    return 0;
#}}}

if __name__ == '__main__' :
    # Check argv
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp();
        sys.exit(1);

    outFile="";
    inFile="";
    iteration=99999; # by default, using the last round
    blastdb="";

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
            elif sys.argv[i] == "-blastdb" or sys.argv[i] == "--blastdb":
                blastdb=sys.argv[i+1];
                i = i + 2;
            elif sys.argv[i] == "-evalue" or sys.argv[i] == "--evalue":
                evalue=float(sys.argv[i+1]);
                i = i + 2;
            elif sys.argv[i] == "-round" or sys.argv[i] == "--round":
                iteration=int(sys.argv[i+1]);
                i = i + 2;
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i]);
                sys.exit(1);
        else:
            inFile=sys.argv[i];
            i+=1;
           

    if inFile == "":
        print >> sys.stderr,"Error! Input file not set.";
    if blastdb == "":
        print >> sys.stderr,"Error! blastdb not set.";

    fpout = sys.stdout;
    if outFile != "":
        fpout = open(outFile,"w");

    try :
        BlastM9toFasta(inFile,blastdb, iteration, fpout);
        if fpout != sys.stdout:
            fpout.close();

    except :
        print >>sys.stderr, "except for the input file: %s" % inFile;
        raise ;

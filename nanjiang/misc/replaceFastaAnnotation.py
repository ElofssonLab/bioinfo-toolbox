#!/usr/bin/env python
# replace the annotation of the fasta file by supplied annotation files
import sys,re,os;
import myfunc;

usage="""
Usage:  replaceFastaAnnotation.py -fasta fastaFile -anno annotationFile 
Note: the annotationFile is one record per line and starts with '>'
Options:
  -o         <file> : outputfile
  -h|--help         : print this help message and exit
Created 2010-08-20, updated 2010-08-20, Nanjiang
"""

def PrintHelp():
    print usage;

def ReplaceFastaAnnotation(fastaFile, annotationFile, fpout):
# this version read in the annotation file first 
    cntAnno=0;
    cntLineFasta=0;
    lineFa ="";
    lineAnno="";

    fpFa = open(fastaFile, "r");
    fpAnno = open(annotationFile, "r");
    lineAnno = fpAnno.readline(); cntAnno+=1;

    while lineAnno:
        lineAnno = lineAnno.rstrip('\n').strip();
        if lineAnno and lineAnno[0] == ">":
            if not (lineFa and lineFa[0]) == ">":
                lineFa = fpFa.readline();cntLineFasta+=1;
                lineFa = lineFa.rstrip('\n').strip();
            if lineFa and lineFa[0] == ">":
                #starting a new sequnce
                idFa=myfunc.GetSeqIDFromAnnotation(lineFa);
                idAnno=myfunc.GetSeqIDFromAnnotation(lineAnno);
                if idFa != idAnno:
                    print >> sys.stderr,"annotation line does not match, annoRecord=%d, fastaLine=%d" %(cntAnno+1, cntLineFasta+1);
                    print >> sys.stderr,"idFa=%s, idAnno=%s" %(idFa, idAnno);
                    sys.exit(1);
                print >> fpout,"%s" % lineAnno;
                lineFa = fpFa.readline();cntLineFasta+=1;
                lineFa = lineFa.rstrip("\n").strip();
                while lineFa and lineFa[0] != ">":
                    print >> fpout,"%s" %lineFa;
                    lineFa = fpFa.readline();cntLineFasta+=1;
                    lineFa = lineFa.rstrip("\n").strip();
        lineAnno = fpAnno.readline(); cntAnno+=1;
    fpFa.close();
    fpAnno.close();
    return 0;



if __name__ == '__main__' :
    # Check argv
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp();
        sys.exit(1);

    outFile="";
    fastaFile="";
    annotationFile="";

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
            elif sys.argv[i] == "-fasta" or sys.argv[i] == "--fasta":
                fastaFile=sys.argv[i+1];
                i = i + 2;
            elif sys.argv[i] == "-anno" or sys.argv[i] == "--anno":
                annotationFile=sys.argv[i+1];
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
        print >> sys.stderr,"Error! fastafile not set.";
    if annotationFile == "":
        print >> sys.stderr,"Error! annotation file not set.";

    fpout = sys.stdout;
    if outFile != "":
        fpout = open(outFile,"w");

    try :
        ReplaceFastaAnnotation(fastaFile,annotationFile, fpout);
        if fpout != sys.stdout:
            fpout.close();

    except :
        print >>sys.stderr, "except for the input file: %s" % fastaFile;
        raise ;

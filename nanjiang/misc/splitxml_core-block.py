#!/usr/bin/env python
# 
import sys,re,os;
DEBUG=0;


usage="""
Usage:  splitxml_core-block.py [Options] [-i] core-block-xml-file
Options:
  -i       <file> : input file
  -outpath <dir>  : output the result to dir, default=./
  -ext     <str>  : set the file extension, default = aa
  -q              : quiet mode, do not write any messages
  -nameseq        : name the splitted files by $rootname_i, i = 0,1,2...
                  : this is useful when it is the filename extracted from the 
                  : annotation line having a duplicated name
  -h|--help       : print this help message and exit

Created 2010-10-22, updated 2010-10-22, Nanjiang

Examples:
    splitfasta.py example.fasta 
    splitfasta.py example.fasta -outpath outdir
"""

def PrintHelp():
    print usage;


def SplitXML_Core_Block(inFile, outpath):#{{{
# The faster version
    rootname=os.path.basename(os.path.splitext(inFile)[0]);
    fpin = open(inFile, "r");
    cntLine = 0;
    line = fpin.readline();
    cntLine +=1;
    if DEBUG: sys.stdout.write("%d: %s"%(cntLine, line));
    cntFile=0;
    while line:
        if line.find (':') >= 0:
            cntFile +=1;
            line = fpin.readline();
            cntLine +=1;
            if DEBUG: sys.stdout.write("%d:%s"%(cntLine, line));
            outFile=line.strip();
            if outFile[0] == '>' or outFile[0] == '<':
                sys.exit(1);
            outFile=outpath + os.sep + outFile;
            if not isQuiet:
                print "outfile: ", outFile;
            line = fpin.readline();
            cntLine +=1;
            if DEBUG: sys.stdout.write("%d:%s"%(cntLine, line));
            fpout = open (outFile, "w");
            line = fpin.readline();
            if DEBUG: print line;
            while line:
                if line[0] != ':':
                    fpout.write("%s"%line);
                else:
                    fpout.close();
                    break;
                line = fpin.readline();
        else:
            line = fpin.readline();
            cntLine +=1;
            if DEBUG: sys.stdout.write("%d: else: %s"%(cntLine, line));


    fpin.close();
        
    return cntFile;

#}}}

if __name__ == '__main__' :
    # Check argv
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp();
        sys.exit(1);

    outpath="./";
    inFile="";
    extension="aa";
    isQuiet=False;
    isNameFileSequentially=False;

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
            elif sys.argv[i] == "-ext" or sys.argv[i] == "--ext":
                extension=(sys.argv[i+1]);
                i = i + 2;
            elif sys.argv[i] == "-outpath" or sys.argv[i] == "--outpath":
                outpath=sys.argv[i+1];
                i = i + 2;
            elif sys.argv[i] == "-q" or sys.argv[i] == "--q" or sys.argv[i] == "--quiet":
                isQuiet=True;
                i = i + 1;
            elif sys.argv[i] == "-nameseq" or sys.argv[i] == "--nameseq":
                isNameFileSequentially=True;
                i = i + 1;
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i]);
                sys.exit(1);
        else:
            inFile=sys.argv[i];
            i+=1;
           

    if inFile == "":
        print >> sys.stderr,"Error! Input file not set.";

    os.system("mkdir -p %s"%outpath);

    try :
        numSeq = SplitXML_Core_Block(inFile,outpath);
        if not isQuiet:
            print >> sys.stdout,"%d sequences are splitted and output to %s"%(numSeq, outpath);

    except :
        print >>sys.stderr, "except for the input file: %s" % inFile;
        raise ;

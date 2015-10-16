#!/usr/bin/env python
# 
import os
import sys
import myfunc
usage="""
usage:  fasta2oneline.py Fasta-file [-o outfile]

convert fasta format sequence file to one line format
each line contains
ID sequence

Options:
  -q             Quiet mode
  -h, --help     Print this help message and exit

Created 2011-10-21, updated 2011-10-21, Nanjiang Shu 
"""
BLOCK_SIZE=100000;

def PrintHelp():
    print usage;

def WriteOnelineSeq(seqWithAnno, idSet, fpout):#{{{
    "Write sequence in oneline, sequences with redundant IDs are ignored"
    seqid = myfunc.GetSeqIDFromAnnotation(seqWithAnno[0:seqWithAnno.find('\n')]);
    aaSeq = seqWithAnno[seqWithAnno.find("\n"):].replace('\n','').replace(' ','');
    if seqid not in idSet:
        seqWithAnno+="\n";
        fpout.write("%s %s\n"%(seqid, aaSeq));
        idSet.add(seqid);
    
#}}}
def Fasta2Oneline(infile,outfile): #{{{
    fpout=sys.stdout;
    if outfile != "":
        try:
            fpout = open(outfile, "wb");
        except IOError:
            fpout = sys.stdout;
            pass


    fpin = None;
    try:
        fpin=open(infile,"rb");
    except IOError:
        print >> sys.stderr, "Failed to open file %s for read"%(infile);
        raise;

    idSet = set([]);

    isFirstSeq=True;
    totalLength=0;
    buff = fpin.read(BLOCK_SIZE);
    brokenSeqWithAnnoLine=""; ##for the annotation line broken by BLOCK read
    while buff:
        beg=0;
        end=0;
        while 1:
            if brokenSeqWithAnnoLine:
                if brokenSeqWithAnnoLine[len(brokenSeqWithAnnoLine)-1] == "\n":
                    end=buff.find(">");
                else:
                    end=buff.find("\n>");
                if end >= 0:
                    seqWithAnno = brokenSeqWithAnnoLine + buff[0:end];
                    WriteOnelineSeq(seqWithAnno,idSet, fpout);
                    brokenSeqWithAnnoLine = "";
                    beg=end;
                else:
                    brokenSeqWithAnnoLine += buff;
                    break;

            beg=buff.find(">",beg);
            end=buff.find("\n>",beg+1);
            if beg >= 0:
                if end >=0:
                    seqWithAnno=buff[beg:end];
                    WriteOnelineSeq(seqWithAnno,idSet, fpout);
                    beg=end;
                else:
                    brokenSeqWithAnnoLine=buff[beg:];
                    break;
            else:
                break;
        buff = fpin.read(BLOCK_SIZE);
    
    if brokenSeqWithAnnoLine:
        seqWithAnno=brokenSeqWithAnnoLine;
        WriteOnelineSeq(seqWithAnno,idSet, fpout);

    fpin.close();   
    if fpout != sys.stdout and fpout != None:
        fpout.close();
    
    return 0;

#}}}

if __name__ == '__main__' :
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        sys.exit()

    isQuiet=False;
    outfile="";
    infile="";

    i = 1;
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            infile=sys.argv[i];
            isNonOptionArg=False;
            i += 1;
        elif sys.argv[i] == "--":
            isNonOptionArg=True;
            i += 1;
        elif sys.argv[i][0] == "-":
            if sys.argv[i] ==  "-h" or  sys.argv[i] == "--help":
                PrintHelp();
                sys.exit();
            elif sys.argv[i] == "-o" or sys.argv[i] == "--o" or sys.argv[i] == "-outfile" or sys.argv[i] == "--outfile":
                outfile=sys.argv[i+1];
                i += 2;
            elif sys.argv[i] == "-q":
                isQuiet=True;
                i += 1;
            else:
                print >> sys.stderr, "Error! Wrong argument:", sys.argv[i];
                sys.exit()
        else:
            infile=sys.argv[i];
            i += 1

    if infile == "" or not os.path.exists(infile):
        print >> sys.stderr, "infile not set";
        sys.exit(1);

    Fasta2Oneline(infile,outfile);


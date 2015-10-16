
#!/usr/bin/env python
# get the length of the sequences in the fasta file
import sys,os;
import myfunc ;

BLOCK_SIZE=100000;
isPrintID=False;
isJustPrintSum=False;

usage="""
Usage:   getseqlen.py [Options] fastafile
Options:
  -i|--printid            : print id as well
  -o               <file> : outputfile
  -just-print-sum         : just print the total number of amino acids
  -bs|--block-size <int>  : size for blocks when reading file, default = 50000
  -h|--help               : print this help message and exit

Created 2010-09-02, updated 2010-10-31, Nanjiang
Examples:
    getseqlen.py test.fa

Note: this python code is about 7 times faster than the perl script getseqlen.pl
"""

def PrintHelp():
    print usage;

def GetFromRawSeq(seqWithAnno):#{{{
#     begseq=seqWithAnno.find("\n");
#     seq=seqWithAnno[begseq:];
#     seq=seq.replace('\n','').replace(' ','');
    length=len(seqWithAnno[seqWithAnno.find("\n"):].replace('\n','').replace(' ',''));
    if not isJustPrintSum:
        if isPrintID:
            seqID=myfunc.GetSeqIDFromAnnotation(seqWithAnno);
            fpout.write("%s\t"%seqID);
        fpout.write("%d\n" % length);
    return length;
#}}}
def Getseqlen(inFile, fpout):#{{{
# The faster version
    isFirstSeq=True;
    totalLength=0;
    fpin = open(inFile, "r");
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
                    length=GetFromRawSeq(seqWithAnno);
                    totalLength +=length;
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
                    length=GetFromRawSeq(seqWithAnno);
                    totalLength +=length;
                    beg=end;
                else:
                    brokenSeqWithAnnoLine=buff[beg:];
                    break;
            else:
                break;

        buff = fpin.read(BLOCK_SIZE);
    fpin.close();
    if brokenSeqWithAnnoLine:
        seqWithAnno=brokenSeqWithAnnoLine;
        length=GetFromRawSeq(seqWithAnno);
        totalLength +=length;

    if isJustPrintSum:
        fpout.write("%d\n"%totalLength);
        
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
            elif sys.argv[i] == "-i" or sys.argv[i] == "--printid":
                isPrintID=True;
                i+=1;
            elif sys.argv[i] == "-just-print-sum" or sys.argv[i] == "--just-print-sum":
                isJustPrintSum=True;
                i+=1;
            elif sys.argv[i] == "-o" or sys.argv[i] == "--outfile":
                outFile=sys.argv[i+1];
                i = i + 2;
            elif sys.argv[i] == "-bs" or sys.argv[i] == "--block-size" or sys.argv[i] == "-block-size":
                BLOCK_SIZE=int(sys.argv[i+1]);
                if BLOCK_SIZE < 0:
                    print >> sys.stderr,"Error! BLOCK_SIZE should >0";
                    sys.exit(1);
                i = i + 2;
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i]);
                sys.exit(1);
        else:
            inFile=sys.argv[i];
            i+=1;
           

    if inFile == "":
        print >> sys.stderr,"Error! Input file not set.";

    fpout = sys.stdout;
    if outFile != "":
        fpout = open(outFile,"w");

    try :
        Getseqlen(inFile,fpout);
        if fpout != sys.stdout:
            fpout.close();

    except :
        print >>sys.stderr, "except for the input file: %s" % inFile;
        raise ;

#!/usr/bin/env python
# convert the fasta format to swissprot format, add the information (e.g.
# secondary structure) if supplied
# predicted secondary structures
import sys,re,os,math;
import myfunc;

DEBUG=False;

p_threshold_C  = 0.4;
p_threshold_HE = 0.4;

usage="""
Usage:  fasta2swissprot.py [Options] fasta-seq-files ...

Note: the output file will be named as $outpath/$rootname.swiss

Options:
  -l         <file>: set the list file for fasta files
  -outpath   <dir> : set the output path for fasta files with gap penalties added, 
                   : default=the same as the fasta-seq-file
  -sspath    <dir> : path for the predicted secondary structures
                   : it must be set, and the secondary structure files is stored as $sspath/$rootname/$rootname_i
  -fss       0|1   : format for the predicted secondary structures, default=0
                   : 0 for psipred output
  -q               : quiet mode
  -h|--help        : print this help message and exit
  
  Other parameters related to the model
  -p-threshold-c <float>
  -p-threshold-he <float>

Created 2010-11-08, updated 2010-11-08, Nanjiang
"""

def PrintHelp():
    print usage;

def AddFileList(fileList, inFile):#{{{
# return fileList
    fpin = open(inFile, "r");
    lines = fpin.readlines();
    fpin.close();
    for line in lines:
        line=line.strip();
        if line:
            fileList.append(line);
    return (fileList);
#}}}

def ReadInPredSS(inFile, ssFormat):#{{{
    fpin = open(inFile, "r");
    lines = fpin.readlines();
    fpin.close();
    predSS={};
    if ssFormat == 0:
        #aaSeq="";
        ssSeq="";
        prob_C=[];
        prob_H=[];
        prob_E=[];
        for line in lines:
            strs=line.split();
            if len(strs) == 6 and strs[0].isdigit():
                #aaSeq+=strs[1];
                ssSeq+=strs[2];
                prob_C.append(float(strs[3]));
                prob_H.append(float(strs[4]));
                prob_E.append(float(strs[5]));
        #predSS['aaSeq']=aaSeq;
        predSS['ss']=ssSeq;
        predSS['C']=prob_C;
        predSS['H']=prob_H;
        predSS['E']=prob_E;
    return predSS;
#}}}

def GetSSEList(ssFile, seqLength, ssFormat):#{{{
    predSS = ReadInPredSS(ssFile, ssFormat);
    if seqLength != len(predSS['ss']):
        print >> sys.stderr, "Error! aaSeqLength = %d, but the ssSeqLength = %d, for file %s" %(seqLength, len(predSS['ss']), ssFile);
        sys.exit(1);
    sseList=[];
    cntSSE=0;
    ssSeq=predSS['ss'];
    i = 0;
    while i < seqLength:
        if ssSeq[i] == 'H' or ssSeq[i] == 'E':
            state = ssSeq[i];
            j = 0;
            pArray = [];
            while  i+j < seqLength and ssSeq[i+j] == state:
                pArray.append(predSS[state][i+j]);
                j += 1;
            if j >= 3:
                avgP =  sum(pArray, 0.0)/len(pArray); 
                if avgP > p_threshold_HE:
                    sseList.append({});
                    if state == 'H':
                        sseList[cntSSE]['type']='helix';
                    else:
                        sseList[cntSSE]['type']='strand';
                    sseList[cntSSE]['start']=i;
                    sseList[cntSSE]['end']=i+j;
                    cntSSE +=1;
            i+=j;
        else:
            i+= 1;
    return sseList;
#}}}
def Fasta2SwissProt(fastaFile, ssPath, outPath,ssFormat):#{{{
# add gap penalties
    rootname = os.path.basename(os.path.splitext(fastaFile)[0]);
    inFilePath=os.path.dirname(fastaFile);
    if inFilePath == "":
        inFilePath="./";
    (annotationList, seqList) = myfunc.ReadFasta_without_annotation(fastaFile);

    if outPath == "":
        localOutPath=inFilePath;
    else:
        localOutPath=outPath;

    outFile="%s/%s.swiss"%(localOutPath,rootname)
    fpout = open(outFile, "w");

    if ssFormat ==0:
        for i in range(len(seqList)) :
            aaSeq=seqList[i];
            seqLength = len(aaSeq);

            if DEBUG: 
                print >> sys.stdout, "%d:%s" %(i, annotationList[i])

            sseList = [];
            if ssPath != "":
                ssFile="%s/%s/%s_%d.ss2"%(ssPath, rootname, rootname, i);
                (sseList) = GetSSEList(ssFile, seqLength, ssFormat);

# write out the result
            fpout.write("%-4s %s\n"%("ID", annotationList[i]));
            for sse in sseList:
                if sse['type']=='helix':
                    fpout.write("%-4s %s "%("FT","HELIX"));
                elif sse['type']=='strand':
                    fpout.write("%-4s %s "%("FT", "STRAND"));
                fpout.write("%d %d\n"% (sse['start'], sse['end']));

            fpout.write("%-4s SEQUENCE %d AA;\n"%("SQ", seqLength));
            j=0;
            cntBlock=0;
            fpout.write("%5s"%(""));
            while j < seqLength:
                fpout.write("%s"% aaSeq[j:j+10]);
                j +=10;
                cntBlock+=1;
                if j >= seqLength:
                    fpout.write("\n");
                    break;
                else:
                    if cntBlock < 6:
                        fpout.write(" ");
                    else: 
                        fpout.write("\n%5s"%(""));
                        cntBlock =0;
            fpout.write("//\n");
            
        fpout.close();
            
    return len(seqList);
#}}}

if __name__ == '__main__' :
    # Check argv
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp();
        sys.exit(1);

    outPath="";
    listFile="";  
    ssPath="";
    fastaFileList=[];
    ssFormat=0; # 0 for psipred
    isQuiet=False;

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
            elif sys.argv[i] ==  "-outpath" or  sys.argv[i] == "--outpath":
                outPath=sys.argv[i+1];
                i+=2;
            elif sys.argv[i] ==  "-sspath" or  sys.argv[i] == "--sspath":
                ssPath=sys.argv[i+1];
                i+=2;
            elif sys.argv[i] ==  "-l" or  sys.argv[i] == "--list":
                listFile=sys.argv[i+1];
                i+=2;
            elif sys.argv[i] ==  "-ssf" or  sys.argv[i] == "--ssf":
                ssFormat=int(sys.argv[i+1]);
                i+=2;
            elif sys.argv[i] ==  "-q" or  sys.argv[i] == "--q" or sys.argv[i] == "--quiet":
                isQuiet=True;
                i+=1;
            elif sys.argv[i] ==  "-p-threshold-c" or  sys.argv[i] == "--p-threshold-c":
                p_threshold_C=float(sys.argv[i+1]);
                i+=2;
            elif sys.argv[i] ==  "-p-threshold-he" or  sys.argv[i] == "--p-threshold-he":
                p_threshold_HE=float(sys.argv[i+1]);
                i+=2;
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i]);
                sys.exit(1);
        else:
            fastaFileList.append(sys.argv[i]);
            i+=1;
           

    if len(fastaFileList) == 0 and listFile == "":
        print >> sys.stderr,"Error! no input is set";
        sys.exit(1);

    try :
        if outPath != "" :
            os.system("mkdir -p %s" % outPath);
        if listFile != "":
            fastaFileList = AddFileList(fastaFileList, listFile);

        cntFile= 0;
        for fastaFile in fastaFileList:
            numseq = Fasta2SwissProt(fastaFile, ssPath, outPath, ssFormat);
            if not isQuiet:
                print >> sys.stdout, "%d\t%d sequences in the file %s have been converted to swissprot format" %(cntFile, numseq, fastaFile);
            cntFile+=1;

    except :
        print >>sys.stderr, "except for the function:%s"%sys._getframe().f_code.co_name ;
        raise ;

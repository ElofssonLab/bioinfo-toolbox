#!/usr/bin/env python
# obtain the multiple sequence alignment of (predicted) seoncdary structures
# according to the sequence alignment
import sys,re,os;

# not finished
BLOCK_SIZE=100000;
CHAR_GAP='-';
LENGTH_ALN=80;

usage="""
Usage:  matchMSAss.py [Options] msaFile  ssFile
Options:
  -l         <file>     : set the list file, one line per record, each line contains "msaFile ssFile"
  -o         <file>     : outputfile
  -msaformat msf|aln|fa : set the format of the msaFile, if not set, it will be detected automatically
  -ssformat  psipred|fa : set the format of the ss file, if not set, it will be detected automatically
  -outformat aln|fa     : set the format of the output file, default=aln
  -nchar        <int>   : set maximal number of residues printed in each line when the outformat is aln, default=80
  -h|--help             : print this help message and exit

Created 2010-10-27, updated 2010-11-16, Nanjiang
"""

def PrintHelp():
    print usage;

def DetectSSFormat(inFile):#{{{
#return fileFormat
    try:
        fileFormat="";
        fpin = open(inFile);
        line = fpin.readline();
        while line:  #read in the first non blank line
            line = line.strip();
            if line:
                break;
        fpin.close();
        if line[0] == '>':
            fileFormat='fa';
        elif line.find("#PSIPRED") == 0:
            fileFormat='psipred'
        else:
            print >> sys.stderr,"unrecognized ss format for the file %s"%inFile;
            sys.exit(1);
        return fileFormat;
    except:
        print >> sys.stderr,"can not open the file %s"%inFile;
        raise;
#}}}
def DetectMSAFormat(inFile):#{{{
#return fileFormat
    try:
        fileFormat="";
        fpin = open(inFile);
        line = fpin.readline();
        while line:  #read in the first non blank line
            line = line.strip();
            if line:
                break;
        fpin.close();
        if line[0] == '>':
            fileFormat='fa';
        elif line.upper().find("CLUSTAL") >= 0:
            fileFormat='aln';
        elif line.find("PileUp") >= 0:
            fileFormat='msf';
        else:
            print >> sys.stderr,"unrecognized msa format for the file %s"%inFile;
            sys.exit(1);
        return fileFormat;
    except:
        print >> sys.stderr,"can not open the file %s"%inFile;
        raise;
    #}}}

def ReadInMSA_msf(lines):#{{{
    msa=[];
    cntAnnoLine=0;
    isAlignRegion=False;
    cntSeq = 0;
    for line in lines:
        line=line.strip();
        if line:
            if not isAlignRegion: 
                if line.find("Name:") == 0:
                    msa.append({});
                    line=line.lstrip("Name:").strip();
                    msa[cntAnnoLine]['annoLine'] = line;
                    msa[cntAnnoLine]['id']=line.split()[0];
                    msa[cntAnnoLine]['seq']="";
                    cntAnnoLine +=1;
                elif line.find("//") == 0:
                    isAlignRegion=True;
                    cntSeq=0;

            elif isAlignRegion:
                strs=line.split();
                if not strs[len(strs)-1].isdigit():
                    if msa[cntSeq]['id'] == strs[0]:
                        for i in range(1,len(strs)):
                            msa[cntSeq]['seq'] += (strs[i].replace('.','-')); #set the gap symbol to '-'
                    else:
                        print >> sys.stderr, "%s: msf format error"%msfFile;

                    cntSeq +=1;
                    if cntSeq == cntAnnoLine:
                        cntSeq = 0;
    return msa;
                #}}}
def ReadInSS_psipred(lines):#{{{
    secStruc=[];
    cntSeq=0;
    i = 0;
#     for line in lines:
#         sys.stdout.write("%s"%line);

    while i < len(lines):
        if lines[i].find("# PSIPRED VFOR") == 0:
# a new record
            i += 2;
            ss = "";
            while i < len(lines) and lines[i][0] != '#':
                if len(lines[i]) > 15:
                    ss+=(lines[i][7]);
                i += 1;
            secStruc.append({}) ;
            secStruc[cntSeq]['ss'] =ss;
            cntSeq+=1;
        i += 1;
    return secStruc;
#}}}
def ReadInSS_fasta(lines):#{{{
    secStruc=[];
    cntSeq=0;
    i = 0;
    while i < len(lines):
        if lines[i][0] == '>':
# a new record
            ss = "";
            i+=1;
            while i < len(lines) and lines[i][0]!= '>':
                ss+=(lines[i].strip());
                i += 1;
            secStruc.append({}) ;
            secStruc[cntSeq]['ss'] =ss;
            cntSeq+=1;
        else:
            i += 1;
    return secStruc;
#}}}

def ReadInMSA(msaFile,msaFormat):#{{{
    msa=[];
    try:
        fpin = open (msaFile,"r");
        lines=fpin.readlines();
        fpin.close();
    except:
        print >> sys.stderr,"can not open the file %s"%msaFile;
        raise;
    if msaFormat=='msf':
        msa = ReadInMSA_msf(lines);
    elif msaFormat == 'aln':
        print >> sys.stderr,"not implemented";
        sys.exit(1);
    elif msaFormat == 'fa':
        print >> sys.stderr,"not implemented";
        sys.exit(1);
    return msa;
#}}}
def ReadInSecondaryStructure(ssFile,ssFormat):#{{{
    secStruc=[];
    try:
        fpin = open (ssFile,"r");
        lines=fpin.readlines();
        fpin.close();
    except:
        print >> sys.stderr,"can not open the file %s"%ssFile;
        raise;
    if ssFormat=='psipred':
        secStruc = ReadInSS_psipred(lines);
    elif ssFormat == 'fa':
        secStruc = ReadInSS_fasta(lines);
    return secStruc;
#}}}
def GetMSA_SS(secStruc,msa):#{{{
    #return msaSS
    msaSS=[];
#     print msaSS;
    msaLength=len(msa[0]['seq']);
    for i in range(len(msa)):
        ss=secStruc[i]['ss'];
        msaseq=msa[i]['seq'];
        ssMSA="";
        cnt=0;
        for j in range(msaLength):
            if msaseq[j]==CHAR_GAP:
                ssMSA+=(CHAR_GAP);
            else: 
                ssMSA+=(ss[cnt]);
                cnt +=1;
        msaSS.append({});
        msaSS[i]['ss_align'] =ssMSA;
        msaSS[i]['id']=msa[i]['id'];
#         print msaSS[i]['ss_align'];
#         print msaSS[i]['id'];
#     print msaSS;
    return msaSS;
#}}}
def WriteMSA_SS(msaSS, fpout, outFormat):#{{{
    numMSA=len(msaSS);
    if outFormat == "fa":
        for i in range(numMSA):
            fpout.write(">%s\n"%msaSS[i]['id']);
            fpout.write("%s\n"%msaSS[i]['ss_align']);
    elif outFormat == "aln":
        beg=0;
        end=LENGTH_ALN;
        alnLength=len(msaSS[0]['ss_align']) ;
        idLengthList=[];
        for i in range(numMSA): idLengthList.append(len(msaSS[i]['id']));
        maxIDLength=max(idLengthList);
        while beg < alnLength:
            for i in range(numMSA):
                fpout.write("%-*s   %s\n"%(maxIDLength, msaSS[i]['id'], msaSS[i]['ss_align'][beg:end]));
            beg=end;
            end+=LENGTH_ALN;
            print >> fpout;
            print >> fpout;
#}}}

def MatchMSASS(msaFile, ssFile, fpout):#{{{
    if not isMSAFormatSet:
        msaFormat = DetectMSAFormat(msaFile);
    if not isSSFormatSet:
        ssFormat = DetectSSFormat(ssFile);
    msa = ReadInMSA(msaFile, msaFormat);
#    print msa;
    secStruc = ReadInSecondaryStructure(ssFile,ssFormat);
    #print secStruc;

    if len(msa) != len(secStruc) :
        print >> sys.stderr,"numseq in msa (%d, %s) and numseq in secStruc (%d, %s) is not equal, exit!"%(len(msa),msaFile, len(secStruc), ssFile);
        sys.exit(1);
    msaSS = GetMSA_SS(secStruc, msa);

    WriteMSA_SS(msaSS,fpout,outFormat );
#}}}

if __name__ == '__main__' :
    # Check argv

    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp();
        sys.exit(1);

    outFile="";
    msaFile="";   # MSA file
    ssFile="";    # secondary structure file

    ssFormat="";  #auto detect by default
    msaFormat=""; #auto detect by default
    outFormat="aln";

    listFile="";

    isSSFormatSet = False;
    isMSAFormatSet = False;

    cntInputFile=0;
    
    i = 1;#{{{ argument parsing
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
            elif sys.argv[i] == "-l" or sys.argv[i] == "--l" or sys.argv[i] == "--list":
                listFile=sys.argv[i+1];
                i = i + 2;
            elif sys.argv[i] == "-msa" or sys.argv[i] == "--msa":
                msaFile=sys.argv[i+1];
                i = i + 2;
            elif sys.argv[i] == "-ss" or sys.argv[i] == "--ss":
                ssFile=sys.argv[i+1];
                i = i + 2;
            elif sys.argv[i] == "-o" or sys.argv[i] == "--outfile":
                outFile=sys.argv[i+1];
                i = i + 2;
            elif sys.argv[i] == "-nchar" or sys.argv[i] == "--nchar":
                LENGTH_ALN=int(sys.argv[i+1]);
                i = i + 2;
            elif sys.argv[i] == "-msaformat" or sys.argv[i] == "--msaformat" or sys.argv[i] == "-mf":
                msaFormat=sys.argv[i+1].lower();
                if not (msaFormat == "aln" or msaFormat == "msf" or msaFormat == "fa"):
                    print >> sys.stderr,"Error! msaFormat should be aln, or msf or fa";
                    sys.exit(1);
                isMSAFormatSet = True;
                i = i + 2;
            elif sys.argv[i] == "-ssformat" or sys.argv[i] == "--ssformat" or sys.argv[i] == "-sf":
                ssFormat=sys.argv[i+1].lower();
                if not (ssFormat == "psipred" or ssFormat == "fa" ):
                    print >> sys.stderr,"Error! ssFormat should be psipred or fa";
                    sys.exit(1);
                isSSFormatSet=True;
                i = i + 2;
            elif sys.argv[i] == "-outformat" or sys.argv[i] == "--outformat" or sys.argv[i] == "-of":
                outFormat=sys.argv[i+1].lower();
                if not (outFormat== "aln" or outFormat == "fa" ):
                    print >> sys.stderr,"Error! outFormat should be aln or fa";
                    sys.exit(1);
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
            if cntInputFile == 0:
                msaFile=sys.argv[i];
            else:
                ssFile=sys.argv[i];
            cntInputFile += 1;
            i+=1;
           #}}}

    pairList=[];
    cntPair = 0;
    if listFile == "" and (msaFile == "" or ssFile == "") :
        print >> sys.stderr,"%s: Error! input not set."%(sys.argv[0]);
        sys.exit(1);
    
    if (msaFile != "" and ssFile != ""):
        pairList.append([]);
        pairList[cntPair].append(msaFile);
        pairList[cntPair].append(ssFile);
        cntPair +=1;
    if listFile != "":
        try:
            fpin=open(listFile);
            lines=fpin.readlines(listFile);
            fpin.close();
            for line in lines:
                line=line.strip();
                if line and line[0] != '#':
                    strs=line.split();
                    if len(strs) == 2:
                        pairList[cntPair].append(str[0]);
                        pairList[cntPair].append(str[1]);
                        cntPair +=1;
        except:
            print >> sys.stderr,"can not open the list file %s"%listFile;
            raise;


    fpout = sys.stdout;
    if outFile != "":
        fpout = open(outFile,"w");

    try :
        for pair in pairList:
            msaFile=pair[0];
            ssFile=pair[1];
            MatchMSASS(msaFile, ssFile, fpout);

        if fpout != sys.stdout:
            fpout.close();

    except :
        print >>sys.stderr, "except in main " ;
        raise ;


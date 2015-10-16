#!/usr/bin/env python
# convert modhmm xml output file into txt file
import sys,re,os;

# Note 2010-08-26, This code is obsolete now, use the xslt still, but split the
# fasta sequence into smaller pieces if it is too big. 
# the xml file should be read in by better languages to avoid errors.
# this code works only with 
# modhmms with -L option on, not maintained anymore
usage="""
Usage:   my_modhmmxml2txt.py [Options] [-i] xmlfile
Options:
  -i         <file> : input file
  -o         <file> : outputfile
  -h|--help         : print this help message and exit
  -m|--method 1|2|3 : select the method, default = 3, the fastest version
Created 2010-08-23, updated 2010-08-23, Nanjiang
"""

def PrintHelp():
    print usage;

def Modhmmxml2txt(inFile,fpout):#{{{
    cntSeq=0;
    fpin = open(inFile, "r");
    line = fpin.readline();
    r = re.compile("[<>]");
    while line:
        line = line.rstrip("\n")
        if line:
            line=line.lstrip("<").rstrip(">");
            strs=r.split(line);
            if strs[0] == "hmm_name":
                hmmname=strs[1];
                print >>fpout,"# Scores for HMM: '%s'"%hmmname;
            elif strs[0] == "pure_seq_name_a":# a new record, ID
                cntSeq+=1;
                seqID=strs[1];
                print >> fpout;
                print >> fpout, "Seq ID: %s" % seqID;
            elif strs[0] == "seqlength":
                length=int(strs[1]);
                print >> fpout,"Seq length: %d"%length;
                fpout.write("Labeling: ");
                for i in range(length):
                    line=fpin.readline().rstrip("\n").lstrip("<").rstrip(">");
                    fpout.write("%c"%r.split(line)[1]);
                fpout.write("\n");
            elif strs[0] == "posteriorprobabilities":
                line=fpin.readline().rstrip("\n").lstrip("<").rstrip(">");
                strs = r.split(line);
                if strs[0] == "labels":
                    line=fpin.readline().rstrip("\n").lstrip("<").rstrip(">");
                    strs = r.split(line);
                    while strs[0] != "/labels":
                        state = strs[1];
                        fpout.write("%-6s"%state);
                        line=fpin.readline().rstrip("\n").lstrip("<").rstrip(">");
                        strs = r.split(line);
                    fpout.write("\n");
            elif strs[0] == "post_prob_label_matrix":
                line=fpin.readline().rstrip("\n").lstrip("<").rstrip(">");
                strs = r.split(line);
                while strs[0] != "/post_prob_label_matrix":
                    if strs[0]=="post-prob":
                        nstr=len(strs)/4;
                        for j in range(nstr):
                            prob_score=strs[1+j*4].strip();
                            fpout.write("%-6s"%prob_score);
                        fpout.write("\n");
                    line=fpin.readline().rstrip("\n").lstrip("<").rstrip(">");
                    strs = r.split(line);
        line = fpin.readline();
    fpin.close();
    return cntSeq;
#}}}
def Modhmmxml2txt2(inFile,fpout):#{{{
#a little faster version
    cntSeq=0;
    fpin = open(inFile, "r");
    line = fpin.readline();
#    r = re.compile("[<>]");
    r=re.compile("[(<)(>)]");
    while line:
        if line[0] != "\n": #if not empty line
            strs=r.split(line);
            if strs[1] == "hmm_name":
                hmmname=strs[2];
                print >>fpout,"# Scores for HMM: '%s'"%hmmname;
            elif strs[1] == "pure_seq_name_a":# a new record, ID
                cntSeq+=1;
                seqID=strs[2];
                print >> fpout;
                print >> fpout, "Seq ID: %s" % seqID;
            elif strs[1] == "seqlength":
                length=int(strs[2]);
                print >> fpout,"Seq length: %d"%length;
                fpout.write("Labeling: ");
                for i in range(length):
                    line=fpin.readline();
                    fpout.write("%c"%r.split(line)[2]);
                fpout.write("\n");
            elif strs[1] == "posteriorprobabilities":
                line=fpin.readline();
                strs = r.split(line);
                if strs[1] == "labels":
                    line=fpin.readline();
                    strs = r.split(line);
                    while strs[1] != "/labels":
                        if strs[1] == "label":
                            state = strs[2];
                            fpout.write("%-6s"%state);
                        line=fpin.readline();
                        strs = r.split(line);
                    fpout.write("\n");
            elif strs[1] == "post_prob_label_matrix":
                line=fpin.readline();
                strs = r.split(line);
                while strs[1] != "/post_prob_label_matrix":
                    if strs[1]=="post-prob":
                        nstr=len(strs)/4;
                        for j in range(nstr):
                            prob_score=strs[2+j*4].strip();
                            fpout.write("%-6s"%prob_score);
                        fpout.write("\n");
                    line=fpin.readline();
                    strs = r.split(line);
        line = fpin.readline();
    fpin.close();
    return cntSeq;
#}}}
def Modhmmxml2txt3(inFile,fpout):#{{{
#much faster version
    cntSeq=0;
    fpin = open(inFile, "r");
    line = fpin.readline();
#    r = re.compile("[<>]");
    r=re.compile("[(<)(>)]");
    while line:
        if line[0] != "\n": #if not empty line
            strs=r.split(line);
            if strs[1] == "hmm_name":
                hmmname=strs[2];
                print >>fpout,"# Scores for HMM: '%s'"%hmmname;
            elif strs[1] == "pure_seq_name_a":# a new record, ID
                cntSeq+=1;
                seqID=strs[2];
                print >> fpout;
                print >> fpout, "Seq ID: %s" % seqID;
            elif strs[1] == "seqlength":
                length=int(strs[2]);
                print >> fpout,"Seq length: %d"%length;
                fpout.write("Labeling: ");
                for i in range(length):
                    line=fpin.readline();
                    fpout.write("%c"%line[7]);
                fpout.write("\n");
            elif strs[1] == "posteriorprobabilities":
                line=fpin.readline();
                strs = r.split(line);
                if strs[1] == "labels":
                    line=fpin.readline();
                    strs = r.split(line);
                    while strs[1] != "/labels":
                        if strs[1] == "label":
                            state = strs[2];
                            fpout.write("%-6s"%state);
                        line=fpin.readline();
                        strs = r.split(line);
                    fpout.write("\n");
            elif strs[1] == "post_prob_label_matrix":
                for i in range(length):
                    line=fpin.readline();
                    line=fpin.readline();
                    nstr=6;
                    for j in range(nstr):
                        prob_score=line[11+j*28:16+j*28];
                        fpout.write("%-6s"%prob_score);
                    fpout.write("\n");
        line = fpin.readline();
    fpin.close();
    return cntSeq;
#}}}
if __name__ == '__main__' :
    # Check argv
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp();
        sys.exit(1);

    outFile="";
    inFile="";
    method=3;

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
            elif sys.argv[i] == "-m" or sys.argv[i] == "--method" or sys.argv[i] == "-method":
                method=int(sys.argv[i+1]);
                if method < 1 or method > 3:
                    print >> sys.stderr,"Error! method should be 1, 2 or 3";
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
        if method == 1:
            numSeq = Modhmmxml2txt(inFile, fpout);
        elif method == 2:
            numSeq = Modhmmxml2txt2(inFile, fpout);
        else: 
            numSeq = Modhmmxml2txt3(inFile, fpout);

        if fpout != sys.stdout:
            fpout.close();

    except :
        print >>sys.stderr, "except for the input file: %s" % inFile;
        raise ;

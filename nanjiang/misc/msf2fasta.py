#!/usr/bin/env python
# convert msf to fasta

import sys, os

usage="""
Usage:  msf2fasta.py msf_files
  convert the msf format to fasta format
Options:
  -o <file> : output to file
  -q              : quiet mode
  -h|--help       : print this help message and exit
Created 2010-10-19, updated 2010-10-19, Nanjiang
"""

def PrintHelp():
    print usage

def MSF2FASTA(msfFile, fpout):
    fpin = open(msfFile, "r");
    cntAnnoLine=0;
    cntSeq=0;
    seqList=[];
    isAlignRegion = False;
    
    line=fpin.readline();
    while line:
        line=line.strip();
        if line:
            if not isAlignRegion: 
                if line.find("Name:") == 0:
                    seqList.append({});
                    line=line.lstrip("Name:").strip();
                    seqList[cntAnnoLine]['annoLine'] = line;
                    seqList[cntAnnoLine]['id']=line.split()[0];
                    seqList[cntAnnoLine]['seq']="";
                    cntAnnoLine +=1;
                elif line.find("//") == 0:
                    isAlignRegion=True;
                    cntSeq=0;

            elif isAlignRegion:
                strs=line.split();
                if not strs[len(strs)-1].isdigit():
                    if seqList[cntSeq]['id'] == strs[0]:
                        for i in range(1,len(strs)):
                            seqList[cntSeq]['seq'] += (strs[i].replace('.',''));
                    else:
                        print >> sys.stderr, "%s: msf format error"%msfFile;

                    cntSeq +=1;
                    if cntSeq == cntAnnoLine:
                        cntSeq = 0;
                
        line=fpin.readline();
    fpin.close();
# write the sequence
    for i in range(cntAnnoLine):
        print >> fpout, ">%s"%seqList[i]['id'];
        print >> fpout, "%s"%seqList[i]['seq'];

numArgv=len(sys.argv)
if numArgv < 2:
    PrintHelp()
    sys.exit()

isQuiet=False
outFile=""
fileList=[]
i = 1
isNonOptionArg=False
while i < numArgv:
    if isNonOptionArg == True:
        idList.append(sys.argv[i])
        isNonOptionArg=False
        i = i + 1
    elif sys.argv[i] == "--":
        isNonOptionArg=True
        i = i + 1
    elif sys.argv[i][0] == "-":
        if sys.argv[i] ==  "-h" or  sys.argv[i] == "--help":
            PrintHelp()
            sys.exit()
        elif sys.argv[i] == "-o" or sys.argv[i] == "--outfile":
            outFile=sys.argv[i+1]
            i = i + 2
        elif sys.argv[i] == "-q":
            isQuiet=True
            i = i + 1
        else:
            print "Error! Wrong argument:", sys.argv[i]
            sys.exit()
    else:
        fileList.append(sys.argv[i])
        i = i + 1

fpout = sys.stdout;
if outFile != "":
    fpout = open(outFile,"w");

try:
    for msfFile in fileList:
        MSF2FASTA(msfFile, fpout);
    if fpout != sys.stdout: 
        fpout.close();
except:
    print >> sys.stderr,"except";
    raise;

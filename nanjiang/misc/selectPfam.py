#!/usr/bin/env python
# select families in Pfam 
import sys,re,os
import math
import tempfile

BLOCK_SIZE = 100000
DEBUG=False

# ChangeLog 2012-11-07
# for output format fastaseq, '-' should also be considered as gap

usage="""
Usage:    selectPfam.py [Options] [-i] pfam-datafile

Options:
  -i         <file>    : input file
  -o         <file>    : outputfile
  -outpath   <dir>     : set the output path, when files are splitted
  -split               : when this is enabled, files will be splitted, the rootname is AC
  -keyword   <str>     : supply the key word to be searched in '#=GF DE'
  -matchopt  0|1       : keyword match option, 0 for AND and 1 for OR, default=0
  -outformat <str>     : set the output format, default=msf
                       : format can be, fasta, fastaseq, msf, gcg, swiss. fastaseq means with gaps removed
  -h|--help            : print this help message and exit
Created 2010-11-10, updated 2012-11-07, Nanjiang

Examples:
    selectPfam.py Pfam-A.seed -keyword membrane -keyword transmembrane -split -outpath Pfam_mem
"""

def PrintHelp():
    print usage

def ReadPfamRecord(pfamRecord):#{{{
#return aln
    aln={}
    aln['alnseq_id']=[]
    aln['alnseq'] = []

    lines=pfamRecord.split("\n")
    for line in lines:
        if not line:
            continue
        if line[0] == "#":
            if line.find("#=GF ID") == 0:
                aln['id']=line.split()[2]
            elif line.find("#=GF AC") == 0:
                aln['ac']=line.split()[2].split('.')[0]
            elif line.find("#=GF DE") == 0:
                aln['def']=line[8:].strip()
            elif line.find("#=GF SQ") == 0:
                aln['numseq']=int(line.split()[2])
        elif line != "//": 
            strs=line.split()
            if len(strs) == 2:
                aln['alnseq_id'].append(strs[0])
                aln['alnseq'].append(strs[1])
            else:
                print >> sys.stderr, "wrong alignment line:%s"%line
    if aln['numseq'] != len(aln['alnseq']):
        print >> sys.stderr, "GF SQ=%d != numseq in the record %d for id %s"%(aln['numseq'], len(aln['alnseq'], aln['id']))
    return aln
    #}}}
def IsMatchedKeyword(defLine, keyWordList, opt):#{{{
    defLine_upper=defLine.upper()
    isMatch=True
    if opt=='AND':
        isMatch = True
        for key in keyWordList:
            if defLine_upper.find(key.upper()) < 0:
                isMatch=False
                break
    elif opt=='OR':
        isMatch=False
        for key in keyWordList:
            if defLine_upper.find(key.upper()) >= 0:
                isMatch=True
                break
    return isMatch
#}}}
def GetDefinitionPfam(pfamRecord):#{{{
    begpos=pfamRecord.find("#=GF DE")
    endpos=pfamRecord.find("\n", begpos+1)
    return pfamRecord[begpos:endpos]
#}}}

def WriteAlnFasta(nameList,seqList,fpout):#{{{
    for i in range(len(nameList)):
        fpout.write(">%s\n"%nameList[i])
        fpout.write("%s\n"%seqList[i])
#}}}
def WriteAlnFastaSeq(nameList,seqList,fpout):#{{{
    for i in range(len(nameList)):
        fpout.write(">%s\n"%nameList[i])
        fpout.write("%s\n"%seqList[i].replace('-','').replace('.',''))
#}}}
def WriteAlnMSF(nameList, seqList, seqType, fpout):#{{{
    fpout.write("PileUp\n")
    fpout.write("\n\n\n")
    alignLength = len(seqList[0])
    fpout.write( '   MSF: %i  Type: %s  Check: 1 ..\n\n' % (alignLength, seqType.upper()) )

    for i in range(len(seqList)):
        seqList[i].replace('-','.')

# write header block
    for name in nameList:
        fpout.write( ' Name: %-32s Len:%6i  Check: %4i  Weight:  1.00\n' % (name, alignLength, 1) )
    fpout.write( '\n//\n\n\n\n' )

    lengthRow=50;   # each row write at most 50 residues
    numBlock=lengthRow/10
    numRow=int(math.ceil(alignLength/float(lengthRow)))

    for row in range(numRow):
        for i in range(len(seqList)):
            cntBlock=0
            j=lengthRow*row
            seq=seqList[i]
            fpout.write('%-32s '%nameList[i][:32])
            while j < alignLength:
                fpout.write("%s"% seq[j:j+10])
                j +=10
                cntBlock+=1
                if j >= alignLength:
                    fpout.write("\n")
                    break
                else:
                    if cntBlock < 5:
                        fpout.write(" ")
                    else: 
                        fpout.write("\n")
                        break
        fpout.write("\n\n")
#}}}

def ReadNextPfamRecord(fpin,oldbuff):#{{{
    pfamRecord=""
    newbuff=""
    endpos=oldbuff.find("\n//")
    while endpos < 0:
        tmpbuff=fpin.read(BLOCK_SIZE)
        if not tmpbuff:
            break
        oldbuff+=tmpbuff
        endpos=oldbuff.find("\n//")
    pfamRecord=oldbuff[:endpos+4]
    newbuff=oldbuff[endpos+4:]
    return (pfamRecord, newbuff)
#}}}
def WritePfamRecord(aln,pfamRecord, outFormat, fpout):#{{{
    if outFormat == 'swiss':
        fpout.write("%s"%pfamRecord)
    elif outFormat == 'msf' or outFormat=='gcg':
        WriteAlnMSF(aln['alnseq_id'], aln['alnseq'],'P', fpout)
    elif outFormat=='fasta':
        WriteAlnFasta(aln['alnseq_id'], aln['alnseq'], fpout)
    elif outFormat=='fastaseq':
        WriteAlnFastaSeq(aln['alnseq_id'], aln['alnseq'], fpout)
#}}}
def SelectPfam(inFile, keyWordList, outFile):#{{{
    try:
        fpout=sys.stdout
        if outFile != "":
            fpout = open(outFile,"w")

        fpin = open(inFile,"r")
        buff = ""
        (pfamRecord, buff)=ReadNextPfamRecord(fpin, buff)
        while pfamRecord:
            #fpout.write("%s"%pfamRecord)
            defLine=GetDefinitionPfam(pfamRecord)
            if DEBUG:
                print pfamRecord[0:200],"\n\ndefLine:", defLine, "\n\n"
            if IsMatchedKeyword(defLine, keyWordList, opt):
                #fpout.write("%s"%pfamRecord)
                aln=ReadPfamRecord(pfamRecord)
                WritePfamRecord(aln, pfamRecord, outFormat, fpout)

            (pfamRecord, buff)=ReadNextPfamRecord(fpin, buff)

        if fpout != sys.stdout:
            fpout.close()
        fpin.close()

        return 0
    except:
        print >>sys.stderr, "except for the input file: %s" % inFile
        raise
#}}}
def SelectPfam_split(inFile, keyWordList, outpath):#{{{
    try:
        fpin = open(inFile,"r")
        buff = ""
        (pfamRecord, buff)=ReadNextPfamRecord(fpin, buff)
        while pfamRecord:
            #fpout.write("%s"%pfamRecord)
            defLine=GetDefinitionPfam(pfamRecord)
            if DEBUG:
                print pfamRecord[0:200],"\n\ndefLine:", defLine, "\n\n"
            if IsMatchedKeyword(defLine, keyWordList, opt):
                #fpout.write("%s"%pfamRecord)
                aln=ReadPfamRecord(pfamRecord)
                if outFormat=='swiss':
                    ext="swiss"
                elif outFormat=='msf' or outFormat=='gcg':
                    ext="msf"
                elif outFormat=="fastaseq":
                    ext="fa"
                elif outFormat =="fasta":
                    ext="aln.fa"
                else:
                    ext=""
                outFile=outpath+os.sep+aln['ac']+'.'+ext
                fpout = open(outFile,"w")
                WritePfamRecord(aln, pfamRecord, outFormat, fpout)
                fpout.close()

            (pfamRecord, buff)=ReadNextPfamRecord(fpin, buff)

        fpin.close()

        return 0
    except:
        print >>sys.stderr, "except for the input file: %s" % inFile
        raise
#}}}
if __name__ == '__main__' :
    # Check argv
    argv = sys.argv
    numArgv=len(argv)
    if numArgv < 2:
        PrintHelp()
        sys.exit(1)

    outFile=""
    outpath=""
    inFile=""
    isSplit=False
    isQuiet=False
    keyWordList=[]
    outFormat='msf'
    opt='AND'

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            isNonOptionArg=False
            i = i + 1
        elif argv[i] == "--":
            isNonOptionArg=True
            i = i + 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                sys.exit(0)
            elif argv[i] in ["-i" , "--infile"]:
                inFile=argv[i+1]
                i = i + 2
            elif argv[i] == "-o" or argv[i] == "--outfile":
                outFile=argv[i+1]
                i = i + 2
            elif argv[i] == "-outpath" or argv[i] == "--outpath":
                outpath=argv[i+1]
                i = i + 2
            elif argv[i] == "-keyword" or argv[i] == "--keyword":
                keyWordList.append(argv[i+1])
                i = i + 2
            elif argv[i] == "-split" or argv[i] == "--split":
                isSplit=True
                i = i + 1
            elif argv[i] == "-outformat" or argv[i] == "--outformat":
                outFormat=argv[i+1].lower()
                i = i + 2
            elif argv[i] == "-matchopt" or argv[i] == "--matchopt":
                tint=int(argv[i+1])
                if tint == 0:
                    opt='AND'
                elif tint == 1:
                    opt = 'OR'
                else:
                    print >> sys.stderr,"-matchopt can be only 0 or 1"
                    sys.exit(1)
                i = i + 2
            elif argv[i] == "-q" or argv[i] == "--q" or argv[i] == "-quiet" or argv[i]=="--quiet":
                isQuiet=True
                i = i + 1
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % argv[i])
                sys.exit(1)
        else:
            inFile=argv[i]
            i+=1

    if inFile == "":
        print >> sys.stderr,"Error! Input file not set."

    try :
        if not isSplit:
            SelectPfam(inFile, keyWordList, outFile)
        else:
            if outpath == "":
                print >> sys.stderr, "when isSplit is enabled, outpath must be set"
                sys.exit(1)
            os.system("mkdir -p %s"%outpath)
            SelectPfam_split(inFile, keyWordList, outpath)

    except :
        print >>sys.stderr, "except for the input file: %s" % inFile
        raise 

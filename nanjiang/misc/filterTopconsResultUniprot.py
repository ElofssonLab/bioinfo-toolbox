#!/usr/bin/env python
# Description:
# Filter the result file of topcons prediction
# filter uniref id by uniprot id and repalce the annotation line
import os
import sys
import myfunc
usage = """
usage:  filterTopconsResultUniprot.py  topcons_resultfile -idwithanno FILE -o OUTFILE
Description:
# Filter the result file by topcons
# filter uniref id by uniprot id and repalce the annotation line

Options:
  -q              Quiet mode
  -h, --help      Print this help message and exit

Created 2012-06-08, updated 2012-06-08, Nanjiang Shu 
"""
BLOCK_SIZE=100000;

def PrintHelp():
    print usage
def WriteRecord(rd, fpout): #{{{
    if not fpout:
        return 1
    else:
        for line in rd['otherLineListBefore']:
            fpout.write("%s\n" %(line))
        fpout.write("%s %s\n" %("Sequence name:", rd['anno']))
        for line in rd['otherLineListAfter']:
            fpout.write("%s\n" %(line))
        return 0
#}}}
def ReadSeqIDWithAnnoFile(infile):#{{{
    try:
        seqid2AnnoDict = {}
        fpin = open(infile, "r")
        line = fpin.readline()
        while line:
            line = line.strip()
            if line and line[0] != "#":
                strs = line.split("\t")
                if len(strs) == 2:
                    seqid2AnnoDict[strs[0]] = strs[1]
            line = fpin.readline()
        fpin.close()
        return seqid2AnnoDict
    except IOError:
        print >> sys.stderr, "Failed to read file", infile
        return {}
#}}}
def ExtractFromTopconsResult(recordContent, seqid2AnnoDict):#{{{
    record = {}
    lines = recordContent.split("\n")
    numLine = len(lines)
    if numLine < 1:
        return {}
    i = 0
    anno = ""
    otherLineListBefore = []
    otherLineListAfter = []
    isSeqNameFound = False
    for line in lines:
        if line.find("Sequence name") == 0:
            anno = line[15:]
            isSeqNameFound = True
        else:
            if not isSeqNameFound:
                otherLineListBefore.append(line)
            else:
                otherLineListAfter.append(line)

    seqid = myfunc.GetSeqIDFromAnnotation(anno)
    if seqid.find("UniRef") == 0:
        seqid = seqid.split("_")[1]
    if seqid in seqid2AnnoDict:
        rd = {}
        rd['anno'] = seqid2AnnoDict[seqid]
        rd['otherLineListBefore'] = otherLineListBefore
        rd['otherLineListAfter'] = otherLineListAfter
        return rd
    else:
        return {}
#}}}

def Read_topcons_result_from_buffer(buff, recordList, seqid2AnnoDict, isEOFreached):#{{{
    if not buff:
        return ""
    unprocessedBuffer="";
    beg=0;
    end=0;
    while 1:
        beg=buff.find("TOPCONS result file",beg);
        if beg >= 0:
            end=buff.find("\nTOPCONS result file",beg+1);
            if end >=0:
                recordContent = buff[beg:end];
                record = ExtractFromTopconsResult(recordContent, seqid2AnnoDict)
                if record != {}:
                    recordList.append(record)
                beg = end;
            else:
                unprocessedBuffer = buff[beg:];
                break;
        else:
            unprocessedBuffer = buff[end:];
            break;
    if isEOFreached and unprocessedBuffer:
        recordContent = unprocessedBuffer
        record = ExtractFromTopconsResult(recordContent, seqid2AnnoDict)
        if record != {}:
            recordList.append(record)
        unprocessedBuffer = ""
    return unprocessedBuffer;
    #}}}
def FilterTopconsResultUniprot(infile, seqid2AnnoDict, fpout):#{{{
    try:
        fpin = open(infile, "r")
        unprocessedBuffer="";
        isEOFreached = False;
        processedTopoIDSet = set([]);
        isBeginOutput = True
        while 1:
            buff = fpin.read(BLOCK_SIZE);
            if len(buff) < BLOCK_SIZE:
                isEOFreached=True;
            buff = unprocessedBuffer + buff;
            recordList = [];
            unprocessedBuffer = Read_topcons_result_from_buffer(
                    buff, recordList, seqid2AnnoDict, isEOFreached);
            if len(recordList) > 0: 
                if isBeginOutput:
                    fpout.write("##############################################################################\n")
                    isBeginOutput = False
                for record in recordList:
                    WriteRecord(record, fpout)
            if isEOFreached == True:
                break;
        fpin.close();
    except IOError:
        print >> sys.stderr, "Failed to open file %s for read"%(infile);
        return 1
#}}}
def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile = ""
    infile = ""
    idwithannofile = ""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            infile = argv[i]
            isNonOptionArg = False
            i += 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-o", "--o"]:
                outfile = argv[i+1]
                i += 2
            elif argv[i] in ["-idwithanno", "--idwithanno"]:
                idwithannofile = argv[i+1]
                i += 2
            elif argv[i] in ["-q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            infile = argv[i]
            i += 1

    if infile == "":
        print >> sys.stderr, "infile not set.exit"
        return 1
    if idwithannofile == "":
        print >> sys.stderr, "idwithannofile not set.exit"
        return 1

    seqid2AnnoDict = ReadSeqIDWithAnnoFile(idwithannofile)
    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    FilterTopconsResultUniprot(infile, seqid2AnnoDict, fpout)
    myfunc.myclose(fpout)
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

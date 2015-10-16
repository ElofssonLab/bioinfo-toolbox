#!/usr/bin/env python
# Description:
# convert phobius (format long) result file to fasta file
import os
import sys
import myfunc
usage = """
usage:  phobius_to_fasta.py phobius_long_result_file -outpath DIR
Description:

Options:
  -outpath DIR    Set ouput path
  -q              Quiet mode
  -h, --help      Print this help message and exit

Created 2012-08-27, updated 2012-08-27, Nanjiang Shu  
"""
BLOCK_SIZE=100000;

def PrintHelp():
    print usage

def PosTM2Topo(posTM, seqLength, NtermState):
    if NtermState == "":
        return ""
    topList = []
    statelist = ["i", "o"]
    idx = 0
    if NtermState == 'i':
        idx = 0
    else:
        idx = 1
        
    state = statelist[idx]
    if len(posTM) < 1:
        topList += [state]*seqLength
    else:
        for j in xrange(len(posTM)):
            state = statelist[idx%2]
            if j == 0:
                seglen = posTM[j][0]
            else:
                seglen = posTM[j][0] - posTM[j-1][1]
            topList += [state]*seglen
            topList += ['M'] * (posTM[j][1]-posTM[j][0])
            idx += 1
        #print posTM, seqLength
        if posTM[len(posTM)-1][1] < seqLength:
            state = statelist[idx%2]
            topList += [state] * (seqLength - posTM[len(posTM)-1][1])
    top = "".join(topList)
    return top

def ExtractFromPhobiusResult(recordContent):#{{{
    record = {}
    lines = recordContent.split("\n")
    numLine = len(lines)
    i = 0
    anno = ""
    predtopo = ""
    seqLength = 0
    posTM = []
    NtermState = ""
    while i < numLine:
#         print i
        if lines[i].find("ID") == 0:
            anno = lines[i][5:]
            i += 1
        elif lines[i].find('FT') == 0:
            strs = lines[i][5:].split()
            if len(strs) >= 3:
                if lines[i].find("NON CYTOPLASMIC") != -1:
                    if NtermState == "":
                        NtermState = "o"
                elif lines[i].find("CYTOPLASMIC") != -1:
                    if NtermState == "":
                        NtermState = "i"
                try:
                    b = int(strs[1])
                except (ValueError, TypeError):
                    print >> sys.stderr, "Error record. <%s>"%lines[i]
                    return {}
                try:
                    e = int(strs[2])
                except (ValueError, TypeError):
                    print >> sys.stderr, "Error record. <%s>"%lines[i]
                    return {}

                seqLength = e

                if strs[0] == "TRANSMEM":
                    posTM.append((b-1, e))
            i += 1
        else:
            i += 1
    if len(posTM) > 0:
        record['anno'] = anno
        record['seqlength'] = seqLength
        record['predtopo'] = PosTM2Topo(posTM, seqLength, NtermState)
        return record
    else:
        return {}
#}}}

def Read_phobius_result_from_buffer(buff, recordList, isEOFreached):#{{{
    if not buff:
        return ""
    unprocessedBuffer="";
    beg=0;
    end=0;
    while 1:
        beg=buff.find("ID",beg);
        if beg >= 0:
            end=buff.find("\n//",beg+1);
            if end >=0:
                recordContent = buff[beg:end];
                record = ExtractFromPhobiusResult(recordContent)
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
        record = ExtractFromPhobiusResult(recordContent)
        if record != {}:
            recordList.append(record)
        unprocessedBuffer = ""
    return unprocessedBuffer;
    #}}}

def Phobius2Fasta(infile, outpath):#{{{
    try:
        rootname=os.path.basename(os.path.splitext(infile)[0])
        outfile = outpath + os.sep + rootname + "_PHOBIUS.topo"
        fpout = open(outfile, "w")
        fpin = open(infile, "r")
        unprocessedBuffer="";
        isEOFreached = False;
        processedTopoIDSet = set([]);
        while 1:
            buff = fpin.read(BLOCK_SIZE);
            if len(buff) < BLOCK_SIZE:
                isEOFreached=True;
            buff = unprocessedBuffer + buff;
            recordList = [];
            unprocessedBuffer = Read_phobius_result_from_buffer(
                    buff, recordList, isEOFreached);
            if len(recordList) > 0: 
                for record in recordList:
                    seqid = myfunc.GetSeqIDFromAnnotation(record['anno'])

                    if record['predtopo'] != "":
                        fpout.write(">%s\n"%(
                            record['anno']))
                        fpout.write("%s\n"%record['predtopo'])
            if isEOFreached == True:
                break;
        fpin.close();
        fpout.close();
        print "Results have been output to"
        print "\t%s"%outfile
    except IOError:
        print >> sys.stderr, "Failed to open file %s for read"%(infile);
        raise
#}}}

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = "./"
    phobius_result_file = ""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            phobius_result_file = argv[i]
            isNonOptionArg = False
            i += 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-outpath", "--outpath"]:
                outpath = argv[i+1]
                i += 2
            elif argv[i] in ["-q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            phobius_result_file = argv[i]
            i += 1

    Phobius2Fasta(phobius_result_file, outpath)
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}

if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

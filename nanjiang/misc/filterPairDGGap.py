#!/usr/bin/env python

# given the data produced by compareMSATopo.py -mode 0
# filter the compared topology by parameters such as DG value, gap fraction
# and topology prediction reliability score.
# output the alignment of the filtered pairs

import os,sys;
import libtopologycmp as lcmp;
import myfunc;

BLOCK_SIZE=100000;
usage="""
usage:   filterPairDGGap.py paircmp-file -aln alignment-file [-o OUTFILE]

Description: Filter the pairwisely compared topologies by parameters such as DG
             values, Gap fraction and topology prediction reliability score.
             output filtered alignments
Options:
  -o OUTFILE   Output the result to outfile
  -q           Quiet mode
  -write-paircmp FILE
               output paircmp data to file
  -h, --help   Print this help message and exit

Selection control options:
  -gap  FLOAT  Select only the TM with gap fraction >= threshold, (default: 0.5)
  -dg   FLOAT  Select only the TM with DG values <= threshold, (default: 1.0)
  -ps   FLOAT  Select only the TM with topology prediction reliability >=
               threshold, (default: 0.5)
  -min-seqidt  FLOAT, (default: 0)
  -max-seqidt  FLOAT, (default: 100)
               Select only pairs with global sequence identity within [minSeqIDT, maxSeqIDT]

Created 2012-02-03, updated 2012-02-03, Nanjiang Shu 
"""

def PrintHelp():
    print usage;

def IsAnaHasInternalVariation(ana):#{{{
    if ana == {}:
        return False;
    if 'internal' in ana and len(ana['internal']) > 0:
        return True;
    else:
        return False;
        #}}}
def IsHasInternalVariation(record):#{{{
    if record == {}:
        return False;
    if record['cmpclass'] != 'DIFF':
        return False;
    if (IsAnaHasInternalVariation(record['ana1']) or
        IsAnaHasInternalVariation(record['ana2'])):
        return True;
    else:
        return False;
#}}}

def FilterPairCmpResult(recordList, parameters):#{{{
    """
    Filter paircmp result by parameters
    return newList;
    """
    newList = [];
    for record in recordList:
        if record == {}:
            continue;
        
        if record['cmpclass'] != 'DIFF':
            continue;
        if (record['seqidt'] < parameters['minSeqIDT'] or record['seqidt'] >
                parameters['maxSeqIDT']):
            continue;

        newRecord = {};
        newana1  = lcmp.SelectAnaDIFF(record['ana1'], parameters);
        newana2  = lcmp.SelectAnaDIFF(record['ana2'], parameters);
        if newana1 != {} or newana2 != {}:
            lcmp.CopyGeneralInfo_pairwise(newRecord, record);
            newRecord['ana1'] = newana1;
            newRecord['ana2'] = newana2;
            newList.append(newRecord);
    return newList;
#}}}

def GetPairTopoAln(pairalnTopoFile):#{{{
    (idList, annoList, seqList) = myfunc.ReadFasta(pairalnTopoFile);
    numPair = len(idList)/2;
    pairTopoAlnDict = {};
    for i in xrange(numPair):
        pair = {};
        pair['id1'] = idList[i*2];
        pair['id2'] = idList[i*2+1];
        pair['anno1'] = annoList[i*2];
        pair['anno2'] = annoList[i*2+1];
        pair['seq1'] = seqList[i*2];
        pair['seq2'] = seqList[i*2+1];
        key = "%s-%s"%(idList[i*2], idList[i*2+1]);
        pairTopoAlnDict[key] = pair;
    return pairTopoAlnDict;
#}}}

def main():#{{{
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1;

    parameters={};
    parameters['minGapFraction'] = 0.5;
    parameters['maxGapFraction'] = 1.0;
    parameters['minDGvalue'] = -999999.0;
    parameters['maxDGvalue'] = 1.0;
    parameters['minSeqIDT'] = 0.0;
    parameters['maxSeqIDT'] = 100.0;

    infile="";
    outfile="";
    outPaircmpfile = "";
    pairalnTopoFile = "";
    isQuiet=False;

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
            elif (sys.argv[i] == '-o' or sys.argv[i] == '--o' or sys.argv[i]
                    == "-outfile" or sys.argv[i] == "--outfile"):
                outfile=sys.argv[i+1];
                i += 2;
            elif sys.argv[i] == "-gap" or sys.argv[i] == "--gap":
                parameters['minGapFraction'] = float(sys.argv[i+1]);
                i += 2;
            elif sys.argv[i] == "-dg" or sys.argv[i] == "--dg":
                parameters['maxDGvalue'] = float(sys.argv[i+1]);
                i += 2;
            elif sys.argv[i] in ["-min-seqidt", "--min-seqidt"]:
                parameters['minSeqIDT'] = float(sys.argv[i+1]);
                i += 2;
            elif sys.argv[i] in ["-max-seqidt", "--max-seqidt"]:
                parameters['maxSeqIDT'] = float(sys.argv[i+1]);
                i += 2;
            elif sys.argv[i] in ["-write-paircmp", "--write-paircmp"]:
                outPaircmpfile = sys.argv[i+1];
                i += 2;
            elif sys.argv[i] in ["-aln", "--aln"]:
                pairalnTopoFile = sys.argv[i+1];
                i += 2;
            elif sys.argv[i] == "-q":
                isQuiet=True;
                i += 1;
            else:
                print >> sys.stderr, "Error! Wrong argument:", sys.argv[i];
                return -1;
        else:
            infile=sys.argv[i];
            i += 1
    if infile == "":
        print >> sys.stderr, "infile not set. Exit.";
        return -1;
    elif not os.path.exists(infile):
        print >> sys.stderr, "infile %s does not exists. Exit."%infile;

    if pairalnTopoFile == "":
        print >> sys.stderr, "pairalnTopoFile not set. Exit.";
        return -1;

    pairTopoAlnDict = GetPairTopoAln(pairalnTopoFile);
    # pairTopoAlnDict[id1-id2]['id1] ['id2'] ['anno1'] ['anno2] ['seq1']
    # ['seq2']

    rootname=os.path.basename(os.path.splitext(infile)[0]);
    fpout = sys.stdout;
    fppaircmp = None;

    if outPaircmpfile != "" :
        fppaircmp = open(outPaircmpfile, "w")


    if outfile != "":
        try:
            fpout = open(outfile,"w");
        except IOError:
            print >>sys.stderr, "Failed to write to file %s."%outfile;
            print >> sys.stderr, "Reset output to sys.stdout.";
            fpout = sys.stdout;
            pass;

    fpin = open (infile, "rb");
    if not fpin:
        print >> sys.stderr, "Failed to open input file %s"%(infile);
        return -1;

    unprocessedBuffer="";
    cntTotalOutputRecord = 0;
    cntTotalReadInRecord = 0;
    isEOFreached = False;
    while 1:
        buff = fpin.read(BLOCK_SIZE);
        if buff == "":
            isEOFreached = True;
        buff = unprocessedBuffer + buff;
        pairCmpRecordList=[];
        unprocessedBuffer = lcmp.ReadPairCmpResultFromBuffer(buff,pairCmpRecordList);
        if len(pairCmpRecordList) > 0: 
            #WritePairCmpRecord(pairCmpRecordList,fpout);
           filteredList =  FilterPairCmpResult(pairCmpRecordList, parameters);
           for record in filteredList:
               if IsHasInternalVariation(record):
                   if fppaircmp != None:
                       li = [];
                       li.append(record);
                       (status, cntTotalOutputRecord ) = lcmp.WritePairCmpRecord(li, cntTotalOutputRecord, fppaircmp);
                   key = "%s-%s"%(record['id1'], record['id2']);
                   pair = pairTopoAlnDict[key];
                   fpout.write(">%s\n"%pair['anno1']);
                   fpout.write("%s\n"%pair['seq1']);
                   fpout.write(">%s\n"%pair['anno2']);
                   fpout.write("%s\n"%pair['seq2']);
           cntTotalReadInRecord += len(pairCmpRecordList);
        if isEOFreached == True:
            break;

    fpin.close();
    print "cntTotalReadInRecord =", cntTotalReadInRecord;
    print "cntTotalOutputRecord =", cntTotalOutputRecord;

    if fpout != None and fpout != sys.stdout:
        fpout.close();
    if fppaircmp != None:
        fppaircmp.close();

    return 0;

#}}}
if __name__ == '__main__' :
    sys.exit(main());

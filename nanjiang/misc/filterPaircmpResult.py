#!/usr/bin/env python

# given the data produced by compareMSATopo.py -mode 0
# filter the compared topology by parameters such as DG value, gap fraction
# and topology prediction reliability score.

import os,sys;
import libtopologycmp as lcmp;

BLOCK_SIZE=100000;
usage="""
usage:   filterPaircmpResult.py paircmp-file [-o OUTFILE]

Description: Filter the pairwisely compared topologies by parameters such as DG
             values, Gap fraction and topology prediction reliability score.
Options:
  -o OUTFILE   Output the result to outfile
  -q           Quiet mode
  -h, --help   Print this help message and exit

Selection control options:
  -gap  FLOAT  Select only the TM with gap fraction >= threshold, (default: 0.5)
  -dg   FLOAT  Select only the TM with DG values <= threshold, (default: 1.0)
  -ps   FLOAT  Select only the TM with topology prediction reliability >=
               threshold, (default: 0.5)

Created 2011-11-07, updated 2011-11-07, Nanjiang Shu 
"""

def PrintHelp():
    print usage;

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
            newList.append(record);
        else:
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

    infile="";
    outfile="";
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

    rootname=os.path.basename(os.path.splitext(infile)[0]);
    fpout = sys.stdout;
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
           (status, cntTotalOutputRecord) = lcmp.WritePairCmpRecord(filteredList, cntTotalOutputRecord, fpout);
           cntTotalReadInRecord += len(pairCmpRecordList);
        if isEOFreached == True:
            break;

    fpin.close();
    print "cntTotalReadInRecord =", cntTotalReadInRecord;
    print "cntTotalOutputRecord =", cntTotalOutputRecord;

    if fpout != None and fpout != sys.stdout:
        fpout.close();

    return 0;

#}}}
if __name__ == '__main__' :
    sys.exit(main());

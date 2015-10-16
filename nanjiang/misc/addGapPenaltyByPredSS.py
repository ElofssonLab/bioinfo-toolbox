#!/usr/bin/env python
# add the gap penalty information to the fasta sequence according to the
# predicted secondary structures
import sys,re,os,math;
import myfunc

DEBUG=False;

p_shift_C      = 0.5;
p_shift_HE     = 0.53;
p_threshold_C  = 0.5;
p_threshold_HE = 0.5;
weight_C       = 0.6;
weight_HE      = 0.4;

usage="""
Usage: addGapPenaltyByPredSS.py [Options] fasta-seq-files ...

Note: the output file will be named as $outpath/$rootname.gp.fa

Options:
  -l         <file>: set the list file for fasta files
  -outpath   <dir> : set the output path for fasta files with gap penalties added, 
                   : default=the same as the fasta-seq-file
  -sspath    <dir> : path for the predicted secondary structures
                   : it must be set, and the secondary structure files is stored as $sspath/$rootname/$rootname_i
  -fss       0|1   : format for the predicted secondary structures, default=0
                   : 0 for psipred output
  -method-g  <int> : method for calculating gap penalties, 0-5, default=5
  -q               : quiet mode
  -h|--help        : print this help message and exit
  
  Other parameters related to the model
  -p-shift-c <float> 
  -p-shift-he <float> 
  -weight-c <float>
  -weight-he <float>
  -p-threshold-c <float>
  -p-threshold-he <float>

Created 2010-10-24, updated 2010-11-08, Nanjiang
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
    predSS = {};
    if ssFormat == 0:
        #aaSeq="";
        ssSeq = ""
        prob_C=[];
        prob_H=[];
        prob_E=[];
        for line in lines:
            strs = line.split();
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
def ReadInPredSS_2(inFile, ssFormat):#{{{
#this procedure is actually slower
    try:
        filesize=os.path.getsize(inFile);
        fpin = open(inFile, "r");
        buff = fpin.read(filesize);
        fpin.close();
        predSS={};
        lines=buff.split("\n");
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
    except:
        print >> sys.stderr,"except to read file %s"%inFile;
        raise;
#}}}

def GetGapPenalty(ssFile, seqLength, ssFormat):#{{{
# return (gpoArray, gpeArray, tgpeArray)
# increase gapopens for secondary structure elements predicted with the
# probability > threshold
    predSS = ReadInPredSS(ssFile, ssFormat);
    if seqLength != len(predSS['ss']):
        print >> sys.stderr, "Error! aaSeqLength = %d, but the ssSeqLength = %d, for file %s" %(seqLength, len(predSS['ss']), ssFile);
        sys.exit(1);
        
    gpoArray=[1]*seqLength;
    gpeArray=[1]*seqLength;
    tgpeArray=[1]*seqLength;
    ssSeq=predSS['ss'];
    if DEBUG:
        print >> sys.stdout, "%s" % ssSeq;
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
                    gpo = 1.0+ (avgP - p_shift_HE)*weight_HE;
                    gpoArray[i:i+j] = [gpo]*j ; 
                    if DEBUG:
                        print >> sys.stdout, "gpo=%g" % gpo; 
            i+= j;
        i+= 1;
    return (gpoArray, gpeArray, tgpeArray);
#}}}
def GetGapPenalty1(ssFile, seqLength, ssFormat):#{{{
# return (gpoArray, gpeArray, tgpeArray)
# increase gapopens for secondary structure elements predicted with the
# probability > threshold. But the gpo is calculated based on the probability
# of individual residue position, no the average probability
    predSS = ReadInPredSS(ssFile, ssFormat);
    if seqLength != len(predSS['ss']):
        print >> sys.stderr, "Error! aaSeqLength = %d, but the ssSeqLength = %d, for file %s" %(seqLength, len(predSS['ss']), ssFile);
        sys.exit(1);
        
    gpoArray=[1]*seqLength;
    gpeArray=[1]*seqLength;
    tgpeArray=[1]*seqLength;
    ssSeq=predSS['ss'];
    if DEBUG:
        print >> sys.stdout, "%s" % ssSeq;
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
                if avgP > p_shift_HE:
                    for jj in range(j):
                        incPercent=weight_HE*(pArray[jj]-p_shift_HE);
                        gpoArray[i+jj] = gpoArray[i+jj]*(1.0+incPercent) ; 
                        gpeArray[i+jj] = gpeArray[i+jj]*(1.0+incPercent) ; 
                        tgpeArray[i+jj] = tgpeArray[i+jj]*(1.0+incPercent) ; 
                        if DEBUG:
                            print >> sys.stdout, "incPercent=%g" % incPercent*100; 
            i+= j;
        i+= 1;
    return (gpoArray, gpeArray, tgpeArray);
#}}}
def GetGapPenalty2(ssFile, seqLength, ssFormat):#{{{
# return (gpoArray, gpeArray, tgpeArray)
# decrease gapopens for random coils predicted with 
# probability > threshold. 
    predSS = ReadInPredSS(ssFile, ssFormat);
    if seqLength != len(predSS['ss']):
        print >> sys.stderr, "Error! aaSeqLength = %d, but the ssSeqLength = %d, for file %s" %(seqLength, len(predSS['ss']), ssFile);
        sys.exit(1);
        
    gpoArray=[1]*seqLength;
    gpeArray=[1]*seqLength;
    tgpeArray=[1]*seqLength;
    ssSeq=predSS['ss'];
    if DEBUG:
        print >> sys.stdout, "%s" % ssSeq;
    i = 0;
    while i < seqLength:
        if ssSeq[i] == 'C':
            state = ssSeq[i];
            j = 0;
            pArray = [];
            while i+j < seqLength and ssSeq[i+j] == state:
                pArray.append(predSS[state][i+j]);
                j += 1;
            if j >= 1:
                avgP =  sum(pArray, 0.0)/len(pArray); 
                if avgP > p_threshold_C:
                    for jj in range(j):
                        decPercent=weight_C*(pArray[jj]-p_shift_C);
                        gpoArray[i+jj] = gpoArray[i+jj]*(1.0-decPercent) ; 
                        gpeArray[i+jj] = gpeArray[i+jj]*(1.0-decPercent) ; 
                        tgpeArray[i+jj] = tgpeArray[i+jj]*(1.0-decPercent) ; 
                        if DEBUG:
                            print >> sys.stdout, "decPercent=%g" % decPercent*100; 
            i+= j;
        i+= 1;
    return (gpoArray, gpeArray, tgpeArray);
#}}}
def GetGapPenalty3(ssFile, seqLength, ssFormat):#{{{
# return (gpoArray, gpeArray, tgpeArray)
# increase gapopens for SSEs and 
# decrease gapopens for random coils predicted with 
# probability > threshold. 
# gpos are estimated based on probabilities on each residue position
#
# This procedure takes >50% of the time

    gpoArray=[1]*seqLength;
    gpeArray=[1]*seqLength;
    tgpeArray=[1]*seqLength;
    #return (gpoArray, gpeArray, tgpeArray);

    predSS = ReadInPredSS(ssFile, ssFormat);
    if seqLength != len(predSS['ss']):
        print >> sys.stderr, "Error! aaSeqLength = %d, but the ssSeqLength = %d, for file %s" %(seqLength, len(predSS['ss']), ssFile);
        sys.exit(1);
        
    # the following procedure takes about 20% of time

    ssSeq=predSS['ss'];
    if DEBUG:
        print >> sys.stdout, "%s" % ssSeq;
    #increase the gap penalties for secondary structure elements
    i = 0;
    while i < seqLength:
        if ssSeq[i] == 'H' or ssSeq[i] == 'E':
            state = ssSeq[i];
            j = 0;
            pArray = [];
            while  i+j < seqLength and ssSeq[i+j] == state:
                pArray.append(predSS[state][i+j]);
                j += 1;
            if j >= 2:
                avgP =  sum(pArray, 0.0)/len(pArray); 
                if avgP > p_threshold_HE:
                    for jj in range(j):
                        incPercent=weight_HE*(pArray[jj]-p_shift_HE);
                        gpoArray[i+jj] = gpoArray[i+jj]*(1.0+incPercent) ; 
                        gpeArray[i+jj] = gpeArray[i+jj]*(1.0+incPercent) ; 
                        tgpeArray[i+jj] = tgpeArray[i+jj]*(1.0+incPercent) ; 
                        if DEBUG:
                            print >> sys.stdout, "incPercent=%g" % incPercent*100; 
            i+= j;
        i+= 1;
    #decrease the gap penalties for random coils
    i = 0;
    while i < seqLength:
        if ssSeq[i] == 'C':
            state = ssSeq[i];
            j = 0;
            pArray = [];
            while i+j < seqLength and ssSeq[i+j] == state:
                pArray.append(predSS[state][i+j]);
                j += 1;
            if j >= 1:
                avgP =  sum(pArray, 0.0)/len(pArray); 
                if avgP > p_threshold_C:
                    for jj in range(j):
                        decPercent=weight_C*(pArray[jj]-p_shift_C);
                        gpoArray[i+jj] = gpoArray[i+jj]*(1.0-decPercent) ; 
                        gpeArray[i+jj] = gpeArray[i+jj]*(1.0-decPercent) ; 
                        tgpeArray[i+jj] = tgpeArray[i+jj]*(1.0-decPercent) ; 
                        if DEBUG:
                            print >> sys.stdout, "decPercent=%g" % decPercent*100; 
            i+= j;
        i+= 1;
    return (gpoArray, gpeArray, tgpeArray);
#}}}
def GetGapPenalty4(ssFile, seqLength, ssFormat):#{{{

# return (gpoArray, gpeArray, tgpeArray)
# increase gapopens for SSEs and 
# decrease gapopens for random coils predicted with 
# probability > threshold. 
# gpos are estimated based on average probabilities

    predSS = ReadInPredSS(ssFile, ssFormat);
    if seqLength != len(predSS['ss']):
        print >> sys.stderr, "Error! aaSeqLength = %d, but the ssSeqLength = %d, for file %s" %(seqLength, len(predSS['ss']), ssFile);
        sys.exit(1);
        
    gpoArray=[1]*seqLength;
    gpeArray=[1]*seqLength;
    tgpeArray=[1]*seqLength;
    ssSeq=predSS['ss'];
    if DEBUG:
        print >> sys.stdout, "%s" % ssSeq;
    #increase the gap penalties for secondary structure elements
    i = 0;
    while i < seqLength:
        if ssSeq[i] == 'H' or ssSeq[i] == 'E':
            state = ssSeq[i];
            j = 0;
            pArray = [];
            while  i+j < seqLength and ssSeq[i+j] == state:
                pArray.append(predSS[state][i+j]);
                j += 1;
            if j >= 2:
                avgP =  sum(pArray, 0.0)/len(pArray); 
                if avgP > p_threshold_HE:
                    incPercent=weight_HE*(avgP-p_shift_HE);
                    for jj in range(j):
                        gpoArray[i+jj] = gpoArray[i+jj]*(1.0+incPercent) ; 
                        gpeArray[i+jj] = gpeArray[i+jj]*(1.0+incPercent) ; 
#                         tgpeArray[i+jj] = tgpeArray[i+jj]*(1.0+incPercent) ; 
                        if DEBUG:
                            print >> sys.stdout, "incPercent=%g" % incPercent*100; 
            i+= j;
        i+= 1;
    #decrease the gap penalties for random coils
    i = 0;
    while i < seqLength:
        if ssSeq[i] == 'C':
            state = ssSeq[i];
            j = 0;
            pArray = [];
            while i+j < seqLength and ssSeq[i+j] == state:
                pArray.append(predSS[state][i+j]);
                j += 1;
            if j >= 1:
                avgP =  sum(pArray, 0.0)/len(pArray); 
                if avgP > p_threshold_C:
                    decPercent=weight_C*(avgP-p_shift_C);
                    for jj in range(j):
                        gpoArray[i+jj] = gpoArray[i+jj]*(1.0-decPercent) ; 
#                         gpeArray[i+jj] = gpeArray[i+jj]*(1.0-decPercent) ; 
#                         tgpeArray[i+jj] = tgpeArray[i+jj]*(1.0-decPercent) ; 
                        if DEBUG:
                            print >> sys.stdout, "decPercent=%g" % decPercent*100; 
            i+= j;
        i+= 1;
    return (gpoArray, gpeArray, tgpeArray);
#}}}
def GetGapPenalty5(ssFile, seqLength, ssFormat):#{{{
# return (gpoArray, gpeArray, tgpeArray)
# increase gapopens for SSEs and 
# decrease gapopens for random coils predicted with 
# probability > threshold. 
# gpos are estimated based on probabilities on each residue position
# 
# set different weights for the middle of the helix/sheet and the two ends of
# the helix/sheet
#
# This procedure takes >50% of the time

#     gpoArray=[54.94941]*seqLength;
#     gpeArray=[8.52492]*seqLength;
#     tgpeArray=[4.42410]*seqLength;
    gpoArray=[1.0]*seqLength;
    gpeArray=[1.0]*seqLength;
    tgpeArray=[1.0]*seqLength;
    #return (gpoArray, gpeArray, tgpeArray);

    predSS = ReadInPredSS(ssFile, ssFormat);
    if seqLength != len(predSS['ss']):
        print >> sys.stderr, "Error! aaSeqLength = %d, but the ssSeqLength = %d, for file %s" %(seqLength, len(predSS['ss']), ssFile);
        sys.exit(1);
        
    # the following procedure takes about 20% of time

    ssSeq=predSS['ss'];
    if DEBUG:
        print >> sys.stdout, "%s" % ssSeq;
    #increase the gap penalties for secondary structure elements
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
                    mid=(j+1)/2;
                    for jj in range(1,j-1):

                        posShift=1.0-abs(jj+1-mid)/float(mid); # posShift (0.0-1.0)
                        if jj == 0 or jj == j-1:
                            weightPos = 0.05;
                        else:
                            weightPos=posShift;
                        #print >> sys.stderr, "%s %d (%d) = %g"%(state, jj,j, weightPos);

                        #weightPos=1;
                        incPercent=weight_HE*(pArray[jj]-p_shift_HE)*weightPos;
                        #incPercent=weight_HE*(avgP-p_shift_HE)*weightPos;
                        gpoArray[i+jj] += incPercent; 
                        gpeArray[i+jj] += incPercent; 
                        if not i+jj <=0 or i+jj >=seqLength -1:
                            tgpeArray[i+jj] += incPercent; 
                        if DEBUG:
                            print >> sys.stdout, "incPercent=%g" % incPercent*100; 
            i+= j;
        i+= 1;
    #decrease the gap penalties for random coils
    i = 0;
    while i < seqLength:
        if ssSeq[i] == 'C':
            state = ssSeq[i];
            j = 0;
            pArray = [];
            while i+j < seqLength and ssSeq[i+j] == state:
                pArray.append(predSS[state][i+j]);
                j += 1;
            if j >= 1:
                avgP =  sum(pArray, 0.0)/len(pArray); 
                if avgP > p_threshold_C:
                    for jj in range(1,j-1):
                        decPercent=weight_C*(pArray[jj]-p_shift_C);
                        #decPercent=weight_C*(avgP-p_shift_C);
                        gpoArray[i+jj] -= decPercent ; 
                        gpeArray[i+jj] -= decPercent ; 
                        if not i+jj <=0 or i+jj >=seqLength -1:
                            tgpeArray[i+jj] -= decPercent ; 
                        if DEBUG:
                            print >> sys.stdout, "decPercent=%g" % decPercent*100; 
            i+= j;
        i+= 1;
    return (gpoArray, gpeArray, tgpeArray);
#}}}
def GetGapPenalty6(ssFile, seqLength, ssFormat):#{{{
# return (gpoArray, gpeArray, tgpeArray)
# increase gapopens for SSEs and 
# decrease gapopens for random coils predicted with 
# probability > threshold. 
# gpos are estimated based on probabilities on each residue position
# 
# set different weights for the middle of the helix/sheet and the two ends of
# the helix/sheet
# set weight_pos on random coils as well
#
# This procedure takes >50% of the time

#     gpoArray=[54.94941]*seqLength;
#     gpeArray=[8.52492]*seqLength;
#     tgpeArray=[4.42410]*seqLength;
    gpoArray=[1.0]*seqLength;
    gpeArray=[1.0]*seqLength;
    tgpeArray=[1.0]*seqLength;
    #return (gpoArray, gpeArray, tgpeArray);

    predSS = ReadInPredSS(ssFile, ssFormat);
    if seqLength != len(predSS['ss']):
        print >> sys.stderr, "Error! aaSeqLength = %d, but the ssSeqLength = %d, for file %s" %(seqLength, len(predSS['ss']), ssFile);
        sys.exit(1);
        
    # the following procedure takes about 20% of time

    ssSeq=predSS['ss'];
    if DEBUG:
        print >> sys.stdout, "%s" % ssSeq;
    #increase the gap penalties for secondary structure elements
    i = 0;
    while i < seqLength:
        if ssSeq[i] == 'H' or ssSeq[i] == 'E':
            state = ssSeq[i];
            j = 0;
            pArray = [];
            while  i+j < seqLength and ssSeq[i+j] == state:
                pArray.append(predSS[state][i+j]);
                j += 1;
            if j >= 5:
                avgP =  sum(pArray, 0.0)/len(pArray); 
                if avgP > p_threshold_HE:
                    mid=(j+1)/2;
                    for jj in range(1,j-1):

                        posShift=1.0-abs(jj+1-mid)/float(mid); # posShift (0.0-1.0)
                        if jj == 0 or jj == j-1:
                            weightPos = 0.05;
                        else:
                            weightPos=posShift;
                        #print >> sys.stderr, "%s %d (%d) = %g"%(state, jj,j, weightPos);

                        #weightPos=1;
                        incPercent=weight_HE*(pArray[jj]-p_shift_HE)*weightPos;
                        #incPercent=weight_HE*(avgP-p_shift_HE)*weightPos;
                        gpoArray[i+jj] += incPercent; 
                        gpeArray[i+jj] += incPercent; 
                        if not i+jj <= 0 or i+jj >=seqLength -1:
                            tgpeArray[i+jj] += incPercent; 
                        if DEBUG:
                            print >> sys.stdout, "incPercent=%g" % incPercent*100; 
            i+= j;
        i+= 1;
    #decrease the gap penalties for random coils
    i = 0;
    while i < seqLength:
        if ssSeq[i] == 'C':
            state = ssSeq[i];
            j = 0;
            pArray = [];
            while i+j < seqLength and ssSeq[i+j] == state:
                pArray.append(predSS[state][i+j]);
                j += 1;
            if j >= 1:
                avgP =  sum(pArray, 0.0)/len(pArray); 
                if avgP > p_threshold_C:
                    mid=(j+1)/2;
                    for jj in range(1,j-1):
                        if j < 5:
                            weightPos=1.0;
                        else:
                            posShift=1.0-abs(jj+1-mid)/float(mid); # posShift (0.0-1.0)
                            if jj == 0 or jj == j-1:
                                weightPos = 0.05;
                            else:
                                weightPos=posShift;
                        decPercent=weight_C*(pArray[jj]-p_shift_C)*weightPos;
                        #decPercent=weight_C*(avgP-p_shift_C);
                        gpoArray[i+jj] -= decPercent ; 
                        gpeArray[i+jj] -= decPercent ; 
                        if not i+jj <=0 or i+jj >=seqLength -1:
                            tgpeArray[i+jj] -= decPercent/2 ; 
                        if DEBUG:
                            print >> sys.stdout, "decPercent=%g" % decPercent*100; 
            i+= j;
        i+= 1;
    return (gpoArray, gpeArray, tgpeArray);
#}}}

def AddGapPenaltyByPredSS(fastaFile, ssPath, outPath,ssFormat):#{{{
# add gap penalties
    rootname = os.path.basename(os.path.splitext(fastaFile)[0]);
    inFilePath=os.path.dirname(fastaFile);
    if inFilePath == "":
        inFilePath='./';
    (annotationList, seqList) = myfunc.ReadFasta_without_id(fastaFile);

    if ssFormat ==0:
        if outPath == "":
            localOutPath=inFilePath;
        else:
            localOutPath=outPath;

        outFile="%s/%s.gp.fa"%(localOutPath,rootname)
        fpout = open(outFile, "w");
        for i in range(len(seqList)) :
            if DEBUG: 
                print >> sys.stdout, "%d:%s" %(i, annotationList[i])
            ssFile="%s/%s/%s_%d.ss2"%(ssPath, rootname, rootname, i);
            seqLength = len(seqList[i]);
#             gpoArray=[0]*seqLength;
#             gpeArray=[0]*seqLength;
#             tgpeArray=[0]*seqLength;
            if method_g == 0:
                (gpoArray, gpeArray, tgpeArray) = GetGapPenalty(ssFile, seqLength, ssFormat);
            elif method_g == 1:
                (gpoArray, gpeArray, tgpeArray) = GetGapPenalty1(ssFile, seqLength, ssFormat);
            elif method_g == 2:
                (gpoArray, gpeArray, tgpeArray) = GetGapPenalty2(ssFile, seqLength, ssFormat);
            elif method_g == 3:
                (gpoArray, gpeArray, tgpeArray) = GetGapPenalty3(ssFile, seqLength, ssFormat);
            elif method_g == 4:
                (gpoArray, gpeArray, tgpeArray) = GetGapPenalty4(ssFile, seqLength, ssFormat);
            elif method_g == 5:
                (gpoArray, gpeArray, tgpeArray) = GetGapPenalty5(ssFile, seqLength, ssFormat);
            elif method_g == 6:
                (gpoArray, gpeArray, tgpeArray) = GetGapPenalty6(ssFile, seqLength, ssFormat);
            else:  
                (gpoArray, gpeArray, tgpeArray) = GetGapPenalty5(ssFile, seqLength, ssFormat);

# write out the result
            fpout.write(">%s\n"% annotationList[i]);
            fpout.write("%s\n"% seqList[i]);
            fpout.write("{gpo: ");
            for j in range(seqLength):
                fpout.write("%g "% gpoArray[j]);
            fpout.write("}\n");
            fpout.write("{gpe: ");
            for j in range(seqLength):
                fpout.write("%g "% gpeArray[j]);
            fpout.write("}\n");
            fpout.write("{tgpe: ");
            for j in range(seqLength):
                fpout.write("%g "% tgpeArray[j]);
            fpout.write("}\n");
            
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
    method_g=5;

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
            elif sys.argv[i] ==  "-method-g" or  sys.argv[i] == "--method-g":
                method_g=int(sys.argv[i+1]);
                i+=2;
            elif sys.argv[i] ==  "-q" or  sys.argv[i] == "--q" or sys.argv[i] == "--quiet" or sys.argv[i] == "-quiet":
                isQuiet=True;
                i+=1;
            elif sys.argv[i] ==  "-p-shift-c" or  sys.argv[i] == "--p-shift-c":
                p_shift_C=float(sys.argv[i+1]);
                i+=2;
            elif sys.argv[i] ==  "-p-shift-he" or  sys.argv[i] == "--p-shift-he":
                p_shift_HE=float(sys.argv[i+1]);
                i+=2;
            elif sys.argv[i] ==  "-p-threshold-c" or  sys.argv[i] == "--p-threshold-c":
                p_threshold_C=float(sys.argv[i+1]);
                i+=2;
            elif sys.argv[i] ==  "-p-threshold-he" or  sys.argv[i] == "--p-threshold-he":
                p_threshold_HE=float(sys.argv[i+1]);
                i+=2;
            elif sys.argv[i] ==  "-weight-c" or  sys.argv[i] == "--weight-c":
                weight_C=float(sys.argv[i+1]);
                i+=2;
            elif sys.argv[i] ==  "-weight-he" or  sys.argv[i] == "--weight-he":
                weight_HE=float(sys.argv[i+1]);
                i+=2;
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i]);
                sys.exit(1);
        else:
            fastaFileList.append(sys.argv[i]);
            i+=1;
           

    if ssPath == "":
        print >> sys.stderr,"Error! ssPath not set.";
        sys.exit(1);
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
            numseq = AddGapPenaltyByPredSS(fastaFile, ssPath, outPath, ssFormat);
            if not isQuiet:
                print >> sys.stdout, "%d\t%d sequences in the file %s have been added gap penalties" %(cntFile, numseq, fastaFile);
            cntFile+=1;

    except :
        print >>sys.stderr, "except for the function:%s"%sys._getframe().f_code.co_name ;
        raise ;

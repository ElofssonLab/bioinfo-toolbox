#!/usr/bin/env python
# given the 
# topoalnwithdgscore file in the format
# #Topology alignment 1
# P23895    UniRef100_B7MWC5      INV     4     4      INV     4     4     DIFF     4     5 98.18  1e-25  UniRef100_B7MWC5   1.2659   6.3975  -0.1359      yes    170
# >sp|P23895|EMRE_ECOLI Multidrug transporter emrE OS=Escherichia coli (strain K12) GN=emrE PE=1 SV=1
# ------------------------------------------------------------iiMMMMMMMMMMMMMMMMMMMMMoooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiMMMMMMMMMMMMMMMMMMMMMoooooMMMMMMMMMMMMMMMMMMMMMiiiii
# /dgscoers/
# >UniRef100_Q1RAK5 evalue=5e-26 pid=98.18% alnLength=110 bitscore=119
# iiiiiiiiMMMMMMMMMMMMMMMMMMMMMoooooooooooooooooooooooooooooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMooooMMMMMMMMMMMMMMMMMMMMMiiiiiMMMMMMMMMMMMMMMMMMMMMooooo
# /dgscores/
# 
# output the selected alignments with DG 
# 
# 
import sys,re,os;
import myfunc ;
import comptopo as ct;

isPrintAllRecord=False;
method_inv=1;
method_diff=1;
threshold_dg_inv=0.0;

usage="""
Usage:    selTopoAlnWithDGScore.py [Options] [-i] topoalnwithdgscoreFile
Options:
  -o         <file>    : outputfile
  -invm      <int>     : method for selecting inverse topology, default = %d
                       : 0. sum(dgscore) < threshold_dg_inv
                       : 1. max(dgscore) < threshold_dg_inv
  -tdginv    <float>   : threshold_dg_inv, default = %g
  -h|--help            : print this help message and exit

Created 2010-09-06, updated 2010-09-07, Nanjiang
"""%(method_inv, threshold_dg_inv)

def PrintHelp():
    print usage;

def ReadTopoAlnWithDGScore(inFile):#{{{
# return topoAlnWithDGScoreList indexID2
# indexID2 is a hash table to quickly locate the index of
# topoAlnWithDGScoreList given the id2
    topoAlnWithDGScoreList=[];
    indexID2={};
    fpin = open(inFile, "r");
    lines = fpin.readlines();
    fpin.close();
    i = 0;
    cntAlign = 0;
    while i < len(lines):
        line = lines[i];
        if line.find("#Alignment") >=0:
            topoAlnWithDGScoreList.append({});
            alnIndex=int(line.split()[1]);
            cmpLine=lines[i+1].lstrip("#").strip();
            annoLine1=lines[i+2].lstrip(">").strip();
            annoLine2=lines[i+5].lstrip(">").strip();
            topo1=lines[i+3].strip();
            topo2=lines[i+6].strip();
            dgscores1=lines[i+4].replace("/","").split();
            dgscores1=[float(dg) for dg in dgscores1];
            dgscores2=lines[i+7].replace("/","").split();
            dgscores2=[float(dg) for dg in dgscores2];

            strs=cmpLine.split();
            id1=strs[0];
            id2=strs[1];
            cmpgapless=strs[2];
            cmplocal=strs[5];
            cmpglobal=strs[8];
            pid=float(strs[11]);
            evalue=float(strs[12]);
            normLogodds=float(strs[14]);
            logodds=float(strs[15]);
            reversi=float(strs[16]);
            isTMPro=strs[17];

            topoAlnWithDGScoreList[cntAlign]['alnindex'] = alnIndex; 
            topoAlnWithDGScoreList[cntAlign]['cmpline'] = cmpLine; 
            topoAlnWithDGScoreList[cntAlign]['annoline1'] = annoLine1; 
            topoAlnWithDGScoreList[cntAlign]['annoline2'] = annoLine2; 
            topoAlnWithDGScoreList[cntAlign]['id1'] = id1; 
            topoAlnWithDGScoreList[cntAlign]['id2'] = id2; 
            topoAlnWithDGScoreList[cntAlign]['topo1'] = topo1; 
            topoAlnWithDGScoreList[cntAlign]['topo2'] = topo2; 
            topoAlnWithDGScoreList[cntAlign]['dgscores1'] = dgscores1; 
            topoAlnWithDGScoreList[cntAlign]['dgscores2'] = dgscores2; 
            topoAlnWithDGScoreList[cntAlign]['cmpgapless'] =cmpgapless; 
            topoAlnWithDGScoreList[cntAlign]['cmplocal'] =cmplocal; 
            topoAlnWithDGScoreList[cntAlign]['cmpglobal'] =cmpglobal; 
            topoAlnWithDGScoreList[cntAlign]['pid'] =pid; 
            topoAlnWithDGScoreList[cntAlign]['evalue'] =evalue; 
            topoAlnWithDGScoreList[cntAlign]['normlogodds'] =normLogodds; 
            topoAlnWithDGScoreList[cntAlign]['logodds'] =logodds; 
            topoAlnWithDGScoreList[cntAlign]['reversi'] =reversi; 
            topoAlnWithDGScoreList[cntAlign]['isTMPro'] =isTMPro; 
            indexID2[id2]=(cntAlign);
            cntAlign+=1;
            i=i+7
        i+=1;
    return (topoAlnWithDGScoreList, indexID2);
#}}}
def PrintTopoAlnWithDGScore(topoAlnWithDGScore, fpout):#{{{
    fpout.write("#Alignment %d\n"%(topoAlnWithDGScore['alnindex']));
    fpout.write("#%s\n"%topoAlnWithDGScore['cmpline']);

    fpout.write(">%s\n"%topoAlnWithDGScore['annoline1']);
    fpout.write("%s\n"%topoAlnWithDGScore['topo1']);
    fpout.write("/");
    for dg in topoAlnWithDGScore['dgscores1']:
        fpout.write("%.3f " % dg);
    fpout.write("/\n");

    fpout.write(">%s\n"%topoAlnWithDGScore['annoline2']);
    fpout.write("%s\n"%topoAlnWithDGScore['topo2']);
    fpout.write("/");
    for dg in topoAlnWithDGScore['dgscores2']:
        fpout.write("%.3f " % dg);
    fpout.write("/\n");

    fpout.write("\n");
#}}}
def coverage(a1,b1,a2,b2):#{{{
    return (min(b1,b2)-max(a1,a2));
#}}}

def IsHasOverlap(start,end,startList,endList):#{{{
#start, end are integers
#startList and endList are integer list
    isOverlap=False;
    for i in range(len(startList)):
        if coverage(start,end,startList[i],endList[i]) >0:
            isOverlap=True;
            break;
    return isOverlap;
#}}}
def SelTopoAlnWithDGScore(topoAlnWithDGScoreList, fpout):#{{{
    cntAlign=0;
    for i in range(len(topoAlnWithDGScoreList)):
        if topoAlnWithDGScoreList[i]['isTMPro'] != 'yes':
            continue;
        if (topoAlnWithDGScoreList[i]['cmpgapless'] == "OK" and topoAlnWithDGScoreList[i]['cmplocal'] == "OK" and topoAlnWithDGScoreList[i]['cmpglobal'] == "OK" ):
            continue;
        # 1. -----!DIFF: in case of the same number of TM segments, applying
        # the first rule, that is, either sumdgscore2 or maxdgscore2 should be < threshold_dg_inv
        if topoAlnWithDGScoreList[i]['cmpglobal'] != "DIFF":
            sumdgscore2 = sum(dg for dg in topoAlnWithDGScoreList[i]['dgscores2']);
            maxdgscore2 = max(topoAlnWithDGScoreList[i]['dgscores2']);
            if method_inv == 0:
                if sumdgscore2 < threshold_dg_inv:
                    PrintTopoAlnWithDGScore(topoAlnWithDGScoreList[i],fpout);
                    cntAlign+=1;
            elif method_inv == 1:
                if maxdgscore2 < threshold_dg_inv:
                    PrintTopoAlnWithDGScore(topoAlnWithDGScoreList[i],fpout);
                    cntAlign+=1;
        # 2. -----DIFF: in case of inserted/deleted TM segments, applying the
        # second rule 
        else:
            #select inserted/deleted TM regions
            topo1=topoAlnWithDGScoreList[i]['topo1'];
            topo2=topoAlnWithDGScoreList[i]['topo2'];
# find iterators for TM segments
# e.g. if the topo is
# topo="--iMo--ooo---ooooMMMMiiiiiMMM--MMM--M-M-M-oooo---i----i"
# it will find 
# M
# MMMM
# MMM--MMM--M-M-M-
            r1=re.finditer("([M-]*M[M-]*)",topo1);
            r2=re.finditer("([M-]*M[M-]*)",topo2);
            start1=[];
            end1=[];
            for jr in r1:
                start1.append(jr.start(0));
                end1.append(jr.end(0));
            start2=[];
            end2=[];
            for jr in r2:
                start2.append(jr.start(0));
                end2.append(jr.end(0));

            isIndelTMWithNegativeDG=True;# whether the inserted/deleted TM segments are having negative DG values, any non negative indelTM segment can fail this boolean variable. In addition, if no indelTM are found, the first rule is applied for selection, that is, either sumdgscore2 or maxdgscore2 should be below threshold_dg_inv

            isHasNoMatchSeg=False; #whether there exists a non matched TM segments in the topology alignment

            for jj in range(len(start1)):
                if not IsHasOverlap(start1[jj],end1[jj], start2,end2):
                    isHasNoMatchSeg=True;
                    if topoAlnWithDGScoreList[i]['dgscores1'][jj] > 0.0:
                        isIndelTMWithNegativeDG = False;

            for jj in range(len(start2)):
                if not IsHasOverlap(start2[jj],end2[jj], start1,end1):
                    isHasNoMatchSeg=True;
                    if topoAlnWithDGScoreList[i]['dgscores2'][jj] > 0.0:
                        isIndelTMWithNegativeDG = False;

            if isHasNoMatchSeg:# if there exists indelTM, the second rule is applied
                if isIndelTMWithNegativeDG:
                    PrintTopoAlnWithDGScore(topoAlnWithDGScoreList[i], fpout);
                    cntAlign+=1;
            else:#if no indelTM is found, the first rule is applied
                if method_inv == 0:
                    sumdgscore2 = sum(dg for dg in topoAlnWithDGScoreList[i]['dgscores2']);
                    if sumdgscore2 < threshold_dg_inv:
                        PrintTopoAlnWithDGScore(topoAlnWithDGScoreList[i],fpout);
                        cntAlign+=1;
                elif method_inv == 1:
                    maxdgscore2 = max(topoAlnWithDGScoreList[i]['dgscores2']);
                    if maxdgscore2 < threshold_dg_inv:
                        PrintTopoAlnWithDGScore(topoAlnWithDGScoreList[i],fpout);
                        cntAlign+=1;


    print >> fpout, "#Number of selected alignment: %d" %cntAlign;
    return 0;
#}}}
if __name__ == '__main__' :
    # Check argv
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp();
        sys.exit(1);

    outFile="";
    topoalnwithdgscoreFile="";  

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
            elif sys.argv[i] ==  "-i" or  sys.argv[i] == "--infile":
                topoalnwithdgscoreFile=sys.argv[i+1];
                i+=2;
            elif sys.argv[i] ==  "-o" or  sys.argv[i] == "--outfile":
                outFile=sys.argv[i+1];
                i+=2;
            elif sys.argv[i] ==  "-invm" or  sys.argv[i] == "--invm":
                method_inv=int(sys.argv[i+1]);
                i+=2;
            elif sys.argv[i] ==  "-tdginv" or  sys.argv[i] == "--tdginv":
                threshold_dg_inv=float(sys.argv[i+1]);
                i+=2;
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i]);
                sys.exit(1);
        else:
            topoalnwithdgscoreFile=sys.argv[i];
            i+=1;
           

    if topoalnwithdgscoreFile == "":
        print >> sys.stderr,"Error! topoalnwithdgscoreFile not set.";
        sys.exit(1);

    fpout = sys.stdout;
    if outFile != "":
        fpout = open(outFile,"w");

    try :
        (topoAlnWithDGScoreList , indexID2)= ReadTopoAlnWithDGScore(topoalnwithdgscoreFile);
        SelTopoAlnWithDGScore(topoAlnWithDGScoreList, fpout);
        if fpout != sys.stdout:
            fpout.close();

    except :
        print >>sys.stderr, "except for the function:%s"%sys._getframe().f_code.co_name ;
        raise ;

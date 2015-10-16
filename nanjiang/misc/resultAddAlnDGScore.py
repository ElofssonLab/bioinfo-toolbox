#!/usr/bin/env python
# given the 
# 1. the topocmp or topocmp_and_hmmscore or result file
#   P23895    UniRef100_B7MWC5      INV     4     4      INV     4     4     DIFF     4     5 98.18  1e-25  UniRef100_B7MWC5   1.2659   6.3975  -0.1359      yes    170
# 2. the topoaln file in the format
# #Topology alignment 1
# >sp|P23895|EMRE_ECOLI Multidrug transporter emrE OS=Escherichia coli (strain K12) GN=emrE PE=1 SV=1
# ------------------------------------------------------------iiMMMMMMMMMMMMMMMMMMMMMoooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiMMMMMMMMMMMMMMMMMMMMMoooooMMMMMMMMMMMMMMMMMMMMMiiiii
# >UniRef100_Q1RAK5 evalue=5e-26 pid=98.18% alnLength=110 bitscore=119
# iiiiiiiiMMMMMMMMMMMMMMMMMMMMMoooooooooooooooooooooooooooooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMooooMMMMMMMMMMMMMMMMMMMMMiiiiiMMMMMMMMMMMMMMMMMMMMMooooo
# 3 and 4. dgscorelist in the format
# P23895 PYIYLGGAILAEVIGTTLMKF 1.464
# P23895 WPSVGTIICYCASFWLLAQTL 1.268
# 
# output the topoalnwithdgscore  file  in the format
# #Topology alignment 1
# P23895    UniRef100_B7MWC5      INV     4     4      INV     4     4     DIFF     4     5 98.18  1e-25  UniRef100_B7MWC5   1.2659   6.3975  -0.1359      yes    170
# >sp|P23895|EMRE_ECOLI Multidrug transporter emrE OS=Escherichia coli (strain K12) GN=emrE PE=1 SV=1
# ------------------------------------------------------------iiMMMMMMMMMMMMMMMMMMMMMoooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiMMMMMMMMMMMMMMMMMMMMMoooooMMMMMMMMMMMMMMMMMMMMMiiiii
# /dgscoers/
# >UniRef100_Q1RAK5 evalue=5e-26 pid=98.18% alnLength=110 bitscore=119
# iiiiiiiiMMMMMMMMMMMMMMMMMMMMMoooooooooooooooooooooooooooooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMooooMMMMMMMMMMMMMMMMMMMMMiiiiiMMMMMMMMMMMMMMMMMMMMMooooo
# /dgscores/
# 
# 
# 
import sys,re,os;
import myfunc ;
import comptopo as ct;

isPrintAllRecord=False;

usage="""
Usage:   resultAddAlnDGScore.py [Options] -topocmp topcmpFile -topoaln  topoAlnFile -qdg queryDGListFile -tdg targetDGListFile
Options:
  -o         <file>    : outputfile
  -all                 : print result for all records in topcmpFile, by default, print only records with different topologies
  -h|--help            : print this help message and exit
Created 2010-09-06, updated 2010-09-06, Nanjiang
"""

def PrintHelp():
    print usage;

def CheckTopoCMPFile(inFile):#{{{
    isCorrectFileFormat=True;
    fpin = open(inFile, "r");
    lines = fpin.readlines();
    fpin.close();
    for line in lines:
        line=line.strip();
        if line[0] != '#':
            strs=line.split();
            if len(strs) < 18 or strs[17] not in ['yes','no']:
                isCorrectFileFormat=False;
                break;
    return isCorrectFileFormat;
#}}}
def ReadTopoCMP(inFile):#{{{
    fpin = open(inFile, "r");
    lines = fpin.readlines();
    fpin.close();

    topoCMPList=[];
    cnt=0;
    for line in lines:
        strs=line.split();
        if len(strs) > 18 and line[0] != "#" :
            topoCMPList.append({});
            topoCMPList[cnt]['id1']=strs[0];
            topoCMPList[cnt]['id2']=strs[1];
            topoCMPList[cnt]['cmpgapless']=strs[2];
            topoCMPList[cnt]['cmplocal']=strs[5];
            topoCMPList[cnt]['cmpglobal']=strs[8];
            topoCMPList[cnt]['isTMPro']=strs[17];
            topoCMPList[cnt]['line']=line.strip();
            cnt +=1;
    return topoCMPList;
#}}}
def ReadPairwiseTopoAlignment(inFile):#{{{
# return (idList1, idList2,annoLineList1, annoLineList2, topoSeqList1,topoSeqList2, indexID2, pidList, evalueList)
# return topoAlnList, indexID2
# indexID2 is a hash table to quickly locate the index of
# topoAlnWithDGScoreList given the id2
    indexID2={};
    topoAlnList=[];
    fpin = open(inFile, "r");
    lines = fpin.readlines();
    fpin.close();
    i = 0;
    numAlign=0;
    cntAlign = 0;
    while i < len(lines):
        line = lines[i];
        if line.find("#Number of alignments") >=0 :
            numAlign=int(line.split(":")[1]);
        elif line.find("#Topology alignment") >=0:
            topoAlnList.append({});
            alnIndex=int(line.split()[2]);
            annoLine1=lines[i+1].lstrip(">").strip();
            annoLine2=lines[i+3].lstrip(">").strip();
            topo1=lines[i+2].strip();
            topo2=lines[i+4].strip();
            id1=ct.GetSeqIDFromAnnotation(annoLine1);
            id2=ct.GetSeqIDFromAnnotation(annoLine2);
            pid=float(ct.GetPIDFromAnnotation(annoLine2));
            evalue=float(ct.GetEvalueFromAnnotation(annoLine2));

            topoAlnList[cntAlign]['alnindex'] = alnIndex;
            topoAlnList[cntAlign]['annoline1'] = annoLine1;
            topoAlnList[cntAlign]['annoline2'] = annoLine2;
            topoAlnList[cntAlign]['id1'] = id1;
            topoAlnList[cntAlign]['id2'] = id2;
            topoAlnList[cntAlign]['topo1'] = topo1;
            topoAlnList[cntAlign]['topo2'] = topo2;
            topoAlnList[cntAlign]['pid'] =pid; 
            topoAlnList[cntAlign]['evalue'] =evalue; 

            indexID2[id2]=(cntAlign);
            cntAlign+=1;
            i=i+5
        i+=1;
    #check
    if numAlign != len(topoAlnList):
        print >> sys.stderr, "The number of alignment read in (%d) does not match the annotation (%d)" % (len(idList1), numAlign);
    return (topoAlnList, indexID2);
#}}}
def ReadDGList(inFile):#{{{
# return dgList (dictionary)
    dgList={};
    fpin = open(inFile, "r");
    lines = fpin.readlines();
    fpin.close();
    for line in lines:
        strs=line.split();
        if len(strs) == 3:
            seqID=strs[0];
            if not seqID in dgList.keys():
                dgList[seqID]=[];
            dgList[seqID].append(strs[2]);
    return dgList;
#}}}
def ResultAddAlnDGScore(topoCMPList, indexID2, topoAlnList, queryDGList, targetDGList, fpout):#{{{
    for i in range(len(topoCMPList)):
        if topoCMPList[i]['isTMPro'] != 'yes':
            continue;
        if isPrintAllRecord or (topoCMPList[i]['cmpgapless'] != "OK" or topoCMPList[i]['cmplocal'] != "OK" or topoCMPList[i]['cmpglobal'] != "OK" ):
            id2=topoCMPList[i]['id2'];
            idx=indexID2[id2];
            fpout.write("#Alignment %d\n"%(topoAlnList[idx]['alnindex']));
            fpout.write("#%s\n"%topoCMPList[i]['line']);

            fpout.write(">%s\n"%topoAlnList[idx]['annoline1']);
            fpout.write("%s\n"%topoAlnList[idx]['topo1']);
            fpout.write("/");
            for dg in queryDGList[topoCMPList[i]['id1']]:
                fpout.write("%s " % dg);
            fpout.write("/\n");

            fpout.write(">%s\n"%topoAlnList[idx]['annoline2']);
            fpout.write("%s\n"%topoAlnList[idx]['topo2']);
            fpout.write("/");
            for dg in targetDGList[topoCMPList[i]['id2']]:
                fpout.write("%s " % dg);
            fpout.write("/\n");

            print >> fpout;
    return 0;
#}}}
if __name__ == '__main__' :
    # Check argv
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp();
        sys.exit(1);

    outFile="";
    targetDGListFile="";
    queryDGListFile ="";
    topoAlnFile="";   # topoalnfile
    topoCMPFile="";   # topoCMPFile

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
            elif sys.argv[i] == "-qdg" or sys.argv[i] == "--qdg":
                queryDGListFile=sys.argv[i+1];
                i = i + 2;
            elif sys.argv[i] == "-tdg" or sys.argv[i] == "--tdg":
                targetDGListFile=sys.argv[i+1];
                i = i + 2;
            elif sys.argv[i] == "-topoaln" or sys.argv[i] == "--topoaln":
                topoAlnFile=sys.argv[i+1];
                i = i + 2;
            elif sys.argv[i] == "-topocmp" or sys.argv[i] == "--topocmp":
                topoCMPFile=sys.argv[i+1];
                i = i + 2;
            elif sys.argv[i] == "-o" or sys.argv[i] == "--outfile":
                outFile=sys.argv[i+1];
                i = i + 2;
            elif sys.argv[i] == "-all" or sys.argv[i] == "--all":
                isPrintAllRecord=True;
                i = i + 1;
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i]);
                sys.exit(1);
        else:
            print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i]);
            sys.exit(1);
           

    if targetDGListFile == "":
        print >> sys.stderr,"Error! targetDGListFile not set.";
        sys.exit(1);
    if topoAlnFile == "":
        print >> sys.stderr,"Error! topoAlnFile not set.";
        sys.exit(1);
    if topoCMPFile == "":
        print >> sys.stderr,"Error! topoCMPFile not set.";
        sys.exit(1);
    if queryDGListFile == "":
        print >> sys.stderr,"Error! queryDGListFile not set.";
        sys.exit(1);

    fpout = sys.stdout;
    if outFile != "":
        fpout = open(outFile,"w");

    try :
        # check the format of input file
        if not CheckTopoCMPFile(topoCMPFile):
            sys.stderr >> "The format of the input topoCMPFile '%s' is not correct " %topoCMPFile;
            sys.exit(1);

        (topoAlnList, indexID2) = ReadPairwiseTopoAlignment(topoAlnFile);
        (queryDGList) = ReadDGList(queryDGListFile);
#         print queryDGList;
#         print;
        (targetDGList) = ReadDGList(targetDGListFile);
        (topoCMPList) = ReadTopoCMP(topoCMPFile);
#         print indexID2;
        ResultAddAlnDGScore(topoCMPList, indexID2, topoAlnList, queryDGList, targetDGList, fpout);
        if fpout != sys.stdout:
            fpout.close();

    except :
        print >>sys.stderr, "except for the function:%s"%sys._getframe().f_code.co_name ;
        raise ;

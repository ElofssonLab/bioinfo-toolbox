#!/usr/bin/env python
# given the data produced by drawmsatopo.sh
# report the result in html table

# ChangeLog 2011-11-17
#    1. definition set to '-' if not found
#    2. table format changed.
#
#    Before:
#    PfamID Def SeqID nTMIDT NtCons NtQuery CtCons CtQuery ICons0 ...  Figure
#    PFxxxx xxxx xxx  10     1 (1)  2 (1,2)  0      0      2(6,7) 0     PNG
#
#    New
#    subtable
#    PfamID Def nTMcons nHit nHit_nr Figure
#    xxxx   xxx 10      50   9        PNG
#    subtable
#    No. hitID nTM  N-Cons N-Hit C-Cons C-Hit I-Cons-0 I-Hit-0 ... TM_Map
#    1   xxx    9   1 (0)   -      -      -      -       -         Cons: 0 9 9
#                                                                   Hit: - 9 9
# ChangeLog 2011-11-21 
#    1. working with the grouped diff.ana result file
#    2. TM map array output as well.
#       

import os;
import sys;
import libtopologycmp as lcmp;

usage="""
Usage:  reportDIFFAnaResult.py  [-l FILE] [-outpath DIR] 
                                [-datapath DIR] [-pfamdef FILE ] [-q]
                                ID [ID ...]

Description: Report the result of MSA topology comparison in HTML format

Options:
  -outpath      DIR  Set ouput path, (default: ./)
  -htmlname     STR  Set the name of the main html file, (default: index)
  -datapath     DIR  Set the data path
  -figtreepath  DIR  Set the path for phylo trees
  -ordermsapath DIR  Set the data path ordered topomsa
  -pfamdef      FILE Pfam definition file, (default:
                     /data3/data/pfam/pfamA.seed.ac-delist)  
  -famstat      FILE Write general statistic to FILE
  -writehtml  y|n    Whether write html to outpath, (default: yes)
  -l       FILE      set the idListFile
  -q                 Quiet mode
  -h, --help         Print this help message and exit

Selection control options:
  -min-nseq  INT    Set the minimal number of sequences in family, (default: 20)
  -max-nseq  INT    Set the maximal number of sequences in family, (default: 2000)
  -min-peridt FLOAT Set the minimal percentage of the IDT group, (default: 30.0)
  -max-peridt FLOAT Set the maximal percentage of the IDT group, (default: 90.0)
  -gap  FLOAT       Set the threshold (minimum) for gap fraction, (default: 0.5)
  -dg   FLOAT       Set the threshold (maximum) for DG value, (default: 1.0)
  -keyword STRING   Output those with Pfam definition matching the keyword,
                    case-insensitive

Created 2011-09-12, updated 2012-03-20, Nanjiang Shu

Examples:
    reportDIFFAnaResult.py PF10002 -datapath figure_mode2_cmp1 -outpath rst
"""

comparisonClassNameListAll=[[]]*3; # 3 is the number of method_comparison;
comparisonClassNameListAll[0]=["OK","SHIFT","INV","DIFF"];  #associated with method_comparison 0
comparisonClassNameListAll[1]=["OK","SHIFT","INV","INV_SHIFT", "DIFF"]; #associated with method_comparison 1 
comparisonClassNameList = comparisonClassNameListAll[1]; 
STR_NON_AVAILABLE="-";

categoryList=[0,999, 1,2,3,4,5];
tableNameDict={999:"tableall", 0:"tableempty", 1:"tabel1", 2:"table2",
        3:"table3", 4:"table4", 5:"table5"};
tableTitleDict={
        999:"All",
        0:"No unmapped TMs", 
        1:"With unmapped TMs at all place", 
        2:"With unmapped TMs at the internal region and one of the terminal", 
        3:"With unmapped TMs only at the internal region", 
        4:"With unmapped TMs at both terminals but not the internal region",
        5:"With unmapped TMs at only one of the terminals",
        };
BLOCK_SIZE = 100000
# typeTopoAnaDict
# 0: started from $id.fa with multiple homologous sequences
# 1: started from $id.fa with a single sequence, and the homologous sequences
# were obtained by blast searching

def PrintHelp():
    print usage;


def countseq(inFile):#{{{
    try:
        isFirstSeq=True;
        isPreviousBuffEndWithNewLine=False;
        fpin = open(inFile, "r");
        buff = fpin.read(BLOCK_SIZE);
        cntSeq = 0;
        while buff:
            if isFirstSeq and buff[0] == '>':
                cntSeq +=1;
                isFirstSeq = False;
            if isPreviousBuffEndWithNewLine and buff[0] == '>':
                cntSeq += 1;
                isPreviousBuffEndWithNewLine = False;
            cntSeq += buff.count("\n>");
            if buff[len(buff)-1] == '\n':
                isPreviousBuffEndWithNewLine = True;
            buff = fpin.read(BLOCK_SIZE);

        fpin.close();
        return cntSeq;
    except IOError:
        print >> sys.stderr,"Fail to open file %s"%(inFile);
        return 0;

#}}}
def ReadPfamDEList(infile):#{{{
    pfamDefDict={};
    try:
        fpin=open(infile);
        lines=fpin.readlines();
        fpin.close()
        for line in lines:
            strs=line.split(':');
            pfamid=strs[0].strip();
            pfamDefDict[pfamid]=strs[1].strip();
        return pfamDefDict;
    except IOError: 
        print >> sys.stderr, "%s: Error reading %s"%(sys.argv[0], infile);
        return 1;
#}}}
def GetMaxNumInternalDiff(anaFamList):#{{{
# get the maximal number of internal diff occurrences for all sequences
    maxNumInternalCons = 0;
    maxNumInternalQuery = 0;
    for anaFam in anaFamList:
        if len(anaFam['paircmp']) > 0:
            for record in anaFam['paircmp']:
                if 'ana1' in record and record['ana1'] != {}:
                    if ('internal' in record['ana1'] and 
                            len(record['ana1']['internal']) > maxNumInternalCons):
                        maxNumInternalCons = len(record['ana1']['internal']); 
                if 'ana2' in record and record['ana2'] != {}:
                    if ('internal' in record['ana2'] and
                            len(record['ana2']['internal']) >
                            maxNumInternalQuery):
                        maxNumInternalQuery = len(record['ana2']['internal']); 
    return (maxNumInternalCons, maxNumInternalQuery);
#}}}

def GetTypeTopoAna(pfamid):#{{{
    difffile1 = g_params['datapath'] + os.sep + pfamid + '.diff.ana';
    difffile2 = g_params['datapath'] + os.sep + pfamid + '.homology.diff.ana';
    if os.path.exists(difffile1):
        return 0;# with multiple sequences
    elif os.path.exists(difffile2):
        return 0;# with single sequences
    else:
        return -1;# file does not exists
#}}}
def CopyPairCmpGeneralInfo(paircmpFrom, paircmpTo ):#{{{
    paircmpTo['id1'] = paircmpFrom['id1'];
    paircmpTo['id2'] = paircmpFrom['id2'];
    paircmpTo['numTM'] = paircmpFrom['numTM'];
    paircmpTo['length'] = paircmpFrom['length'];
    return paircmpTo;
#}}}
def CopyAnaFamGeneralInfo(anaFamFrom, anaFamTo ):#{{{
    anaFamTo['pfamid']      = anaFamFrom['pfamid'];
    anaFamTo['numseq']      = anaFamFrom['numseq'];
    anaFamTo['numTM_IDT']   = anaFamFrom['numTM_IDT'];
    anaFamTo['numTMcons']   = anaFamFrom['numTMcons'];
    anaFamTo['numIDTgroup'] = anaFamFrom['numIDTgroup'];
    anaFamTo['perIDTgroup'] = anaFamFrom['perIDTgroup'];
    anaFamTo['numCMPclass'] = [];
    anaFamTo['perCMPclass'] = [];
    for i in range(len(comparisonClassNameList)):
        anaFamTo['numCMPclass'].append( anaFamFrom['numCMPclass'][i]) ;
        anaFamTo['perCMPclass'].append( anaFamFrom['perCMPclass'][i]) ;
    return anaFamTo;
#}}}
def WriteAnaFam(anaFam, fpout):#{{{
    if anaFam != {} and fpout != None:
        fpout.write("//Begin CMPMSA\n");
        lcmp.WriteOverallInfo_msa(comparisonClassNameList, anaFam['pfamid'],
                anaFam['numseq'], anaFam['numTM_IDT'], anaFam['numTMcons'],
                anaFam['numIDTgroup'], anaFam['numCMPclass'], fpout);
# Write detailed topology variation info
        cntRecord = 0;
        for record in anaFam['paircmp']:
            print >> fpout, "//Begin record", cntRecord + 1;
            lcmp.WriteOverallInfo_pairwise(record['id1'], record['id2'],
                    record['seqidt'], record['cmpclass'],
                    record['numTM1'], record['numTM2'], record['seqLength1'],
                    record['seqLength2'], fpout);
            if 'member' in record and len(record['member']) > 0:
                maxSizeID = max([len(record['member'][i]['id']) for i in
                    len(record['member'])]);
                for i in xrange(len(record['member'])):
                    mbr = record['member'][i];
                    fpout.write("Member-%d: " % i + "%*s %4d " % (maxSizeID,
                        mbr['id'], mbr['seqLength']));
                fpout.write(" dgscore (%d): "  % len(mbr['dgscore']));
                for j in range(len(mbr['dgscore'])):
                    fpout.write("%6.3f " % mbr['dgscore'][j]);
            PrintMappedArray(record['mapArray1'], record['mapArray2'],
                record['id1'], record['id2'], fpout);
            print >> fpout, "NtermTopo1", record['NtermTopo1'];
            print >> fpout, "NtermTopo2", record['NtermTopo2'];
            if 'ana1' in record and record['ana1'] != {}:
                lcmp.WriteAna(record['ana1'], fpout, "1") ;
            if 'ana2' in record and record['ana2'] != {}:
                lcmp.WriteAna(record['ana2'], fpout, "2") ;
            print >> fpout, "//End record", cntRecord+1;
            cntRecord += 1;
        fpout.write("//End CMPMSA\n");
#}}}

def WriteFamStatistic(outfile):
    try:
        fpout = open(outfile, "w")
        fpout.write("#%-3s %4s %8s %9s %9s %9s %8s %10s %10s %10s\n"%(
            "DG", "Gap", 
            "nFam_sel",
            "nFam_Nterm", "nFam_Cterm", "nFam_inter", "nFamAll",
            "nFam_numseq", "nFam_peridtgroup", "nFam_1and2"))
        fpout.write("%4.1f %4.1f %8d %9d %9d %9d %8d %10d %10d %10d\n"%(
            g_params['maxDGvalue'],
            g_params['minGapFraction'],
            g_params['numSelectedAnaFam'],
            g_params['numFamWithNterm'],
            g_params['numFamWithCterm'],
            g_params['numFamWithInter'],
            g_params['numAnaFam'],
            g_params['numFam_meet_numseq_requirement'],
            g_params['numFam_meet_peridtgroup_requirement'],
            g_params['numFam_meet_numseq_and_peridtgroup_requirement'],
            ))
        fpout.close()
    except IOError:
        print >> sys.stderr, "Failed to write to file %s"%outfile
        return 1

def GetCategory(record):#{{{      0
    category = 0;    #non assigned
    if record == {}:
        category = 0; #empty 
    else:
        isHasNtermCons=False;
        isHasCtermCons=False;
        isHasNtermQuery=False;
        isHasCtermQuery=False;
        isHasInterCons=False;
        isHasInterQuery=False;
        if record['ana1'] != {} and record['ana1']['Nterm'] != {}:
            isHasNtermCons=True;
        if record['ana1'] != {} and record['ana1']['Cterm'] != {}:
            isHasCtermCons=True;
        if record['ana2'] != {} and record['ana2']['Nterm'] != {}:
            isHasNtermQuery=True;
        if record['ana2'] != {} and record['ana2']['Cterm'] != {}:
            isHasCtermQuery=True;
        if record['ana1'] != {} and len(record['ana1']['internal'])>0:
            isHasInterCons=True;
        if record['ana2'] != {} and len(record['ana2']['internal'])>0:
            isHasInterQuery=True;
        isHasNterm = isHasNtermCons | isHasNtermQuery;
        isHasCterm = isHasCtermCons | isHasCtermQuery;
        isHasInter = isHasInterCons | isHasInterQuery;

        biNum=isHasNterm*4+isHasCterm*2+isHasInter;
        if biNum == 7:
            category=1;  # 111
        elif biNum == 3 or biNum==5:
            category=2; #  101 or 011
        elif biNum == 1:
            category=3; #  001
        elif biNum == 6:
            category=4;  # 110
        elif biNum == 2 or biNum==4:
            category=5; #  010 or 100
        else:
            category=0; # 000 empty
    return category;
#}}}
def IsHasNterm(anaFam):#{{{
    if anaFam == {} or 'paircmp' not in anaFam or len(anaFam['paircmp']) < 1:
        return False;
    else:
        for record in anaFam['paircmp']:
            if 'ana1' in record and record['ana1'] != {} and record['ana1']['Nterm'] != {}:
                return True;
            if 'ana2' in record and record['ana2'] != {} and record['ana2']['Nterm'] != {}:
                return True;
        return False;
#}}}
def IsHasCterm(anaFam):#{{{
    if anaFam == {} or 'paircmp' not in anaFam or len(anaFam['paircmp']) < 1:
        return False;
    else:
        for record in anaFam['paircmp']:
            if 'ana1' in record and record['ana1'] != {} and record['ana1']['Cterm'] != {}:
                return True;
            if 'ana2' in record and record['ana2'] != {} and record['ana2']['Cterm'] != {}:
                return True;
        return False;
#}}}
def IsHasInter(anaFam):#{{{
    if anaFam == {} or 'paircmp' not in anaFam or len(anaFam['paircmp']) < 1:
        return False;
    else:
        for record in anaFam['paircmp']:
            if 'ana1' in record and record['ana1'] != {} and len(record['ana1']['internal']) > 0:
                return True;
            if 'ana2' in record and record['ana2'] != {} and len(record['ana2']['internal']) > 0:
                return True;
        return False;
#}}}

def WriteUnmappedRecordHTMLCell(subAna,fpout):#{{{
    print >> fpout, "<td>";
    print >> fpout, "<b>%d</b><br>"%(subAna['numTMunmapped']);
    strIndex="";
    strGapPercentage="";
    strDGvalue="";
    for i in range(len(subAna['index'])):
        strIndex += "%d "%subAna['index'][i] ;
        strGapPercentage += "%.1f "%subAna['gapFraction'][i] ;
        strDGvalue += "%.2f "%subAna['DGvalue'][i] ;
    fpout.write("<font size=2>index(%s)<br>gap(%s)<br>DG(%s)</font>"%(strIndex.strip(), strGapPercentage.strip(), strDGvalue.strip()));
    print >> fpout, "</td>";
#}}}
def WriteHTMLHeader(title, fpout):#{{{
    print >> fpout, "<HTML>";
    print >> fpout, "<head>";
    print >> fpout, "<title>%s</title>"%(title);
    print >> fpout, "<script src=\"sorttable.js\"></script>"
    print >> fpout, "<style type=\"text/css\" media=\"all\"> @import \"layout.css\";";
    print >> fpout, "<!--";
    #print >> fpout, "td {font-family: \"SansSerif\", \"SansSerif\", mono; font-size: 10px; text-align: center; }";
    print >> fpout, "table.sortable thead {"
    print >> fpout, "    background-color:#eee;"
    print >> fpout, "    color:#666666;"
    print >> fpout, "    font-weight: bold;"
    print >> fpout, "    cursor: default;"
    print >> fpout, "}"

    print >> fpout, "table.maininfo";
    print >> fpout, "{";
    print >> fpout, "   border: 1px solid black;";
    print >> fpout, "   margin: 10px 10px 10px 10px;";
    print >> fpout, "   padding: 2px 10px 2px 10px;";
    print >> fpout, "}";
    print >> fpout, "table.maininfo td";
    print >> fpout, "{";
    print >> fpout, "   border: 0px solid black;";
    print >> fpout, "   padding: 2px 20px 2px 20px;";
    print >> fpout, "}";
    print >> fpout, "table.content";
    print >> fpout, "{";
    print >> fpout, "   border: 0px solid black;";
    print >> fpout, "   margin: 5px 5px 5px 5px;";
    print >> fpout, "   padding: 5px 5px 5px 5px;";
    print >> fpout, "}";
    print >> fpout, "table.content td";
    print >> fpout, "{";
    print >> fpout, "   font-family:Arial; font-size: 14px;";
    print >> fpout, "   vertical-align: top;";
    print >> fpout, "   text-align: center;";
    print >> fpout, "   border: 3px solid black;";
    print >> fpout, "   padding: 15px 5px 15px 5px;";
    print >> fpout, "}";
    print >> fpout, "table.record";
    print >> fpout, "{";
    print >> fpout, "   border: 0px solid black;";
    print >> fpout, "   padding: 2px 10px 2px 10px;";
    print >> fpout, "}";
    print >> fpout, "table.record td, th";
    print >> fpout, "{";
    print >> fpout, "   font-family:Arial; font-size: 14px;";
    print >> fpout, "   vertical-align: top;";
    print >> fpout, "   text-align: center;";
    print >> fpout, "   border: 0px solid black;";
    print >> fpout, "   padding: 5px 2px 5px 2px;";
    print >> fpout, "}";
    print >> fpout, "table.entry";
    print >> fpout, "{";
    print >> fpout, "   border: 1px solid black;";
    print >> fpout, "   border-collapse: collapse;";
    print >> fpout, "   padding: 2px 10px 2px 10px;";
    print >> fpout, "}";
    print >> fpout, "table.entry td, th";
    print >> fpout, "{";
    print >> fpout, "   font-family:Arial; font-size: 14px;";
    print >> fpout, "   vertical-align: middle;";
    print >> fpout, "   text-align: center;";
    print >> fpout, "   max-width: 200px;";
    print >> fpout, "   border: 1px solid black;";
    print >> fpout, "   padding: 2px 2px 2px 2px;";
    print >> fpout, "}";
    print >> fpout, "table.hits";
    print >> fpout, "{";
    print >> fpout, "   border: 1px solid black;";
    print >> fpout, "   border-collapse: collapse;";
    print >> fpout, "   padding: 2px 10px 2px 10px;";
    print >> fpout, "}";
    print >> fpout, "table.hits td, th";
    print >> fpout, "{";
    print >> fpout, "   font-family:Arial; font-size: 14px;";
    print >> fpout, "   vertical-align: top;";
    print >> fpout, "   text-align: center;";
    print >> fpout, "   max-width: 150px;";
    print >> fpout, "   border: 1px solid black;";
    print >> fpout, "   padding: 2px 2px 2px 2px;";
    print >> fpout, "}";
    print >> fpout, "table.subtable";
    print >> fpout, "{";
    print >> fpout, "   border: 0px solid black;";
    print >> fpout, "   border-collapse: collapse;";
    print >> fpout, "   padding: 2px 10px 2px 10px;";
    print >> fpout, "}";
    print >> fpout, "table.subtable td, th";
    print >> fpout, "{";
    print >> fpout, "   font-family:Arial; font-size: 11px;";
    print >> fpout, "   vertical-align: top;";
    print >> fpout, "   text-align: center;";
    print >> fpout, "   border: 0px solid black;";
    print >> fpout, "   padding: 2px 2px 2px 2px;";
    print >> fpout, "}";
    print >> fpout, "-->";
    print >> fpout, "</style>";
    print >> fpout, "</head>";
    print >> fpout, "<BODY>";
#}}}
def WriteHTMLTail(fpout):#{{{
    print >> fpout, "</BODY>";
    print >> fpout, "</HTML>";
#}}}
def WriteHTMLGeneralInfo(fpout, parameter):#{{{
    (minNumSeq, maxNumSeq, minPerIDTgroup, maxPerIDTgroup, minGapFraction,
            maxGapFraction, minDGvalue, maxDGvalue, numAnaFam,
            numSelectedAnaFam, numBadCons,
            numFamWithNterm, numFamWithCterm, numFamWithInter) = parameter;
    print >> fpout, "<h2>Filtering criteria</h2>";
    print >> fpout, "<table class=maininfo>"
    print >> fpout, "<tr>" ;
    print >> fpout, "<td>Number of sequences in multiple alignment</td>";
    print >> fpout, "<td>%d - %d</td>"%(minNumSeq, maxNumSeq);
    print >> fpout, "<td>%d families</td>"%(g_params['numFam_meet_numseq_requirement']);
    print >> fpout, "</tr>" ;
    print >> fpout, "<tr>" ;
    print >> fpout, "<td>Percentage of largest identical group</td>";
    print >> fpout, "<td>%.1f - %.1f</td>"%(minPerIDTgroup, maxPerIDTgroup );
    print >> fpout, "<td>%d families</td>"%(g_params['numFam_meet_peridtgroup_requirement']);
    print >> fpout, "</tr>" ;
    print >> fpout, "<tr>" ;
    print >> fpout, "<td>Gap percentage</td>";
    print >> fpout, "<td>&ge; %.1f</td>"%(minGapFraction*100);
    print >> fpout, "</tr>" ;
    print >> fpout, "<tr>" ;
    print >> fpout, "<td>DG value</td>";
    print >> fpout, "<td>&le; %.1f</td>"%(maxDGvalue);
    print >> fpout, "</tr>" ;
    print >> fpout, "</table>" ;

    print >> fpout, "<h2>Data information</h2>";
    print >> fpout, "<table class=maininfo>"

    print >> fpout, "<tr>" ;
    print >> fpout, "<td>Number of input families</td>";
    print >> fpout, "<td>%d</td>"%(numAnaFam);
    print >> fpout, "</tr>" ;

    print >> fpout, "<tr>" ;
    print >> fpout, "<td>Number of families satisfying numseq and perIDTgroup restriction</td>";
    print >> fpout, "<td>%d</td>"%(g_params['numFam_meet_numseq_and_peridtgroup_requirement']);
    print >> fpout, "</tr>" ;

#     print >> fpout, "<tr>" ;
#     print >> fpout, "<td>Number of dropped families with bad consensus topology</td>";
#     print >> fpout, "<td>%d</td>"%(numBadCons);
#     print >> fpout, "</tr>" ;

    print >> fpout, "<tr>" ;
    print >> fpout, "<td>Number of selected families</td>";
    print >> fpout, "<td>%d</td>"%(numSelectedAnaFam);
    print >> fpout, "</tr>" ;

    print >> fpout, "<tr>" ;
    print >> fpout, "<td>Number of selected families with Unmapped TMs at N-terminal</td>";
    print >> fpout, "<td>%d</td>"%(numFamWithNterm);
    print >> fpout, "</tr>" ;

    print >> fpout, "<tr>" ;
    print >> fpout, "<td>Number of selected families with Unmapped TMs at C-terminal</td>";
    print >> fpout, "<td>%d</td>"%(numFamWithCterm);
    print >> fpout, "</tr>" ;

    print >> fpout, "<tr>" ;
    print >> fpout, "<td>Number of selected families with Unmapped TMs at Internal regions</td>";
    print >> fpout, "<td>%d</td>"%(numFamWithInter);
    print >> fpout, "</tr>" ;

    print >> fpout, "</table>" ;
#}}}
def WriteTableHeader(maxNumInternalCons, maxNumInternalQuery, fpout):#{{{
    headerItemList=[];
    headerItemList.append("PfamID");
    headerItemList.append("Definition");
    headerItemList.append("SeqID");
    headerItemList.append("nTMIDT");
    headerItemList.append("N-Cons");
    headerItemList.append("N-Query");
    headerItemList.append("C-Cons");
    headerItemList.append("C-Query");
    for i in range(maxNumInternalCons):
        headerItemList.append("I-Cons-%d"%i);
    for i in range(maxNumInternalQuery):
        headerItemList.append("I-Query-%d"%i);
    headerItemList.append("Image TopoMSA");

    print >> fpout, "<tr>";
    for item in headerItemList:
        print >> fpout, "<th>";
        print >> fpout, item;
        print >> fpout, "</th>";
    print >> fpout, "</tr>";
#}}}
def WriteNonAvailableHTMLCell(fpout):#{{{
    print >> fpout, "<td>%s</td>"%STR_NON_AVAILABLE;
#}}}
def WriteImageTopoMSAHTMLCell(pfamid, fpout):#{{{
    htmlname = g_params['htmlname'];
    outpath = g_params['outpath'];
    datapath = g_params['datapath'];
    typeTopoAnaDict = g_params['typeTopoAnaDict'];
    sourceImageFile="";
    sourceThumb_imageFile = "";
    figuredir = outpath + os.sep + htmlname;
    if typeTopoAnaDict[pfamid] == 0:
        sourceImageFile = (datapath + os.sep +
                pfamid + '.sorted.orig.topomsa.png'); 
        sourceThumb_imageFile = (datapath + os.sep +'thumb.'+
                pfamid+'.sorted.orig.topomsa.png'); 
        imageFile = (figuredir + os.sep + pfamid + '.sorted.orig.topomsa.png');
        thumb_imageFile = (figuredir + os.sep + 'thumb.'+ pfamid +
                '.sorted.orig.topomsa.png');
    else:
        sourceImageFile = (datapath + os.sep +
                pfamid+'.homology.sorted.orig.topomsa.png'); 
        sourceThumb_imageFile = (datapath + os.sep +'thumb.'+
                pfamid+'.homology.sorted.orig.topomsa.png'); 
        imageFile = (figuredir + os.sep + pfamid +
                '.homology.sorted.orig.topomsa.png');
        thumb_imageFile = (figuredir + os.sep + 'thumb.'+ pfamid +
                '.homology.sorted.orig.topomsa.png');

    if not os.path.exists(imageFile):
        os.system("/bin/cp -uf %s %s"%(sourceImageFile, imageFile));
    if not os.path.exists(thumb_imageFile):
        os.system("/bin/cp -uf %s %s"%(sourceThumb_imageFile, thumb_imageFile));

    print >> fpout, "<td>";
#    print >> fpout, "Image";
    print >> fpout, ("<a href=\"%s\"target=\"_blank\">"
            % (htmlname + os.sep + os.path.basename(imageFile)));
    print >> fpout, ("<img src=\"%s\">" % (htmlname + os.sep +
        os.path.basename(thumb_imageFile)));
    print >> fpout, "</a>";
    print >> fpout, "</td>";
#}}}
def WriteOrderedTopoMSAHTMLCell(pfamid, fpout):#{{{
    htmlname = g_params['htmlname'];
    outpath = g_params['outpath'];
    ordermsapath = g_params['ordermsapath'];
    sourceImageFile = "";
    sourceThumb_imageFile = "";
    figuredir = outpath + os.sep + htmlname;
    ext = '.reordered.topomsa.png'
    imageSourceFile = (ordermsapath + os.sep + pfamid + ext); 
    imageTargetFile = (outpath + os.sep + htmlname + os.sep + pfamid + ext); 
    thumbImageSourceFile = (ordermsapath + os.sep + 'thumb.' + pfamid + ext); 
    thumbImageTargetFile = (outpath + os.sep + htmlname + os.sep + 'thumb.' + pfamid + ext); 
    if os.path.exists(imageSourceFile):
        os.system("/bin/cp -uf %s %s"%(imageSourceFile, imageTargetFile));
    if os.path.exists(thumbImageSourceFile):
        os.system("/bin/cp -uf %s %s"%(thumbImageSourceFile, thumbImageTargetFile));

    print >> fpout, "<td>";
    print >> fpout, ("<a href=\"%s\"target=\"_blank\">"
            % (htmlname + os.sep + os.path.basename(imageTargetFile)));
    print >> fpout, ("<img src=\"%s\">" % (htmlname + os.sep +
        os.path.basename(thumbImageTargetFile)));
    print >> fpout, "</a>";
    print >> fpout, "</td>";
#}}}
def WriteFigTreeMSAHTMLCell(pfamid, fpout):#{{{
    htmlname = g_params['htmlname'];
    outpath = g_params['outpath'];
    figtreepath = g_params['figtreepath'];
    sourceImageFile = "";
    sourceThumb_imageFile = "";
    figuredir = outpath + os.sep + htmlname;
    ext = '-itol.jpg'
    imageSourceFile = (figtreepath + os.sep + pfamid + ext); 
    imageTargetFile = (outpath + os.sep + htmlname + os.sep + pfamid + ext); 
    thumbImageSourceFile = (figtreepath + os.sep + 'thumb.' + pfamid + ext); 
    thumbImageTargetFile = (outpath + os.sep + htmlname + os.sep + 'thumb.' + pfamid + ext); 
    if os.path.exists(imageSourceFile):
        os.system("/bin/cp -uf %s %s"%(imageSourceFile, imageTargetFile));
    if os.path.exists(thumbImageSourceFile):
        os.system("/bin/cp -uf %s %s"%(thumbImageSourceFile, thumbImageTargetFile));

    print >> fpout, "<td>";
    print >> fpout, ("<a href=\"%s\"target=\"_blank\">"
            % (htmlname + os.sep + os.path.basename(imageTargetFile)));
    print >> fpout, ("<img src=\"%s\">" % (htmlname + os.sep +
        os.path.basename(thumbImageTargetFile)));
    print >> fpout, "</a>";
    print >> fpout, "</td>";
#}}}
def WriteTMMapHTMLCell(MapTMLine, fpout):#{{{
    text = "";
    text += "<td>\n";
    text += "<table class=subtable>\n";
    for i in xrange(2):
        mapList = MapTMLine[i].split(':')[1].split();
        text += "<tr>\n";
        for s in mapList:
            text += "<td>%s</td>"%s;
        text += "\n</tr>\n";
    text += "</table>\n";
    text += "</td>";
    print >> fpout, text;
#}}}

def WriteTableContent(anaFamList, maxNumInternalCons, maxNumInternalQuery,category, fpout):#{{{
    pfamDefDict = g_params['pfamDefDict'];
    for anaFam in anaFamList:
        if len(anaFam['paircmp'])>0:
            for anadiff in anaFam['paircmp']:
                if GetCategory(anadiff) != category :
                    continue;
                if g_params['keyword'] != "" and pfamDefDict[anaFam['pfamid']].lower().find(g_params['keyword']) == -1:
                    continue;
                print >> fpout, "<tr>";
                itemList=[];
                itemList.append(anaFam['pfamid']);
                if anaFam['pfamid'] in pfamDefDict:
                    itemList.append(pfamDefDict[anaFam['pfamid']]);
                else:
                    itemList.append(STR_NON_AVAILABLE);
                seqid = anadiff['seqid']
                seqidURL = 'http://www.ncbi.nlm.nih.gov/protein/' + seqid
                itemList.append('<a href=\"%s\" target=\"_blank\">%s</a>'%(seqid, seqidURL))
                itemList.append("%d"%anaFam['numTM_IDT']);
                for item in itemList:
                    print >> fpout, "<td>";
                    print >> fpout, item;
                    print >> fpout, "</td>";

                if anadiff['anaCons'] != {} and anadiff['anaCons']['Nterm'] != {}:
                    WriteUnmappedRecordHTMLCell(anadiff['anaCons']['Nterm'], fpout);
                else:
                    WriteNonAvailableHTMLCell(fpout);

                if anadiff['anaQuery'] != {} and anadiff['anaQuery']['Nterm'] != {}:
                    WriteUnmappedRecordHTMLCell(anadiff['anaQuery']['Nterm'], fpout);
                else:
                    WriteNonAvailableHTMLCell(fpout);

                if anadiff['anaCons'] != {} and anadiff['anaCons']['Cterm'] != {}:
                    WriteUnmappedRecordHTMLCell(anadiff['anaCons']['Cterm'], fpout);
                else:
                    WriteNonAvailableHTMLCell(fpout);


                if anadiff['anaQuery'] != {} and anadiff['anaQuery']['Cterm'] != {}:
                    WriteUnmappedRecordHTMLCell(anadiff['anaQuery']['Cterm'], fpout);
                else:
                    WriteNonAvailableHTMLCell(fpout);

                if anadiff['anaCons'] != {} and len(anadiff['anaCons']['internal'])>0:
                    for i in range(maxNumInternalCons):
                        if i < len(anadiff['anaCons']['internal']):
                            WriteUnmappedRecordHTMLCell(
                                    anadiff['anaCons']['internal'][i], fpout);
                        else:
                            WriteNonAvailableHTMLCell(fpout);
                else:
                    for i in range(maxNumInternalCons):
                        WriteNonAvailableHTMLCell(fpout);

                if anadiff['anaQuery'] != {} and len(anadiff['anaQuery']['internal'])>0:
                    for i in range(maxNumInternalQuery):
                        if i < len(anadiff['anaQuery']['internal']):
                            WriteUnmappedRecordHTMLCell(
                                    anadiff['anaQuery']['internal'][i], fpout);
                        else:
                            WriteNonAvailableHTMLCell(fpout);
                else:
                    for i in range(maxNumInternalQuery):
                        WriteNonAvailableHTMLCell(fpout);


                WriteImageTopoMSAHTMLCell(anaFam['pfamid'], fpout);

                print >> fpout, "</tr>";
#}}}

def WriteTableHeaderMain(fpout):#{{{
    headerItemList=[];
    headerItemList.append("EntryID");
    headerItemList.append("Definition");
    headerItemList.append("nTMCons");
    headerItemList.append("numSeq_IDT");
    headerItemList.append("numSeq");
    headerItemList.append("%IDT");
    if g_params['figtreepath'] != "":
        headerItemList.append("Phylo Tree");
    if g_params['ordermsapath'] != "":
        headerItemList.append("Topology MSA ordered according to phylo tree");
    headerItemList.append("Topology MSA grouped by topology comparison");
    print >> fpout, "<tr>";
    for item in headerItemList:
        print >> fpout, "<th>%s</th>"%item;
    print >> fpout, "</tr>";
#}}}
def WriteTableContentMain(anaFam, fpout):#{{{
    pfamDefDict = g_params['pfamDefDict'];
    pfamid = anaFam['pfamid']
    pfamURL = 'http://pfam.sanger.ac.uk/family/' + pfamid
    itemList=[];
    itemList.append('<a href=\"%s\" target=\"_blank\">%s</a>'%(pfamURL, pfamid))
    if anaFam['pfamid'] in pfamDefDict:
        itemList.append(pfamDefDict[anaFam['pfamid']]);
    else:
        itemList.append(STR_NON_AVAILABLE);
    itemList.append(anaFam['numTMcons']);
    itemList.append(anaFam['numIDTgroup']);
    itemList.append(anaFam['numseq']);
    itemList.append("%.1f"%(anaFam['numIDTgroup']/float(anaFam['numseq'])*100));
    print >> fpout, "<tr>";
    for item in itemList:
        print >> fpout, "<td>";
        print >> fpout, item;
        print >> fpout, "</td>";

    if g_params['figtreepath'] != "":
        WriteFigTreeMSAHTMLCell(anaFam['pfamid'], fpout)
    if g_params['ordermsapath'] != "":
        WriteOrderedTopoMSAHTMLCell(anaFam['pfamid'], fpout)
    WriteImageTopoMSAHTMLCell(anaFam['pfamid'], fpout);
    print >> fpout, "</tr>";
#}}}
def WriteSubtableMain(anaFam, fpout):#{{{
    print >> fpout, "<table class=entry>";
    WriteTableHeaderMain(fpout);
    WriteTableContentMain(anaFam, fpout);
    print >> fpout, "</table>";
#}}}

def WriteTableHeaderHits(maxNumInternalCons, maxNumInternalQuery, fpout):#{{{
    headerItemList=[];
    headerItemList.append("No.");
    headerItemList.append("HitID");
    headerItemList.append("Source");
    headerItemList.append("nTM");
    headerItemList.append("Length");
    headerItemList.append("N-Cons");
    headerItemList.append("N-Hit");
    headerItemList.append("C-Cons");
    headerItemList.append("C-Hit");
    for i in range(maxNumInternalCons):
        headerItemList.append("I-Cons-%d"%i);
    for i in range(maxNumInternalQuery):
        headerItemList.append("I-Hit-%d"%i);
    headerItemList.append("TM mapping");
    print >> fpout, "<tr>";
    for item in headerItemList:
        print >> fpout, "<th>%s</th>"%item;
    print >> fpout, "</tr>";
#}}}
def WriteTableContentHits(record, cnt, maxNumInternalCons, #{{{
        maxNumInternalQuery, fpout):
    hitsDefDict = g_params['hitsDefDict'];
    print >> fpout, "<tr>";
    itemList = [];
    itemList.append(cnt);
    if 'member' in record and len(record['member']) > 1:
        idURL= '<a href=\"http://www.ncbi.nlm.nih.gov/protein/%s\"target=\"_blank\">%s</a>'%(record['id2'],record['id2'])
        item = "%s<br>\n"%idURL
        item += "<br><br>"
        item += "(Including %d proteins)\n"%len(record['member']);
        item += "<table class=subtable>\n"
        item += "<tr>\n";
        item += "<th>%s</th>\n"%"ID";
        item += "<th>%s</th>\n"%"Length";
        item += "</tr>\n";
        for mbr in record['member']:
            item += "<tr>";
            idURL= '<a href=\"http://www.ncbi.nlm.nih.gov/protein/%s\"target=\"_blank\">%s</a>'%(mbr['id'], mbr['id'])
            item += "<td>%s</td>"%idURL;
            item += "<td>%d</td>"%mbr['seqLength'];
            item += "</tr>\n";
        item += "</table>\n";
        itemList.append(item);
    else:
        itemList.append('<a href=\"http://www.ncbi.nlm.nih.gov/protein/%s\"target=\"_blank\">%s</a>'%(record['id2'], record['id2']));

    if record['id2'] in hitsDefDict:
        itemList.append(hitsDefDict[record['id2']]);
    else:
        itemList.append(STR_NON_AVAILABLE);

    itemList.append(record['numTM2']);
    itemList.append(record['seqLength2']);
    for item in itemList:
        print >> fpout, "<td>";
        print >> fpout, item;
        print >> fpout, "</td>";

    for termstr in ['Nterm', 'Cterm']:
        for anastr in ['ana1', 'ana2']:
            if (anastr in record and record[anastr] != {} and termstr in
                    record[anastr] and record[anastr][termstr] != {}):
                WriteUnmappedRecordHTMLCell(record[anastr][termstr], fpout);
            else:
                WriteNonAvailableHTMLCell(fpout);

    maxNumInternal = (maxNumInternalCons, maxNumInternalQuery);
    anastrlist = ['ana1', 'ana2'];
    for i in range(len(anastrlist)):
        anastr = anastrlist[i];
        if (anastr in record and record[anastr] != {} and 'internal' in
                record[anastr] and len(record[anastr]['internal']) > 0):
            for j in range(maxNumInternal[i]):
                if j < len(record[anastr]['internal']):
                    WriteUnmappedRecordHTMLCell(record[anastr]['internal'][j],
                            fpout);
                else:
                    WriteNonAvailableHTMLCell(fpout);
        else:
            for j in range(maxNumInternal[i]):
                WriteNonAvailableHTMLCell(fpout);
    if 'mapTMline' in record and len(record['mapTMline']) == 2:
        WriteTMMapHTMLCell(record['mapTMline'], fpout);
#         print record['mapTMline'];
    else:
        WriteNonAvailableHTMLCell(fpout);
    print >> fpout, "</tr>";
#}}}
# def WriteSubtableHits(record, maxNumInternalCons, maxNumInternalQuery, fpout):#{{{
#     print >> fpout, "<table class=\"sortable\" border=1>";
#     WriteTableHeaderHits(maxNumInternalCons, maxNumInternalQuery, fpout);
#     WriteTableContentHits(record, 1, maxNumInternalCons,
#             maxNumInternalQuery,fpout);
#     print >> fpout, "</table>";
# #}}}

def WriteHTMLTable(tablename, tabletitle, maxNumInternalCons, maxNumInternalQuery, selectedAnaFamList, category, fpout):#{{{
    print >> fpout, "<a name=\"%s\"></a><h2>%s</h2>"%(tablename,tabletitle);
    print >> fpout, "<table class=content>";
    pfamDefDict = g_params['pfamDefDict'];
    for anaFam in selectedAnaFamList:
        if g_params['keyword'] != "" and pfamDefDict[anaFam['pfamid']].lower().find(g_params['keyword']) == -1:
            continue;
        if ('paircmp' in anaFam and len(anaFam['paircmp']) > 0):
            cntValidAnadiff = 0;
            for record in anaFam['paircmp']:
                if (record['cmpclass'] == 'DIFF' and category == 999 or
                        GetCategory(record) == category):
                    cntValidAnadiff += 1;
            if cntValidAnadiff > 0:
                print >> fpout, "<tr><td>";
                print >> fpout, "<table class=record>"

                print >> fpout, "<tr><td>";
                WriteSubtableMain(anaFam, fpout);
                print >> fpout, "</td></tr>";
                print >> fpout, "<tr><td>";
                print >> fpout, "<table class=hits>";
                WriteTableHeaderHits(maxNumInternalCons, maxNumInternalQuery,
                        fpout);
                cnt = 0;
                for record in anaFam['paircmp']:
                    if (category == 999 or GetCategory(record) == category):
                        WriteTableContentHits(record, cnt+1, maxNumInternalCons,
                                maxNumInternalQuery,fpout);
                        cnt += 1;
                print >> fpout, "</table>";
                print >> fpout, "</td></tr>";

                print >> fpout, "</table>";
                print >> fpout, "</td></tr>";

    print >> fpout, "</table>";
#}}}
def WriteHTML(selectedAnaFamList, htmlname, outpath):#{{{
# Write selectedAnaFamList in HTML format to the outpath
# The name of the main html file is called index.html if no name is give
    htmlfilename = outpath + os.sep + htmlname + ".html";
    figuredir = outpath + os.sep + htmlname;
    (maxNumInternalCons, 
            maxNumInternalQuery) = GetMaxNumInternalDiff(selectedAnaFamList);

    local_javalibpath = os.environ['HOME'] + os.sep + 'libjavascript'
    jsfileList = []
    jsfileList.append(local_javalibpath + os.sep + 'sorttable.js')
# copy jsfile to outpath 
    for f in jsfileList:
        if os.path.exists(f):
            os.system("cp -f %s %s"%(f, outpath))

    try:
        os.system("mkdir -p %s"%(figuredir)) ;
        fpout = open(htmlfilename,"w");
        title="List of unmapped TM helices";
        os.system("cp %s %s"%("/data3/doc/html/layout.css", outpath))
        WriteHTMLHeader(title, fpout);

#Write header line
        print >> fpout, "<dir id=\"Header\">";
        print >> fpout, "<h1>Analysis of topology variation in protein family</h2>";
        print >> fpout, "</dir>";

#write left panel
        print >> fpout, "<dir id=\"Menu\" style=\"line-height:24px\">";
        for i in categoryList[1:]:
            if i == 999: 
                print >> fpout, "<li><a class=\"menulist\" href=\"index.html#%s\">All</a></li>" % (tableNameDict[i]);
            else:
                print >> fpout, "<li><a class=\"menulist\" href=\"index.html#%s\">Category %d</a></li>" % (tableNameDict[i], i);

        print >> fpout, "<br>";
        print >> fpout, "<hr shade=\"noshade\" size=\"2\" width=\"100%\">";
        print >> fpout, "</dir>";

#write right panel        
        print >> fpout, "<dir id=\"Content\">";
        para1 = (g_params['minNumSeq'], g_params['maxNumSeq'],
                g_params['minPerIDTgroup'], g_params['maxPerIDTgroup'],
                g_params['minGapFraction'],
                g_params['maxGapFraction'],g_params['minDGvalue'],
                g_params['maxDGvalue'], g_params['numAnaFam'],
                g_params['numSelectedAnaFam'], g_params['numBadCons'],
                g_params['numFamWithNterm'], g_params['numFamWithCterm'],
                g_params['numFamWithInter']) ;
        WriteHTMLGeneralInfo(fpout, para1);

        for i in categoryList[1:]:
            WriteHTMLTable(tableNameDict[i], tableTitleDict[i], maxNumInternalCons, maxNumInternalQuery, selectedAnaFamList,i, fpout);

        print >> fpout, "</dir>";

        WriteHTMLTail(fpout);
        fpout.close();
        if not g_params['isQuiet']:
            print  "Result has been output to %s"%htmlfilename;

    except IOError:
        print >> sys.stderr, "%s: Error to write HTML file %s to %s. Exit." %(sys.argv[0], htmlfilename, figuredir);
        raise;
        return 1;
#}}}

def SelectAnaTermInternal(subAna, parameter):#{{{
    (minGapFraction, maxGapFraction, minDGvalue, maxDGvalue) = parameter;
    newSubAna = {};
    if subAna != {} and subAna['numTMunmapped'] > 0:
        cntTM = 0;
        for j in range(len(subAna['index'])):
            if subAna['gapFraction'][j] >= minGapFraction and subAna['gapFraction'][j] <= maxGapFraction and subAna['DGvalue'][j] >= minDGvalue and subAna['DGvalue'][j] <= maxDGvalue:
                if cntTM == 0:
                    newSubAna['index']=[];
                    newSubAna['gapFraction']=[];
                    newSubAna['DGvalue']=[];
                newSubAna['index'].append(subAna['index'][j]);
                newSubAna['gapFraction'].append(subAna['gapFraction'][j]);
                newSubAna['DGvalue'].append(subAna['DGvalue'][j]);
                cntTM += 1;
        if cntTM > 0:
            newSubAna['numTMunmapped']  = cntTM;
    return newSubAna;

#}}}
def SelectAnaDIFF(ana, parameter):#{{{
    newAna={};
    newAna['Nterm'] = {};
    newAna['Cterm'] = {};
    newAna['internal'] = [];
    
    if 'Nterm' in ana:
        newAna['Nterm'] = SelectAnaTermInternal(ana['Nterm'], parameter);
    if 'Cterm' in ana:
        newAna['Cterm'] = SelectAnaTermInternal(ana['Cterm'], parameter);

    if 'internal' in ana and len(ana['internal']) > 0:
        for anaInter in ana['internal']:
            selAnaInter = SelectAnaTermInternal(anaInter, parameter); 
            if selAnaInter != {}:
                newAna['internal'].append(selAnaInter);
    if (newAna['Nterm'] == {} and newAna['Cterm'] == {} and
            len(newAna['internal']) == 0):
        newAna = {};
    return newAna;
#}}}
def SelectAnaFam(anaFamList, parameter):#{{{
    (minNumSeq, maxNumSeq, minNumTM_IDT, maxNumTM_IDT, minPerIDTgroup,
            maxPerIDTgroup, minGapFraction, maxGapFraction, minDGvalue,
            maxDGvalue) = parameter;
    para1 = (minGapFraction, maxGapFraction, minDGvalue, maxDGvalue);
    newList = [];
    cnt=0;
    cntNumBadCons=0;
    for anaFam in anaFamList:
        if anaFam == {}:
            continue;
        if anaFam['numTM_IDT'] != anaFam['numTMcons']:
            cntNumBadCons += 1;
            continue;
        if anaFam['numseq'] < minNumSeq or anaFam['numseq'] > maxNumSeq:
            continue;
        if (anaFam['perIDTgroup'] < minPerIDTgroup or anaFam['perIDTgroup'] >
                maxPerIDTgroup):
            continue;
        if anaFam['numTM_IDT'] < minNumTM_IDT or anaFam['numTM_IDT'] > maxNumTM_IDT:
            continue;

        newAnaFam = anaFam.copy();
        newAnaFam['paircmp'] = [];
        cntPaircmp = 0;
        for record in anaFam['paircmp']:
#             print record;
            if record['cmpclass'] != 'DIFF':
                continue;
            ana1  = SelectAnaDIFF(record['ana1'], para1);
            ana2 = SelectAnaDIFF(record['ana2'], para1);
            if ana1 != {} or ana2 != {}:
                newAnaFam['paircmp'].append(record);
                newAnaFam['paircmp'][cntPaircmp]['ana1'] = ana1;
                newAnaFam['paircmp'][cntPaircmp]['ana2'] = ana2;
                cntPaircmp +=1;
        if cntPaircmp > 0:
#            print "newAnaFam", cnt, newAnaFam;
            newList.append(newAnaFam);
#            print "newList-%d"%(cnt), newList;
            cnt +=1;
#    print "newList=", newList;
    return (newList, cntNumBadCons);
#}}}

def ScanfOverallInfo_msa(line, anaFam):#{{{
    strs = line.split();
    anaFam['pfamid'] = strs[1];
    anaFam['numseq'] = int(strs[2]);
    anaFam['numTM_IDT'] = int(strs[3]);
    anaFam['numTMcons'] = int(strs[4]);
    anaFam['numIDTgroup'] = int(strs[5]);

    anaFam['numCMPclass'] = [];
    for i in range(6,11):
        anaFam['numCMPclass'].append(int(strs[i]));

    anaFam['perIDTgroup'] = float(strs[11]);

    anaFam['perCMPclass'] = [];
    for i in range(12,17):
        anaFam['perCMPclass'].append(float(strs[i]));

    anaFam['paircmp'] = [];
    return 0;
#}}}

def ReadInAnaData(pfamid, datapath, typeTopoAnaDict):#{{{
    anaFile = "";
    anaFam = {};
    if typeTopoAnaDict[pfamid] == 0:
        anaFile = datapath + os.sep + pfamid + '.grouped.diff.ana';
        if not os.path.exists(anaFile):
            anaFile = datapath + os.sep + pfamid + '.diff.ana';
        if not os.path.exists(anaFile):
            return anaFam;
    else:
        anaFile = datapath + os.sep+pfamid + '.homology.grouped.diff.ana';
        if not os.path.exists(anaFile):
            anaFile = datapath + os.sep+pfamid + '.homology.diff.ana';
        if not os.path.exists(anaFile):
            return anaFam;
    try: 
        fpin=open(anaFile,"rb");
        buff = fpin.read();
        fpin.close();
        beg = buff.find("//Begin CMPMSA");
        end = buff.find("//End CMPMSA");
        if beg >= 0 and end >= 0:
            recordBuff = buff[beg:end];
        else:
            print >> sys.stderr, "record not found in file %s" %anaFile;
            return anaFam;

        pos = beg;
# 1. Get the overall information
        tpb = recordBuff.find("\nMultipleComparison", pos);
        if tpb > 0:
            tpe = recordBuff.find('\n',tpb+1);
            if tpe > 0:
                ScanfOverallInfo_msa(recordBuff[tpb+1:tpe], anaFam);
                pos = tpe;
# 2. get all pairwise comparison records (compared to consensus)
        lcmp.ReadPairCmpResultFromBuffer(recordBuff[pos:], anaFam['paircmp']);
        return anaFam;
    except IOError:
        print >> sys.stderr, ("%s: file %s does not exist. Return." %
                (sys.argv[0], anaFile));
        #raise;
        return anaFam;
#}}}

def main(g_params):#{{{
    numArgv=len(sys.argv);
    if numArgv < 2:
        PrintHelp()
        return 1;

    idList = [];
    idListFile = "";
    i = 1;
    isNonOptionArg=False
    while i < numArgv:#{{{
        if isNonOptionArg == True:
            idList.append(sys.argv[i]);
            isNonOptionArg=False;
            i += 1;
        elif sys.argv[i] == "--":
            isNonOptionArg=True;
            i += 1;
        elif sys.argv[i][0] == "-":
            if sys.argv[i] ==  "-h" or  sys.argv[i] == "--help":
                PrintHelp();
                return 1;
            elif sys.argv[i] == "-outpath" or sys.argv[i] == "--outpath":
                g_params['outpath'] = sys.argv[i+1];
                i += 2;
            elif sys.argv[i] == "-datapath" or sys.argv[i] == "--datapath":
                g_params['datapath'] = sys.argv[i+1];
                i += 2;
            elif sys.argv[i] in ["-figtreepath", "--figtreepath"]:
                g_params['figtreepath'] = sys.argv[i+1];
                i += 2;
            elif sys.argv[i] in ["-ordermsapath", "--ordermsapath"]:
                g_params['ordermsapath'] = sys.argv[i+1];
                i += 2;
            elif sys.argv[i] == "-htmlname" or sys.argv[i] == "--htmlname":
                g_params['htmlname'] = sys.argv[i+1];
                i += 2;
            elif sys.argv[i] == "-pfamdef" or sys.argv[i] == "--pfamdef":
                g_params['pfamACDEListFile'] = sys.argv[i+1];
                i += 2;
            elif sys.argv[i] in [ "-famstat", "--famstat"]:
                g_params['outFamStatisticFile'] = sys.argv[i+1];
                i += 2;
            elif (sys.argv[i] in ["-writehtml", "--writehtml"]):
                if (sys.argv[i+1].lower())[0] == "y": 
                    g_params['isWriteHTML'] = True;
                else:
                    g_params['isWriteHTML'] = False;
                i = i + 2;
            elif sys.argv[i] == "-l" :
                idListFile=sys.argv[i+1];
                i += 2;
            elif sys.argv[i] == "-min-nseq" or sys.argv[i] == "--min-nseq":
                g_params['minNumSeq'] = int(sys.argv[i+1]);
                i += 2;
            elif sys.argv[i] == "-max-nseq" or sys.argv[i] == "--max-nseq":
                g_params['maxNumSeq'] = int(sys.argv[i+1]);
                i += 2;
            elif sys.argv[i] == "-min-peridt" or sys.argv[i] == "--min-peridt":
                g_params['minPerIDTgroup'] = float(sys.argv[i+1]);
                i += 2;
            elif sys.argv[i] == "-max-peridt" or sys.argv[i] == "--max-peridt":
                g_params['maxPerIDTgroup'] = float(sys.argv[i+1]);
                i += 2;
            elif sys.argv[i] == "-gap" or sys.argv[i] == "--gap":
                g_params['minGapFraction'] = float(sys.argv[i+1]);
                i += 2;
            elif sys.argv[i] == "-dg" or sys.argv[i] == "--dg":
                g_params['maxDGvalue'] = float(sys.argv[i+1]);
                i += 2;
            elif sys.argv[i] == "-keyword" or sys.argv[i] == "--keyword":
                g_params['keyword'] = sys.argv[i+1];
                i += 2;
            elif sys.argv[i] == "-q":
                g_params['isQuiet'] = True;
                i += 1;
            else:
                print >> sys.stderr, "Error! Wrong argument:", sys.argv[i];
                return 1;
        else:
            idList.append(sys.argv[i])
            i += 1
#}}}
    if idListFile != "":
        try:
            fp=open(idListFile,"r");
            idList+=fp.read().split();
            fp.close();
        except IOError:        
            print >> sys.stderr, "%s: file %s does not exist."%(sys.argv[0]);
            pass;

    if not os.path.exists(g_params['datapath']):
        print >> sys.stderr, "%s: datapath %s does not exist. Exit." %(sys.argv[0]);
        return 1;

    anaFamList = [];
#    print idList;
    
    for idd in idList:
        g_params['typeTopoAnaDict'][idd] = GetTypeTopoAna(idd);

    for pfamid in idList:
        tmpAnaFam = ReadInAnaData(pfamid, g_params['datapath'],
            g_params['typeTopoAnaDict'])
        if tmpAnaFam != {}:
            anaFamList.append(tmpAnaFam)
    g_params['numAnaFam'] = len(anaFamList);
    g_params['pfamDefDict'] = ReadPfamDEList(g_params['pfamACDEListFile']);
    g_params['hitsDefDict'] = ReadPfamDEList(g_params['hitsDefFile']);

# Get number of families within [minNumSeq, maxNumSeq]
# Get number of families within [minPerIDTgroup, maxPerIDTgroup]
    tc1 = 0
    tc2 = 0
    tc1and2 = 0
    for anaFam in anaFamList:
        is1 = False
        is2 = False
        if ('numseq' in anaFam and anaFam['numseq'] >= g_params['minNumSeq']
                and anaFam['numseq'] <= g_params['maxNumSeq']):
            is1 = True;
        if ('perIDTgroup' in anaFam and anaFam['perIDTgroup'] >= g_params['minPerIDTgroup']
                and anaFam['perIDTgroup'] <= g_params['maxPerIDTgroup']):
            is2 = True;
        tc1 += is1
        tc2 += is2
        tc1and2 += (is1 & is2)
    g_params['numFam_meet_numseq_requirement'] = tc1;
    g_params['numFam_meet_peridtgroup_requirement'] = tc2;
    g_params['numFam_meet_numseq_and_peridtgroup_requirement'] = tc1and2;

#checking code, 
# checked 2011-09-13 12:28:00 Tuesday Week 37: no error
#    print "Original ana File"
#    for anaFam in anaFamList:
#        WriteAnaFam(anaFam, sys.stdout);
#checking code
    parameter = (g_params['minNumSeq'], g_params['maxNumSeq'],
            g_params['minNumTM_IDT'], g_params['maxNumTM_IDT'],
            g_params['minPerIDTgroup'], g_params['maxPerIDTgroup'],
            g_params['minGapFraction'], g_params['maxGapFraction'],
            g_params['minDGvalue'], g_params['maxDGvalue']);
    (selectedAnaFamList, g_params['numBadCons'] ) = SelectAnaFam(anaFamList,
            parameter);
    g_params['numSelectedAnaFam'] = len(selectedAnaFamList);

    for anaFam in selectedAnaFamList:
        if IsHasNterm(anaFam):
            g_params['famWithNtermList'].append(anaFam['pfamid']);
        if IsHasCterm(anaFam):
            g_params['famWithCtermList'].append(anaFam['pfamid']);
        if IsHasInter(anaFam):
            g_params['famWithInterList'].append(anaFam['pfamid']);
    g_params['numFamWithNterm'] = len(g_params['famWithNtermList']);
    g_params['numFamWithCterm'] = len(g_params['famWithCtermList']);
    g_params['numFamWithInter'] = len(g_params['famWithInterList']);


    if g_params['outFamStatisticFile'] != "":
        print "Write statistics to %s"%g_params['outFamStatisticFile']
        WriteFamStatistic(g_params['outFamStatisticFile'])
#checking code, 
# checked 2011-09-13 14:55:33 Tuesday Week 37, no error
#    print "Selected ana File"
#    for anaFam in selectedAnaFamList:
#        WriteAnaFam(anaFam, sys.stdout);
#checking code
    if g_params['isWriteHTML']:
        print "Write html to %s"%g_params['outpath']
        WriteHTML(selectedAnaFamList, g_params['htmlname'], g_params['outpath']);
    return 0;
#}}}

if __name__ == '__main__' :
    g_params = {};
    g_params['numBadCons'] = 0;
    g_params['minNumSeq'] = 20;
    g_params['maxNumSeq'] = 2000;
    g_params['minPerIDTgroup'] = 30.0;
    g_params['maxPerIDTgroup'] = 90.0;
    g_params['minGapFraction'] = 0.5;
    g_params['maxGapFraction'] = 1.0;
    g_params['minDGvalue'] = -999999.0;
    g_params['maxDGvalue'] = 1.0;

    g_params['numFamWithNterm'] = 0;
    g_params['numFamWithCterm'] = 0;
    g_params['numFamWithInter'] = 0;
    g_params['famWithNtermList'] = [];
    g_params['famWithCtermList'] = [];
    g_params['famWithInterList'] = [];

    g_params['outFamStatisticFile'] = ""
#optional
    g_params['minNumTM_IDT'] = 1;
    g_params['maxNumTM_IDT'] = 100;
    g_params['htmlname'] = "index";
    g_params['isWriteHTML'] = True
    g_params['pfamACDEListFile'] = "/data3/data/pfam/pfamA.seed.ac-delist";
    g_params['isQuiet'] = False;
    g_params['outpath'] = "./";
    g_params['datapath'] = "";
    g_params['ordermsapath'] = "";
    g_params['figtreepath'] = "";
    g_params['keyword'] = "";
    g_params['numAnaFam'] = 0;
    g_params['numSelectedAnaFam'] = 0;
    g_params['pfamDefDict'] = {};
    g_params['hitsDefDict'] = {};
    g_params['typeTopoAnaDict'] = {};
    g_params['hitsDefFile'] = "/data3/wk/MPTopo/transporter_ASBT/hits.tax.def";
    sys.exit(main(g_params));

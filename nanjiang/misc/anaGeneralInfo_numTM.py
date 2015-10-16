#!/usr/bin/env python
# 
import os,sys
import tempfile
import random
usage="""
usage:   anaGeneralInfo_numTM.py [-outpath DIR] [-q] FILE

Options:
  -outpath DIR    set ouput path
  -q              quiet mode
  -h, --help      print this help message and exit

Created 2011-09-20, updated 2011-09-20, Nanjiang Shu 
"""
#"  -pfamc   FILE   set pfamc file, (default: /data3/data/pfam/Pfam-C)"

MAX_NTM_TO_SHOW=15
def PrintHelp():
    print usage;

def ReadPfamC(pfamcFile):#{{{
# @params
# pfamC_MBSet       a set of all MB (pfamid)
# pfamC_ACSet       a set of all AC (pfam clan accession number)
# PfamC_Dict_MB     a dictionary, PfamC_Dict_MB[pfamid] = [clan-accession, clan-definition]
# PfamC_Dict_AC     a dictionary, PfamC_Dict_AC[clan-ac] 
#                   = {definition:  "clan-def"; MB: [MB-list]}
    pfamC_MBSet=set([]);
    pfamC_ACSet=set([]);
    pfamC_Dict_AC={}; # each item [clan-accession, clan-definition]
    pfamC_Dict_MB={}; # each item {definition: "clan-def", MB: [MB-list]}
    try:
        fpin=open(pfamcFile, "r");
        lines=fpin.readlines();
        fpin.close();
        accession="";
        definition="";
        for line in lines:
#            print line;
            tag = line[:7];
            if tag == "#=GF AC":
                strs=line.split();
                accession = strs[2].split(".")[0];
                if not accession in pfamC_ACSet:
                    pfamC_ACSet.add(accession);
                    pfamC_Dict_AC[accession]={};
                    pfamC_Dict_AC[accession]['MB']=[];
                    pfamC_Dict_AC[accession]['definition']="";
            elif tag == "#=GF DE":
                definition = line[7:].strip();
                pfamC_Dict_AC[accession]['definition']=definition;
            elif tag == "#=GF MB":
                pfamid=line.split()[2].split(".")[0].strip(";");
                if not pfamid in pfamC_MBSet:
                    pfamC_Dict_MB[pfamid]=[];
                pfamC_Dict_MB[pfamid].append(accession);
                pfamC_Dict_MB[pfamid].append(definition);
                pfamC_MBSet.add(pfamid);
                pfamC_Dict_AC[accession]['MB'].append(pfamid);
        return (pfamC_MBSet, pfamC_ACSet, pfamC_Dict_AC, pfamC_Dict_MB);
    except IOError:
        print >> sys.stderr,"%s: Failed to read the pfamcFile %s. Exit."%(sys.argv[0],pfamcFile);
        raise;
        sys.exit();
#}}}
def WriteHTMLTail(fpout):#{{{
    print >> fpout, "</BODY>";
    print >> fpout, "</HTML>";
#}}}
def CountLine(infile):
    try:
        fpin=open(infile,"r");
        lines=fpin.readlines();
        num=len(lines);
        fpin.close();
        return num;
    except IOError:
        print >> sys.stderr,"Fail to open file %s"%(inFile);
        raise;
        return -1;

def WriteHTMLHeader(title, fpout):#{{{
    print >> fpout, "<HTML>";
    print >> fpout, "<title>%s</title>"%(title);
    print >> fpout, "<!--";
#     print >> fpout, "td {font-family: \"SansSerif\", \"SansSerif\", mono; font-size: 10px; text-align: center; }";
    print >> fpout, "table, th, td";
    print >> fpout, "{";
    print >> fpout, "   border: 1px solid black;";
    print >> fpout, "   margin: 10px 10px 20px 20px;";
    print >> fpout, "   padding: 1px 1px 2px 2px;";
    print >> fpout, "}";
    print >> fpout, "-->";
    print >> fpout, "</style>";
    print >> fpout, "<BODY>";
#}}}
def WriteHTMLGeneralInfo(fpout):#{{{
    print >> fpout, "<h2>Color scheme 1</h2>";
    print >> fpout, "<table border=\"0\">"
    print >> fpout, "<tr>" ;
    print >> fpout, "<td>Largest identical group (IDT)</td>";
    print >> fpout, "<td><font color=\"green\">%s</font></td>"%("Green");
    print >> fpout, "</tr>" ;

    print >> fpout, "<tr>" ;
    print >> fpout, "<td>Same to the consensus (OK)</td>";
    print >> fpout, "<td><font color=\"red\">%s</font></td>"%("Red");
    print >> fpout, "</tr>" ;

    print >> fpout, "<tr>" ;
    print >> fpout, "<td>Shifted to the consensus (SHIFT)</td>";
    print >> fpout, "<td><font color=\"pink\">%s</font></td>"%("Pink");
    print >> fpout, "</tr>" ;


    print >> fpout, "<tr>" ;
    print >> fpout, "<td>Inverse to the consensus (INV)</td>";
    print >> fpout, "<td><font color=\"blue\">%s</font></td>"%("Blue");
    print >> fpout, "</tr>" ;

    print >> fpout, "<tr>" ;
    print >> fpout, "<td>Inverse and shifted to the consensus (INV_SHIFT)</td>";
    print >> fpout, "<td><font color=\"lightgreen\">%s</font></td>"%("Lightgreen");
    print >> fpout, "</tr>" ;

    print >> fpout, "<tr>" ;
    print >> fpout, "<td>Different to the consensus (DIFF)</td>";
    print >> fpout, "<td><font color=\"black\"> %s</font></td>"%("Black");
    print >> fpout, "</tr>" ;

    print >> fpout, "</table>" ;

    print >> fpout, "<h2>Color scheme 2</h2>";
    print >> fpout, "<table border=\"0\">"
    print >> fpout, "<tr>" ;
    print >> fpout, "<td>Largest identical group (IDT)</td>";
    print >> fpout, "<td><font color=\"green\">%s</font></td>"%("Green");
    print >> fpout, "</tr>" ;
    print >> fpout, "<tr>" ;
    print >> fpout, "<td>Same Nterm topology</td>";
    print >> fpout, "<td><font color=\"pink\">%s</font></td>"%("Pink");
    print >> fpout, "</tr>" ;
    print >> fpout, "<tr>" ;
    print >> fpout, "<td>Different Nterm topology</td>";
    print >> fpout, "<td><font color=\"blue\">%s</font></td>"%("Blue");
    print >> fpout, "</tr>" ;
    print >> fpout, "</table>" ;

#}}}
def WriteTableHeader(fpout):#{{{
    headerItemList=[];
    headerItemList.append("numTM");
    for i in range(1,MAX_NTM_TO_SHOW):
        headerItemList.append("%d"%i);

    print >> fpout, "<tr>";
    for item in headerItemList:
        print >> fpout, "<th>";
        print >> fpout, item;
        print >> fpout, "</th>";
    print >> fpout, "</tr>";
#}}}

def WriteImageHTMLCell(pngfile, fpout):#{{{
    print >> fpout, "<td align=\"center\">";
#    print >> fpout, "Image";
    print >> fpout, "<a href=\"%s\"target=\"_blank\">" %(pngfile);
    print >> fpout, "<img src=\"%s\" width=\"100\" >" %(pngfile);
    print >> fpout, "</a>";
    print >> fpout, "</td>";
#}}}
def WriteTableContent(fpout):#{{{
    nseqs=[20,20,50,50,100,100]
    ntype=['new','cmp11']*3;
    
    for m in range(len(ntype)):
        print >> fpout, "<tr>";
        print >> fpout, "<td align=\"center\">";
        print >> fpout, "numSeq >= %d"%nseqs[m];
        print >> fpout, "</td>";
        for i  in range(1,MAX_NTM_TO_SHOW):
            pngfile="all.sorted.generalinfo.%s.nseq%d-%d.nTM%d-%d.png"%(ntype[m],nseqs[m],10000,i,i );
            WriteImageHTMLCell(pngfile, fpout);
        print >> fpout, "</tr>";
#}}}
def WriteHTMLTable(tablename, tabletitle,  fpout):#{{{
    print >> fpout, "<a name=\"%s\"></a><h2>%s</h2>"%(tablename,tabletitle);
    print >> fpout, "<table border=\"1\">";
    WriteTableHeader(fpout);
    WriteTableContent(fpout);
    print >> fpout, "</table>";
#}}}
def WriteHTML():#{{{
    htmlname="index"
    htmlfilename=outpath+os.sep+htmlname+".html";
    figuredir="./";
    try:
        fpout=open(htmlfilename,"w");
        title="Fraction of topology comparison classes in Pfam families";
        WriteHTMLHeader(title, fpout);
        WriteHTMLGeneralInfo(fpout);
        tablename="table1" ;
        tabletitle="Fraction of topology comparison classes in Pfam families with different numbers of TMs";
        WriteHTMLTable(tablename, tabletitle, fpout);
        WriteHTMLTail(fpout);
        fpout.close();
    except IOError:
        print >> sys.stderr, "%s: Error to write HTML file %s to %s. Exit." %(sys.argv[0], htmlfilename, figuredir);
        raise;
        sys.exit(1);

#}}}

if __name__ == '__main__' :
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        sys.exit()

    isQuiet=False;
    outpath="./";
    infile="";
    htmlname="index";
    pfamcFile="/data3/data/pfam/Pfam-C";

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
            elif sys.argv[i] == "-outpath" or sys.argv[i] == "--outpath":
                outpath=sys.argv[i+1];
                i += 2;
            elif sys.argv[i] == "-pfamc" or sys.argv[i] == "--pfamc":
                pfamcFile=sys.argv[i+1];
                i += 2;
            elif sys.argv[i] == "-q":
                isQuiet=True;
                i += 1;
            else:
                print >> sys.stderr, "Error! Wrong argument:", sys.argv[i];
                sys.exit()
        else:
            infile=sys.argv[i];
            i += 1

    try:
        os.system("mkdir -p %s"%(outpath)) 
        WriteHTML();
    except:
        print >> sys.stderr, "%s: Exit with failure." %(sys.argv[0]);
        raise;



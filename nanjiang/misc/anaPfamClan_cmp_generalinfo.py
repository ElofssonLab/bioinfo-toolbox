#!/usr/bin/env python
# 
import os,sys
import tempfile
import random
usage="""
usage:  anaPfamClan_cmp_generalinfo.py [-outpath DIR] [-q] FILE

Options:
  -outpath DIR    set ouput path
  -pfamc   FILE   set pfamc file, (default: /data3/data/pfam/Pfam-C)
  -q              quiet mode
  -h, --help      print this help message and exit

Created 2011-09-20, updated 2011-09-20, Nanjiang Shu 
"""

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
    print >> fpout, "<h2>Color scheme</h2>";
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
#}}}
def WriteTableHeader(fpout):#{{{
    headerItemList=[];
    headerItemList.append("ID");
    headerItemList.append("Definition");
    headerItemList.append("No. of famlies");
    headerItemList.append("Image");

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
    print >> fpout, "<img src=\"%s\" width=\"228\" height=\"164\">" %(pngfile);
    print >> fpout, "</a>";
    print >> fpout, "</td>";
#}}}
def WriteTableContent(pngfilelist,pfamC_Dict_AC, fpout):#{{{
    recordList=[];
    cnt=0;
    for f in pngfilelist:
        recordList.append({});
        datafile=os.path.splitext(f)[0]+".txt";
        numline=CountLine(datafile);
        recordList[cnt]['numline']=numline;
        recordList[cnt]['datafile']=datafile;
        recordList[cnt]['pngfile']=f;
        if f.find("pfamclan") != -1:
            acClan=os.path.basename(f).split(".")[0];
            definition=pfamC_Dict_AC[acClan]['definition'];
            recordList[cnt]['id']="Pfam Clan %s"%acClan;
            recordList[cnt]['definition']=definition;
        elif f.find("allwithoutclan") != -1:
            recordList[cnt]['id']="AllWithoutClan";
            recordList[cnt]['definition']="All families of TM proteins but without those in the Pfam Clan";
        elif f.find("allwithoneexample") != -1:
            recordList[cnt]['id']="AllWithOne";
            recordList[cnt]['definition']="All families of TM proteins and for those in the Pfam Clan, take only one example randomly";
        elif f.find("all.general") != -1:
            recordList[cnt]['id']="All";
            recordList[cnt]['definition']="All families of TM proteins";
        cnt +=1;
            
    # sort the list
    newRecordList= sorted(recordList, key=lambda k: k['numline'], reverse=True);
    
    for record in newRecordList:
        #print record;
        itemList=[];
        itemList.append(record['id']);
        itemList.append(record['definition']);
        itemList.append(record['numline']);
        print >> fpout, "<tr>";
        for item in itemList:
            print >> fpout, "<td align=\"center\">";
            print >> fpout, item;
            print >> fpout, "</td>";
        
        WriteImageHTMLCell("data"+os.sep+os.path.basename(record['pngfile']), fpout);
        print >> fpout, "</tr>";
#}}}
def WriteHTMLTable(tablename, tabletitle, pngfilelist, pfamC_Dict_AC, fpout):#{{{
    print >> fpout, "<a name=\"%s\"></a><h2>%s</h2>"%(tablename,tabletitle);
    print >> fpout, "<table border=\"1\">";
    WriteTableHeader(fpout);
    WriteTableContent(pngfilelist, pfamC_Dict_AC,fpout);
    print >> fpout, "</table>";
#}}}
def WriteHTML(pngfilelist,pfamC_Dict_AC, outpath):
    htmlfilename=outpath+os.sep+htmlname+".html";
    figuredir=outpath+os.sep+"data";
    try:
        fpout=open(htmlfilename,"w");
        title="Fraction of topology comparison classes in Pfam families";
        WriteHTMLHeader(title, fpout);
        WriteHTMLGeneralInfo(fpout);
        tablename="table1" ;
        tabletitle="Fraction of different topology comparison classes in Pfam families";
        WriteHTMLTable(tablename, tabletitle, pngfilelist, pfamC_Dict_AC, fpout);
        WriteHTMLTail(fpout);
        fpout.close();
    except IOError:
        print >> sys.stderr, "%s: Error to write HTML file %s to %s. Exit." %(sys.argv[0], htmlfilename, figuredir);
        raise;
        sys.exit(1);

#Write header line
        print >> fpout, "<dir id=\"Header\">";
        print >> fpout, "<h1>Analysis of topology variation in protein family</h2>";
        print >> fpout, "</dir>";


def AnaPfamClanCMPGeneralInfo(infile, pfamC_MBSet, pfamC_ACSet,pfamC_Dict_AC , pfamC_Dict_MB, outpath):#{{{
    try:
        fpin=open(infile, "r");
        lines=fpin.readlines();
        fpin.close();
        rootname=os.path.basename(os.path.splitext(infile)[0]);

        infoClanDict={};
        for acClan in pfamC_ACSet:
            infoClanDict[acClan]={};
            infoClanDict[acClan]['definition']=pfamC_Dict_AC[acClan]['definition'];
            infoClanDict[acClan]['numinfo']=0;
            infoClanDict[acClan]['info']=[];
        outdatapath=outpath+os.sep+"data";
        os.system("mkdir -p %s"%(outdatapath)) 

        allWithoutClanFile        = outdatapath + os.sep + "allwithoutclan.generalinfo.txt";
        allWithOneExampleClanFile = outdatapath + os.sep + "allwithoneexample.generalinfo.txt";
        sortedInfile = outdatapath + os.sep + rootname + ".sorted.txt";

        tmpfile1=tempfile.mktemp();
        tmpfile2=tempfile.mktemp();
        fpout1=open(tmpfile1,"w");
        fpout2=open(tmpfile2,"w");

        for line in lines:
            tag=line[0:5];
            if tag == 'Compa':
                strs=line.split();
                pfamid=strs[1];
                #print "pfamid=",pfamid
                if pfamid in pfamC_MBSet:
                    acClan=pfamC_Dict_MB[pfamid][0];
                    #print acClan
                    infoClanDict[acClan]['info'].append(line.strip("\n"));
                    infoClanDict[acClan]['numinfo']+=1;
                else:
                    fpout1.write("%s"%line);
                    fpout2.write("%s"%line);

        #write out
        pngfilelist=[];
        for acClan in pfamC_ACSet:
            if infoClanDict[acClan]['numinfo'] > 0:
                randidx=random.randint(0, infoClanDict[acClan]['numinfo']-1);
                fpout2.write("%s\n"%(infoClanDict[acClan]['info'][randidx]));
                outfile=outdatapath+os.sep+acClan+".generalinfo.pfamclan.txt";
                tmpfile=tempfile.mktemp();
                fpout=open(tmpfile,"w");
                for info in infoClanDict[acClan]['info']:
                    fpout.write("%s\n"%info);
                fpout.close();
                # sort 
                os.system("sort -k12,12rg -k13,13rg -k14,14rg -k15,15rg -k16,16rg -k17,17rg %s > %s"%(tmpfile, outfile));
                if not isQuiet:
                    print "%s output"%(outfile);
                os.remove(tmpfile);
                #create image
                os.system("../bin/drawcmp1allstate.py -outpath %s %s"%(outdatapath, outfile)); #DEBUG
                pngfilelist.append(os.path.splitext(outfile)[0]+".png");

        fpout1.close();
        fpout2.close();
        os.system("sort -k12,12rg -k13,13rg -k14,14rg -k15,15rg -k16,16rg -k17,17rg %s > %s"%(tmpfile1, allWithoutClanFile));
        os.system("sort -k12,12rg -k13,13rg -k14,14rg -k15,15rg -k16,16rg -k17,17rg %s > %s"%(tmpfile2, allWithOneExampleClanFile));
        os.system("sort -k12,12rg -k13,13rg -k14,14rg -k15,15rg -k16,16rg -k17,17rg %s > %s"%(infile, sortedInfile));
        os.system("../bin/drawcmp1allstate.py -outpath %s %s"%(outdatapath,sortedInfile));
        os.system("../bin/drawcmp1allstate.py -outpath %s %s"%(outdatapath, allWithoutClanFile));
        os.system("../bin/drawcmp1allstate.py -outpath %s %s"%(outdatapath, allWithOneExampleClanFile));
        pngfilelist.append(os.path.splitext(sortedInfile)[0]+".png");
        pngfilelist.append(os.path.splitext(allWithoutClanFile)[0]+".png");
        pngfilelist.append(os.path.splitext(allWithOneExampleClanFile)[0]+".png");
        os.remove(tmpfile1);
        os.remove(tmpfile2);

        #create HTML file
        WriteHTML(pngfilelist,pfamC_Dict_AC, outpath);
    except IOError:
        print >> sys.stderr,"%s: Failed to read %s"%(sys.argv[0],infile);
        raise;
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

    if infile != "":
        try:
            os.system("mkdir -p %s"%(outpath)) 
            (pfamC_MBSet, pfamC_ACSet, pfamC_Dict_AC, pfamC_Dict_MB) = ReadPfamC(pfamcFile);
#            print pfamC_ACSet;
            AnaPfamClanCMPGeneralInfo(infile, pfamC_MBSet, pfamC_ACSet,pfamC_Dict_AC , pfamC_Dict_MB, outpath);
        except:
            print >> sys.stderr, "%s: Exit with failure." %(sys.argv[0]);
            raise;

    else:
        print >> sys.stderr,"infile not set. Exit";


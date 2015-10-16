#!/usr/bin/env python
# create html file for figtree
# e.g. for figtree_kalignp_nr95

# ChangeLog 2012-03-15 
#  create sortable table using sorttable.js
# ChangeLog 2012-03-20
#  Add reordered topomsa (in the same order as phylo tree)
# ChangeLog 2012-06-11
#  -treepath -datapath 
# ChangeLog 2012-09-13
#  pdf trees are also copied

import os
import sys
import libtopologycmp as lcmp
import myfunc

usage="""
Usage:  reportFigTree.py   -l pfamidlistfile [ID [ID...]]
                           -datapath DIR
                           -topomsapath DIR

Description: Report the result of Pfam tree

Options:
  -treepath DIR      Set path for tree images
  -datapath DIR      Set path for aligment files
  -topomsapath DIR   Set path for topology MSA
  -ordermsapath DIR  Set path for topology MSA reordered according to phylo tree
  -l    LISTFILE     Set pfamidlistfile
  -outpath DIR       Set ouput path, (default: ./)
  -htmlname STR      Set the name of the main html file, (default: index)
  -pfamdef FILE      Pfam definition file, (default:
                     /data3/data/pfam/pfamA.seed.ac-delist)  
  -q                 Quiet mode
  -h, --help         Print this help message and exit

Selection control options:
  -min-numseq INT Set the minimal number of sequences for family
                  (default: 20)

Created 2012-03-14, updated 2012-09-13, Nanjiang Shu

Examples:
     reportFigTree.py PF00854 -datapath figtree_kalignp_nr95 -outpath htmlTreeKalignP_nr95
"""

BLOCK_SIZE = myfunc.BLOCK_SIZE
# typeTopoAnaDict
# 0: started from $id.fa with multiple homologous sequences
# 1: started from $id.fa with a single sequence, and the homologous sequences
# were obtained by blast searching

def PrintHelp():
    print usage;

def GetNumCluster(infile):
    try:
        fpin = open(infile)
        line = fpin.readline()
        lastAnnoLine = ""
        while line:
            if line and line[0] == ">":
                lastAnnoLine = line
            line = fpin.readline()
        fpin.close()
        clusternum = myfunc.GetClusterNoFromAnnotation(lastAnnoLine)
        if clusternum != None:
            return clusternum
        else:
            return -1
    except IOError:
        print >> sys.stderr, "Failed to read infile %s"%infile
        return -1

def ReadIDList(infile):
    try:
        fpin = open(infile,"r")
        li = fpin.read().split()
        fpin.close()
        return li
    except IOError:
        print "Failed to read listfile ", infile

def CountNonBlankLine(infile):
    try:
        fpin = open(infile,"r")
        lines = fpin.readlines()
        fpin.close()
        cnt = 0
        for line in lines:
            line = line.strip()
            if line:
                cnt += 1

        return cnt
    except IOError:
        print "Failed to read file ", infile
        return -1


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
        filetype=0
        if infile.find("tsv") != -1:
            filetype=1 #tsv file
        
        fpin=open(infile);
        lines=fpin.readlines();
        fpin.close()
        for line in lines:
            if line and line[0] != "#":
                strs = []
                if filetype == 0:
                    strs = line.split(':');
                    pfamid = strs[0].strip();
                    pfamdef = strs[1].strip()
                else:
                    strs = line.split("\t")
                    pfamid = strs[0].strip()
                    pfamdef = strs[4].strip()
                pfamDefDict[pfamid]= pfamdef
        return pfamDefDict;
    except IOError: 
        print >> sys.stderr, "%s: Error reading %s"%(sys.argv[0], infile);
        return 1;
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
    print >> fpout, "   border: 1px solid black;";
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
def WriteNonAvailableHTMLCell(fpout):#{{{
    print >> fpout, "<td>%s</td>"%STR_NON_AVAILABLE;
#}}}


def WriteSubtableHits(record, maxNumInternalCons, maxNumInternalQuery, fpout):#{{{
    print >> fpout, "<table class=hits>";
    WriteTableHeaderHits(maxNumInternalCons, maxNumInternalQuery, fpout);
    WriteTableContentHits(record, 1, maxNumInternalCons,
            maxNumInternalQuery,fpout);
    print >> fpout, "</table>";
#}}}

def WriteHTMLTable(tablename, tabletitle, idList, pfamDefDict, #{{{
        datapath, topomsapath, ordermsapath, htmlname, outpath, fpout): 
    numInputID = len (idList)

    ordermsapath = g_params['ordermsapath']
    datapath = g_params['datapath']
    topomsapath = g_params['topomsapath']
    treepath = g_params['treepath']

    print >> fpout, "<a name=\"%s\"></a><h4>%s</h4>"%(tablename,tabletitle);
    print >> fpout, "<table class=\"sortable\" border=1>";
    cntOutputID = 0

    headerItemList=[];
    headerItemList.append("No.");
    headerItemList.append("PfamID");
    headerItemList.append("Definition");
    headerItemList.append("numSeq");
    headerItemList.append("numCluster");
    headerItemList.append("Phylo Tree");
    if ordermsapath != "" and os.path.exists(ordermsapath):
        headerItemList.append("Topology MSA ordered according to phylo tree");
    if topomsapath != "" and os.path.exists(topomsapath):
        headerItemList.append("Topology MSA grouped by topology comparison");


    print >> fpout, "<tr>";
    for item in headerItemList:
        print >> fpout, "<th>";
        print >> fpout, item;
        print >> fpout, "</th>";
    print >> fpout, "</tr>";

    for i in xrange(numInputID):
        pfamid = idList[i]
        pfamURL = 'http://pfam.sanger.ac.uk/family/' + pfamid
        if pfamid in pfamDefDict:
            pfamDef = pfamDefDict[pfamid]
        else:
            pfamDef = '-'
        
        topomsafile = datapath + os.sep + pfamid + '.sorted.orig.topomsa.fa'
        if os.path.exists(topomsafile):
            numSeq = myfunc.CountFastaSeq(topomsafile)
        else:
            numSeq = -1
        if numSeq < g_params['MIN_NUMSEQ']:
            continue
        cntOutputID += 1
        
        print >> fpout, "<tr>"
#---------------------------
        print >> fpout, '<td>'
        print >> fpout, '%d'%(cntOutputID)
        print >> fpout, '</td>'
#---------------------------
        print >> fpout, '<td>'
        print >> fpout, '<a href=\"%s\" target=\"_blank\">%s</a>'%(pfamURL, pfamid)
        print >> fpout, '</td>'
#---------------------------
        print >> fpout, '<td>'
        print >> fpout, '%s'%pfamDef
        print >> fpout, '</td>'
#---------------------------
        print >> fpout, '<td>'
        print >> fpout, '%d'%numSeq
        print >> fpout, '</td>'
#---------------------------
        numCluster = -1
        clusteredmsafile = datapath + os.sep + pfamid + '.clustered.orig.topomsa.fa'
        numCluster = GetNumCluster(clusteredmsafile)
        print >> fpout, '<td>'
        if numCluster == -1:
            print >> fpout, '%s'%"-"
        else:
            print >> fpout, '%d'%numCluster

        print >> fpout, '</td>'
#---------------------------
        ext = '-itol.jpg'
        extpdf = '-itol.pdf'
        print >> fpout, '<td>'
        imageSourceFile = g_params['treepath'] + os.sep + pfamid + ext
        imageSourceFilePDF = g_params['treepath'] + os.sep + pfamid + extpdf
        imageTargetFile = outpath + os.sep + htmlname + os.sep + pfamid + ext
        imageTargetFilePDF = outpath + os.sep + htmlname + os.sep + pfamid + extpdf
        thumbImageSourceFile = g_params['treepath'] + os.sep + 'thumb.' + pfamid + ext
        thumbImageTargetFile = outpath + os.sep + htmlname + os.sep + 'thumb.' + pfamid + ext
        if os.path.exists(imageSourceFile):
            os.system("%s %s %s"%(g_params['CP_EXE'], imageSourceFile,
                imageTargetFile));
        if os.path.exists(imageSourceFilePDF):
            os.system("%s %s %s"%(g_params['CP_EXE'], imageSourceFilePDF,
                imageTargetFilePDF));
        if os.path.exists(thumbImageSourceFile):
            os.system("%s %s %s"%(g_params['CP_EXE'], thumbImageSourceFile,
                thumbImageTargetFile));
        print >> fpout, ("<a href=\"%s\"target=\"_blank\">"
                % (htmlname + os.sep + os.path.basename(imageTargetFile)));
        print >> fpout, ("<img src=\"%s\">" % (htmlname + os.sep +
            os.path.basename(thumbImageTargetFile)));
        print >> fpout, "</a>";
        print >> fpout, '</td>'
#---------------------------
        if ordermsapath != "" and os.path.exists(ordermsapath):
            print >> fpout, '<td>'
            ext = '.reordered.topomsa.png'
            imageSourceFile = ordermsapath + os.sep + pfamid + ext
            imageTargetFile = outpath + os.sep + htmlname + os.sep + pfamid + ext
            thumbImageSourceFile = ordermsapath + os.sep + 'thumb.' + pfamid + ext
            thumbImageTargetFile = outpath + os.sep + htmlname + os.sep + 'thumb.' + pfamid + ext
            if os.path.exists(imageSourceFile):
                os.system("%s %s %s"%(g_params['CP_EXE'], imageSourceFile,
                    imageTargetFile));
            if os.path.exists(thumbImageSourceFile):
                os.system("%s %s %s"%(g_params['CP_EXE'], thumbImageSourceFile,
                    thumbImageTargetFile));
            print >> fpout, ("<a href=\"%s\"target=\"_blank\">"
                    % (htmlname + os.sep + os.path.basename(imageTargetFile)));
            print >> fpout, ("<img src=\"%s\">" % (htmlname + os.sep +
                os.path.basename(thumbImageTargetFile)));
            print >> fpout, "</a>";
            print >> fpout, '</td>'
#---------------------------
        if topomsapath != "" and os.path.exists(topomsapath):
            ext = '.sorted.orig.topomsa.png'
            print >> fpout, '<td>'
            imageSourceFile = topomsapath + os.sep + pfamid + ext
            imageTargetFile = outpath + os.sep + htmlname + os.sep + pfamid + ext
            thumbImageSourceFile = topomsapath + os.sep + 'thumb.' + pfamid + ext
            thumbImageTargetFile = outpath + os.sep + htmlname + os.sep + 'thumb.' + pfamid + ext
            if os.path.exists(imageSourceFile):
                os.system("%s %s %s"%(g_params['CP_EXE'], imageSourceFile,
                    imageTargetFile));
            if os.path.exists(thumbImageSourceFile):
                os.system("%s %s %s"%(g_params['CP_EXE'], thumbImageSourceFile,
                    thumbImageTargetFile));
            print >> fpout, ("<a href=\"%s\"target=\"_blank\">"
                    % (htmlname + os.sep + os.path.basename(imageTargetFile)));
            print >> fpout, ("<img src=\"%s\">" % (htmlname + os.sep +
                os.path.basename(thumbImageTargetFile)));
            print >> fpout, "</a>";
            print >> fpout, '</td>'
#---------------------------
        print >> fpout, "</tr>"
    print >> fpout, "</table>";
#}}}
def WriteHTML(idList, pfamDefDict, datapath, topomsapath, ordermsapath, #{{{
        htmlname, outpath): 
# The name of the main html file is called index.html if no name is give
    htmlfilename = outpath + os.sep + htmlname + ".html";
    figuredir = outpath + os.sep + htmlname;
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
        title="Tree of Pfam families";
        WriteHTMLHeader(title, fpout);

#Write header line
        print >> fpout, "<dir id=\"Header\">";
        print >> fpout, "<h2>Phylogenetic tree of protein family</h2>";
        print >> fpout, "<p><u>Click header labels to sort</u></p>";
        print >> fpout, "<p>"
        print >> fpout, "Description of the phylogenetic tree<br>"
        print >> fpout, "<ul>"
        print >> fpout, "<li>Outer ring: number of TM helices</li>"
        print >> fpout, "<li>Middle ring: N-terminal state, inside (<font color=\"red\">red</font>) or outside (<font color=\"blue\">blue</font>)</li>"
        print >> fpout, "<li>"
        print >> fpout, "Inner ring: topology comparison classes"
        print >> fpout, "<ul>"
        print >> fpout, "<li>Largest identical group (IDT): <font color=\"green\">Green</font></li>"
        print >> fpout, "<li>Shifted to the consensus (SHIFT): <font color=\"pink\">Pink</font></li>"
        print >> fpout, "<li>Inverse to the consensus (INV): <font color=\"blue\">Blue</font></li>"
        print >> fpout, "<li>Inverse and shifted to the consensus (INV_SHIFT): <font color=\"lightgreen\">Lightgreen</font></li>"
        print >> fpout, "<li>Different to the consensus (DIFF): <td><font color=\"black\"> Black</font></li>"
        print >> fpout, "</ul>"
        print >> fpout, "</li>"
        print >> fpout, "<li>Minimal number of sequences in family: %d</li>"%(g_params['MIN_NUMSEQ'])
        print >> fpout, "</ul>"

        print >> fpout, "</dir>";

#write right panel        
        print >> fpout, "<dir id=\"Content\">";
        tablename = 'table1'

        tabletitle = ""
        #tabletitle = "Phylogenetic trees of Pfam families<br>Outer ring: number of TM helices<br>Middle ring: N-terminal state, inside (<font color=\"red\">red</font>) or outside (<font color=\"blue\">blue</font>)<br>Inner ring: topology comparison classes<br>Minimal number of sequences in family: %d"%(g_params['MIN_NUMSEQ'])
        WriteHTMLTable(tablename, tabletitle, idList, pfamDefDict,
                datapath, topomsapath, ordermsapath, htmlname, outpath, fpout);

        print >> fpout, "</dir>";

        WriteHTMLTail(fpout);
        fpout.close();
        print  "Result has been output to %s"%htmlfilename;

    except IOError:
        print >> sys.stderr, "%s: Error to write HTML file %s to %s. Exit." %(sys.argv[0], htmlfilename, figuredir);
        raise;
        return 1;
#}}}


def main(g_params):#{{{
    numArgv=len(sys.argv);
    if numArgv < 2:
        PrintHelp()
        return 1;

    i = 1;
    idList = []
    g_params['MIN_NUMSEQ'] = 20
    idListFile = ''
    isNonOptionArg=False
    g_params['datapath'] = "";
    g_params['topomsapath'] = "";
    g_params['ordermsapath'] = "";
    g_params['treepath'] = ""
    pfamACDEListFile = '/data3/data/pfam/pfamA.seed.ac-delist';
    outpath = ""
    htmlname = 'index'
    while i < numArgv:#{{{
        if isNonOptionArg == True:
            idList.append(sys.argv[i])
            isNonOptionArg=False;
            i += 1;
        elif sys.argv[i] == "--":
            isNonOptionArg=True;
            i += 1;
        elif sys.argv[i][0] == "-":
            if sys.argv[i] ==  "-h" or  sys.argv[i] == "--help":
                PrintHelp();
                return 1;
            elif sys.argv[i] == "-l" or sys.argv[i] == "--l":
                idListFile = sys.argv[i+1];
                i += 2;
            elif sys.argv[i] == "-datapath" or sys.argv[i] == "--datapath":
                g_params['datapath'] = sys.argv[i+1];
                i += 2;
            elif sys.argv[i] in ["-topomsapath", "--topomsapath"]:
                g_params['topomsapath'] = sys.argv[i+1];
                i += 2;
            elif sys.argv[i] in ["-treepath", "--treepath"]:
                g_params['treepath'] = sys.argv[i+1];
                i += 2;
            elif sys.argv[i] in ["-ordermsapath", "--ordermsapath"]:
                g_params['ordermsapath'] = sys.argv[i+1];
                i += 2;
            elif sys.argv[i] == "-outpath" or sys.argv[i] == "--outpath":
                outpath = sys.argv[i+1];
                i += 2;
            elif sys.argv[i] == "-htmlname" or sys.argv[i] == "--htmlname":
                htmlname = sys.argv[i+1];
                i += 2;
            elif sys.argv[i] == "-min-numseq" or sys.argv[i] == "--min-numseq":
                g_params['MIN_NUMSEQ'] = int(sys.argv[i+1]);
                i += 2;
            elif sys.argv[i] == "-pfamdef" or sys.argv[i] == "--pfamdef":
                pfamACDEListFile = sys.argv[i+1];
                i += 2;
            elif sys.argv[i] == "-q":
                isQuiet = True;
                i += 1;
            else:
                print >> sys.stderr, "Error! Wrong argument:", sys.argv[i];
                return 1;
        else:
            idList.append(sys.argv[i])
            i += 1
#}}}

    if os.path.exists(idListFile):
        idList += ReadIDList(idListFile)
    numInputID = len(idList)
    if (numInputID) <= 0:
        print >> sys.stderr, "Error. No ID set!"
        return 1

    #Read In pfamDefList
    if not os.path.exists(pfamACDEListFile):
        print >> sys.stderr, "Error! file pfamACDEListFile (%s) does not exist." %pfamACDEListFile;
        return 1;
    pfamDefDict = ReadPfamDEList(pfamACDEListFile);
#     print 'len(pfamDefDict)=', len(pfamDefDict);
    #Read in seqinfoList

    if not os.path.exists(outpath):
        os.system("mkdir -p %s"%outpath)

    g_params['OS'] = os.uname()[0];
    if g_params['OS'].find('Linux') != -1:
        g_params['CP_EXE'] = "/bin/cp -uf"
    else:
        g_params['CP_EXE'] = "/bin/cp -f"

    datapath = g_params['datapath']
    topomsapath = g_params['topomsapath']
    treepath = g_params['treepath']
    ordermsapath = g_params['ordermsapath']


    WriteHTML(idList, pfamDefDict, datapath, topomsapath, ordermsapath, htmlname, outpath);
    return 0;
#}}}

if __name__ == '__main__' :
    g_params = {};
    sys.exit(main(g_params));

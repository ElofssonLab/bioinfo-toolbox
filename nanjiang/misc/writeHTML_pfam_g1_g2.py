#!/usr/bin/python
# Description:
import os
import sys
import myfunc
import subprocess

progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

usage_short="""
Usage: %s -l pairlistfile -datapath DIR -outpath DIR
"""%(progname)

usage_ext="""
Description:
    Write HTML webpage for selected pairwise alignment between the largest two
    different groups

Format of the pairlistfile
#seqid1 seqid2 seqidt   famid         pfamdef numSeqCls1 numSeqCls2 numSeq nTM1 nTM2  isSP isPDB

OPTIONS:
  -l LISTFILE       Set the listfile
  -treepath DIR     Set path for tree files
  -msapath  DIR     Set path for msa files
  -maxperfamily INT Set the maximum number of pairs to output for each family, (default: 10)
  -swissprot INT    Whether output only when the pair are both from swissprot
  -seqlen  FILE     Seqlen file
  -shortid2fullid   FILE 
                    Mapfile from uniprot ac to full uniprot id
  -q                Quiet mode
  -h, --help        Print this help message and exit

Created 2013-12-20, updated 2014-03-10, Nanjiang Shu 
"""
usage_exp="""
Examples:
    %s -l cluster.prokar.nr0.9.txt.selected.pairlistwithfamid -datapath sel_prokar -outpath HTML_sel_prokar
"""%(progname)

def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}

def WriteHTMLHeader(title, fpout):#{{{
    print >> fpout, "<HTML>"
    print >> fpout, "<head>"
    print >> fpout, "<title>%s</title>"%(title)
    print >> fpout, "<script src=\"sorttable.js\"></script>"
    print >> fpout, "<style type=\"text/css\" media=\"all\"> @import \"layout.css\";"
    print >> fpout, "td"
    print >> fpout, "{"
    print >> fpout, "   font-family:Arial; font-size: 12px;"
    print >> fpout, "   vertical-align: top;"
    print >> fpout, "   text-align: center;"
    print >> fpout, "   border: 1px solid black;"
    print >> fpout, "   padding: 5px 5px 5px 5px;"
    print >> fpout, "}"
    print >> fpout, "<!--"
    #print >> fpout, "td {font-family: \"SansSerif\", \"SansSerif\", mono; font-size: 10px; text-align: center; }"

    print >> fpout, "table.sortable thead {"
    print >> fpout, "    background-color:#eee;"
    print >> fpout, "    color:#666666;"
    print >> fpout, "    font-weight: bold;"
    print >> fpout, "    cursor: default;"
    print >> fpout, "}"

    print >> fpout, "table.maininfo"
    print >> fpout, "{"
    print >> fpout, "   border: 1px solid black;"
    print >> fpout, "   margin: 10px 10px 10px 10px;"
    print >> fpout, "   padding: 2px 10px 2px 10px;"
    print >> fpout, "}"
    print >> fpout, "table.maininfo td"
    print >> fpout, "{"
    print >> fpout, "   border: 0px solid black;"
    print >> fpout, "   padding: 2px 20px 2px 20px;"
    print >> fpout, "}"
    print >> fpout, "table.content"
    print >> fpout, "{"
    print >> fpout, "   border: 1px solid black;"
    print >> fpout, "   margin: 5px 5px 5px 5px;"
    print >> fpout, "   padding: 5px 5px 5px 5px;"
    print >> fpout, "}"
    print >> fpout, "table.content td"
    print >> fpout, "{"
    print >> fpout, "   font-family:Arial; font-size: 14px;"
    print >> fpout, "   vertical-align: top;"
    print >> fpout, "   text-align: center;"
    print >> fpout, "   border: 3px solid black;"
    print >> fpout, "   padding: 15px 5px 15px 5px;"
    print >> fpout, "}"
    print >> fpout, "table.record"
    print >> fpout, "{"
    print >> fpout, "   border: 0px solid black;"
    print >> fpout, "   padding: 2px 10px 2px 10px;"
    print >> fpout, "}"
    print >> fpout, "table.record td, th"
    print >> fpout, "{"
    print >> fpout, "   font-family:Arial; font-size: 14px;"
    print >> fpout, "   vertical-align: top;"
    print >> fpout, "   text-align: center;"
    print >> fpout, "   border: 0px solid black;"
    print >> fpout, "   padding: 5px 2px 5px 2px;"
    print >> fpout, "}"
    print >> fpout, "table.entry"
    print >> fpout, "{"
    print >> fpout, "   border: 1px solid black;"
    print >> fpout, "   border-collapse: collapse;"
    print >> fpout, "   padding: 2px 10px 2px 10px;"
    print >> fpout, "}"
    print >> fpout, "table.entry td, th"
    print >> fpout, "{"
    print >> fpout, "   font-family:Arial; font-size: 14px;"
    print >> fpout, "   vertical-align: middle;"
    print >> fpout, "   text-align: center;"
    print >> fpout, "   max-width: 200px;"
    print >> fpout, "   border: 1px solid black;"
    print >> fpout, "   padding: 2px 2px 2px 2px;"
    print >> fpout, "}"
    print >> fpout, "table.hits"
    print >> fpout, "{"
    print >> fpout, "   border: 1px solid black;"
    print >> fpout, "   border-collapse: collapse;"
    print >> fpout, "   padding: 2px 10px 2px 10px;"
    print >> fpout, "}"
    print >> fpout, "table.hits td, th"
    print >> fpout, "{"
    print >> fpout, "   font-family:Arial; font-size: 14px;"
    print >> fpout, "   vertical-align: top;"
    print >> fpout, "   text-align: center;"
    print >> fpout, "   max-width: 150px;"
    print >> fpout, "   border: 1px solid black;"
    print >> fpout, "   padding: 2px 2px 2px 2px;"
    print >> fpout, "}"
    print >> fpout, "table.subtable"
    print >> fpout, "{"
    print >> fpout, "   border: 0px solid black;"
    print >> fpout, "   border-collapse: collapse;"
    print >> fpout, "   padding: 2px 10px 2px 10px;"
    print >> fpout, "}"
    print >> fpout, "table.subtable td, th"
    print >> fpout, "{"
    print >> fpout, "   font-family:Arial; font-size: 11px;"
    print >> fpout, "   vertical-align: top;"
    print >> fpout, "   text-align: center;"
    print >> fpout, "   border: 0px solid black;"
    print >> fpout, "   padding: 2px 2px 2px 2px;"
    print >> fpout, "}"
    print >> fpout, "-->"
    print >> fpout, "</style>"
    print >> fpout, "</head>"
    print >> fpout, "<BODY>"
#}}}
def WriteHTMLTail(fpout):#{{{
    print >> fpout, "</BODY>"
    print >> fpout, "</HTML>"
#}}}

def WriteSubTable(pairlist, datapath, outpath, fpout):
    seqLenDict = g_params['seqLenDict']
    uniprotAC2FullSeqIDMap = g_params['uniprotAC2FullSeqIDMap']
    numInputID = len (pairlist)
    print >> fpout, "<table class=\"sortable\" border=1>"

    targetpath = outpath + os.sep + "data"
    if not os.path.exists(targetpath):
        os.system("mkdir -p %s"%(targetpath))

    cntOutputID = 0

    headerItemList=[]
    headerItemList.append("No.")
    headerItemList.append("ID1")
    headerItemList.append("ID2")
    headerItemList.append("SeqIDT")
    headerItemList.append("Len1")
    headerItemList.append("Len2")
    headerItemList.append("isSwissProt")
    headerItemList.append("isPDB")
    headerItemList.append("Sequence alignment")
    headerItemList.append("Figure pairwise alignment")

# pairlist [ seqid1, seqid2, seqidt, isSP, isPDB]
    print >> fpout, "<tr>"
    for item in headerItemList:
        print >> fpout, "<th>"
        print >> fpout, item
        print >> fpout, "</th>"
    print >> fpout, "</tr>"

    cntpair_SP = 0
    cntpair_All = 0
    for i in xrange(numInputID):
        li = pairlist[i]
        seqid1 = li[0]
        seqid2 = li[1]
        seqidt = li[2]
        isSP = li[3]
        isPDB = li[4]

        if isSP < g_params['isSP_threshold']:
            continue

        if cntpair_All >= g_params['max_num_output_per_family']:
            continue

        try:
            seqLen1 = seqLenDict[seqid1]
        except KeyError:
            seqLen1 = -1
        try:
            seqLen2 = seqLenDict[seqid2]
        except KeyError:
            seqLen2 = -1

        try:
            full_seqid1 = uniprotAC2FullSeqIDMap[seqid1]
        except KeyError:
            full_seqid1 = seqid1
        try:
            full_seqid2 = uniprotAC2FullSeqIDMap[seqid2]
        except KeyError:
            full_seqid2 = seqid2


        cntOutputID += 1
        print >> fpout, "<tr>"
# 1. No ---------------------------
        print >> fpout, '<td>'
        print >> fpout, '%d'%(cntOutputID)
        print >> fpout, '</td>'
# 2. ID1 ---------------------------
        print >> fpout, '<td>'
        print >> fpout, '<a href=\"http://www.uniprot.org/uniprot/%s\"target=\"_blank\">%s</a>'%(seqid1, full_seqid1)
        print >> fpout, '</td>'
# 3. ID2 ---------------------------
        print >> fpout, '<td>'
        print >> fpout, '<a href=\"http://www.uniprot.org/uniprot/%s\"target=\"_blank\">%s</a>'%(seqid2, full_seqid2)
        print >> fpout, '</td>'
# 4. SeqIDT ---------------------------
        print >> fpout, '<td>'
        print >> fpout, '%.1f'%(seqidt)
        print >> fpout, '</td>'
# 5. Len1 ---------------------------
        print >> fpout, '<td>'
        print >> fpout, '%d'%(seqLen1)
        print >> fpout, '</td>'
# 6. Len2 ---------------------------
        print >> fpout, '<td>'
        print >> fpout, '%d'%(seqLen2)
        print >> fpout, '</td>'
# 7 isSP ---------------------------
        print >> fpout, '<td>'
        print >> fpout, '%d'%(isSP)
        print >> fpout, '</td>'
# 8 isPDB ---------------------------
        print >> fpout, '<td>'
        print >> fpout, '%d'%(isPDB)
        print >> fpout, '</td>'
# 9 sequence alignment  ---------------------------
        sourcefile = (datapath + os.sep + "%s_%s"%(seqid1,
            seqid2) + ".topoaln.krbias.html")
        cmd = ["/bin/cp", "-f", sourcefile, targetpath]
        try:
            subprocess.check_call(cmd)
        except subprocess.CalledProcessError, e:
            print e
        print >> fpout, '<td>'
        print >> fpout, ("<a href=\"%s\" target=\"_blank\">"
                % ("data" + os.sep +
                    os.path.basename(sourcefile)))
        print >> fpout, os.path.basename(sourcefile)
        print >> fpout, "</a>"
        print >> fpout, '</td>'

# 10 Figure pairwise alignment -----------------------
        sourcefile = (datapath + os.sep + "%s_%s"%(seqid1,
            seqid2) + ".topoaln.krbias-crop.png")
        cmd = ["/bin/cp", "-f", sourcefile, targetpath]
        try:
            subprocess.check_call(cmd)
        except subprocess.CalledProcessError, e:
            print e

        sourcefile_thumb = (datapath + os.sep + "%s_%s"%(seqid1,
            seqid2) + ".topoaln.krbias-crop.thumb.png")
        cmd = ["/bin/cp", "-f", sourcefile_thumb, targetpath]
        try:
            subprocess.check_call(cmd)
        except subprocess.CalledProcessError, e:
            print e

        print >> fpout, '<td>'
        print >> fpout, ("<a href=\"%s\" target=\"_blank\">"
                % ( "data" + os.sep +
                    os.path.basename(sourcefile)))
        print >> fpout, ("<img src=\"%s\" target=\"_blank\" height=\"100\" width=\"200\">"
                % ( "data" + os.sep +
                    os.path.basename(sourcefile_thumb)))
        print >> fpout, "</a>"
        print >> fpout, '</td>'
# Finish writing the row =====================================
        print >> fpout, "</tr>"


        if isSP >= g_params['isSP_threshold']:
            cntpair_SP += 1
        cntpair_All += 1

    print >> fpout, "</table>"
#}}}
def WriteHTMLTable(tablename, tabletitle, pairInfoList, datapath, htmlname, #{{{
        outpath, fpout):
    numInputID = len (pairInfoList)
    print >> fpout, "<a name=\"%s\"></a><h4>%s</h4>"%(tablename,tabletitle)
    print >> fpout, "<table class=\"sortable\" border=1>"

    targetpath = outpath + os.sep + "data"
    if not os.path.exists(targetpath):
        os.system("mkdir -p %s"%(targetpath))

    cntOutputID = 0

    headerItemList=[]
    headerItemList.append("No.")
    headerItemList.append("PfamID")
    headerItemList.append("PfamDef")
    headerItemList.append("NumSeq")
    headerItemList.append("#Group1")
    headerItemList.append("#Group2")
    headerItemList.append("nTM_Group1")
    headerItemList.append("nTM_Group2")
    headerItemList.append("Figure MSA")
    headerItemList.append("Figure Tree")
    headerItemList.append("Example pairwise alignment between group 1 and group 2")

    print >> fpout, "<tr>"
    for item in headerItemList:
        print >> fpout, "<th>"
        print >> fpout, item
        print >> fpout, "</th>"
    print >> fpout, "</tr>"

    for i in xrange(numInputID):
        info = pairInfoList[i]
        pfamid = info[0]
        dt = info[1]
        #checking whether it has pairs with sequences both from swissprot
        pairlist = dt['pairlist']
        numpair = len(pairlist)
        isQualify_isSP = False
        for jj in xrange(numpair):
            isSP = pairlist[jj][3]
            if isSP >= g_params['isSP_threshold']:
                isQualify_isSP = True
                break

        if not isQualify_isSP:
            continue


        cntOutputID += 1
        print >> fpout, "<tr>"
# 1. No ---------------------------
        print >> fpout, '<td>'
        print >> fpout, '%d'%(cntOutputID)
        print >> fpout, '</td>'
# 2. PfamID ---------------------------
        pfamURL = 'http://pfam.sanger.ac.uk/family/' + pfamid
        print >> fpout, '<td>'
        print >> fpout, '<a href=\"%s\" target=\"_blank\">%s</a>'%(pfamURL, pfamid)
        print >> fpout, '</td>'
# 3 PfamDef ---------------------------
        pfamdef = dt['pfamdef']
        print >> fpout, '<td>'
        print >> fpout, '%s'%(pfamdef)
        print >> fpout, '</td>'
# 4 NumSeq ---------------------------
        numseq = dt['numseq']
        print >> fpout, '<td>'
        print >> fpout, '%d'%(numseq)
        print >> fpout, '</td>'
# 5 #Group1---------------------------
        print >> fpout, '<td>'
        print >> fpout, '%d'%(dt['numSeqCls1'])
        print >> fpout, '</td>'
# 6 #Group2---------------------------
        print >> fpout, '<td>'
        print >> fpout, '%d'%(dt['numSeqCls2'])
        print >> fpout, '</td>'
# 7 #nTM_Group1---------------------------
        print >> fpout, '<td>'
        print >> fpout, '%d'%(dt['nTM_Group1'])
        print >> fpout, '</td>'
# 8 #nTM_Group2---------------------------
        print >> fpout, '<td>'
        print >> fpout, '%d'%(dt['nTM_Group2'])
        print >> fpout, '</td>'

# 9 #Figure MSA---------------------------
        if g_params['msapath'] != "" and os.path.exists(g_params['msapath']):
            ext = '.sorted.orig.topomsa.png'
            print >> fpout, '<td>'
            imageSourceFile = g_params['msapath'] + os.sep + pfamid + ext
            imageTargetFile = outpath + os.sep + htmlname + os.sep + pfamid + ext
            thumbImageSourceFile = g_params['msapath'] + os.sep + 'thumb.' + pfamid + ext
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

# 10 #Figure Tree---------------------------
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

# 11 #pairwise alignment ---------------------------
        pairlist = dt['pairlist']
        # seqid1, seqid2, seqidt, isSP, isPDB
        print >> fpout, '<td>'
        WriteSubTable(pairlist, datapath, outpath, fpout)
        print >> fpout, '</td>'

# Finish writing the row =====================================
        print >> fpout, "</tr>"
    print >> fpout, "</table>"
#}}}
def ReadPairInfo(infile):#{{{
    """
    Format of the pairlistfile
    #seqid1 seqid2 seqidt   famid         pfamdef numSeqCls1 numSeqCls2 numSeq nTM1 nTM2  isSP isPDB

    Output:
        pairInfoDict  {pfamid: {'':, ''}}
    """
    hdl = myfunc.ReadLineByBlock(infile)
    if hdl.failure:
        return {}
    lines = hdl.readlines()
    dt = {}
    while lines != None:
        for line in lines:
            if not line or line[0] == "#":
                continue
            strs = line.split()
            if len(strs) >= 12:
                pfamid = strs[3]
                if not pfamid in dt:
                    dt[pfamid] = {}
                    dt[pfamid]['pfamdef'] = strs[4]
                    dt[pfamid]['numSeqCls1'] = int(strs[5])
                    dt[pfamid]['numSeqCls2'] = int(strs[6])
                    dt[pfamid]['numseq'] = int(strs[7])
                    dt[pfamid]['nTM_Group1'] = int(strs[8])
                    dt[pfamid]['nTM_Group2'] = int(strs[9])
                    dt[pfamid]['pairlist'] = []
                seqid1 = strs[0]
                seqid2 = strs[1]
                seqidt = float(strs[2])
                isSP = int(strs[10])
                isPDB = int(strs[11])
                dt[pfamid]['pairlist'].append((seqid1, seqid2, seqidt, isSP,
                    isPDB))
        lines = hdl.readlines()
    hdl.close()
    return dt
#}}}



def WriteHTML(pairInfoList, datapath, outpath):#{{{
    """
    Write HTML
    """
    htmlname = g_params['htmlname']
# The name of the main html file is called index.html if no name is give
    htmlfilename = outpath + os.sep + htmlname + ".html"
    figuredir = outpath + os.sep + htmlname
    local_javalibpath = os.environ['HOME'] + os.sep + 'libjavascript'
    jsfileList = []
    jsfileList.append(local_javalibpath + os.sep + 'sorttable.js')
# copy jsfile to outpath 
    for f in jsfileList:
        if os.path.exists(f):
            os.system("cp -f %s %s"%(f, outpath))

    try:
        os.system("mkdir -p %s"%(figuredir)) 
        fpout = open(htmlfilename,"w")
        title="Pairs with different topologies but within the same protein family"
        WriteHTMLHeader(title, fpout)

        print >> fpout, "<dir id=\"Content\">"
        tablename = 'table1'
        tabletitle = ""
        WriteHTMLTable(tablename, tabletitle, pairInfoList, datapath, htmlname,
                outpath, fpout)
        print >> fpout, "</dir>"
        WriteHTMLTail(fpout)
        fpout.close()
        print  "Result has been output to %s"%(htmlfilename)
    except IOError:
        raise

#}}}
def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = ""
    datapath = ""
    pairListFile = ""
    seqlenFile = ""
    shortid2fullidFile = ""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            isNonOptionArg = False
            i += 1
            return 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-outpath", "--outpath"]:
                (outpath, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-datapath", "--datapath"]:
                (datapath, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-treepath", "--treepath"]:
                (g_params['treepath'], i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-maxperfamily", "--maxperfamily"]:
                (g_params['max_num_output_per_family'], i) = myfunc.my_getopt_int(argv, i)
            elif argv[i] in ["-swissprot", "--swissprot"]:
                (g_params['isSP_threshold'], i) = myfunc.my_getopt_int(argv, i)
            elif argv[i] in ["-msapath", "--msapath"]:
                (g_params['msapath'], i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-seqlen", "--seqlen"]:
                (seqlenFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-shortid2fullid", "--shortid2fullid"]:
                (shortid2fullidFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-l", "--l"] :
                (pairListFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True; i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            print >> sys.stderr, "Error! Wrong argument:", argv[i]
            return 1


    if pairListFile != "":
        pairInfoDict = ReadPairInfo(pairListFile)


    if len(pairInfoDict) < 1:
        print >> sys.stderr, "no pair info read in. exit"
        return 1

    if datapath == "":
        print >> sys.stderr, "datapath not set"
        return 1
    elif not os.path.exists(datapath):
        print >> sys.stderr, "datapath %s does not exist"%(datapath)
        return 1

    if outpath == "":
        print >> sys.stderr, "outpath not set"
        return 1
    elif not os.path.exists(outpath):
        cmd = ["mkdir", "-p", outpath]
        subprocess.check_call(cmd)


    g_params['OS'] = os.uname()[0];
    if g_params['OS'].find('Linux') != -1:
        g_params['CP_EXE'] = "/bin/cp -uf"
    else:
        g_params['CP_EXE'] = "/bin/cp -f"


    # read seqlen file
    if seqlenFile != "":
        g_params['seqLenDict'] = myfunc.ReadSeqLengthDict(seqlenFile)

    if shortid2fullidFile != "":
        g_params['uniprotAC2FullSeqIDMap'] = myfunc.ReadID2IDMap(shortid2fullidFile)
    #Sort by descending order of numseq
    pairInfoList = sorted(pairInfoDict.items(), key=lambda x:x[1]['numseq'],
            reverse=True)

    #pairInfoList ['pfamid', {'numseq': 123, 'pairlist':[]}]
    WriteHTML(pairInfoList, datapath, outpath)


#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['htmlname'] = "index"
    g_params['treepath'] = ""
    g_params['msapath'] = ""
    g_params['seqLenDict'] = {}
    g_params['uniprotAC2FullSeqIDMap']  = {}
    g_params['max_num_output_per_family'] = 10 #set the maximum number of pairs to output for each family
    g_params['isSP_threshold'] = 0 # set the threshold for isSP, when the threshold ==3, it means both sequences should be from swissprot
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

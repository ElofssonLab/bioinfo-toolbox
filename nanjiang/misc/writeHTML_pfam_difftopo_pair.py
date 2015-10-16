#!/usr/bin/python
# Description:
import os
import sys
import myfunc
import subprocess
import tempfile
import libtopologycmp as lcmp

progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))
rundir = os.path.dirname(sys.argv[0])

usage_short="""
Usage: %s -datapath DIR -outpath DIR -basename STR
"""%(progname)

usage_ext="""
Description:
    Write HTML webpage for different topology pairs
    datapath are e.g. topcons_single_agreement_44/hhalign/rstpair_topcons_cov0
    basename must be set, and data files are named as
    ${bname}.paircmp                                paircmp file
    ${bname}.pairaln                                Sequence pairwise alignment file
    ${bname}.topoaln.fa                             Topology pairwise alignment file
    ${bname}_.FULL_ALIGNED.TM2SP.pairinfo.txt       pairinfo file


OPTIONS:
  -description STR  Set description text
  -pfamseqpath DIR  Set path for sequence files of each Pfam family
  -pdb2sp     FILE  Set file pdb2sp table
  -topodb      STR  Database for the predicted topology
  -seqdb       STR  Database of sequence
  -seq2pfam FILE    Seqid to Pfam id map file
  -pfam2seq FILE    Pfamid to seqid map file
  -seqmsapath FILE  Sequence msa path, must be set
  -pfamdef FILE     pfam definition file, (default:
                    /data3/data/pfam/pfam27.0/Pfam-A.clans.tsv)
  -shortid2fullid   FILE 
                    Mapfile from uniprot ac to full uniprot id

OPTIONAL arguments
  -msapath  DIR      Set path for msa files
  -treepath DIR      Set path for tree files
  -pairalnpath DIR   Set path for tree files
  -alignrange STR    Select alignment with different alignment ranges
                     all, full, part, (default: full)
  -min-seqidt FLOAT  (default: 0)
  -max-seqidt FLOAT  (default: 100)
                      Select only pairs with global sequence identity within
                      [minSeqIDT, maxSeqIDT]

  -q                Quiet mode
  -h, --help        Print this help message and exit


Created 2014-09-29, updated 2014-09-29, Nanjiang Shu 
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
    print >> fpout, "   font-family:Arial; font-size: 11px;"
    print >> fpout, "   vertical-align: top;"
    print >> fpout, "   text-align: center;"
    print >> fpout, "   word-wrap: break-all;"
    print >> fpout, "   border: 1px solid black;"
    print >> fpout, "   padding: 5px 5px 5px 5px;"
    print >> fpout, "}"
    print >> fpout, ".wrapword td th"
    print >> fpout, "{"
    print >> fpout, "   font-family:Arial; font-size: 11px;"
    print >> fpout, "   vertical-align: top;"
    print >> fpout, "   text-align: center;"
    print >> fpout, "   width: 5px;"
    print >> fpout, "   word-wrap: break-word;"
    print >> fpout, "   border: 0px solid black;"
    print >> fpout, "   padding: 1px 1px 1px 1px;"
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

def WriteSubTable(dataTable, pfamid,  outpath, fpout):#{{{
    uniprotAC2FullSeqIDMap = g_params['uniprotAC2FullSeqIDMap']

    numDiffPair = 0
    numPair = 0
    for cmpclass in g_params['cmpClassList_mp3_cmpdup'][0:]:
        try:
            numPair_thisclass =  len(dataTable[pfamid]['difftopopair'][cmpclass])
        except KeyError:
            numPair_thisclass = 0
        numPair += numPair_thisclass
        if cmpclass != "IDT":
            numDiffPair += numPair_thisclass

    print >> fpout, "<table class=\"sortable\" border=1>"

    targetpath = outpath + os.sep + "data"
    if not os.path.exists(targetpath):
        os.system("mkdir -p %s"%(targetpath))

    cntOutputID = 0

    headerItemList=[]
    headerItemList.append("No.")
    headerItemList.append("CmpClass")
    headerItemList.append("ID1")
    headerItemList.append("ID2")
    headerItemList.append("SeqIDT")
    headerItemList.append("nTM1")
    headerItemList.append("nTM2")
    headerItemList.append("Len1")
    headerItemList.append("Len2")
    headerItemList.append("Sequence alignment")
    headerItemList.append("Figure pairwise alignment")

# pairlist [ seqid1, seqid2, seqidt, isSP, isPDB]
    print >> fpout, "<tr>"
    for item in headerItemList:
        if item in ["ID1", "ID2"]:
            th_class = "th class=\"wrapword\""
        else:
            th_class = "th"
        print >> fpout, "<%s>"%(th_class)
        print >> fpout, item
        print >> fpout, "</th>"
    print >> fpout, "</tr>"

    for cmpclass in g_params['cmpClassList_mp3_cmpdup'][1:]:
        try:
            pairList = dataTable[pfamid]['difftopopair'][cmpclass]
        except KeyError:
            pairList = []

        if len(pairList) <= 0:
            continue

        pairList = sorted(pairList, key=lambda x:x[8], reverse=True) # sorted by sequence identity

        for tup in pairList:
            seqid1 = tup[0]
            seqid2 = tup[1]
            nTM1 = tup[4]
            nTM2 = tup[5]
            seqLen1 = tup[6]
            seqLen2 = tup[7]
            seqidt = tup[8]
            cntOutputID += 1
            print >> fpout, "<tr>"
# 1. No ---------------------------
            print >> fpout, '<td>'
            print >> fpout, '%d'%(cntOutputID)
            print >> fpout, '</td>'
# 2. Cmpclass ---------------------------
            print >> fpout, '<td>'
            print >> fpout, '%s'%(cmpclass)
            print >> fpout, '</td>'
# 3. ID1 ---------------------------
            try:
                full_seqid1 = uniprotAC2FullSeqIDMap[seqid1]
                strs = full_seqid1.split("|")
                full_seqid1 = "<br>".join(strs[1:])
            except KeyError:
                full_seqid1 = seqid1

            try:
                pdblist = g_params['uniprot2pdbMap'][seqid1]
            except KeyError:
                pdblist = []
            print >> fpout, '<td class=\"wrapword\">'
            print >> fpout, '<a href=\"http://www.uniprot.org/uniprot/%s\"target=\"_blank\">%s</a>'%(seqid1, full_seqid1)
            if len(pdblist) > 0:
                print >> fpout, "<br><br>PDB: "
                for pdbid in pdblist: 
                    print >> fpout, "<a href=\"http://www.rcsb.org/pdb/explore.do?structureId=%s\"target=\"_blank\">%s </a>"%(pdbid, pdbid)
            print >> fpout, '</td>'
# 4. ID2 ---------------------------
            try:
                full_seqid2 = uniprotAC2FullSeqIDMap[seqid2]
                strs = full_seqid2.split("|")
                full_seqid2 = "<br>".join(strs[1:])
            except KeyError:
                full_seqid2 = seqid2
            try:
                pdblist = g_params['uniprot2pdbMap'][seqid2]
            except KeyError:
                pdblist = []
            print >> fpout, '<td class=\"wrapword\">'
            print >> fpout, '<a href=\"http://www.uniprot.org/uniprot/%s\"target=\"_blank\">%s</a>'%(seqid2, full_seqid2)
            if len(pdblist) > 0:
                print >> fpout, "<br><br>PDB: "
                for pdbid in pdblist: 
                    print >> fpout, "<a href=\"http://www.rcsb.org/pdb/explore.do?structureId=%s\"target=\"_blank\">%s </a>"%(pdbid, pdbid)
            print >> fpout, '</td>'
            print >> fpout, '</td>'
# 5. SeqIDT ---------------------------
            print >> fpout, '<td>'
            print >> fpout, '%.1f'%(seqidt)
            print >> fpout, '</td>'
# 6. nTM1 ---------------------------
            print >> fpout, '<td>'
            print >> fpout, '%d'%(nTM1)
            print >> fpout, '</td>'
# 7. nTM1 ---------------------------
            print >> fpout, '<td>'
            print >> fpout, '%d'%(nTM2)
            print >> fpout, '</td>'
# 8. Len1 ---------------------------
            print >> fpout, '<td>'
            print >> fpout, '%d'%(seqLen1)
            print >> fpout, '</td>'
# 9. Len2 ---------------------------
            print >> fpout, '<td>'
            print >> fpout, '%d'%(seqLen2)
            print >> fpout, '</td>'
# 10 sequence alignment  ---------------------------
            sourcefile = (g_params['pairalnpath'] + os.sep + "%s_%s"%(seqid1,
                seqid2) + ".topoaln.krbias.html")
            if not os.path.exists(sourcefile):
                PrepareDataForPairaln(seqid1, seqid2,  g_params['pairalnpath'])

            if os.path.exists(sourcefile):
                cmd = ["/bin/cp", "-f", sourcefile, targetpath]
                try:
                    subprocess.check_call(cmd)
                except subprocess.CalledProcessError, e:
                    print e
            print >> fpout, '<td>'
            print >> fpout, ("<a href=\"%s\" target=\"_blank\">"
                    % ("data" + os.sep +
                        os.path.basename(sourcefile)))
            #print >> fpout, os.path.basename(sourcefile)
            print >> fpout, "seq_align"
            print >> fpout, "</a>"
            print >> fpout, '</td>'

# 11 Figure pairwise alignment -----------------------
            sourcefile = (g_params['pairalnpath'] + os.sep + "%s_%s"%(seqid1,
                seqid2) + ".topoaln.krbias-crop.png")
            if os.path.exists(sourcefile):
                cmd = ["/bin/cp", "-f", sourcefile, targetpath]
                try:
                    subprocess.check_call(cmd)
                except subprocess.CalledProcessError, e:
                    print e

            sourcefile_thumb = (g_params['pairalnpath'] + os.sep + "%s_%s"%(seqid1,
                seqid2) + ".topoaln.krbias-crop.thumb.png")
            if os.path.exists(sourcefile_thumb):
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

    print >> fpout, "</table>"
#}}}
def PrepareDataForTopoanaTMPro(pfamid, seqwithtopo_idlist, outpath):#{{{
    seqwithtopo_set = set(seqwithtopo_idlist)
    # 1 pfamid.fa
    seqfile = "%s/%s.fa"%(outpath, pfamid)
    fpout = open(seqfile, "w")
    for seqid in seqwithtopo_set:
        record = g_params['hdl_seqdb'].GetRecord(seqid)
        if record:
            fpout.write(record)
    fpout.close()

    # 2. pfamid.topo
    topofile = "%s/%s.topo"%(outpath, pfamid)
    fpout = open(topofile, "w")
    for seqid in seqwithtopo_set:
        record = g_params['hdl_topodb'].GetRecord(seqid)
        if record:
            fpout.write(record)
    fpout.close()


    # 3. get msafile
    source_msafile = "%s/%s.clustalo10.mfa"%(g_params['seqmsapath'], pfamid)
    msafile = "%s/%s.clustalo10.mfa"%(outpath, pfamid)
    if os.path.exists(source_msafile):
        cmd = ["/bin/cp", "-f", source_msafile, msafile]
        try:
            subprocess.check_call(cmd)
        except subprocess.CalledProcessError, e:
            print e


    # 4. produce figures
    cmd = ["%s/topoana_TMpro.sh"%(rundir), "-anamode", "2", "-outpath", g_params['msapath'], "-datapath", g_params['msapath'], "-maxseq", "10000", "-extfa", ".fa" ,"-exttopo", ".topo", "-extmsa", ".clustalo10.mfa", pfamid]
    if g_params['isDEBUG']:
        print " ".join(cmd)
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError, e:
        print e
#}}}
def PrepareDataForPairaln(seqid1, seqid2, outpath):#{{{

    topoalnfile = "%s/%s.topoaln.fa"%(g_params['datapath'], g_params['basename'])
    if myfunc.checkfile(topoalnfile, "topoalnfile") != 0:
        return 1
    # 1 seqid1-seqid2.topoaln.fa
    ext_topoaln = ".topoaln.fa"
    pair_topoalnfile = "%s/%s_%s.topoaln.fa"%(outpath, seqid1, seqid2)
    cmd = ["%s/selectPairaln.py"%(rundir), "-pairaln", topoalnfile, 
            "-outpath", outpath, "-ext", ext_topoaln, "-split", seqid1, seqid2]
    try:
        subprocess.check_output(cmd)
    except subprocess.CalledProcessError, e:
        print e

    if myfunc.checkfile(pair_topoalnfile, "pair_topoalnfile") != 0:
        return 1

    # 2. seqid1_seqid2.fa
    seqfile = "%s/%s_%s.fa"%(outpath, seqid1,seqid2)
    fpout = open(seqfile, "w")
    for seqid in [seqid1,seqid2]:
        record = g_params['hdl_seqdb'].GetRecord(seqid)
        if record:
            fpout.write(record)
    fpout.close()

    if myfunc.checkfile(seqfile, "seqfile") != 0:
        return 1
    # 3 draw pairaln
    method_shrink = "2"
    method_plot = "mat"
    shrinkrateTM = "3"
    maxHoldLoop = "4"
    cmd = ["%s/drawMSATopo.py"%(rundir), pair_topoalnfile,  "-aaseq",
            seqfile, "-pdg", "yes", "-shrink", "yes", "-m-shrink",
            method_shrink, "-method", method_plot, "-text", "yes",
            "-krbias", "-shrinkrateTM", shrinkrateTM, "-max-hold-loop",
            maxHoldLoop]
    try:
        subprocess.check_output(cmd)
    except subprocess.CalledProcessError, e:
        print "[Command failed]", " ".join(cmd)
        print e
#}}}

def WriteHTMLTable(tablename, tabletitle, dataTable,  htmlname, #{{{
        outpath, fpout):
    numInputID = len (dataTable)
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
    headerItemList.append("NumSeq<br>All")
    headerItemList.append("NumSeq<br>TPS_44")
    headerItemList.append("NumSeq<br>Used")
    headerItemList.append("Figure MSA")
    headerItemList.append("Figure Tree")
    headerItemList.append("Pair")
    headerItemList.append("INV(#,%)")
    headerItemList.append("Dup")
    headerItemList.append("TM2GAP")
    headerItemList.append("Mixed")
    headerItemList.append("TM2SEQ")
    headerItemList.append("TM2SP")
    headerItemList.append("Pairwise alignment with different topology")

    print >> fpout, "<tr>"
    for item in headerItemList:
        print >> fpout, "<th>"
        print >> fpout, item
        print >> fpout, "</th>"
    print >> fpout, "</tr>"

    inputIDList = dataTable.keys()

    for pfamid in inputIDList:
        numDiffPair = 0
        numPair = 0
        for cmpclass in g_params['cmpClassList_mp3_cmpdup'][0:]:
            try:
                numPair_thisclass =  len(dataTable[pfamid]['difftopopair'][cmpclass])
            except KeyError:
                numPair_thisclass = 0
            numPair += numPair_thisclass
            if cmpclass != "IDT":
                numDiffPair += numPair_thisclass

        if numDiffPair <= 0:
            continue

        info = dataTable[pfamid]

        try:
            numSeqAll = len(g_params['pfamid2seqidDict'][pfamid])
        except KeyError:
            numSeqAll = 0

        
        try:
            seqwithtopo_idlist = list(set(g_params['pfamid2seqidDict'][pfamid]) & set(g_params['hdl_topodb'].indexedIDList))
        except KeyError:
            seqwithtopo_idlist = []
        numSeqWithTopology = len(seqwithtopo_idlist)




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
        try:
            pfamdef = g_params['pfamidDefDict'][pfamid]
        except KeyError:
            pfamdef = ""
        print >> fpout, '<td>'
        print >> fpout, '%s'%(pfamdef)
        print >> fpout, '</td>'
# 4 NumSeqAll ---------------------------
        print >> fpout, '<td>'
        print >> fpout, '%d'%(numSeqAll)
        print >> fpout, '</td>'
# 5 NumSeq with topology ---------------------------
        print >> fpout, '<td>'
        print >> fpout, '%d'%(numSeqWithTopology)
        print >> fpout, '</td>'
# 6 NumSeqUsed ---------------------------
        try:
            numSeqUsed = len(dataTable[pfamid]['set_seqid'])
        except KeyError:
            numSeqUsed = 0
        print >> fpout, '<td>'
        print >> fpout, '%d'%(numSeqUsed)
        print >> fpout, '</td>'

# 7 Figure MSA---------------------------
        if 1:
            ext = '.reordered.topomsa.png'
            print >> fpout, '<td>'
            imageSourceFile = g_params['msapath'] + os.sep + pfamid + ext

            if not os.path.exists(imageSourceFile):
                PrepareDataForTopoanaTMPro(pfamid, seqwithtopo_idlist, g_params['msapath'])

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

# 8 Figure Tree---------------------------
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

# 9 numPair ---------------------------
        print >> fpout, '<td>'
        print >> fpout, '%d'%(numPair)
        print >> fpout, '</td>'
# 10-15 INV, DUP, TM2GAP, Mixed, TM2SEQ, TM2SP ---------------------------
        for cmpclass in g_params['cmpClassList_mp3_cmpdup'][1:]:
            try:
                nn = len(dataTable[pfamid]['difftopopair'][cmpclass])
            except KeyError:
                nn = 0
            print >> fpout, '<td>'
            if nn > 0:
                ss1 = "%d, %.1f%%"%(nn, myfunc.FloatDivision(nn, numPair))
            else:
                ss1 = "-"
            print >> fpout, '%s'%(ss1)
            print >> fpout, '</td>'

# 11 #pairwise alignment ---------------------------
        print >> fpout, '<td>'
        WriteSubTable(dataTable, pfamid, outpath, fpout)
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



def WriteHTML(dataTable, outpath):#{{{
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

        print >> fpout, "<p><h2>%s</h2></p>"%(g_params['description'])
        print >> fpout, "<dir id=\"Content\">"
        tablename = 'table1'
        tabletitle = ""
        WriteHTMLTable(tablename, tabletitle, dataTable,  htmlname,
                outpath,  fpout)
        print >> fpout, "</dir>"
        WriteHTMLTail(fpout)
        fpout.close()
        print  "Result has been output to %s"%(htmlfilename)
    except IOError:
        raise

#}}}
def FilterPairCmpResult(recordList):#{{{
    """
    Filter paircmp result by g_params
    return newList
    """
    newList = []
    pairListSet = set([])
    numInputRecord = len(recordList)
    for record in recordList:
        if record == {}:
            continue
        key = "%s-%s"%(record['id1'], record['id2'])
        cmpclass = record['cmpclass']
        if (cmpclass.find('UNALIGNED') == 0 
                or cmpclass.find('AMBIGUOUS') == 0):
            continue
        if record['alignrange'] != g_params['alignrange']:
            continue
        newList.append(record)

    numOutputRecord = len(newList)
    if g_params['isDEBUG']:
        if numOutputRecord < numInputRecord:
            print "%d pairs dropped" % (numInputRecord-numOutputRecord)
    return newList
#}}}
def AddAllSeqInPairCmp(dataTable, pairCmpRecordList, seqid2pfamidDict):#{{{
    """
    Add all proteins used in paircmp to each pfam group
    """
    for record in pairCmpRecordList:
        id1 = record['id1']
        id2 = record['id2']
        try:
            pfamidlist1 = seqid2pfamidDict[id1]
            pfamidlist2 = seqid2pfamidDict[id2]
            common_pfamidlist = list(set(pfamidlist1) & set(pfamidlist2))
        except KeyError:
            common_pfamidlist = []
        for pfamid in common_pfamidlist:
            if not pfamid in dataTable:
                dataTable[pfamid] = {}
                dataTable[pfamid]['set_seqid'] = set([])
                dataTable[pfamid]['difftopopair'] = {}
            dataTable[pfamid]['set_seqid'].add(id1)
            dataTable[pfamid]['set_seqid'].add(id2)
#}}}

def ReadPairInfo_cmpclass(infile):#{{{
    """
    Read files e.g.  Pfam2.mp3.cmpdup.max10000.maxpair100.hhalign_.cmpdup.FULL_ALIGNED.DUP.pairinfo.txt
    """
    try:
        fpin = open(infile, "r")
        lines = fpin.readlines()
        fpin.close()
        li = []
        for line in lines:
            strs = line.split()
            if len(strs) >= 9:
#                            id1    id2    NtermState1 NtermState2  nTM1  nTM2
                tup =  (strs[0],strs[1],strs[2],strs[3], int(strs[4]), int(strs[5]), 
#                       len1       len2          seqidt1         pfamidlist
                    int(strs[6]), int(strs[7]), float(strs[8]), strs[9:])
                #print >> sys.stderr, "line=", line, "\n", "tup=", tup
                li.append(tup)
        return li
    except IOError:
        print >> sys.stderr, "Failed to read infile %s"%(infile)
        return []
#}}}
def AddPairInfo(dataTable, pairinfoList, cmpclass):#{{{
    """
    Add pairinfoList to 'difftopopair':{'INV':[(id1,id2)]},'TM2GAP':{},{}}
    pairinfoList = [(id1, id2, NtermState1, NtermState2, nTM1, nTM2, len1, len2, seqidt1), ()]
    """
    for tup in pairinfoList:
        pfamidlist = tup[9]
        seqidt = tup[8]
        if not (seqidt >= g_params['minSeqIDT'] and seqidt < g_params['maxSeqIDT']):
            continue
        #print >> sys.stderr, tup, pfamidlist
        for pfamid in pfamidlist:
            if cmpclass not in dataTable[pfamid]['difftopopair']:
                dataTable[pfamid]['difftopopair'][cmpclass] = []
            dataTable[pfamid]['difftopopair'][cmpclass].append(tup[0:9])


#}}}
def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = ""
    pairListFile = ""
    seqlenFile = ""
    shortid2fullidFile = ""
    seqid2pfamidMapFile = ""
    pfamDefFile = '/data3/data/pfam/pfam27.0/Pfam-A.clans.tsv'
    topodb = ""
    seqdb = ""
    pdb2spFile = ""

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
            elif argv[i] in ["-topodb", "--topodb"]:
                (topodb, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-pdb2sp", "-pdb2sp", "-pdbtosp","--pdbtosp"]:
                (pdb2spFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-seqdb", "--seqdb"]:
                (seqdb, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-seqmsapath", "--seqmsapath"]:
                (g_params['seqmsapath'], i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-datapath", "--datapath"]:
                (g_params['datapath'], i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-seq2pfam", "--seq2pfam"]:
                (seqid2pfamidMapFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-pfam2seq", "--pfam2seq"]:
                (pfamid2seqidMapFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-description", "--description"]:
                (g_params['description'], i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-pfamdef", "--pfamdef"]:
                (pfamDefFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-alignrange", "--alignrange"]:
                g_params['alignrange'],i  = myfunc.my_getopt_str(argv,i)
                if not g_params['alignrange'] in ['all', 'full', 'part']:
                    print >> sys.stderr, "alignrange must be one of [all, full, part]"
                    return 1
                else:
                    if g_params['alignrange'] == 'full':
                        g_params['alignrange'] = 'FULL_ALIGNED'
                    elif g_params['alignrange'] == 'part':
                        g_params['alignrange'] = 'PART_ALIGNED'
            elif argv[i] in ["-basename", "--basename"]:
                (g_params['basename'], i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-treepath", "--treepath"]:
                (g_params['treepath'], i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-pairalnpath", "--pairalnpath"]:
                (g_params['pairalnpath'], i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-maxperfamily", "--maxperfamily"]:
                (g_params['max_num_output_per_family'], i) = myfunc.my_getopt_int(argv, i)
            elif argv[i] in ["-min-seqidt", "--min-seqidt"]:
                g_params['minSeqIDT'], i  = myfunc.my_getopt_float(argv, i)
            elif argv[i] in ["-max-seqidt", "--max-seqidt"]:
                g_params['maxSeqIDT'], i  = myfunc.my_getopt_float(argv, i)
            elif argv[i] in ["-shortid2fullid", "--shortid2fullid"]:
                (shortid2fullidFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-debug", "--debug"]:
                if argv[i+1][0].lower() == 'y':
                    g_params['isDEBUG'] = True
                else:
                    g_params['isDEBUG'] = False
                i += 2
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True; i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            print >> sys.stderr, "Error! Wrong argument:", argv[i]
            return 1


    if g_params['basename'] == "":
        print >> sys.stderr,"basename not set. exit"
        return 1
    if myfunc.checkfile(g_params['datapath'], "datapath") != 0:
        return 1
    if myfunc.checkfile(seqid2pfamidMapFile, "seqid2pfamidMapFile") != 0:
        return 1
    if myfunc.checkfile(pfamid2seqidMapFile, "pfamid2seqidMapFile") != 0:
        return 1

    if myfunc.checkfile(topodb+"0.db", "topodb") != 0:
        return 1
    if myfunc.checkfile(seqdb+"0.db", "seqdb") != 0:
        return 1
    if myfunc.checkfile(g_params['seqmsapath'], "seqmsapath") != 0:
        return 1

    if pdb2spFile != "":
        (g_params['pdb2uniprotMap'], g_params['uniprot2pdbMap']) = myfunc.ReadPDBTOSP(pdb2spFile)

    if g_params['datapath'] == "":
        print >> sys.stderr, "datapath not set"
        return 1
    elif not os.path.exists(g_params['datapath']):
        print >> sys.stderr, "datapath %s does not exist"%(g_params['datapath'])
        return 1

    if outpath == "":
        print >> sys.stderr, "outpath not set"
        return 1
    elif not os.path.exists(outpath):
        cmd = ["mkdir", "-p", outpath]
        subprocess.check_call(cmd)


    paircmpfile = "%s/%s.paircmp"%(g_params['datapath'], g_params['basename'])
    if myfunc.checkfile(paircmpfile, "paircmpfile") != 0:
        return 1

    (g_params['pfamidDefDict'], g_params['clanidDefDict']) = lcmp.ReadPfamDefFile(pfamDefFile)

    g_params['seqid2pfamidDict'] = myfunc.ReadFam2SeqidMap(seqid2pfamidMapFile)
    g_params['pfamid2seqidDict'] = myfunc.ReadFam2SeqidMap(pfamid2seqidMapFile)

    tmpdir = tempfile.mkdtemp()
    if g_params['msapath'] == "":
        g_params['msapath'] = tmpdir
    if g_params['treepath'] == "":
        g_params['treepath'] = tmpdir
    if g_params['pairalnpath'] == "":
        g_params['pairalnpath'] = tmpdir


    pairCmpRecordList = []
    unprocessedBuffer=""
    cntTotalReadInRecord = 0
    cntTotalOutputRecord = 0
    isEOFreached = False
    try:
        fpin = open(paircmpfile, "r")
    except IOError:
        print >> sys.stderr, "Failed to open input file %s"%(paircmpfile)
        return 1
    while 1:
        buff = fpin.read(myfunc.BLOCK_SIZE)
        if buff == "":
            isEOFreached = True
        buff = unprocessedBuffer + buff
        rdList=[]
        unprocessedBuffer = lcmp.ReadPairCmpResultFromBuffer(buff,rdList)
        rdList = FilterPairCmpResult(rdList)
        cntTotalReadInRecord += len(rdList)
        pairCmpRecordList += rdList
        if isEOFreached == True:
            break
    fpin.close()

    print "cntTotalReadInRecord =", cntTotalReadInRecord

    g_params['hdl_seqdb'] = myfunc.MyDB(seqdb)
    g_params['hdl_topodb'] = myfunc.MyDB(topodb)

    g_params['OS'] = os.uname()[0];
    if g_params['OS'].find('Linux') != -1:
        g_params['CP_EXE'] = "/bin/cp -uf"
    else:
        g_params['CP_EXE'] = "/bin/cp -f"


    if shortid2fullidFile != "":
        g_params['uniprotAC2FullSeqIDMap'] = myfunc.ReadID2IDMap(shortid2fullidFile)

    addname = ""
    if g_params['alignrange'] != 'all':
        addname += ".%s"%(g_params['alignrange'])


    dataTable = {}
# structure of dataTable
#    dataTable[pfamid] = {'set_seqid':set(), 'difftopopair':[{'INV':[(id1,id2)]},{'TM2GAP':},{}}

# first read in pairCmpRecordList
    AddAllSeqInPairCmp(dataTable, pairCmpRecordList, g_params['seqid2pfamidDict'])

    pairInfoFileList = []
    for cmpclass in g_params['cmpClassList_mp3_cmpdup'][0:]:
        ss = "%s/%s_.cmpdup.FULL_ALIGNED.%s.pairinfo.txt" %(g_params['datapath'], g_params['basename'], cmpclass)
        pairInfoFileList.append(ss)
        pairinfoList = ReadPairInfo_cmpclass(ss)
        AddPairInfo(dataTable, pairinfoList, cmpclass)
#     print "\n".join(pairInfoFileList)
    if g_params['isDEBUG']:#{{{
        for pfamid in dataTable:
            print pfamid
            print "\tset_seqid"
            print dataTable[pfamid]['set_seqid']
            print "\tdifftopopair"
            for cls in dataTable[pfamid]['difftopopair']:
                print "\t\t",cls
                for tup in  dataTable[pfamid]['difftopopair'][cls]:
                    print "\t\t\t", tup#}}}

    WriteHTML(dataTable,  outpath)

    os.system("rm -rf %s"%(tmpdir))

#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['minSeqIDT'] = 0.0
    g_params['maxSeqIDT'] = 100.0
    g_params['pdb2uniprotMap'] = {}
    g_params['uniprot2pdbMap'] = {}
    g_params['isDEBUG'] = False
    g_params['htmlname'] = "index"
    g_params['treepath'] = ""
    g_params['msapath'] = ""
    g_params['datapath'] = ""
    g_params['basename'] = ""
    g_params['pairalnpath'] = ""
    g_params['description'] = ""
    g_params['seqmsapath'] = ""
    g_params['seqLenDict'] = {}
    g_params['pfamid2seqidDict'] = {}
    g_params['seqid2pfamidDict'] = {}
    g_params['pfamidDefDict'] = {}
    g_params['clanidDefDict'] = {}
    g_params['hdl_seqdb'] = None
    g_params['hdl_topodb'] = None
    g_params['uniprotAC2FullSeqIDMap']  = {}
    g_params['max_num_output_per_family'] = 10 #set the maximum number of pairs to output for each family
    g_params['cmpClassList_mp3_cmpdup'] = ["IDT","INV","DUP", "TM2GAP", "TM2SEQ", "TM2SP", "Mixed"]
    g_params['isSP_threshold'] = 0 # set the threshold for isSP, when the threshold ==3, it means both sequences should be from swissprot
    g_params['alignrange'] = 'full'
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

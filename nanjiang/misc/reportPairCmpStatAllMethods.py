#!/usr/bin/env python
# create html file for figtree
# e.g. for figtree_kalignp_nr95

# ChangeLog 2012-03-15 
#  create sortable table using sorttable.js
# ChangeLog 2012-03-20
#  Add reordered topomsa (in the same order as phylo tree)

import os
import sys
import libtopologycmp as lcmp

usage="""
Usage:  reportPairCmpStatAllMethods.py

Description: Report the result of Pfam tree

Options:
  -datapath DIR      Set path for tree images
  -outpath DIR       Set ouput path, (default: ./)
  -htmlname STR      Set the name of the main html file, (default: index)
  -q                 Quiet mode
  -h, --help         Print this help message and exit

Created 2012-05-04, updated 2012-05-04, Nanjiang Shu  

Examples:
"""

BLOCK_SIZE = myfunc.BLOCK_SIZE
# typeTopoAnaDict
# 0: started from $id.fa with multiple homologous sequences
# 1: started from $id.fa with a single sequence, and the homologous sequences
# were obtained by blast searching

def PrintHelp():
    print usage


def ReadPfamDEList(infile):#{{{
    pfamDefDict={}
    try:
        fpin=open(infile)
        lines=fpin.readlines()
        fpin.close()
        for line in lines:
            strs=line.split(':')
            pfamid=strs[0].strip()
            pfamDefDict[pfamid]=strs[1].strip()
        return pfamDefDict
    except IOError: 
        print >> sys.stderr, "%s: Error reading %s"%(sys.argv[0], infile)
        return 1
#}}}


def WriteHTMLHeader(title, fpout):#{{{
    print >> fpout, "<HTML>"
    print >> fpout, "<head>"
    print >> fpout, "<title>%s</title>"%(title)
    print >> fpout, "<script src=\"sorttable.js\"></script>"
    print >> fpout, "<style type=\"text/css\" media=\"all\"> @import \"layout.css\";"
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
def WriteNonAvailableHTMLCell(fpout):#{{{
    print >> fpout, "<td>%s</td>"%STR_NON_AVAILABLE
#}}}


def WriteSubtableHits(record, maxNumInternalCons, maxNumInternalQuery, fpout):#{{{
    print >> fpout, "<table class=hits>"
    WriteTableHeaderHits(maxNumInternalCons, maxNumInternalQuery, fpout)
    WriteTableContentHits(record, 1, maxNumInternalCons,
            maxNumInternalQuery,fpout)
    print >> fpout, "</table>"
#}}}

def WriteHTMLTable(tablename, tabletitle, method_list, datapath, #{{{
        htmlname, outpath, fpout):

    print >> fpout, "<a name=\"%s\"></a><h4>%s</h4>"%(tablename,tabletitle)
    print >> fpout, "<table class=\"sortable\" border=1>"
    cntOutputID = 0

    headerItemList=[]
    headerItemList.append("Methods")
    headerItemList.append("cmpclass")
    headerItemList.append("NCtermInter")
    headerItemList.append("NumTMInter")
    print >> fpout, "<tr>"
    for item in headerItemList:
        print >> fpout, "<th>"
        print >> fpout, item
        print >> fpout, "</th>"
    print >> fpout, "</tr>"

    for i in xrange(len(method_list)):
        method = method_list[i]
        
        print >> fpout, "<tr>"
#---------------------------
        print >> fpout, '<td>'
        print >> fpout, '%s'%(method)
        print >> fpout, '</td>'
#---------------------------
        print >> fpout, '<td>'
        orig_image_file = (datapath + os.sep + "new_rstpair_" + method +
            os.sep +
            "dump.pairaln_type_all_dg10.0_gap0.0_ps0.0-cmpclass.table.txt.eps")

        imageSourceFile = orig_image_file.rstrip(".eps") + ".png"
        thumbImageSourceFile =  orig_image_file.rstrip(".eps") + ".thumb.png" 

        if not os.path.exists(imageSourceFile) or g_params['FORCE_OVERWRITE']:
            os.system("convert -density 150  %s %s"%(orig_image_file, imageSourceFile))
        if not os.path.exists(thumbImageSourceFile) or g_params['FORCE_OVERWRITE']:
            os.system("convert -thumbnail 200  %s %s"%(imageSourceFile, thumbImageSourceFile))
        imageTargetFile = (outpath + os.sep + htmlname + os.sep + method + "_" + 
            os.path.basename(imageSourceFile))
        thumbImageTargetFile = (outpath + os.sep + htmlname + os.sep +  method + "_" +
            os.path.basename(thumbImageSourceFile))
        os.system("%s %s %s"%(g_params['CP_EXE'], imageSourceFile,
            imageTargetFile))
        os.system("%s %s %s"%(g_params['CP_EXE'], thumbImageSourceFile,
            thumbImageTargetFile))
        print >> fpout, ("<a href=\"%s\"target=\"_blank\">"
                % (htmlname + os.sep + os.path.basename(imageTargetFile)))
        print >> fpout, ("<img src=\"%s\">" % (htmlname + os.sep +
            os.path.basename(thumbImageTargetFile)))
        print >> fpout, "</a>"
        print >> fpout, '</td>'
#---------------------------
        print >> fpout, '<td>'
        orig_image_file = (datapath + os.sep + "new_rstpair_" + method +
            os.sep +
            "dump.pairaln_type_all_dg10.0_gap0.0_ps0.0-NCtermInter.table2.c2.txt.eps")

        imageSourceFile = orig_image_file.rstrip(".eps") + ".png"
        thumbImageSourceFile =  orig_image_file.rstrip(".eps") + ".thumb.png" 

        if not os.path.exists(imageSourceFile) or g_params['FORCE_OVERWRITE']:
            os.system("convert -density 150  %s %s"%(orig_image_file, imageSourceFile))
        if not os.path.exists(thumbImageSourceFile) or g_params['FORCE_OVERWRITE']:
            os.system("convert -thumbnail 200  %s %s"%(imageSourceFile, thumbImageSourceFile))
        imageTargetFile = (outpath + os.sep + htmlname + os.sep +  method + "_" +
            os.path.basename(imageSourceFile))
        thumbImageTargetFile = (outpath + os.sep + htmlname + os.sep + method + "_" +
            os.path.basename(thumbImageSourceFile))
        os.system("%s %s %s"%(g_params['CP_EXE'], imageSourceFile,
            imageTargetFile))
        os.system("%s %s %s"%(g_params['CP_EXE'], thumbImageSourceFile,
            thumbImageTargetFile))
        print >> fpout, ("<a href=\"%s\"target=\"_blank\">"
                % (htmlname + os.sep + os.path.basename(imageTargetFile)))
        print >> fpout, ("<img src=\"%s\">" % (htmlname + os.sep +
            os.path.basename(thumbImageTargetFile)))
        print >> fpout, "</a>"
        print >> fpout, '</td>'
#---------------------------
        print >> fpout, '<td>'
        orig_image_file = (datapath + os.sep + "new_rstpair_" + method +
            os.sep +
            "dump.pairaln_type_all_dg10.0_gap0.0_ps0.0-NumTMDistri.seqidt30-100.txt.eps")

        imageSourceFile = orig_image_file.rstrip(".eps") + ".png"
        thumbImageSourceFile =  orig_image_file.rstrip(".eps") + ".thumb.png" 

        if not os.path.exists(imageSourceFile) or g_params['FORCE_OVERWRITE']:
            os.system("convert -density 150  %s %s"%(orig_image_file, imageSourceFile))
        if not os.path.exists(thumbImageSourceFile) or g_params['FORCE_OVERWRITE']:
            os.system("convert -thumbnail 200  %s %s"%(imageSourceFile, thumbImageSourceFile))
        imageTargetFile = (outpath + os.sep + htmlname + os.sep +  method + "_" +
            os.path.basename(imageSourceFile))
        thumbImageTargetFile = (outpath + os.sep + htmlname + os.sep + method + "_" +
            os.path.basename(thumbImageSourceFile))
        os.system("%s %s %s"%(g_params['CP_EXE'], imageSourceFile,
            imageTargetFile))
        os.system("%s %s %s"%(g_params['CP_EXE'], thumbImageSourceFile,
            thumbImageTargetFile))
        print >> fpout, ("<a href=\"%s\"target=\"_blank\">"
                % (htmlname + os.sep + os.path.basename(imageTargetFile)))
        print >> fpout, ("<img src=\"%s\">" % (htmlname + os.sep +
            os.path.basename(thumbImageTargetFile)))
        print >> fpout, "</a>"
        print >> fpout, '</td>'
#---------------------------
        print >> fpout, "</tr>"
    print >> fpout, "</table>"
#}}}
def WriteHTMLTable2(tablename, tabletitle, method_list, datapath, #{{{
        htmlname, outpath, fpout):

    print >> fpout, "<a name=\"%s\"></a><h4>%s</h4>"%(tablename,tabletitle)
    print >> fpout, "<table class=\"sortable\" border=1>"
    cntOutputID = 0

    headerItemList=[]
    headerItemList.append("Methods")
    headerItemList.append("cmpclass")
    headerItemList.append("NCtermInter")
    headerItemList.append("NumTMInter")
    print >> fpout, "<tr>"
    for item in headerItemList:
        print >> fpout, "<th>"
        print >> fpout, item
        print >> fpout, "</th>"
    print >> fpout, "</tr>"

    for i in xrange(len(method_list)):
        method = method_list[i]
        
        print >> fpout, "<tr>"
#---------------------------
        print >> fpout, '<td>'
        print >> fpout, '%s'%(method)
        print >> fpout, '</td>'
#---------------------------
        print >> fpout, '<td>'
        orig_image_file = (datapath + os.sep + "new_rstpair_" + method +
            os.sep +
            "dump.pairaln_type_all_dg0.0_gap0.8_ps0.0-cmpclass.table.txt.eps")

        imageSourceFile = orig_image_file.rstrip(".eps") + ".png"
        thumbImageSourceFile =  orig_image_file.rstrip(".eps") + ".thumb.png" 

        if not os.path.exists(imageSourceFile) or g_params['FORCE_OVERWRITE']:
            os.system("convert -density 150  %s %s"%(orig_image_file, imageSourceFile))
        if not os.path.exists(thumbImageSourceFile) or g_params['FORCE_OVERWRITE']:
            os.system("convert -thumbnail 200  %s %s"%(imageSourceFile, thumbImageSourceFile))
        imageTargetFile = (outpath + os.sep + htmlname + os.sep + method + "_" + 
            os.path.basename(imageSourceFile))
        thumbImageTargetFile = (outpath + os.sep + htmlname + os.sep +  method + "_" +
            os.path.basename(thumbImageSourceFile))
        os.system("%s %s %s"%(g_params['CP_EXE'], imageSourceFile,
            imageTargetFile))
        os.system("%s %s %s"%(g_params['CP_EXE'], thumbImageSourceFile,
            thumbImageTargetFile))
        print >> fpout, ("<a href=\"%s\"target=\"_blank\">"
                % (htmlname + os.sep + os.path.basename(imageTargetFile)))
        print >> fpout, ("<img src=\"%s\">" % (htmlname + os.sep +
            os.path.basename(thumbImageTargetFile)))
        print >> fpout, "</a>"
        print >> fpout, '</td>'
#---------------------------
        print >> fpout, '<td>'
        orig_image_file = (datapath + os.sep + "new_rstpair_" + method +
            os.sep +
            "dump.pairaln_type_all_dg0.0_gap0.8_ps0.0-NCtermInter.table2.c2.txt.eps")

        imageSourceFile = orig_image_file.rstrip(".eps") + ".png"
        thumbImageSourceFile =  orig_image_file.rstrip(".eps") + ".thumb.png" 

        if not os.path.exists(imageSourceFile) or g_params['FORCE_OVERWRITE']:
            os.system("convert -density 150  %s %s"%(orig_image_file, imageSourceFile))
        if not os.path.exists(thumbImageSourceFile) or g_params['FORCE_OVERWRITE']:
            os.system("convert -thumbnail 200  %s %s"%(imageSourceFile, thumbImageSourceFile))
        imageTargetFile = (outpath + os.sep + htmlname + os.sep +  method + "_" +
            os.path.basename(imageSourceFile))
        thumbImageTargetFile = (outpath + os.sep + htmlname + os.sep + method + "_" +
            os.path.basename(thumbImageSourceFile))
        os.system("%s %s %s"%(g_params['CP_EXE'], imageSourceFile,
            imageTargetFile))
        os.system("%s %s %s"%(g_params['CP_EXE'], thumbImageSourceFile,
            thumbImageTargetFile))
        print >> fpout, ("<a href=\"%s\"target=\"_blank\">"
                % (htmlname + os.sep + os.path.basename(imageTargetFile)))
        print >> fpout, ("<img src=\"%s\">" % (htmlname + os.sep +
            os.path.basename(thumbImageTargetFile)))
        print >> fpout, "</a>"
        print >> fpout, '</td>'
#---------------------------
        print >> fpout, '<td>'
        orig_image_file = (datapath + os.sep + "new_rstpair_" + method +
            os.sep +
            "dump.pairaln_type_all_dg0.0_gap0.8_ps0.0-NumTMDistri.seqidt30-100.txt.eps")

        imageSourceFile = orig_image_file.rstrip(".eps") + ".png"
        thumbImageSourceFile =  orig_image_file.rstrip(".eps") + ".thumb.png" 

        if not os.path.exists(imageSourceFile) or g_params['FORCE_OVERWRITE']:
            os.system("convert -density 150  %s %s"%(orig_image_file, imageSourceFile))
        if not os.path.exists(thumbImageSourceFile) or g_params['FORCE_OVERWRITE']:
            os.system("convert -thumbnail 200  %s %s"%(imageSourceFile, thumbImageSourceFile))
        imageTargetFile = (outpath + os.sep + htmlname + os.sep +  method + "_" +
            os.path.basename(imageSourceFile))
        thumbImageTargetFile = (outpath + os.sep + htmlname + os.sep + method + "_" +
            os.path.basename(thumbImageSourceFile))
        os.system("%s %s %s"%(g_params['CP_EXE'], imageSourceFile,
            imageTargetFile))
        os.system("%s %s %s"%(g_params['CP_EXE'], thumbImageSourceFile,
            thumbImageTargetFile))
        print >> fpout, ("<a href=\"%s\"target=\"_blank\">"
                % (htmlname + os.sep + os.path.basename(imageTargetFile)))
        print >> fpout, ("<img src=\"%s\">" % (htmlname + os.sep +
            os.path.basename(thumbImageTargetFile)))
        print >> fpout, "</a>"
        print >> fpout, '</td>'
#---------------------------
        print >> fpout, "</tr>"
    print >> fpout, "</table>"
#}}}
def WriteHTMLTable3(tablename, tabletitle, method_list, datapath, #{{{
        htmlname, outpath, fpout):

    print >> fpout, "<a name=\"%s\"></a><h4>%s</h4>"%(tablename,tabletitle)
    print >> fpout, "<table class=\"sortable\" border=1>"
    cntOutputID = 0

    headerItemList=[]
    headerItemList.append("Methods")
    headerItemList.append("cmpclass")
    headerItemList.append("NCtermInter")
    headerItemList.append("NumTMInter")
    print >> fpout, "<tr>"
    for item in headerItemList:
        print >> fpout, "<th>"
        print >> fpout, item
        print >> fpout, "</th>"
    print >> fpout, "</tr>"

    for i in xrange(len(method_list)):
        method = method_list[i]
        
        print >> fpout, "<tr>"
#---------------------------
        print >> fpout, '<td>'
        print >> fpout, '%s'%(method)
        print >> fpout, '</td>'
#---------------------------
        print >> fpout, '<td>'
        orig_image_file = (datapath + os.sep + "new_rstpair_" + method +
            os.sep +
            "dump.pairaln_type_all_dg10.0_gap0.8_ps0.0-cmpclass.table.txt.eps")

        imageSourceFile = orig_image_file.rstrip(".eps") + ".png"
        thumbImageSourceFile =  orig_image_file.rstrip(".eps") + ".thumb.png" 

        if not os.path.exists(imageSourceFile) or g_params['FORCE_OVERWRITE']:
            os.system("convert -density 150  %s %s"%(orig_image_file, imageSourceFile))
        if not os.path.exists(thumbImageSourceFile) or g_params['FORCE_OVERWRITE']:
            os.system("convert -thumbnail 200  %s %s"%(imageSourceFile, thumbImageSourceFile))
        imageTargetFile = (outpath + os.sep + htmlname + os.sep + method + "_" + 
            os.path.basename(imageSourceFile))
        thumbImageTargetFile = (outpath + os.sep + htmlname + os.sep +  method + "_" +
            os.path.basename(thumbImageSourceFile))
        os.system("%s %s %s"%(g_params['CP_EXE'], imageSourceFile,
            imageTargetFile))
        os.system("%s %s %s"%(g_params['CP_EXE'], thumbImageSourceFile,
            thumbImageTargetFile))
        print >> fpout, ("<a href=\"%s\"target=\"_blank\">"
                % (htmlname + os.sep + os.path.basename(imageTargetFile)))
        print >> fpout, ("<img src=\"%s\">" % (htmlname + os.sep +
            os.path.basename(thumbImageTargetFile)))
        print >> fpout, "</a>"
        print >> fpout, '</td>'
#---------------------------
        print >> fpout, '<td>'
        orig_image_file = (datapath + os.sep + "new_rstpair_" + method +
            os.sep +
            "dump.pairaln_type_all_dg10.0_gap0.8_ps0.0-NCtermInter.table2.c2.txt.eps")

        imageSourceFile = orig_image_file.rstrip(".eps") + ".png"
        thumbImageSourceFile =  orig_image_file.rstrip(".eps") + ".thumb.png" 

        if not os.path.exists(imageSourceFile) or g_params['FORCE_OVERWRITE']:
            os.system("convert -density 150  %s %s"%(orig_image_file, imageSourceFile))
        if not os.path.exists(thumbImageSourceFile) or g_params['FORCE_OVERWRITE']:
            os.system("convert -thumbnail 200  %s %s"%(imageSourceFile, thumbImageSourceFile))
        imageTargetFile = (outpath + os.sep + htmlname + os.sep +  method + "_" +
            os.path.basename(imageSourceFile))
        thumbImageTargetFile = (outpath + os.sep + htmlname + os.sep + method + "_" +
            os.path.basename(thumbImageSourceFile))
        os.system("%s %s %s"%(g_params['CP_EXE'], imageSourceFile,
            imageTargetFile))
        os.system("%s %s %s"%(g_params['CP_EXE'], thumbImageSourceFile,
            thumbImageTargetFile))
        print >> fpout, ("<a href=\"%s\"target=\"_blank\">"
                % (htmlname + os.sep + os.path.basename(imageTargetFile)))
        print >> fpout, ("<img src=\"%s\">" % (htmlname + os.sep +
            os.path.basename(thumbImageTargetFile)))
        print >> fpout, "</a>"
        print >> fpout, '</td>'
#---------------------------
        print >> fpout, '<td>'
        orig_image_file = (datapath + os.sep + "new_rstpair_" + method +
            os.sep +
            "dump.pairaln_type_all_dg10.0_gap0.8_ps0.0-NumTMDistri.seqidt30-100.txt.eps")

        imageSourceFile = orig_image_file.rstrip(".eps") + ".png"
        thumbImageSourceFile =  orig_image_file.rstrip(".eps") + ".thumb.png" 

        if not os.path.exists(imageSourceFile) or g_params['FORCE_OVERWRITE']:
            os.system("convert -density 150  %s %s"%(orig_image_file, imageSourceFile))
        if not os.path.exists(thumbImageSourceFile) or g_params['FORCE_OVERWRITE']:
            os.system("convert -thumbnail 200  %s %s"%(imageSourceFile, thumbImageSourceFile))
        imageTargetFile = (outpath + os.sep + htmlname + os.sep +  method + "_" +
            os.path.basename(imageSourceFile))
        thumbImageTargetFile = (outpath + os.sep + htmlname + os.sep + method + "_" +
            os.path.basename(thumbImageSourceFile))
        os.system("%s %s %s"%(g_params['CP_EXE'], imageSourceFile,
            imageTargetFile))
        os.system("%s %s %s"%(g_params['CP_EXE'], thumbImageSourceFile,
            thumbImageTargetFile))
        print >> fpout, ("<a href=\"%s\"target=\"_blank\">"
                % (htmlname + os.sep + os.path.basename(imageTargetFile)))
        print >> fpout, ("<img src=\"%s\">" % (htmlname + os.sep +
            os.path.basename(thumbImageTargetFile)))
        print >> fpout, "</a>"
        print >> fpout, '</td>'
#---------------------------
        print >> fpout, "</tr>"
    print >> fpout, "</table>"
#}}}

def WriteHTML(method_list, datapath, htmlname, outpath): #{{{
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
        title="topology variation analysis with topologies predicted by different methods"
        WriteHTMLHeader(title, fpout)

#Write header line

#write right panel        
        print >> fpout, "<dir id=\"Content\">"

#=========================================
        tablename = 'table1'
        tabletitle = "Non filtered"
        WriteHTMLTable(tablename, tabletitle, method_list, datapath, htmlname,
                outpath, fpout)
#=========================================
        tablename = 'table2'
        tabletitle = "DG <= 0.0, Gap fraction >= 0.8"
        WriteHTMLTable2(tablename, tabletitle, method_list, datapath, htmlname,
                outpath, fpout)
#=========================================
        tablename = 'table3'
        tabletitle = "Gap fraction >= 0.8"
        WriteHTMLTable3(tablename, tabletitle, method_list, datapath, htmlname,
                outpath, fpout)
#=========================================
        print >> fpout, "</dir>"

        WriteHTMLTail(fpout)
        fpout.close()
        print  "Result has been output to %s"%htmlfilename

    except IOError:
        print >> sys.stderr, "%s: Error to write HTML file %s to %s. Exit." %(sys.argv[0], htmlfilename, figuredir)
        raise
        return 1
#}}}


def main(g_params):#{{{
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    i = 1
    isNonOptionArg=False
    datapath = ""
    outpath = ""
    htmlname = 'index'
    while i < numArgv:#{{{
        if isNonOptionArg == True:
            idList.append(sys.argv[i])
            isNonOptionArg=False
            i += 1
        elif sys.argv[i] == "--":
            isNonOptionArg=True
            i += 1
        elif sys.argv[i][0] == "-":
            if sys.argv[i] ==  "-h" or  sys.argv[i] == "--help":
                PrintHelp()
                return 1
            elif sys.argv[i] == "-l" or sys.argv[i] == "--l":
                idListFile = sys.argv[i+1]
                i += 2
            elif sys.argv[i] == "-datapath" or sys.argv[i] == "--datapath":
                datapath = sys.argv[i+1]
                i += 2
            elif sys.argv[i] in ["-topomsapath", "--topomsapath"]:
                topomsapath = sys.argv[i+1]
                i += 2
            elif sys.argv[i] in ["-ordermsapath", "--ordermsapath"]:
                ordermsapath = sys.argv[i+1]
                i += 2
            elif sys.argv[i] == "-outpath" or sys.argv[i] == "--outpath":
                outpath = sys.argv[i+1]
                i += 2
            elif sys.argv[i] == "-htmlname" or sys.argv[i] == "--htmlname":
                htmlname = sys.argv[i+1]
                i += 2
            elif sys.argv[i] == "-min-numseq" or sys.argv[i] == "--min-numseq":
                g_params['MIN_NUMSEQ'] = int(sys.argv[i+1])
                i += 2
            elif sys.argv[i] == "-pfamdef" or sys.argv[i] == "--pfamdef":
                pfamACDEListFile = sys.argv[i+1]
                i += 2
            elif sys.argv[i] == "-q":
                isQuiet = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", sys.argv[i]
                return 1
        else:
            idList.append(sys.argv[i])
            i += 1
#}}}

    datapath = "."
    method_list = [
            "scampi_single",
            "hmmtop",
            "stmhmm",
            "memsat",
            "topcons_single_method4",
            "scampi_msa",
            "octopus",
            "prodiv",
            "pro",
            "topcons"
            ]

    g_params['OS'] = os.uname()[0]
    if g_params['OS'].find('Linux') != -1:
        g_params['CP_EXE'] = "/bin/cp -uf"
    else:
        g_params['CP_EXE'] = "/bin/cp -f"
    
    if not os.path.exists(outpath):
        os.system("mkdir -p %s"%outpath)
    WriteHTML(method_list, datapath, htmlname, outpath)
    return 0
#}}}

if __name__ == '__main__' :
    g_params = {}
    g_params['FORCE_OVERWRITE'] = True
    sys.exit(main(g_params))

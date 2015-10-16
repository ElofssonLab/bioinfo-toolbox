#!/usr/bin/env python
# given paircmp
# report the result in html table

# table format
# pfamid pfam_description
# sequence identity
# refseq_id_1     sequence_description
# refseq_id_2     sequence_description
# TM mapping
# Internal gap percentage DG value
# alignement figure shrinked
# alignment figure non shrinked   

# definition of paircmp record
#     record['id1']
#     record['id2']
#     record['seqidt']
#     record['seqidt1']
#     record['seqidt2']
#     record['seqdef1']
#     record['seqdef2']
#     record['cmpclass']
#     record['numTM1']
#     record['numTM2']
#     record['seqLength1']
#     record['seqLength2']
#     record['mapTMline']=[]
#     record['mapArray'] = []
#     record['ana1'] = {}
#     record['ana2'] = {}
#     record['member'] = []


import os
import sys
import re
import libtopologycmp as lcmp

rundir = os.path.dirname(sys.argv[0])
binpath = rundir

progname =  os.path.basename(sys.argv[0])
usage="""
Usage:  reportPaircmpAlignment.py paircmpfile -alnfigpath DIR

Description: Report the result of MSA topology comparison in HTML format

Options: 
  -mode 0|1          report pairalignment mode, (default: 1)
                     0: select 
                     1: supply pairs
  -tableformat INT   table format, 0, 1 (default: 1)
                     0: each pair is a subtable
                     1: in a sortable single table
  -topomsapath DIR   Set path for topology MSA
  -ordermsapath DIR  Set path for topology MSA reordered according to phylo tree
  -treepath     DIR  Set path for phylo tree
  -outpath DIR       set ouput path, (default: ./)
  -htmlname STR      set the name of the main html file, (default: index)
  -pfamdef FILE      pfam definition file, (default:
                     /data3/data/pfam/pfamA.seed.ac-delist)  
  -seqdef  FILE      sequence definition file
                     (default: /data3/wk/MPTopo/pfamAna/pfam2-selTM-giid-refseqid-pfamid-description.txt)
  -pickone           Pick only one example from each Pfam family
  -msapath           multiple sequence alignment path
  -msapath2          multiple sequence alignment path 2
  -type STR          selecting type, nterm, cterm, internal, all
                     (default: all)
  -tableinfo FILE    Set pairwise alignment table info, get more pairwise
                     statistics
  -filter-predseq    y | n 
                     Filter predicted sequences
  -q                 quiet mode
  -h, --help         print this help message and exit

Selection control options:
  -gap  FLOAT  Select only the TM with gap fraction >= threshold, (default: 0.5)
  -dg   FLOAT  Select only the TM with DG values <= threshold, (default: 1.0)
  -ps   FLOAT  Select only the TM with topology prediction reliability >=
               threshold, (default: 0.5)
  -min-seqidt  FLOAT, (default: 0)
  -max-seqidt  FLOAT, (default: 100)
               Select only pairs with global sequence identity within [minSeqIDT, maxSeqIDT]

Created 2012-02-09, updated 2013-04-19, Nanjiang Shu  

Examples:
    reportPaircmpAlignment t1.paircmp -alnfigpath figAlnPair -outpath selectedpair
"""%(progname)

BLOCK_SIZE = myfunc.BLOCK_SIZE
# typeTopoAnaDict
# 0: started from $id.fa with multiple homologous sequences
# 1: started from $id.fa with a single sequence, and the homologous sequences
# were obtained by blast searching

def PrintHelp():
    print usage


def GeneratePairAlnFigure(outpath, pairlistfile):#{{{
    pairalndbfile = "Pfam-A-full.perTM75_nseq20.nr100.filter.fragmented.withinclan.max30000.kalignP_filterNonTM.pairaln" 
    out_pairalnfile = outpath + os.sep + "result.pairaln"
    cmd = "%s/selectPairaln.py -l %s -pairaln %s > %s"%(binpath, pairlistfile, pairalndbfile, out_pairalnfile)
    os.system(cmd)
# /data3/wk/MPTopo/src/selectPairaln.py -l data/highseqidt_TM2GAP.pairlist -pairaln Pfam-A-full.perTM75_nseq20.nr100.filter.fragmented.withinclan.max30000.kalignP_filterNonTM.pairaln > data/highseqidt_TM2GAP.pairaln

    alnfigpath = outpath + os.sep + "alnfigpath"
    g_params['alnFigPath'] = alnfigpath
    topofile = "/data3/wk/MPTopo/pfamAna_refpro/pred_topcons/refpro20120604-celluar.selmaxlength-m1.topcons.result_TOPCONS.topo"
    cmd = "%s/drawpairwisetopo.py %s -outpath %s -topofile %s "%(binpath, out_pairalnfile, alnfigpath, topofile)
    os.system(cmd)
#}}}
#/data3/wk/MPTopo/src/drawpairwisetopo.py  data/highseqidt_TM2GAP.pairaln -outpath data/figpairaln -topofile /data3/wk/MPTopo/pfamAna_refpro/pred_topcons/refpro20120604-celluar.selmaxlength-m1.topcons.result_TOPCONS.topo -seqdb Pfam-A-full.perTM75_nseq20.nr100.filter.fragmented.uniq.fasta
                                                              
def AppendRLTYInfo(recordList, rltyDict):#{{{
    for record in recordList:
        id1 = record['id1']
        id2 = record['id2']
        try: 
            rlty = min(rltyDict[id1], rltyDict[id2])
            record['rlty'] = rlty
        except KeyError:
#            print id1, id2, "not in rltyDict"
            record['rlty'] = -100
            pass#}}}
def AppendFamInfo(recordList, pair2famidDict, famDefDict):#{{{
    for record in recordList:
        id1 = record['id1']
        id2 = record['id2']
        try: 
            famid = pair2famidDict[(id1,id2)]
            record['famid'] = famid
            try:
                famdef = famDefDict[famid]
                record['famdef'] =famdef
            except KeyError:
                record['famdef'] = ""
                pass
        except KeyError:
            record['famid'] = ""
            record['famdef'] = ""
            pass#}}}
def AppendSeqAnnoInfo(recordList,seqannoDict):#{{{
    for record in recordList:
        id1 = record['id1']
        id2 = record['id2']
        try: 
            record['seqanno1'] = seqannoDict[id1]
            record['seqanno2'] = seqannoDict[id2]
        except KeyError:
            record['seqanno1'] = ""
            record['seqanno2'] = ""
            pass
#}}}
def ReadInRLTYInfo(infile):#{{{
    try:
        fpin = open(infile, "r")
        lines = fpin.readlines()
        fpin.close()
        rltyDict = {}
        for line in lines:
            strs = line.split()
            if len(strs) == 2:
                seqid = strs[0]
                rlty = float(strs[1])
                rltyDict[seqid] = rlty
        return rltyDict
    except IOError:
        print >> sys.stderr, "Failed to read rltyinfoFile %s"%infile
        return {}
#}}}
def ReadFamDefFile(infile):#{{{
    try:
        fpin = open(infile, "r")
        lines = fpin.readlines()
        fpin.close()
        dt = {}
        for line in lines:
            strs = line.split("\t")
            try:
               pfamid = strs[0]
               clanid = strs[1]
               if clanid == "\N":
                   dt[pfamid] = strs[3]
               else:
                   dt[clanid] = strs[3]
            except IndexError:
                pass
        return dt
    except IOError:
        print >> sys.stderr, "Failed to read file %s"%infile
        return {}
#}}}
def ReadPair2FamIDInfo(infile):#{{{
    try:
        fpin = open(infile, "r")
        lines = fpin.readlines()
        fpin.close()
        dt = {}
        for line in lines:
            strs = line.split()
            if len(strs) == 3:
                id1 = strs[0]
                id2 = strs[1]
                famid = strs[2]
                dt[(id1,id2)] = famid
        return dt
    except IOError:
        print >> sys.stderr, "Failed to read file %s"%infile
        return {}#}}}
def ReadIDWithAnnoInfo(infile):#{{{
    try:
        fpin = open(infile, "r")
        lines = fpin.readlines()
        fpin.close()
        dt = {}
        for line in lines:
            strs = line.split("\t")
            if len(strs) == 2:
                seqid = strs[0]
                anno = strs[1].strip()
                dt[seqid] = anno
        return dt
    except IOError:
        print >> sys.stderr, "Failed to read file %s"%infile
        return {}#}}}
def countseq(inFile):#{{{
    try:
        isFirstSeq=True
        isPreviousBuffEndWithNewLine=False
        fpin = open(inFile, "r")
        buff = fpin.read(BLOCK_SIZE)
        cntSeq = 0
        while buff:
            if isFirstSeq and buff[0] == '>':
                cntSeq +=1
                isFirstSeq = False
            if isPreviousBuffEndWithNewLine and buff[0] == '>':
                cntSeq += 1
                isPreviousBuffEndWithNewLine = False
            cntSeq += buff.count("\n>")
            if buff[len(buff)-1] == '\n':
                isPreviousBuffEndWithNewLine = True
            buff = fpin.read(BLOCK_SIZE)

        fpin.close()
        return cntSeq
    except IOError:
        print >> sys.stderr,"Fail to open file %s"%(inFile)
        return 0

#}}}
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
        print >> sys.stderr, "%s: Error reading %s"%(argv[0], infile)
        return 1
#}}}

def ReadInTableInfo(infile):#{{{
    try:
        fpin = open(infile, "r")
        lines = fpin.readlines()
        fpin.close()
        pairalnStat = {}
        for line in lines:
            if line[0] != "#":
                strs = line.split()
                if len(strs) == 13:
                    id1 = strs[0]
                    id2 = strs[1]
                    seqidt = float(strs[2])
                    alignLen = float(strs[4])
                    seqlen1 = int(strs[5])
                    seqlen2 = int(strs[6])
                    seqidt1 = float(strs[11])
                    seqidt2 = float(strs[12])
                    pairid = id1+'-'+id2
                    pairalnStat[pairid] = {}
                    pairalnStat[pairid]['seqidt'] = seqidt
                    pairalnStat[pairid]['seqidt1'] = seqidt1
                    pairalnStat[pairid]['seqidt2'] = seqidt1
                    pairalnStat[pairid]['seqLength1'] = seqlen1
                    pairalnStat[pairid]['seqLength2'] = seqlen2
        return pairalnStat
    except IOError:
        print >> sys.stderr, "Failed to open tableinfo file %s"%infile
        return {}
#}}}
def AddTableInfo(recordList, pairalnStat):#{{{
    if pairalnStat != {}:
        for record in recordList:
            pairid = record['id1'] + '-' + record['id2']
            try: 
                dd = pairalnStat[pairid]
                record['seqidt'] = dd['seqidt']
                record['seqidt1'] = dd['seqidt1']
                record['seqidt2'] = dd['seqidt2']
            except KeyError:
                pass
#}}}

def AddSeqDefInfo(recordList, seqInfoDict):#{{{
    if seqInfoDict != {}:
        for record in recordList:
            id1 =  record['id1']
            id2 =  record['id2']
            if id1 in seqInfoDict:
                record['seqdef1'] = seqInfoDict[id1]['seqdef']
            if id2 in seqInfoDict:
                record['seqdef2'] = seqInfoDict[id2]['seqdef']
#}}}
def WriteUnmappedRecordHTMLCell(subAna,fpout):#{{{
    print >> fpout, "<td>"
    print >> fpout, "<b>%d</b><br>"%(subAna['numTMunmapped'])
    strIndex=""
    strGapPercentage=""
    strDGvalue=""
    for i in range(len(subAna['index'])):
        strIndex += "%d "%subAna['index'][i] 
        strGapPercentage += "%.1f "%subAna['gapFraction'][i] 
        strDGvalue += "%.2f "%subAna['DGvalue'][i] 
    fpout.write("<font size=2>index(%s)<br>gap(%s)<br>DG(%s)</font>"%(strIndex.strip(), strGapPercentage.strip(), strDGvalue.strip()))
    print >> fpout, "</td>"
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
    print >> fpout, "table.sortable td"
    print >> fpout, "{"
    print >> fpout, "   font-family:Arial; font-size: 11px;"
    print >> fpout, "   vertical-align: middle;"
    print >> fpout, "   text-align: center;"
    print >> fpout, "   max-width: 200px;"
    print >> fpout, "   border: 1px solid black;"
    print >> fpout, "   padding: 2px 2px 2px 2px;"
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
    print >> fpout, "   border: 0px solid black;"
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
    print >> fpout, "   font-family:Arial; font-size: 9px;"
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
def WriteTMMapHTMLCell(MapTMLine, fpout):#{{{
    text = ""
    text += "<td>\n"
    text += "<table class=subtable>\n"
    for i in xrange(2):
        mapList = MapTMLine[i].split(':')[1].split()
        text += "<tr>\n"
        for s in mapList:
            text += "<td>%s</td>"%s
        text += "\n</tr>\n"
    text += "</table>\n"
    text += "</td>"
    print >> fpout, text
#}}}


def WriteSubtableHits(record, maxNumInternalCons, maxNumInternalQuery, fpout):#{{{
    print >> fpout, "<table class=hits>"
    WriteTableHeaderHits(maxNumInternalCons, maxNumInternalQuery, fpout)
    WriteTableContentHits(record, 1, maxNumInternalCons,
            maxNumInternalQuery,fpout)
    print >> fpout, "</table>"
#}}}

def WriteHTMLTable(tablename, tabletitle,pairCmpRecordList, numPair, numFam, #{{{
        outpath, htmlname, fpout): 

    alnFigPath = g_params['alnFigPath']
    MSAPath = g_params['MSAPath']
    MSAPath2 = g_params['MSAPath2']

    print >> fpout, "<a name=\"%s\"></a><h2>%s</h2>"%(tablename,tabletitle)
    print >> fpout, "<table class=content>"
    for i in xrange(numPair):
        record = pairCmpRecordList[i]
        id1 = record['id1']
        id2 = record['id2']
        seqidt = record['seqidt']
        pfamid = record['famid']
        
        print >> fpout, "<tr><td>"

        print >> fpout, "<table class=record>"
        print >> fpout, "<tr><td>"
        print >> fpout, "<b>Pair %d</b>" %(i+1)
        print >> fpout, "</td></tr>"
        

#---------------------------
        print >> fpout, "<tr>"

        print >> fpout, '<td>'
        print >> fpout, 'Sequence 1'
        print >> fpout, '</td>'

        print >> fpout, '<td>'
        print >> fpout, '<a href=\"http://www.uniprot.org/uniprot/%s\"target=\"_blank\">%s</a>'%(id1, id1)
        print >> fpout, '</td>'


        print >> fpout, '<td>'
        print >> fpout, '%s'%record['seqanno1']
        print >> fpout, '</td>'

        print >> fpout, "</tr>"

#---------------------------
        print >> fpout, "<tr>"

        print >> fpout, '<td>'
        print >> fpout, 'Sequence 2'
        print >> fpout, '</td>'

        print >> fpout, '<td>'
        print >> fpout, '<a href=\"http://www.uniprot.org/uniprot/%s\"target=\"_blank\">%s</a>'%(id2, id2)
        print >> fpout, '</td>'

        print >> fpout, '<td>'
        print >> fpout, '%s'%record['seqanno2']
        print >> fpout, '</td>'

        print >> fpout, "</tr>"
#---------------------------
        print >> fpout, "<tr>"

        print >> fpout, '<td>'
        print >> fpout, 'Sequence identity'
        print >> fpout, '</td>'

        print >> fpout, '<td>'
        print >> fpout, '%s'%record['seqidt']
        print >> fpout, '</td>'

        print >> fpout, "</tr>"
#---------------------------
        print >> fpout, "<tr>"

        print >> fpout, '<td>'
        print >> fpout, 'PfamID and definition'
        print >> fpout, '</td>'

        print >> fpout, '<td>'
        if pfamid.find("PF") != -1:
            pfamURL = 'http://pfam.sanger.ac.uk/family/' + pfamid
        else:
            pfamURL = 'http://pfam.sanger.ac.uk/clan/' + pfamid

        print >> fpout, '<a href=\"%s\"target=\"_blank\">%s</a>'%(pfamURL,
                pfamid)
        print >> fpout, '</td>'

        print >> fpout, '<td>'
        print >> fpout, '%s'%record['famdef']
        print >> fpout, '</td>'

        print >> fpout, "</tr>"
#---------------------------
        print >> fpout, "<tr>"

        print >> fpout, '<td>'
        print >> fpout, 'TM mapping'
        print >> fpout, '</td>'

        WriteTMMapHTMLCell(record['mapTMline'], fpout)

        print >> fpout, "</tr>"
#---------------------------
        print >> fpout, "<tr>"

        print >> fpout, '<td>'
        print >> fpout, 'Alignment (shrinked, showing ioM)'
        print >> fpout, '</td>'

        print >> fpout, '<td>'
        bsname = '%s-%s'%(id1,id2)
        ext = '.topoaln.shrinked.png'
        imageSourceFile = alnFigPath + os.sep + bsname + ext
        imageTargetFile = outpath + os.sep + htmlname + os.sep + bsname + ext
        thumb_imageSourceFile = alnFigPath + os.sep + 'thumb.' + bsname + ext
        thumb_imageTargetFile = outpath + os.sep + htmlname + os.sep + 'thumb.' + bsname + ext

        if os.path.exists(imageSourceFile):
            os.system("/bin/cp -uf %s %s"%(imageSourceFile, imageTargetFile))
        if os.path.exists(thumb_imageSourceFile):
            os.system("/bin/cp -uf %s %s"%(thumb_imageSourceFile, thumb_imageTargetFile))
        print >> fpout, ("<a href=\"%s\"target=\"_blank\">"
                % (htmlname + os.sep + os.path.basename(imageTargetFile)))
        print >> fpout, ("<img src=\"%s\">" % (htmlname + os.sep +
            os.path.basename(thumb_imageTargetFile)))
        print >> fpout, "</a>"
        print >> fpout, '</td>'

        print >> fpout, "</tr>"
#---------------------------
        print >> fpout, "<tr>"

        print >> fpout, '<td>'
        print >> fpout, 'Alignment (non shrinked, showing amino acids)'
        print >> fpout, '</td>'

        print >> fpout, '<td>'
        bsname = '%s-%s'%(id1,id2)
        ext = '.topoaln.nonshrinked.png'
        imageSourceFile = alnFigPath + os.sep + bsname + ext
        imageTargetFile = outpath + os.sep + htmlname + os.sep + bsname + ext
        thumb_imageSourceFile = alnFigPath + os.sep + 'thumb.' + bsname + ext
        thumb_imageTargetFile = outpath + os.sep + htmlname + os.sep + 'thumb.' + bsname + ext

        if os.path.exists(imageSourceFile):
            os.system("/bin/cp -uf %s %s"%(imageSourceFile, imageTargetFile))
        if os.path.exists(thumb_imageSourceFile):
            os.system("/bin/cp -uf %s %s"%(thumb_imageSourceFile, thumb_imageTargetFile))
        print >> fpout, ("<a href=\"%s\"target=\"_blank\">"
                % (htmlname + os.sep + os.path.basename(imageTargetFile)))
        print >> fpout, ("<img src=\"%s\">" % (htmlname + os.sep +
            os.path.basename(thumb_imageTargetFile)))
        print >> fpout, "</a>"
        print >> fpout, '</td>'

        print >> fpout, "</tr>"
#---------------------------
        if os.path.exists(MSAPath):
            print >> fpout, "<tr>"

            print >> fpout, '<td>'
            print >> fpout, 'TopoMSA of Pfam %s'%pfamid
            print >> fpout, '</td>'

            print >> fpout, '<td>'
            print >> fpout, 'By KalignP<br>'
            imgdir1 = outpath + os.sep + 'image1'
            os.system("mkdir -p %s"%imgdir1)
            imageSourceFile = MSAPath + os.sep + pfamid + '.grouped.sorted.orig.topomsa.png'
            imageTargetFile = imgdir1 + os.sep + os.path.basename(imageSourceFile)
            thumbSourceFile = MSAPath + os.sep + 'thumb.' + pfamid + '.grouped.sorted.orig.topomsa.png'
            thumbTargetFile = imgdir1 + os.sep + os.path.basename(thumbSourceFile)
            if os.path.exists(imageSourceFile):
                os.system("/bin/cp -uf %s %s"%(imageSourceFile, imageTargetFile))
            if os.path.exists(thumb_imageSourceFile):
                os.system("/bin/cp -uf %s %s"%(thumbSourceFile, thumbTargetFile))
            print >> fpout, ("<a href=\"%s\"target=\"_blank\">"
                    % ('image1' + os.sep + os.path.basename(imageTargetFile)))
            print >> fpout, ("<img src=\"%s\">" % ('image1' + os.sep +
                os.path.basename(thumbTargetFile)))
            print >> fpout, "</a>"
            print >> fpout, '</td>'


            if os.path.exists(MSAPath2):
                print >> fpout, '<td>'
                print >> fpout, 'By ClustalO<br>'
                imgdir2 = outpath + os.sep + 'image2'
                os.system("mkdir -p %s"%imgdir2)
                imageSourceFile = MSAPath2 + os.sep + pfamid + '.grouped.sorted.orig.topomsa.png'
                imageTargetFile = imgdir2 + os.sep + os.path.basename(imageSourceFile)
                thumbSourceFile = MSAPath2 + os.sep + 'thumb.' + pfamid + '.grouped.sorted.orig.topomsa.png'
                thumbTargetFile = imgdir2 + os.sep + os.path.basename(thumbSourceFile)
                if os.path.exists(imageSourceFile):
                    os.system("/bin/cp -uf %s %s"%(imageSourceFile, imageTargetFile))
                if os.path.exists(thumbSourceFile):
                    os.system("/bin/cp -uf %s %s"%(thumbSourceFile, thumbTargetFile))
                print >> fpout, ("<a href=\"%s\"target=\"_blank\">"
                        % ('image2' + os.sep + os.path.basename(imageTargetFile)))
                print >> fpout, ("<img src=\"%s\">" % ('image2' + os.sep +
                    os.path.basename(thumbTargetFile)))
                print >> fpout, "</a>"
                print >> fpout, '</td>'

            print >> fpout, "</tr>"
#---------------------------

        print >> fpout, "</table>"

        print >> fpout, "</td></tr>"

    print >> fpout, "</table>"
#}}}
def WriteHTMLTable2(tablename, tabletitle,pairCmpRecordList, numPair, numFam, #{{{ 
        outpath, htmlname, fpout): 
    """
    Write in a sortable single table
    """
    if len(pairCmpRecordList) == 0:
        return 1

    alnFigPath = g_params['alnFigPath']
    MSAPath = g_params['MSAPath']
    MSAPath2 = g_params['MSAPath2']
    topomsapath = g_params['topomsapath']
    ordermsapath = g_params['ordermsapath']
    treepath = g_params['treepath']

    print >> fpout, "<a name=\"%s\"></a><h2>%s</h2>"%(tablename,tabletitle)
    print >> fpout, "<table class=\"sortable\" border=1>"


    if 'seqidt1' in pairCmpRecordList[0]:
        isWriteSeqIDT1 = True
    else:
        isWriteSeqIDT1 = False
# write html table header
    headerItemList=[]
    headerItemList.append("No.")
    headerItemList.append("seq1")
    headerItemList.append("seq2")
    headerItemList.append("seqIDT")
    if isWriteSeqIDT1:
        headerItemList.append("seqIDT1")
    headerItemList.append("seqLen1")
    headerItemList.append("seqLen2")
    headerItemList.append("pfamID")
    headerItemList.append("pfamDef")
    headerItemList.append("TMmap")
    headerItemList.append("Alignment (shrinked)")
    headerItemList.append("Alignment (Non shrinked)")
    if g_params['treepath'] != "" and os.path.exists(g_params['treepath']):
        headerItemList.append("Phylo Tree")
    if g_params['ordermsapath'] != "" and os.path.exists(g_params['ordermsapath']):
        headerItemList.append("Topology MSA ordered according to phylo tree")
    if g_params['topomsapath'] != "" and os.path.exists(g_params['topomsapath']):
        headerItemList.append("Topology MSA grouped by topology comparison")

    print >> fpout, "<tr>"
    for item in headerItemList:
        print >> fpout, "<th>"
        print >> fpout, item
        print >> fpout, "</th>"
    print >> fpout, "</tr>"

# write html table content
    if g_params['isShowProgress']:
        print 'Output HTML table ...'
    for i in xrange(numPair):
        record = pairCmpRecordList[i]
        id1 = record['id1']
        id2 = record['id2']
        seqLength1 = record['seqLength1']
        seqLength2 = record['seqLength2']
        seqidt = record['seqidt']
        pfamid = record['famid']

        if g_params['isShowProgress']:
            if i%10 == 0:
                print i, "..."
        
        print >> fpout, "<tr>"
#---------------------------
        print >> fpout, "<td>"
        print >> fpout, "%d"%(i+1)
        print >> fpout, "</td>"
#---------------------------
        print >> fpout, '<td>'
        print >> fpout, '<a href=\"http://www.uniprot.org/uniprot/%s\"target=\"_blank\">%s</a>'%(id1, id1)
        print >> fpout, "<br>%s"%record['seqanno1']
        print >> fpout, '</td>'
#---------------------------
        print >> fpout, '<td>'
        print >> fpout, '<a href=\"http://www.uniprot.org/uniprot/%s\"target=\"_blank\">%s</a>'%(id2, id2)
        print >> fpout, "<br>%s"%record['seqanno2']
        print >> fpout, '</td>'
#---------------------------
        print >> fpout, '<td>'
        print >> fpout, '%.1f'%seqidt
        print >> fpout, '</td>'
#---------------------------
        if isWriteSeqIDT1:
            print >> fpout, '<td>'
            if 'seqidt1' in record:
                print >> fpout, '%.1f'%record['seqidt1']
            else:
                print >> fpout, '-'
            print >> fpout, '</td>'
#---------------------------
        print >> fpout, '<td>'
        print >> fpout, '%d'%seqLength1
        print >> fpout, '</td>'
#---------------------------
        print >> fpout, '<td>'
        print >> fpout, '%d'%seqLength2
        print >> fpout, '</td>'
#---------------------------
        print >> fpout, '<td>'
        if pfamid.find("PF") != -1:
            pfamURL = 'http://pfam.sanger.ac.uk/family/' + pfamid
        else:
            pfamURL = 'http://pfam.sanger.ac.uk/clan/' + pfamid
        
        print >> fpout, '<a href=\"%s\"target=\"_blank\">%s</a>'%(pfamURL,
                pfamid)
        print >> fpout, '</td>'
#---------------------------
        print >> fpout, '<td>'
        print >> fpout, '%s'%record['famdef']
        print >> fpout, '</td>'
#---------------------------
        #print >> fpout, '<td>'
        WriteTMMapHTMLCell(record['mapTMline'], fpout)
        #print >> fpout, '</td>'
#---------------------------
        print >> fpout, '<td>'
        bsname = '%s-%s'%(id1,id2)
        ext = '.topoaln.shrinked.png'
        imageSourceFile = alnFigPath + os.sep + bsname + ext
        imageTargetFile = outpath + os.sep + htmlname + os.sep + bsname + ext
        thumb_imageSourceFile = alnFigPath + os.sep + 'thumb.' + bsname + ext
        thumb_imageTargetFile = outpath + os.sep + htmlname + os.sep + 'thumb.' + bsname + ext

        if os.path.exists(imageSourceFile):
            os.system("/bin/cp -uf %s %s"%(imageSourceFile, imageTargetFile))
        if os.path.exists(thumb_imageSourceFile):
            os.system("/bin/cp -uf %s %s"%(thumb_imageSourceFile, thumb_imageTargetFile))
        print >> fpout, ("<a href=\"%s\"target=\"_blank\">"
                % (htmlname + os.sep + os.path.basename(imageTargetFile)))
        print >> fpout, ("<img src=\"%s\">" % (htmlname + os.sep +
            os.path.basename(thumb_imageTargetFile)))
        print >> fpout, "</a>"
        print >> fpout, '</td>'
#---------------------------
        print >> fpout, '<td>'
        bsname = '%s-%s'%(id1,id2)
        ext = '.topoaln.nonshrinked.png'
        imageSourceFile = alnFigPath + os.sep + bsname + ext
        imageTargetFile = outpath + os.sep + htmlname + os.sep + bsname + ext
        thumb_imageSourceFile = alnFigPath + os.sep + 'thumb.' + bsname + ext
        thumb_imageTargetFile = outpath + os.sep + htmlname + os.sep + 'thumb.' + bsname + ext

        if os.path.exists(imageSourceFile):
            os.system("/bin/cp -uf %s %s"%(imageSourceFile, imageTargetFile))
        if os.path.exists(thumb_imageSourceFile):
            os.system("/bin/cp -uf %s %s"%(thumb_imageSourceFile, thumb_imageTargetFile))
        print >> fpout, ("<a href=\"%s\"target=\"_blank\">"
                % (htmlname + os.sep + os.path.basename(imageTargetFile)))
        print >> fpout, ("<img src=\"%s\">" % (htmlname + os.sep +
            os.path.basename(thumb_imageTargetFile)))
        print >> fpout, "</a>"
        print >> fpout, '</td>'
#---------------------------
        if g_params['treepath'] != "" and os.path.exists(g_params['treepath']):
            bsname = pfamid
            ext = '-itol.jpg'
            imageSourceFile = treepath + os.sep + bsname + ext
            imageTargetFile = outpath + os.sep + htmlname + os.sep + bsname + ext
            thumb_imageSourceFile = treepath + os.sep + 'thumb.' + bsname + ext
            thumb_imageTargetFile = outpath + os.sep + htmlname + os.sep + 'thumb.' + bsname + ext

            if os.path.exists(imageSourceFile):
                os.system("/bin/cp -uf %s %s"%(imageSourceFile, imageTargetFile))
            if os.path.exists(thumb_imageSourceFile):
                os.system("/bin/cp -uf %s %s"%(thumb_imageSourceFile, thumb_imageTargetFile))
            print >> fpout, '<td>'
            print >> fpout, ("<a href=\"%s\"target=\"_blank\">"
                    % (htmlname + os.sep + os.path.basename(imageTargetFile)))
            print >> fpout, ("<img src=\"%s\">" % (htmlname + os.sep +
                os.path.basename(thumb_imageTargetFile)))
            print >> fpout, "</a>"
            print >> fpout, '</td>'
#---------------------------
        if g_params['ordermsapath'] != "" and os.path.exists(g_params['ordermsapath']):
            bsname = pfamid
            ext = '.reordered.topomsa.png'
            imageSourceFile = ordermsapath + os.sep + bsname + ext
            imageTargetFile = outpath + os.sep + htmlname + os.sep + bsname + ext
            thumb_imageSourceFile = ordermsapath + os.sep + 'thumb.' + bsname + ext
            thumb_imageTargetFile = outpath + os.sep + htmlname + os.sep + 'thumb.' + bsname + ext

            if os.path.exists(imageSourceFile):
                os.system("/bin/cp -uf %s %s"%(imageSourceFile, imageTargetFile))
            if os.path.exists(thumb_imageSourceFile):
                os.system("/bin/cp -uf %s %s"%(thumb_imageSourceFile, thumb_imageTargetFile))
            print >> fpout, '<td>'
            print >> fpout, ("<a href=\"%s\"target=\"_blank\">"
                    % (htmlname + os.sep + os.path.basename(imageTargetFile)))
            print >> fpout, ("<img src=\"%s\">" % (htmlname + os.sep +
                os.path.basename(thumb_imageTargetFile)))
            print >> fpout, "</a>"
            print >> fpout, '</td>'
#---------------------------
        if g_params['topomsapath'] != "" and os.path.exists(g_params['topomsapath']):
            bsname = pfamid
            ext = '.sorted.orig.topomsa.png' 
            imageSourceFile = topomsapath + os.sep + bsname + ext
            imageTargetFile = outpath + os.sep + htmlname + os.sep + bsname + ext
            thumb_imageSourceFile = topomsapath + os.sep + 'thumb.' + bsname + ext
            thumb_imageTargetFile = outpath + os.sep + htmlname + os.sep + 'thumb.' + bsname + ext

            if os.path.exists(imageSourceFile):
                os.system("/bin/cp -uf %s %s"%(imageSourceFile, imageTargetFile))
            if os.path.exists(thumb_imageSourceFile):
                os.system("/bin/cp -uf %s %s"%(thumb_imageSourceFile, thumb_imageTargetFile))
            print >> fpout, '<td>'
            print >> fpout, ("<a href=\"%s\"target=\"_blank\">"
                    % (htmlname + os.sep + os.path.basename(imageTargetFile)))
            print >> fpout, ("<img src=\"%s\">" % (htmlname + os.sep +
                os.path.basename(thumb_imageTargetFile)))
            print >> fpout, "</a>"
            print >> fpout, '</td>'
#---------------------------
        print >> fpout, "</tr>"
#---------------------------
    if g_params['isShowProgress']:
        print "Finished"

    print >> fpout, "</table>"
#}}}
def WriteHTML(pairCmpRecordList, htmlname, outpath): #{{{
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
        print "Write html to %s"%outpath
        os.system("mkdir -p %s"%(figuredir)) 
        fpout = open(htmlfilename,"w")
        title="Selected pairwise topology comparison with internal topology variation"
        WriteHTMLHeader(title, fpout)

        numPair = len (pairCmpRecordList)
        uniqFamIDSet = set([])
        for record in pairCmpRecordList:
            famid = record['famid']
            uniqFamIDSet.add(famid)
        numFam = len (uniqFamIDSet)
#Write header line
        print >> fpout, "<dir id=\"Header\">"
        print >> fpout, "<h1>Analysis of topology variation in protein family</h2>"

        print >> fpout, "<p><u>Click header labels to sort</u></p>"
        print >> fpout, "<p>"
        print >> fpout, "<br>"
        print >> fpout, "List of pairwise topology comparisons with topology variation" 
        print >> fpout, "<ul>"
        print >> fpout, "<li>"
        print >> fpout, "topology variation type: %s"% g_params['selecttype']
        print >> fpout, "</li>"
        print >> fpout, "<li>"
        print >> fpout, "DG value <= %.1f" %g_params['maxDGvalue']
        print >> fpout, "</li>"
        print >> fpout, "<li>"
        print >> fpout, "Gap percentage of segment aligned to TM helix >= %.1f" %g_params['minGapFraction']
        print >> fpout, "</li>"
        print >> fpout, "<li>"
        print >> fpout, "sequence identity [%.1f, %.1f]" % (
                g_params['minSeqIDT'], g_params['maxSeqIDT'])
        print >> fpout, "</li>"
        print >> fpout, "</ul>"
        print >> fpout, "</p>"
        print >> fpout, "<p>"
        print >> fpout, "%d pairs from %d Pfam families" % (numPair, numFam)
        print >> fpout, "</p>"
        print >> fpout, "<p>"
        print >> fpout, "<ul>"
        print >> fpout, "<li>"
        print >> fpout, "SeqIDT: numIDTres / alignment length"
        print >> fpout, "</li>"
        print >> fpout, "<li>"
        print >> fpout, "SeqIDT1: numIDTres / min(seqlength1, seqlength2)"
        print >> fpout, "</li>"
        print >> fpout, "</ul>"
        print >> fpout, "</p>"
        print >> fpout, "</dir>"

#write right panel        
        print >> fpout, "<dir id=\"Content\">"
        tablename = 'table1'
        tabletitle = ''
        if g_params['htmltableformat'] == 0:
            WriteHTMLTable(tablename, tabletitle,pairCmpRecordList, numPair,
                    numFam, outpath, htmlname, fpout)
        elif g_params['htmltableformat'] == 1:
            WriteHTMLTable2(tablename, tabletitle, pairCmpRecordList, numPair,
                    numFam, outpath, htmlname, fpout)

        print >> fpout, "</dir>"

        WriteHTMLTail(fpout)
        fpout.close()
        print  "Result has been output to %s"%htmlfilename

    except IOError:
        print >> sys.stderr, "%s: Error to write HTML file %s to %s. Exit." %(argv[0], htmlfilename, figuredir)
        return 1
#}}}
def PickOnlyOneForEachPfam(pairCmpRecordList, tup_sorted_by_pfamid): #{{{
    from itertools import groupby
# return the new recordList
# first find the groups with the same pfamid
    pfamidList = [x[1] for x in tup_sorted_by_pfamid]
    uniqList = [ key for key,_ in groupby(pfamidList)]
    numRecordForEachPfamIDList = [pfamidList.count(x) for x in uniqList]

    start = 0
    newList = []
    for iGrp in xrange(len(uniqList)):
        numRd = numRecordForEachPfamIDList[iGrp]
        tupList = []
        for i in xrange(start, start+numRd,1):
            tupList.append((i, pairCmpRecordList[i]['seqidt']))
        sorted_by_seqidt = sorted(tupList, key=lambda tup: tup[1], reverse=True)
        pivot = sorted_by_seqidt[0][0]
        newList.append(pairCmpRecordList[pivot])
        start += numRd
    return newList
#}}}
def IsAnaHasCtermVariation(ana):#{{{
    if ana == {}:
        return False
    if 'Cterm' in ana and ana['Cterm'] != {} > 0:
        return True
    else:
        return False
        #}}}
def IsHasCtermVariation(record):#{{{
    if record == {}:
        return False
    if record['cmpclass'] != 'DIFF':
        return False
    if (IsAnaHasCtermVariation(record['ana1']) or
        IsAnaHasCtermVariation(record['ana2'])):
        return True
    else:
        return False
#}}}
def IsAnaHasNtermVariation(ana):#{{{
    if ana == {}:
        return False
    if 'Nterm' in ana and ana['Nterm'] != {} > 0:
        return True
    else:
        return False
        #}}}
def IsHasNtermVariation(record):#{{{
    if record == {}:
        return False
    if record['cmpclass'] != 'DIFF':
        return False
    if (IsAnaHasNtermVariation(record['ana1']) or
        IsAnaHasNtermVariation(record['ana2'])):
        return True
    else:
        return False
#}}}

def IsAnaHasInternalVariation(ana):#{{{
    if ana == {}:
        return False
    if 'internal' in ana and len(ana['internal']) > 0:
        return True
    else:
        return False
        #}}}
def IsHasInternalVariation(record):#{{{
    if record == {}:
        return False
    if record['cmpclass'] != 'DIFF':
        return False
    if (IsAnaHasInternalVariation(record['ana1']) or
        IsAnaHasInternalVariation(record['ana2'])):
        return True
    else:
        return False
#}}}
def IsIdenticalOrIsoformProtein(seqdef1, seqdef2): 
    seqdef1.lstrip("PREDICTED: ")
    seqdef2.lstrip("PREDICTED: ")
    if seqdef1 == seqdef2:
        return True
    else:
        if (seqdef1.find('isoform') != -1 and seqdef2.find('isoform') != -1 and
            seqdef1[0:10] == seqdef2[0:10]):
            return True
        else:
            return False
def GetSegPos(string, keyC):#{{{
    """
    Get segment of a continue keyC state
    e.g. given a string "0001100022000111100"
    and keyC = '1'
    return [(3,5), (13,17)]
    """
    posList = []
    ex = "(%s+)"%(keyC)
    m = re.finditer(ex,string)
    for i in m:
        posList.append((i.start(0), i.end(0)))
    return posList
#}}}
def FilterSegPos(posList, string, neighbour_char):#{{{
### return only list of "0110"
    newList = []
    N = len(string)
    for (b,e) in posList:
        try:
            if (b>0 and e<N and string[b-1] == neighbour_char and string[e] ==
                    neighbour_char and (e-b)%2 == 0): newList.append((b,e))
        except IndexError:
            pass
#     for (b,e) in posList:
#         if b>0 and string[b-1] != neighbour_char:
#             continue
#         if e < N-1 and string[e] != neighbour_char:
#             continue
#         if (e-b)%2 != 0: # must be even number
#             continue
#         newList.append((b,e))
    return newList
#}}}
def IsHavingInternalTM2GAP(mapArrayList):
    for mapArray in mapArrayList:
        str_maparray_list = ["%d"%x for x in mapArray]
        str_maparray = "".join(str_maparray_list)
        neighbour_char = "0"
        for st in [1]: # TM2GAP
            posContList = GetSegPos(str_maparray, "%d"%st)
            posContList = FilterSegPos(posContList, str_maparray, neighbour_char)
            if len(posContList) > 0:
                return True
    return False

def FilterPairCmpResult(recordList):#{{{
    """
    Filter paircmp result by g_params
    return newList
    """
    newList = []
    for record in recordList:
        if record == {}:
            continue

        cmpclass = record['cmpclass']
        if not cmpclass in  ['TM2GAP', 'TM2GAP_AND_TM2SEQ']:
            continue

        seqidt = lcmp.GetSeqIDT(record, g_params['seqidttype'])
        if (seqidt < g_params['minSeqIDT'] or seqidt > g_params['maxSeqIDT']):
            continue
        #print record['rlty']
        if (record['rlty'] < g_params['min_rlty'] ):
            continue

        id1 = record['id1']
        id2 = record['id2']
        #print  id1, id2, record['mapArray'][0], record['mapArray'][1]
        if IsHavingInternalTM2GAP(record['mapArray']):
            newList.append(record)
    return newList
#}}}
def WritePairList(recordList, outfile):
    try:
        fpout = open(outfile, "w")
        for i in xrange(len(recordList)):
            try: 
                rd = recordList[i]
                fpout.write("%s %s\n"%(rd['id1'], rd['id2']))
            except (IOError, KeyError):
                pass
        fpout.close()
        print "Result output to %s"%outfile
    except IOError:
        print >> sys.stderr, "Failed to write to file %s"%outfile
        return 1

def WriteResult(recordList, outfile):
    try:
        fpout = open(outfile, "w")
        for i in xrange(len(recordList)):
            try: 
                rd = recordList[i]
                fpout.write("\nRecord %d\n"%(i+1))
                fpout.write("%s: %s - %s %6.1f %6.1f\n"%(rd['famid'], rd['id1'],
                    rd['id2'], rd['seqidt'], rd['rlty']))
                for j in xrange(2):
                    fpout.write("%s\n"%(rd['seqanno%d'%(j+1)]))
                fpout.write("SeqLength: %d - %d \n"%(rd['seqLength1'], rd['seqLength2']))
                for j in xrange(2):
                    ss = "".join(["%d"%x for x in rd['mapArray'][j]])
                    fpout.write("MapArray%d: %s \n"%(j+1, ss))
            except (IOError, KeyError, IndexError):
                pass
        fpout.close()
        print "Result output to %s"%outfile
    except IOError:
        print >> sys.stderr, "Failed to write to file %s"%outfile
        return 1

def main(g_params):#{{{
    argv = sys.argv
    numArgv=len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1
    i = 1
    isNonOptionArg=False
    isPickOne = False
    paircmpFile = ""
    pfamACDEListFile = '/data3/data/pfam/pfamA.seed.ac-delist'
    pfamACDEListFile = '/data3/data/pfam/pfam26.0/Pfam-A.clans.tsv'
    seqDefFile = '/data3/wk/MPTopo/pfamAna/pfam2-selTM-giid-refseqid-pfamid-description.txt'
    rltyFile = "/data3/wk/MPTopo/pfamAna_refpro/pred_topcons/refpro20120604-celluar.selmaxlength-m1.topcons.result_TOPCONS.rlty"

    pairlistwithfamidFile = "/data3/wk/MPTopo/pfamAna_refpro/cellular_filter_fragment/pairwise/withinClan/Pfam-A-full.perTM75_nseq20.nr100.filter.fragmented.withinclan.max30000.kalignP_filterNonTM.pairlistwithclanid"
    idwithannoFile = "/data3/wk/MPTopo/pfamAna_refpro/cellular_filter_fragment/pairwise/withinClan/Pfam-A-full.perTM75_nseq20.nr100.filter.fragmented.uniq.idwithanno"

    outpath = ""
    htmlname = 'index'
    tableinfoFile = "Pfam-A-full.perTM75_nseq20.nr100.filter.fragmented.withinclan.max30000.kalignP_filterNonTM.tableinfo"
    while i < numArgv:#{{{
        if isNonOptionArg == True:
            paircmpFile = argv[i]
            isNonOptionArg=False
            i += 1
        elif argv[i] == "--":
            isNonOptionArg=True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] ==  "-h" or  argv[i] == "--help":
                PrintHelp()
                return 1
            elif argv[i] == "-outpath" or argv[i] == "--outpath":
                outpath = argv[i+1]
                i += 2
            elif argv[i] == "-htmlname" or argv[i] == "--htmlname":
                htmlname = argv[i+1]
                i += 2
            elif argv[i] == "-alnfigpath" or argv[i] == "--alnfigpath":
                g_params['alnFigPath'] = argv[i+1]
                i += 2
            elif argv[i] == "-pfamdef" or argv[i] == "--pfamdef":
                pfamACDEListFile = argv[i+1]
                i += 2
            elif argv[i] == "-msapath" or argv[i] == "--msapath":
                g_params['MSAPath'] = argv[i+1]
                i += 2
            elif argv[i] == "-msapath2" or argv[i] == "--msapath2":
                g_params['MSAPath2'] = argv[i+1]
                i += 2
            elif argv[i] == "-tableinfo" or argv[i] == "--tableinfo":
                tableinfoFile = argv[i+1]
                i += 2
            elif argv[i] == "-seqdef" or argv[i] == "--seqdef":
                seqDefFile = argv[i+1]
                i += 2
            elif argv[i] == "-gap" or argv[i] == "--gap":
                g_params['minGapFraction'] = float(argv[i+1])
                i += 2
            elif argv[i] == "-dg" or argv[i] == "--dg":
                g_params['maxDGvalue'] = float(argv[i+1])
                i += 2
            elif argv[i] in ["-min-seqidt", "--min-seqidt"]:
                g_params['minSeqIDT'] = float(argv[i+1])
                i += 2
            elif argv[i] in ["-max-seqidt", "--max-seqidt"]:
                g_params['maxSeqIDT'] = float(argv[i+1])
                i += 2
            elif argv[i] in ["-tableformat", "--tableformat"]:
                g_params['htmltableformat'] = int(argv[i+1])
                i += 2
            elif argv[i] in ["-treepath", "--treepath"]:
                g_params['treepath'] = argv[i+1]
                i += 2
            elif argv[i] in ["-ordermsapath", "--ordermsapath"]:
                g_params['ordermsapath'] = argv[i+1]
                i += 2
            elif argv[i] in ["-topomsapath", "--topomsapath"]:
                g_params['topomsapath'] = argv[i+1]
                i += 2
            elif argv[i] in ["-type", "--type"]:
                g_params['selecttype'] = argv[i+1]
                i += 2
            elif argv[i] in ["-filter-predseq", "--filter-predseq"]:
                if argv[i+1].lower()[0] == 'y':
                    g_params['isFilterPredictedSeq'] = True
                else:
                    g_params['isFilterPredictedSeq'] = False
                i += 2
            elif argv[i] == "-q":
                isQuiet = True
                i += 1
            elif argv[i] in ["-pickone", "--pickone"]:
                isPickOne = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            paircmpFile = argv[i]
            i += 1
#}}}
    g_params['outpath'] = outpath 
    if outpath != "" and not os.path.exists(outpath):
        os.system("mkdir -p %s"%outpath)

    rltyDict = ReadInRLTYInfo(rltyFile)
    pair2famidDict = ReadPair2FamIDInfo(pairlistwithfamidFile)
    seqannoDict = ReadIDWithAnnoInfo(idwithannoFile)

# read paircmprecordlist
    fpin = open(paircmpFile,'r')
    buff = fpin.read()
    fpin.close()
    recordList = []
    unprocessedBuffer = lcmp.ReadPairCmpResultFromBuffer(buff, recordList)
    print "len(recordList) =", len(recordList)

    seqIDListSet = set([])
    for record in recordList:
        seqIDListSet.add(record['id1'])
        seqIDListSet.add(record['id2'])
    print "len(seqIDListSet) =", len(seqIDListSet)

# Read In pairwise alignment info
    pairalnStat = {}
    if tableinfoFile != "" and os.path.exists(tableinfoFile):
        pairalnStat = ReadInTableInfo(tableinfoFile)
        print "len(pairalnStat) =", len(pairalnStat)

    famDefDict = ReadFamDefFile(pfamACDEListFile)
    #print len(rltyDict)
    #print rltyDict['Q33784']
    AppendRLTYInfo(recordList, rltyDict)
    AppendFamInfo(recordList, pair2famidDict, famDefDict)
    AppendSeqAnnoInfo(recordList,seqannoDict )
    AddTableInfo(recordList, pairalnStat)



# #Read In pfamDefList
#     if not os.path.exists(pfamACDEListFile):
#         print >> sys.stderr, "Error! file pfamACDEListFile (%s) does not exist." %pfamACDEListFile
#         return 1
#     pfamDefDict = ReadPfamDEList(pfamACDEListFile)
#     print 'len(pfamDefDict)=', len(pfamDefDict)

# #Read in seqinfoList
#     seqInfoDict = {}
#     if not os.path.exists(seqDefFile):
#         print >> sys.stderr, "Error! file seqDefFile (%s) does not exist." %seqDefFile
#         return 1
#     fpin = open(seqDefFile,"r")
#     line = fpin.readline()
#     line = fpin.readline()
#     while line:
#         strs = line.split('|')
#         if len(strs) == 4:
#             gid = strs[0].strip()
#             if gid in seqIDListSet:
#                 refseqid = strs[1].strip()
#                 pfamid = strs[2].strip()
#                 seqdef = strs[3].strip()
#                 seqInfoDict[gid] = {}
#                 seqInfoDict[gid]['pfamid'] = pfamid 
#                 seqInfoDict[gid]['refseqid'] = refseqid 
#                 seqInfoDict[gid]['seqdef'] = seqdef 
#                 seqInfoDict[gid]['pfamdef'] = pfamDefDict[pfamid] 
#         line = fpin.readline()
#     fpin.close()

#     print 'len(seqInfoDict)=', len(seqInfoDict)

# # add tableinfo to record list
#     print "Add pairwise alignment table info to record..."
#     AddTableInfo(recordList, pairalnStat)
#     print "Add seqdef to record ..."
#     AddSeqDefInfo(recordList, seqInfoDict)

    recordList = FilterPairCmpResult(recordList)

    numFilteredRecordList = len(recordList)
    print "numFilteredRecordList = %d"%numFilteredRecordList 

#     # reorder list according to pfamid
#     numPair = len(filteredRecordList)
#     tupList = [] # list of (index - pfamid)
#     for i in xrange(numPair):
#         thisPfamid = seqInfoDict[filteredRecordList[i]['id1']]['pfamid']
#         tupList.append((i,thisPfamid))
#     sorted_by_pfamid = sorted(tupList, key=lambda tup: tup[1])
#     
#     pairCmpRecordList = []
#     for i in xrange(numPair):
#         pairCmpRecordList.append(filteredRecordList[sorted_by_pfamid[i][0]])
# 
#     if isPickOne:
#         pairCmpRecordList = PickOnlyOneForEachPfam(pairCmpRecordList, sorted_by_pfamid)
# 
#     print "len(pairCmpRecordList) = ", len(pairCmpRecordList)
# 
    out_pairlistfile = outpath + os.sep + "result.pairlist"
    text_outfile = outpath + os.sep + "result.txt"

    WritePairList(recordList, out_pairlistfile)
    WriteResult(recordList, text_outfile)
    GeneratePairAlnFigure(outpath, out_pairlistfile)

    WriteHTML(recordList, htmlname, outpath)

    return 0
#}}}

if __name__ == '__main__' :
    g_params = {}
    g_params['selecttype'] = 'all'
    g_params['htmltableformat'] = 1
    g_params['topomsapath'] = ""
    g_params['treepath'] = ""
    g_params['ordermsapath'] = ""
    g_params['MSAPath2'] = ""
    g_params['alnFigPath'] = ""
    g_params['MSAPath'] = ""
    g_params['outpath'] = ""

    g_params['min_rlty'] = 85.0
    g_params['minGapFraction'] = 0.5
    g_params['maxGapFraction'] = 1.0
    g_params['minDGvalue'] = -999999.0
    g_params['maxDGvalue'] = 0.5
    g_params['minSeqIDT'] = 50.0
    g_params['maxSeqIDT'] = 95.0
    g_params['isShowProgress'] = True
    g_params['isFilterPredictedSeq'] = False
    g_params['seqidttype'] = 1
    sys.exit(main(g_params))

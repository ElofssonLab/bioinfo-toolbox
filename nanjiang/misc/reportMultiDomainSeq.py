#!/usr/bin/env python
import os
import sys
import re
import myfunc
import libtopologycmp as lcmp

rundir = os.path.dirname(sys.argv[0])
binpath = rundir

usage="""
Usage:  %s seqid2fam_map_file [-outpath DIR] [-o OUTFILE]

Description: Report multi-domain sequence information

Options: 
  -outpath DIR       set ouput path, (default: ./)
  -htmlname STR      set the name of the main html file, (default: index)
  -pfamdef FILE      pfam definition file, (default:
                     /data3/data/pfam/pfamA.seed.ac-delist)  
  -seqdef  FILE      sequence definition file
                     (default: /data3/wk/MPTopo/pfamAna/pfam2-selTM-giid-refseqid-pfamid-description.txt)
  -seqlen FILE       sequence length file
  -topodb FILE       set topology database
  -seqdb  FILE       set sequence database
  -pfamscan  FILE    set pfamscan file
  -q                 quiet mode
  -h, --help         print this help message and exit

Created 2013-03-12, updated 2013-03-12, Nanjiang Shu

Examples:
"""

BLOCK_SIZE = myfunc.BLOCK_SIZE

def PrintHelp():
    print usage

def ExtendCoverage(bAll, eAll, b, e):
    return (min(bAll, b), max(eAll, e))
def GetCoveredTM(segment, posTM):#{{{
    newPosTM = []
    indexList = []
    cnt = 0
    for (b,e) in posTM:
        x = myfunc.coverage(b,e,segment[0], segment[1])
        if x/float(e-b) > 0.75:
            newPosTM.append((b,e))
            indexList.append(cnt)
        cnt += 1
    return (newPosTM, indexList)
#}}}
def GroupPfamScanDict(pfamscanDict):#{{{
    """
    Group PfamScan Dict, the new structure is
    dt[seqid] = {'pfamid':data, clanid:data}
    data = {'def': bla, 'alnBeg': 15, 'alnEnd': 150}
    """
    newDict = {}
    for seqid in pfamscanDict:
        hitList = pfamscanDict[seqid]
        newDict[seqid] = {}
        ndt = newDict[seqid]
        for hit in hitList:
            pfamid = hit['pfamid']
            clanid = hit['clanid']
            for item in [pfamid, clanid]:
                if not item in ndt:
                    ndt[item] = {}
                    ndt[item]['def'] = hit['pfamname']
                    ndt[item]['alnBeg'] = hit['alnBeg']
                    ndt[item]['alnEnd'] = hit['alnEnd']
                ndt[item]['alnBeg'], ndt[item]['alnEnd'] = ExtendCoverage(
                        ndt[item]['alnBeg'], ndt[item]['alnEnd'], 
                        hit['alnBeg'], hit['alnEnd'])
    return newDict
#}}}
def GetTopoDict(topoDB, idList):
    hdl = myfunc.MyDB(topoDB)
    if hdl.failure:
        return {}
    dt = {}
    for seqid in idList:
        data = hdl.GetRecord(seqid)
        if data:
            (tmp_id, tmp_anno, tmp_seq) = myfunc.ExtractFromSeqWithAnno(data)
            dt[seqid] = tmp_seq
    hdl.close()
    return dt


def ReadIDWithAnnoInfo(infile):#{{{
    hdl = myfunc.ReadLineByBlock(infile)
    if hdl.failure:
        return {}
    dt = {}
    lines = hdl.readlines()
    while lines != None:
        for line in lines:
            if not line or line[0] == "#":
                continue
            strs = line.split("\t")
            if len(strs) == 2:
                seqid = strs[0]
                anno = strs[1].strip()
                dt[seqid] = anno
        lines = hdl.readlines()
    hdl.close()
    return dt
#}}}
def ReadSeqLengthDict(infile):#{{{
    hdl = myfunc.ReadLineByBlock(infile)
    if hdl.failure:
        return {}
    dt = {}
    lines = hdl.readlines()
    while lines != None:
        for line in lines:
            if not line or line[0] == "#":
                continue
            strs = line.split()
            if len(strs) == 2:
                seqid = strs[0]
                length = int(strs[1])
                dt[seqid] = length
        lines = hdl.readlines()
    hdl.close()
    return dt
#}}}
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

def WriteFastaFile(idList, seqDict, outfile):
    try:
        fpout = open(outfile, "w")
        for seqid in idList:
            try:
                seq = seqDict[seqid]
                fpout.write(">%s\n"%(seqid))
                fpout.write("%s\n"%(seq))
            except KeyError:
                print >> sys.stderr, "%s not found in seqDict"%(seqid)
                pass
        fpout.close()
        return 0
    except IOError:
        print >> sys.stderr, "Failed to write to file %s"%(outfile)
        return 1

def MakeAlignment(groupList, seqDict, topoDict, outpath):#{{{
    cnt = 0
    datapath = "%s%s%s"%(outpath, os.sep, "data")
    if not os.path.exists(datapath):
        os.makedirs(datapath)
    for tup in groupList:
        ss = tup[0]
        famidlist = ss.split("\t")
        ss = "_".join(famidlist)
        seqidlist = tup[2]
        aaseqfile = "%s%s%s.fa"%(datapath, os.sep, ss)
        status = WriteFastaFile(seqidlist, seqDict, aaseqfile)
        if status != 0:
            continue
        mfafile = "%s%s%s.mfa"%(datapath, os.sep, ss)
        cmd = "%s/program/run_kalignP/run_kalignP.sh %s -f fasta -q -o %s" %(
                os.environ['DATADIR3'], aaseqfile, mfafile)
        if not os.path.exists(mfafile) or os.path.getsize(mfafile) <= 0:
            os.system(cmd)
        if not os.path.exists(mfafile):
            continue

        topofile = "%s%s%s.topo"%(datapath, os.sep, ss)
        status = WriteFastaFile(seqidlist, topoDict, topofile)
        if status != 0:
            continue
# match topo
        topomsafile = "%s%s%s.topomsa"%(datapath, os.sep, ss)
        cmd = "%s/wk/MPTopo/bin/matchMSAtopo -msa %s -topo %s -o %s" %(
                os.environ['DATADIR3'], mfafile, topofile, topomsafile)

        if not os.path.exists(topomsafile) or os.path.getsize(topomsafile) <= 0:
            os.system(cmd)

        fasttreefile = "%s%s%s.kalignp.fasttree"%(datapath, os.sep, ss)
        if os.path.exists(mfafile):
            cmd = "FastTree %s > %s"%(mfafile, fasttreefile)
            if not os.path.exists(fasttreefile) or os.path.getsize(fasttreefile) <= 0:
                os.system(cmd)


        orderlistfile="%s%s%s-listorder.txt"%(datapath, os.sep, ss)
        if os.path.exists(fasttreefile):
            cmd = "%s/wk/MPTopo/src/itol_get_tree_listorder.py -datapath %s -outpath %s %s" %(
                    os.environ['DATADIR3'], datapath, datapath, ss)
            if  not os.path.exists(orderlistfile) or os.path.getsize(orderlistfile) <= 0:
                os.system(cmd)

        reorderedtopomsafile = "%s%s%s.reordered.mfa"%(datapath, os.sep, ss)
        if os.path.exists(orderlistfile):
            cmd = "%s/wk/MPTopo/src/reordermsa.py -msafile %s -orderlist %s -o %s" %(
                    os.environ['DATADIR3'], topomsafile, orderlistfile, reorderedtopomsafile)
            if  not os.path.exists(reorderedtopomsafile) or os.path.getsize(reorderedtopomsafile) <= 0:
                os.system(cmd)

        if os.path.exists(reorderedtopomsafile):
            cmd = "%s/wk/MPTopo/src/drawMSATopo.py -text y %s -outpath %s" %(
                    os.environ['DATADIR3'], reorderedtopomsafile, datapath)
            os.system(cmd)

#}}}
def WriteInfo(groupList, seqlenDict, seqannoDict, pfamidDefDict,#{{{
        clanidDefDict, topoDict, groupedPfamScanDict, htmlname, fpout):
    cnt = 0
    for tup in groupList:
        try:
            ss = tup[0]
            seqidlist = tup[2]
            famidlist = ss.split("\t")
            fpout.write("Group %d: %d seqs, %d domains "%(cnt+1, tup[1],
                len(famidlist)))
            for famid in famidlist:
                if famid[0] == 'P':
                    famdef = pfamidDefDict[famid]
                else:
                    famdef = clanidDefDict[famid]
                fpout.write(" %s (%s)"%(famid, famdef))
            fpout.write("\n")

            fpout.write("#%-3s %10s %4s %3s %15s %5s\n"%("No",
                "SeqID", "Len", "nTM", "DomainCoverage", "nTM_within"))
            cntseq = 0
            for seqid in seqidlist:
                try:
                    seqlen = seqlenDict[seqid]
                except KeyError:
                    seqlen = -1
                    pass

                fpout.write("%-4d %10s %4d"%(cntseq+1, seqid, seqlen))

                try:
                    topo = topoDict[seqid]
                    posTM = myfunc.GetTMPosition(topo)
                except KeyError:
                    print >> sys.stderr, "topo not found for %s"%seqid
                    fpout.write("\n")
                    continue

                fpout.write(" %3d"%(len(posTM)))

                pfamscan_hit = groupedPfamScanDict[seqid]
                for famid in famidlist:
                    try:
                        b1 = pfamscan_hit[famid]['alnBeg']
                        e1 = pfamscan_hit[famid]['alnEnd']
                        (posTM_covered, indexList_covered) = GetCoveredTM((b1,e1), posTM)
                        fpout.write("%15s %5s %4s"%("(%d,%d)"%(b1,e1), 
                            "%d TM"%(len(posTM_covered)),
                            "%d-%d"%(indexList_covered[0]+1,
                                indexList_covered[len(indexList_covered)-1]+1)
                            ))
                    except (KeyError):
                        print >> sys.stderr, "%s not in pfamscan_hit"%(famid)
                        pass
                fpout.write("\n")
                cntseq += 1
        except (KeyError, IndexError):
            print >> sys.stderr, "Error for %s"%(tup[0])
            pass
        cnt += 1
        fpout.write("\n")
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
def WriteHTML(seqid2pfamidDict, seqlenDict, seqannoDict, pfamidDefDict, #{{{
        clanidDefDict, topoDict, htmlname, outpath):
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
        title="Proteins with Multiple TM Domains"
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

        seqidt = GetSeqIDT(record, g_params['seqidttype'])
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
    infile = ""
    pfamDefFile = '/data3/data/pfam/pfam26.0/Pfam-A.clans.tsv'
    seqDefFile = '/data3/wk/MPTopo/pfamAna/pfam2-selTM-giid-refseqid-pfamid-description.txt'
    idwithannoFile = "/data3/wk/MPTopo/pfamAna_refpro/pfammap_from_uniprot/Pfam-A-full.perTM75_nseq20.nr100.filter.fragmented.uniq.idwithanno"
    seqLengthFile = "/data3/wk/MPTopo/pfamAna_refpro/pfammap_from_uniprot/Pfam-A-full.perTM75_nseq20.nr100.filter.fragmented.uniq.seqlen"
    topoDB = "/data3/wk/MPTopo/pfamAna_refpro/pred_topcons/refpro20120604-celluar.selmaxlength-m1.topcons.result_TOPCONS.topo"
    seqDB = "/data3/wk/MPTopo/pfamAna_refpro/cellular_filter_fragment/Pfam-A-full.perTM75_nseq20.nr100.filter.fragmented.uniq"
    pfamscanFile = "/data3/wk/MPTopo/pfamAna_refpro/result_pfamscan/Pfam-A-full.perTM75_nseq20.nr100.filter.fragmented.pfamscan"

    outpath = ""
    outfile = ""
    htmlname = 'index'
    while i < numArgv:#{{{
        if isNonOptionArg == True:
            infile = argv[i]
            isNonOptionArg=False
            i += 1
        elif argv[i] == "--":
            isNonOptionArg=True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] ==  "-h" or  argv[i] == "--help":
                PrintHelp()
                return 1
            elif argv[i] in ["-outpath", "--outpath"]:
                outpath = argv[i+1]
                i += 2
            elif argv[i] in ["-o", "--o"]:
                (outfile, i) = myfunc.my_getopt_str(argv,i)
            elif argv[i] in ["-htmlname", "--htmlname"]:
                htmlname = argv[i+1]
                i += 2
            elif argv[i] in ["-seqlen", "--seqlen"]:
                seqLengthFile = argv[i+1]
                i += 2
            elif argv[i] in ["-topodb", "--topodb"]:
                topoDB = argv[i+1]
                i += 2
            elif argv[i] in ["-seqdb", "--seqdb"]:
                seqDB = argv[i+1]
                i += 2
            elif argv[i] in [ "-pfamdef", "--pfamdef"]:
                pfamDefFile = argv[i+1]
                i += 2
            elif argv[i] in ["-pfamscan", "--pfamscan"]:
                pfamscanFile = argv[i+1]
                i += 2
            elif argv[i] in ["-seqdef", "--seqdef"]:
                seqDefFile = argv[i+1]
                i += 2
            elif argv[i] in ["-q"]:
                isQuiet = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            infile = argv[i]
            i += 1
#}}}
    g_params['outpath'] = outpath 
    if outpath == "":
        print >> sys.stderr, "outpath not set"
        return 1
    elif outpath != "" and not os.path.exists(outpath):
        os.makedirs(outpath)

    if infile == "":
        print >> sys.stderr, "infile not set"
        return 1

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)

    seqid2pfamidDict = myfunc.ReadFam2SeqidMap(infile)
# group sequences by domains
    groupDict = {}
    for seqid in seqid2pfamidDict:
        famlist = seqid2pfamidDict[seqid]
        ss = "\t".join(famlist)
        if not ss in groupDict:
            groupDict[ss] = []
        groupDict[ss].append(seqid)
    groupList = []
    for ss in groupDict:
        groupList.append((ss, len(groupDict[ss]), groupDict[ss]))
    groupList = sorted(groupList, key=lambda x:x[1], reverse=True)

    seqlenDict = ReadSeqLengthDict(seqLengthFile)
    seqannoDict = ReadIDWithAnnoInfo(idwithannoFile)
    (pfamidDefDict, clanidDefDict) = lcmp.ReadPfamDefFile(pfamDefFile)

    topoDict = GetTopoDict(topoDB, seqid2pfamidDict.keys())
    seqDict = GetTopoDict(seqDB, seqid2pfamidDict.keys())

    pfamScanDict = myfunc.ReadPfamScan(pfamscanFile)
    groupedPfamScanDict = GroupPfamScanDict(pfamScanDict)

#     WriteHTML(seqid2pfamidDict, seqlenDict, seqannoDict, pfamidDefDict,
#             clanidDefDict, topoDict, htmlname, outpath)
    WriteInfo(groupList, seqlenDict, seqannoDict, pfamidDefDict,
            clanidDefDict, topoDict, groupedPfamScanDict, htmlname, fpout)
    MakeAlignment(groupList, seqDict, topoDict, outpath)
    myfunc.myclose(fpout)
    return 0
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['htmltableformat'] = 1
    g_params['topoDB'] = ""
    g_params['outpath'] = ""
    g_params['isShowProgress'] = True
    return g_params
#}}}

if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

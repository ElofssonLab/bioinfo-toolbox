#!/usr/bin/env python
# Description:
import os
import sys
import myfunc


usage="""

USAGE: sel-gram+-.py tableinfofile gram+-_deffile FORMAT

FORMAT, set the format of the tableinfofile (default: 0)
        0: tableinfo output by uniprottrembldata2table.py
        1: seqid[TAB]orgname

Created 2013-10-08, updated 2014-10-08, Nanjiang Shu
"""
def ReadListFromUniprotTableInfo(infile):#{{{
    try:
        recordList = []
        fpin = open(infile, "r")
        line = fpin.readline()
        while line:
            if line != "" and line[0] != "#":
                strs = line.split("\t")
                if len(strs) >= 6:
                    seqidList = strs[0].split(";")
                    seqidList = filter(None, seqidList)
                    seqid = seqidList[0]
                    if seqid == "":
                        print >> sys.stderr, "Bad record, line =", line
                        continue
                    rd = {}
                    rd['seqid'] = seqid
                    organism = strs[3].strip()
                    organism = organism.split("(")[0] #added 2014-09-10
                    organism = organism.strip()
                    rd['organism'] = organism
                    recordList.append(rd)
            line = fpin.readline()
        fpin.close()
        return recordList
    except IOError:
        print >> sys.stderr, "Failed to read file %s" %infile
        return []
#}}}
def ReadSeqID_ORGNAME(infile):#{{{
    try:
        recordList = []
        fpin = open(infile, "r")
        lines = fpin.readlines()
        fpin.close()
        for line in lines:
            if line != "" and line[0] != "#":
                strs = line.split("\t")
                if len(strs) >= 2:
                    seqid = strs[0].strip()
                    orgname = strs[1].strip()
                    rd = {}
                    rd['seqid'] = seqid
                    rd['organism'] = orgname
                    recordList.append(rd)
        return recordList
    except IOError:
        print >> sys.stderr, "Failed to read infile %s"%(infile)
        return []
#}}}

def ReadGramPNDef(infile):#{{{
    try:
        gramPNDict = {}
        fpin = open(infile, "r")
        line = fpin.readline()
        while line:
            if line != "" and line[0] != "#":
                strs = line.split("\t")
                if len(strs) >= 2:
                    organism = strs[0].strip().upper()
                    grampn = strs[1].strip()
                    gramPNDict[organism] = grampn
            line = fpin.readline()
        fpin.close()
        return gramPNDict
    except IOError:
        print >> sys.stderr, "Failed to read file %s" %infile
        return {}
#}}}
#grampn_deffile = "/data3/wk/MPTopo/pfamAna_refpro/bacterialist.csv" 
#grampn_deffile = "/data3/wk/MPTopo/pfamAna_swissprot/bacteria_gram+-_list_from_kostas.txt"
#tabinfofile = "/data3/data/uniprot/reference_proteome/refpro20120604-Bacteria.tableinfo"

try:
    tabinfofile = sys.argv[1]
except IndexError:
    print usage
    sys.exit(1)

try:
    grampn_deffile = sys.argv[2]
except IndexError:
    print usage
    sys.exit(1)

try:
    table_format = int(sys.argv[3])
except IndexError, TypeError:
    table_format = 0



grampnDefDict = ReadGramPNDef(grampn_deffile)

# for dt in grampnDefDict:
#     print dt, grampnDefDict[dt]
#print grampnDefDict
if table_format == 0:
    recordList = ReadListFromUniprotTableInfo(tabinfofile)
elif table_format == 1:
    recordList = ReadSeqID_ORGNAME(tabinfofile)
else:
    print >> sys.stderr, "Bad table_format = ", table_format
    sys.exit(1)
#print tableInfoDict


rootname = os.path.basename(os.path.splitext(tabinfofile)[0]) 
outpath = os.path.dirname(tabinfofile)
if outpath == "":
    outpath = "."
idListFile_gram_p = outpath + os.sep + rootname + ".gram+.idlist"
idListFile_gram_n = outpath + os.sep + rootname + ".gram-.idlist"
idListFile_gram_unknown = outpath + os.sep + rootname + ".gramunknown.idlist"

fp_gram_p = open(idListFile_gram_p, "w")
fp_gram_n = open(idListFile_gram_n, "w")
fp_gram_unknown = open(idListFile_gram_unknown, "w")

for rd in recordList:
    organism = rd['organism']
    grampn = ""
    if organism.upper() in grampnDefDict:
        grampn = grampnDefDict[organism]
#         sys.stdout.write("<%s>"%grampn)
    else:
        sys.stdout.write("<%s>\n"%organism)

    if grampn == "n":
        fp_gram_n.write("%s\n"%rd['seqid'])
    elif grampn == "p":
        fp_gram_p.write("%s\n"%rd['seqid'])
    else:
        fp_gram_unknown.write("%s\n"%rd['seqid'])

fp_gram_p.close()
fp_gram_n.close()
fp_gram_unknown.close()

print "Result output to"
print "\t", idListFile_gram_p
print "\t", idListFile_gram_n
print "\t", idListFile_gram_unknown


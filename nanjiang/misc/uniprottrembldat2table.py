#!/usr/bin/env python
# Description:  Extract data from uniprot_trembl.dat
# ChangeLog 2014-08-28 
#   add options
#       -keep_no_genename     Keep proteins without gene name
#       -keep_non_refpro      Keep protein not in reference proteome
#       -keep_isoform         Keep proteins in isoforms
# ChangeLog 2014-08-29
#    taxonomic_class "Bacteria." is considered as "Bacteria", for e.g. O34002,
#    P07472, Q54506
# 
import os
import sys
import myfunc
usage = """
usage:  uniprottrembldat2table.py -i uniprot_trembl.dat [-o OUTFILE] 
Description: Extract data from uniprot_trembl.dat
             and output table with the format
             # AC Length GN OS OC
             tab delimited
Options:
  -q              Quiet mode
  -oc STR         Restrict to taxonomic_class
  -keep_no_genename     Keep proteins without gene name
  -keep_non_refpro      Keep protein not in reference proteome
  -keep_isoform         Keep proteins in isoforms
  -h, --help      Print this help message and exit

Created 2012-06-01, updated 2014-08-29, Nanjiang Shu
"""

def PrintHelp():
    print usage

def FilterRecord(recordList, restrictOCList):#{{{
# Filter records
# 1. without gene name
# 2. without keyword "reference proteome" (KW)
# 3. filter out isoforms (those with ID-names ended with -INT)
    if len(recordList) <= 0:
        return []

    newList = []
    isRestrictOC = False
    if len(restrictOCList) > 0 and restrictOCList[0].lower() != "all":
        isRestrictOC = True

    for rd in recordList:
        if  g_params['filter_no_genename'] and rd['genename'] == "" :
            print >> sys.stderr, "%s NO genaname, ignored" % (rd['accession'])
            continue
        if g_params['filter_non_refpro'] and not rd['isRefPro']:
            print >> sys.stderr, "%s not reference proteome, ignored" % (rd['accession'])
            continue
        if g_params['filter_isoform'] and rd['accession'].find("-") != -1:
            print >> sys.stderr, "%s with isoforms, ignored" % (rd['accession'])
            continue
        if isRestrictOC and (rd['taxonomic_class'] not in restrictOCList):
            print >> sys.stderr, (rd['accession'], "not in", restrictOCList,
                    "ignored.")
            continue
        newList.append(rd)
    return newList
#}}}

def WriteRecord(recordList, fpout):#{{{
    for rd in recordList:
        fpout.write("%s\t%d\t%s\t%s\t%s\t%s\t%d\n"%(
            rd['accession'],
            rd['length'],
            rd['genename'],
            rd['organism'],
            rd['taxonomic_class'],
            ";".join(rd['pfamidList']),
            rd['isRefPro']
            ))
#}}}
def ExtractFromUniprotTremblRecord(recordContent):#{{{
    record = {}
    lines = recordContent.split("\n")
    numLine = len(lines)
    i = 0
    # AC can be multiple lines
    str_accession = ""     # AC
    str_genename = ""      # GN
    str_organism = ""      # OS
    pfamidList = [] # 
    str_keyword = ""
    length = 0        # from ID record
    str_taxonomic_class = "" # OC, e.g. Archaes, Becteria
    for line in lines:
        if len(line) > 2:
            tag = line[0:2]
            if tag == "ID":
                strs = line[5:].split()
                nstrs = len(strs)
                length = int (strs[nstrs-2])
            elif tag == "AC":
                str_accession += line[5:]
            elif tag == "GN":
                str_genename += line[5:]
            elif tag == "OS":
                str_organism += line[5:]
            elif tag == "OC":
                str_taxonomic_class += line[5:]
            elif tag == "KW":
                str_keyword += line[5:]
            elif tag == "DR":
                if line[5:].find("Pfam") == 0:
                    strs = line[5:].split(";")
                    pfamidList.append(strs[1].strip())
            elif tag == "SQ":
                break

    #accession
    accessionList = str_accession.split(";")
    accessionList = filter(None, accessionList)
    accessionList = [x.strip() for x in accessionList]
    accession = ";".join(accessionList)
    # genename:
    strs = str_genename.split(";")
    strs = filter(None, strs)
    li = []
    for ss in strs:
        sp1 = ss.split("=")
        if len(sp1) == 1:
            ac = sp1[0].strip()
        else:
            ac = sp1[1].strip()
        li.append(ac)
    genename = ";".join(li)
    # organism
    organism = str_organism.rstrip(".")
    # taxonomic_class
    taxonomic_class = str_taxonomic_class.split(";")[0]
    taxonomic_class = taxonomic_class.strip(".") # added 2014-08-29, this solved Bacteria. for P07472. 
    isRefPro = False
    if str_keyword.find("Reference proteome") != -1:
        isRefPro = True
    else:
        isRefPro = False

    if accession != "":
        record['accession'] = accession
        record['length'] = length
        record['genename'] = genename
        record['organism'] = organism
        record['taxonomic_class'] = taxonomic_class
        record['isRefPro'] = isRefPro
        record['pfamidList'] = pfamidList
        return record
    else:
        return {}
#}}}
def Read_UniprotTremblData_from_buffer(buff, recordList, isEOFreached):#{{{
    if not buff:
        return ""
    unprocessedBuffer = ""
    beg = 0
    end = 0
    while 1:
        beg=buff.find("ID ",beg)
        if beg >= 0:
            end=buff.find("\n//",beg+1)
            if end >= 0:
                recordContent = buff[beg:end]
                record = ExtractFromUniprotTremblRecord(recordContent)
                if record != {}:
                    recordList.append(record)
                beg = end
            else:
                unprocessedBuffer = buff[beg:]
                break
        else:
            unprocessedBuffer = buff[end:]
            break
    if isEOFreached and unprocessedBuffer:
        recordContent = unprocessedBuffer
        record = ExtractFromUniprotTremblRecord(recordContent)
        if record != {}:
            recordList.append(record)
        unprocessedBuffer = ""
    return unprocessedBuffer
#}}}
def UniprotTremblData2Table(datafile, restrictOCList, fpout):#{{{
    try: 
        fpout.write("#AC\tLength\tGN\tOS\tOC\tPfamID\tisRefPro\n")
        fpin = open(datafile, "r")
        unprocessedBuffer=""
        isEOFreached = False
        while 1:
            buff = fpin.read(BLOCK_SIZE)
            if len(buff) < BLOCK_SIZE:
                isEOFreached = True
            buff = unprocessedBuffer + buff
            recordList = []
            unprocessedBuffer = Read_UniprotTremblData_from_buffer(
                    buff, recordList, isEOFreached)
            if len(recordList) > 0: 
                filteredRecordList = FilterRecord(recordList, restrictOCList)
                if len(filteredRecordList) > 0:
                    WriteRecord(filteredRecordList, fpout)
            if isEOFreached == True:
                break
        fpin.close()

    except IOError:
        print >> sys.stderr, "Failed to read datafile ", datafile
        return 1

#}}}
def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile = ""
    datafile = ""
    restrictOCList = []

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            datafile = argv[i]
            isNonOptionArg = False
            i += 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-o", "--o", "-outfile", "--outfile"]:
                outfile = argv[i+1]
                i += 2
            elif argv[i] in ["-i", "--i"] :
                datafile = argv[i+1]
                i += 2
            elif argv[i] in ["-keep_isoform", "--keep_isoform"] :
                g_params['filter_isoform'] = False
                i += 1
            elif argv[i] in ["-keep_non_refpro", "--keep_non_refpro"] :
                g_params['filter_non_refpro'] = False
                i += 1
            elif argv[i] in ["-keep_no_genename", "--keep_no_genename"] :
                g_params['filter_no_genename'] = False
                i += 1
            elif argv[i] in ["-oc", "--oc"] :
                restrictOCList.append(argv[i+1])
                i += 2
            elif argv[i] in ["-q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            datafile = argv[i]
            i += 1

    if not os.path.exists(datafile):
        print >> sys.stderr, "datafile %s not set or not exists. Exit" %(datafile)
        return 1

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    UniprotTremblData2Table(datafile, restrictOCList, fpout)
    myfunc.myclose(fpout)
    return 0
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['filter_no_genename'] = True
    g_params['filter_non_refpro'] = True
    g_params['filter_isoform'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    BLOCK_SIZE=100000
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

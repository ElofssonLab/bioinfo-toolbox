#!/usr/bin/env python
# Description:  Extract TaxID from uniprot_trembl.dat
# 
import os
import sys
import myfunc
usage = """
USAGE:  uniprottrembldat2taxid.py -i uniprot_trembl.dat [-o OUTFILE] 

Description: Extract data from uniprot_trembl.dat
             # AC TaxID
             tab delimited
Options:
  -q              Quiet mode
  -h, --help      Print this help message and exit

Created 2014-09-08, updated 2014-09-08, Nanjiang Shu 
"""

def PrintHelp():
    print usage


def WriteRecord(recordList, fpout):#{{{
    for rd in recordList:
        fpout.write("%s\t%s\n"%(
            rd['accession'],
            rd['ncbi_taxid'])
            )
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
    str_ncbi_taxid = "" # OX, e.g. OX   NCBI_TaxID=654924;
    for line in lines:
        if len(line) > 2:
            tag = line[0:2]
            if tag == "ID":
                strs = line[5:].split()
                nstrs = len(strs)
                length = int (strs[nstrs-2])
            elif tag == "AC":
                str_accession += line[5:]
            elif tag == "OX":
                str_ncbi_taxid += line[5:]
            elif tag == "SQ":
                break

    #accession
    accessionList = str_accession.split(";")
    accessionList = filter(None, accessionList)
    accessionList = [x.strip() for x in accessionList]
    accession = ";".join(accessionList)

    #taxid
    strs = str_ncbi_taxid.split("=")
    ncbi_taxid = ""
    if len(strs) >= 2:
        ncbi_taxid = strs[1].rstrip(";")

    if accession != "":
        record['accession'] = accession
        record['ncbi_taxid'] = ncbi_taxid
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
        fpout.write("#AC\tNCBI_TaxID\n")
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
                WriteRecord(recordList, fpout)
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

#!/usr/bin/env python
# Description:  Select the the longest seq with the same "gene name" in the
# same organism
# 
# ChangeLog 2014-08-29
#   proteins without gene name are kept. before it was filtered
import os
import sys
import myfunc
usage = """
usage:   selectSeq_by_uniprotinfotable.py FILE [-o OUTFILE]
Description: Select the longest seq with the same "gene name" in the same
             organism
             # input format
             # uniref_seqid uniprot_AC Length GN OS OC
             tab delimited
Options:
  -outsorted  FILE  Write sorted table to file
  -q                Quiet mode
  -method 0|1       genename merge method, (default: 1)
                    0: merge only when gene name match exactly
                    1: merge when one of the merge name match
  -h, --help        Print this help message and exit

Created 2012-06-01, updated 2014-08-29, Nanjiang Shu

Example: selectSeq_by_uniprotinfotable.py  uniref90_refpro_20120521.random.sorted.uniprottableinfo
"""

def PrintHelp():
    print usage

def ReadUniprotInfoTable(infile):#{{{
    try:
        uniprotInfoDict = {}
        fpin = open(infile, "r")
        method_merge_genename = g_params['method_merge_genename']
        for line in fpin:
            line = line.rstrip("\n")
            if line and line[0] != "#":
                strs = line.split("\t")
                accession = strs[0]
                length = int(strs[1])
                genename = strs[2]
                organism = strs[3]
                taxonomic_class = strs[4]

                if taxonomic_class == "":
                    taxonomic_class = "NO_TAXO"
                    print >> sys.stderr, accession, "no taxonomic_class, ignored."
                    continue
                if organism == "":
                    organism = "NO_ORGANISM"
                    print >> sys.stderr, accession, "no organism, ignored."
                    continue
                if genename == "":
                    genename = "NO_GENENAME" + "." + accession.split(";")[0]
#                    print >> sys.stderr, accession, "no genename, ignored."
#                    continue
                genenameSet = set(genename.split(";"))
                genenameSet = set(filter(None, genenameSet))

                rd = {}
                rd['length'] = length
                rd['line'] = line

                if method_merge_genename == 0:
                    if taxonomic_class in uniprotInfoDict:
                        taxoDict = uniprotInfoDict[taxonomic_class]
                        if organism in taxoDict:
                            osDict = taxoDict[organism]
                            if genename in osDict:
                                osDict[genename]['rdList'].append(rd)
                            else:
                                osDict[genename] = {}
                                osDict[genename]['rdList'] = []
                                osDict[genename]['rdList'].append(rd)
                        else:
                            taxoDict[organism] = {}
                            osDict = taxoDict[organism]
                            osDict[genename] = {}
                            osDict[genename]['rdList'] = []
                            osDict[genename]['rdList'].append(rd)
                    else:
                        uniprotInfoDict[taxonomic_class] = {}
                        taxoDict = uniprotInfoDict[taxonomic_class]
                        taxoDict[organism] = {}
                        osDict = taxoDict[organism]
                        osDict[genename] = {}
                        osDict[genename]['rdList'] = []
                        osDict[genename]['rdList'].append(rd)
                elif method_merge_genename == 1:
                    if taxonomic_class in uniprotInfoDict:
                        taxoDict =  uniprotInfoDict[taxonomic_class]
                        if organism in taxoDict:
                            osDict = taxoDict[organism]
                            isGeneNameFound = False
                            for gn in osDict:
                                gnDict = osDict[gn]
                                s1 = genenameSet
                                s2 = gnDict['genenameSet']
                                if len(s1 & s2) > 0:
                                    gnDict['rdList'].append(rd)
                                    isGeneNameFound = True
                                    gnDict['genenameSet'] = s1 | s2
                                    break
                            if not isGeneNameFound:
                                osDict[genename] = {}
                                gnDict = osDict[genename]
                                gnDict['rdList'] = []
                                gnDict['genenameSet'] = genenameSet
                                gnDict['rdList'].append(rd)
                        else:
                            taxoDict[organism] = {}
                            osDict = taxoDict[organism]
                            osDict[genename] = {}
                            gnDict = osDict[genename]
                            gnDict['rdList'] = []
                            gnDict['genenameSet'] = genenameSet
                            gnDict['rdList'].append(rd)
                    else:
                        uniprotInfoDict[taxonomic_class] = {}
                        taxoDict = uniprotInfoDict[taxonomic_class] 
                        taxoDict[organism] = {}
                        osDict = taxoDict[organism] 
                        osDict[genename] = {}
                        gnDict = osDict[genename]
                        gnDict['rdList'] = []
                        gnDict['genenameSet'] = genenameSet
                        gnDict['rdList'].append(rd)
        fpin.close()
        return uniprotInfoDict
    except IOError:
        print >> sys.stderr, "Failed to read file ",infile
        return {}
#}}}

def WriteSortedTable(outsortedtablefile, uniprotInfoDict):#{{{
    try: 
        fpout = open (outsortedtablefile, "w")
        for taxonomic_class in uniprotInfoDict:
            taxoDict = uniprotInfoDict[taxonomic_class]
            if taxoDict == {}:
                continue
            for organism in taxoDict:
                osDict = taxoDict[organism]
                if osDict == {}:
                    continue
                for genename in osDict:
                    gnDict = osDict[genename]
                    rdList = gnDict['rdList']
                    rdList = sorted(rdList, key=lambda rd:rd['length'], reverse=True)
                    for rd in rdList:
                        fpout.write("%s\n"%(rd['line']))
        fpout.close()
        return 0
    except IOError:
        print >> sys.stderr, "Failed to write to file ", outsortedtablefile
        return 1
#}}}
def SelectSeq_by_uniprotinfotable(uniprotInfoDict, fpout):#{{{
    for taxonomic_class in uniprotInfoDict:
        taxoDict = uniprotInfoDict[taxonomic_class]
        if taxoDict == {}:
            continue
        for organism in taxoDict:
            osDict = taxoDict[organism]
            if osDict == {}:
                continue
            for genename in osDict:
                gnDict = osDict[genename]
                rdList = gnDict['rdList']
                if genename == "NO_GENENAME":
                    for rd in rdList:
                        fpout.write("%s\n"%(rd['line']))
                        print >> sys.stderr, "NO_GENENAME", rd['line']
                else:
                    lengthList = []
                    for rd in rdList:
                        lengthList.append(rd['length'])

                    maxLength = max(lengthList)
                    idxMaxLength = lengthList.index(maxLength)
                    rd = rdList[idxMaxLength]
                    fpout.write("%s\n"%(rd['line']))
    return 0
#}}}

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile = ""
    infile = ""
    outsortedtablefile = ""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            infile = argv[i]
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
            elif argv[i] in ["-outsorted", "--outsorted"]:
                outsortedtablefile = argv[i+1]
                i += 2
            elif argv[i] in ["-i", "--i"] :
                infile = argv[i+1]
                i += 2
            elif argv[i] in ["-method", "--method"] :
                g_params['method_merge_genename'] = int(argv[i+1])
                i += 2
            elif argv[i] in ["-q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            infile = argv[i]
            i += 1


    if not os.path.exists(infile):
        print >> sys.stderr, "infile %s not set or not exists. Exit" %(infile)
        return 1
    if not g_params['method_merge_genename'] in [0,1]:
        print >> sys.stderr, "method_merge_genename not set as 0 or 1. Exit"
        return 1

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)

    uniprotInfoDict = ReadUniprotInfoTable(infile)

    SelectSeq_by_uniprotinfotable(uniprotInfoDict, fpout)

    if outsortedtablefile != "":
        WriteSortedTable(outsortedtablefile, uniprotInfoDict)

    myfunc.myclose(fpout)
    return 0
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['method_merge_genename'] = 1
    return g_params
#}}}
if __name__ == '__main__' :
    BLOCK_SIZE=100000
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

#!/usr/bin/env python
# Description:
import os
import sys
import myfunc
usage = """
usage:  sortedTopoMSA2taxcolor.py  infile [-o OUTFILE]
                                    
Description:

Options:
  -tableinfo  FILE  Set taxonomy definition file
  -h, --help         Print this help message and exit

Created 2012-10-19, updated 2012-10-19, Nanjiang Shu 
"""
nodename = os.uname()[1]
if nodename.find("uppmax") != -1:
    datadir = "/bubo/home/h3/nanjiang/glob/data"
elif nodename.find("illergard") != -1:
    datadir = "/data3/data"
else:
    datadir = "/data"

def ReadUniprotInfoTable(infile):#{{{
    try:
        seqid2TaxoDict = {}
        fpin = open(infile, "r")
        line = fpin.readline()
        while line:
            line = line.rstrip("\n")
            if line[0] != "#":
                strs = line.split("\t")
                accession = strs[0]
                length = int(strs[1])
                genename = strs[2]
                organism = strs[3]
                taxonomic_class = strs[4]

                idlist = accession.split(";")

                if taxonomic_class == "":
                    taxonomic_class = "NO_TAXO"
                    print >> sys.stderr, accession, "no taxonomic_class, ignored."
                    continue
                if organism == "":
                    organism = "NO_ORGANISM"
                if genename == "":
                    genename = "NO_GENENAME"
                for idd in idlist:
                    if not idd in seqid2TaxoDict:
                        seqid2TaxoDict[idd] = {}
                    seqid2TaxoDict[idd]['taxonomic_class'] = taxonomic_class
                    seqid2TaxoDict[idd]['length'] = length
                    seqid2TaxoDict[idd]['organism'] = organism
                    seqid2TaxoDict[idd]['genename'] = genename
            line = fpin.readline()
        fpin.close()
        return seqid2TaxoDict
    except IOError:
        print >> sys.stderr, "Failed to read file ",infile
        return {}
#}}}

black = "#000000"
red = "#FF0000"
blue = "#0000FF"
green = "#008000"
yellow = "#FFFF00"
purple = "#800080"
# 

colorpalette = [
"#008000", # Green
"#C0C0C0", #Silver
"#00FF00", #Lime
"#808080", #Gray
"#808000", #Olive
"#FFFFFF", #White
"#FFFF00", #Yellow
"#800000", #Maroon
"#000080", #Navy
"#FF0000", #Red
"#0000FF", #Blue
"#800080", #Purple
"#008080", #Teal
"#FF00FF", #Fuchsia
"#00FFFF" #Aqua
]

taxoColorDict = {}
taxoColorDict['Eukaryota']  = purple
taxoColorDict['Bacteria']  = yellow
taxoColorDict['Archaea']  = blue

#print colorpalette
def PrintHelp():
    print usage
def WriteTaxoColor(idList, seqid2taxoDict, fpout):
    for idd in idList:
        if idd in seqid2taxoDict:
            taxo = seqid2taxoDict[idd]['taxonomic_class']
            if taxo in taxoColorDict:
                color = taxoColorDict[taxo]
            else:
                color = black
            #fpout.write("%s,%s,%s,%s\n"%(idd, "range", color, pfamid))
            fpout.write("%s,%s\n"%(idd, color))
        else:
            #fpout.write("%s,%s,%s,%s\n"%(idd, "range", black, "Null"))
            fpout.write("%s,%s\n"%(idd, black ))
    fpout.write("\n")

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1
    infile = ""
    outfile = ""
    tableinfoFile = datadir + os.sep + "uniprot/reference_proteome/refpro20120604-celluar.selmaxlength-m1.nr100.tableinfo"

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
            elif argv[i] in ["-tableinfo", "--tableinfo"]:
                tableinfoFile = argv[i+1]
                i += 2
            elif argv[i] in ["-o", "--o"]:
                outfile = argv[i+1]
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
    if infile == "" or not os.path.exists(infile):
        print >> sys.stderr, "Error. Infile not set. exit"
        return 1
    if tableinfoFile == "" or not os.path.exists(tableinfoFile):
        print >> sys.stderr, "Error. tableinfoFile %s does not exist. exit" %(tableinfoFile)
        return 1
    seqid2TaxoDict = ReadUniprotInfoTable(tableinfoFile)
    if tableinfoFile == {}:
        print >> sys.stderr, "Read tableinfoFile failed."
        return 1

    (idList, topoList) = myfunc.ReadFasta_without_annotation(infile)
    if len(idList) < 1:
        print >> sys.stderr, "Read infile failed."
        return 1


    idList.remove("Consensus")
    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)

    WriteTaxoColor(idList, seqid2TaxoDict, fpout)

    myfunc.myclose(fpout)

    return 0

#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

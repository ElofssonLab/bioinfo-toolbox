#!/usr/bin/env python
# Description:
import os
import sys
import myfunc
usage = """
usage:  sortedTopoMSA2pfamcolor.py  infile [-o OUTFILE]
                                    
Description:

Options:
  -seqid2pfam  FILE  Set mapfile
  -h, --help         Print this help message and exit

Created 2012-10-19, updated 2012-10-22, Nanjiang Shu 
"""
nodename = os.uname()[1]
if nodename.find("uppmax") != -1:
    datadir3 = "/bubo/home/h3/nanjiang/glob"
else:
    datadir3 = "/data3"


black = "#000000"
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

#print colorpalette
def PrintHelp():
    print usage
def WritePfamColorDef(idList, seqid2pfamidDict, fpout):
    pfamidList = []
    for idd in idList:
        if idd in seqid2pfamidDict:
            pfamidList.append(seqid2pfamidDict[idd][0])
    pfamidList = list(set(pfamidList))

#    print "len(pfamidList)", len(pfamidList)

    numColor = len(colorpalette)
    numfam = len(pfamidList)
    colorDict = {}
    for i in xrange(numfam):
        if i < numColor:
            colorDict[pfamidList[i]] = colorpalette[i]
            print colorpalette[i], pfamidList[i]
        else:
            print black, pfamidList[i]

#    print "len(colorDict)", len(colorDict)
#    print "len(idList)", len(idList)

    for idd in idList:
        if idd in seqid2pfamidDict:
            pfamid = seqid2pfamidDict[idd][0]
            if pfamid in colorDict:
                color = colorDict[pfamid]
            else:
                color = black
            fpout.write("%s\t%s\t%s\t%s\n"%(idd, "range", color, pfamid))
            #fpout.write("%s,%s\n"%(idd, color))
        else:
            fpout.write("%s\t%s\t%s\t%s\n"%(idd, "range", black, "Others"))
            #fpout.write("%s,%s\n"%(idd, black ))
# For color definition file, it must be demilited by tab and the ending
# newline is not needed
#    fpout.write("\n")

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1
    infile = ""
    outfile = ""
    seqid2pfamidFile = datadir3 + os.sep + "wk/MPTopo/pfamAna_refpro/pfammap_from_uniprot/refpro20120604-celluar.selmaxlength-m1.nr100.seqid2pfamid"

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
            elif argv[i] in ["-seqid2pfamid", "--seqid2pfamid"]:
                seqid2pfamidFile = argv[i+1]
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
    if seqid2pfamidFile == "" or not os.path.exists(seqid2pfamidFile):
        print >> sys.stderr, "Error. seqid2pfamidFile does not exist. exit"
        return 1
    seqid2pfamidDict = myfunc.ReadFam2SeqidMap(seqid2pfamidFile)
    if seqid2pfamidFile == {}:
        print >> sys.stderr, "Read seqid2pfamidFile failed."
        return 1

    (idList, topoList) = myfunc.ReadFasta_without_annotation(infile)
    if len(idList) < 1:
        print >> sys.stderr, "Read infile failed."
        return 1


    idList.remove("Consensus")
    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)

    WritePfamColorDef(idList, seqid2pfamidDict, fpout)

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

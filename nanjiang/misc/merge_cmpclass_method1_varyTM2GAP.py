#!/usr/bin/env python
# Description:
import os
import sys
import myfunc
usage = """
usage:  %s -datapath DIR [-evodist] [-tm2tm STR]

Description: merge cmpclass table for various TM2GAP thresholds

Options:
  -outpath DIR    Set output path
  -datapath DIR   Set datapath
  -evodist        Use evodist
  -thlist  STR    Threshold list for TM2GAP, default: 0.4 0.5 0.6 0.7 0.8
  -tm2tm   STR    default=0.33
  -q              Quiet mode
  -h, --help      Print this help message and exit

Created 2009-06-08, updated 2011-11-07, Nanjiang Shu
"""%(os.path.basename(sys.argv[0]))

def PrintHelp():
    print usage
def WriteData(dataList, outpath, th_TM2GAPList, th_TM2TM, isEvodist):
    if len(dataList) < 1:
        print >> sys.stderr, "No data read in. Exit"
        return 1

    itemList = dataList[0]['itemList']
    for item in itemList:
        if isEvodist:
            outfile = (outpath + os.sep + "method1.varyTM2GAP.TM2TM%s" % (
                th_TM2TM) + ".evodist%s.cmpclass.txt"%(item))
        else:
            outfile = (outpath + os.sep + "method1.varyTM2GAP.TM2TM%s" % (
                th_TM2TM) + ".seqidt%s.cmpclass.txt"%(item))
        try:
            fpout = open(outfile, "w")
            cnt = 0
            fpout.write("#%3s %8s %9s %6s %6s %6s %6s %18s\n"%("Idx", "th_TM2TM", "th_TM2GAP",
                "IDT", "INV", "TM2GAP", "TM2SEQ", "TM2GAP_AND_TM2SEQ"))
            for dt in dataList:
                #print item, dt[item]
                fpout.write("%-4d %8s %9s %6s %6s %6s %6s %18s\n"%(
                    cnt, dt['th_TM2TM'], dt['th_TM2GAP'],
                    dt[item][0], dt[item][1], dt[item][2],
                    dt[item][3], dt[item][4]))
                cnt += 1
            fpout.close()
#plotfigure
            cmd = "%s/plotCmpClass_mp1_varyTM2GAP.sh %s"%(rundir, outfile)
            os.system(cmd)
            print "%s output"%outfile
        except IOError:
            print >> sys.stderr, "Failed to write to file %s"%(outfile)

def ReadDataFile(infile):#{{{
    try:
        dataDict = {}
        fpin = open(infile,"r")
        lines = fpin.readlines()
        fpin.close()
        itemList = []
        for line in lines:
            if line and line[0] != "#":
                strs = line.split()
                if len(strs) == 9:
                    dataDict[strs[1]] = strs[2:]
                    itemList.append(strs[1])
        if len(itemList) > 0:
            dataDict['itemList'] = itemList
        return dataDict
    except IOError:
        print >> sys.stderr, "Failed to read datafile %s"%(infile)
        return {}
#}}}
def GetDataList(datapath, th_TM2GAPList, th_TM2TM, isEvodist):#{{{
    dataList = []
    for th_TM2GAP in th_TM2GAPList:
        if isEvodist:
            datafile = (datapath + os.sep +
                    "method1.TM2TM%s_TM2GAP%s"%(th_TM2TM, th_TM2GAP)
                    +"_type_all_dg100.0_gap0.0_ps0.0.evodist-cmpclass.table.txt")
        else:
            datafile = (datapath + os.sep +
                    "method1.TM2TM%s_TM2GAP%s"%(th_TM2TM, th_TM2GAP)
                    +"_type_all_dg100.0_gap0.0_ps0.0-cmpclass.table.txt")
        if not os.path.exists(datafile):
            print >> sys.stderr, "datafile %s does not exist. Ignore" %datafile
        dataDict = ReadDataFile(datafile)
        if dataDict != {}:
            dataDict['th_TM2TM'] = th_TM2TM
            dataDict['th_TM2GAP'] = th_TM2GAP
            dataList.append(dataDict)
    return dataList
#}}}

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    
    outpath = ""
    datapath = ""
    isEvodist = False
    th_TM2GAPList = ["0.4", "0.5", "0.6", "0.7", "0.8" ]
    th_TM2TM = "0.33"
    
    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            isNonOptionArg = False
            i += 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-outpath", "--outpath"]:
                outpath = argv[i+1]
                i += 2
            elif argv[i] in ["-datapath", "--datapath"]:
                datapath = argv[i+1]
                i += 2
            elif argv[i] in ["-thlist", "--thlist"]:
                th_TM2GAPList = [argv[i+1].split()]
                i += 2
            elif argv[i] in ["-tm2tm", "--tm2tm"]:
                th_TM2TM = argv[i+1]
                i += 2
            elif argv[i] in ["-evodist"]:
                isEvodist = True
                i += 1
            elif argv[i] in ["-q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            idList.append(argv[i])
            i += 1
    if datapath == "":
        print >> sys.stderr, "datapath not set. exit"
        return 1

    if outpath == "":
        outpath = datapath
    if not os.path.exists(outpath):
        os.system("mkdir -p %s"%outpath)

    dataList = GetDataList(datapath, th_TM2GAPList, th_TM2TM, isEvodist )
    WriteData(dataList, outpath, th_TM2GAPList, th_TM2TM, isEvodist)


#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['isEvodist'] = False
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    rundir = os.path.dirname(sys.argv[0])
    sys.exit(main(g_params))

#!/usr/bin/env python
# Description:
import os
import sys
import myfunc
usage = """
usage:  %s -datapath DIR [-evodist] [-th STR]

Description: merge helixcmpclass table for various TM2GAP thresholds

Options:
  -mode 0|1       Set mode, (default: 0)
                  0: TM2TM
                  1: TM2GAP
  -outpath DIR    Set output path
  -datapath DIR   Set datapath
  -evodist        Use evodist
  -thlist  STR    Threshold list for TM2GAP, 
                  default: 0.1 0.2 0.3 0.33 0.4 0.5 0.6 for mode 0
                  default: 0.4 0.5 0.6 0.7 0.8 for mode 1
  -th      STR    default: 0.5 for mode 0 and 0.33 for mode 1
  -thhigh  FLOAT  Threshold to divided sequence identity or evolutionary
                  distance into high or low, default
                  30 for sequence identity
                  1.5 for evolutionary distance
  -q              Quiet mode
  -h, --help      Print this help message and exit

Created 2012-09-11, updated 2012-09-11, Nanjiang Shu 
"""%(os.path.basename(sys.argv[0]))

def PrintHelp():
    print usage
def WriteData(dataList, outpath, thfixed, thHigh, isEvodist, mode):
    if len(dataList) < 1:
        print >> sys.stderr, "No data read in. Exit"
        return 1

    itemList = ["low", "high", "all"]
    itemStrList = [".lt%g"%thHigh, ".ge%g"%thHigh, ".All"]
    addname = ""
    if mode  == 0:
        addname = "method1.varyTM2TM.TM2GAP%s"%thfixed
    else:
        addname = "method1.varyTM2GAP.TM2TM%s"%thfixed

    for i in xrange(len(itemList)):
        item = itemList[i]
        itemStr = itemStrList[i]
        if isEvodist:
            outfile = (outpath + os.sep + addname +
                    ".evodist%s.helixcmpclass.txt"%(itemStr))
            xlabel = "Evolutionary distance"
        else:
            outfile = (outpath + os.sep + addname +
                    ".seqidt%s.helixcmpclass.txt"%(itemStr))
            xlabel = "Sequence identity"
        try:
            fpout = open(outfile, "w")
            cnt = 0
            fpout.write("#%3s %8s %9s %6s %6s %6s %10s\n"%("Idx", "th_TM2TM", "th_TM2GAP",
                 "TM2TM", "TM2GAP", "TM2SEQ", "Count"))
            for dt in dataList:
                if len(dt[item]) > 0:
#                     print item, dt[item]
                    fpout.write("%-4d %8s %9s %6.3f %6.3f %6.3f %10.0f\n"%(
                        cnt, dt['th_TM2TM'], dt['th_TM2GAP'],
                        dt[item][0], dt[item][1], dt[item][2], dt[item][4]))
                    cnt += 1
            fpout.close()
#plotfigure
            cmd = "%s/plotHelixCmpClass_mp1_vary.sh %s -mode %d"%(rundir,
                    outfile, mode)
            os.system(cmd)
            print "%s output"%outfile
        except IOError:
            print >> sys.stderr, "Failed to write to file %s"%(outfile)


def ReadDataFile(infile, thHigh):#{{{
    try:
        dataDict = {}
        fpin = open(infile,"r")
        lines = fpin.readlines()
        fpin.close()
        lowList = []
        highList = []
        for line in lines:
            if line and line[0] != "#":
                strs = line.split()
                if len(strs) >= 7:
                    [th1, th2] = [float(x) for x in strs[1].split("-")]
                    li1 = [float(strs[j]) for j in range(2,7)]
                    #print "th1=",th1, "thhigh",thHigh
                    if th1 >= thHigh:
                        highList.append(li1)
                    else:
                        lowList.append(li1)
                        #print li1
        
        avgLowList = myfunc.AverageOfFraction(lowList)
        avgHighList = myfunc.AverageOfFraction(highList)
        avgAllList = myfunc.AverageOfFraction(lowList+highList)
        dataDict['low'] = avgLowList
        dataDict['high'] = avgHighList
        dataDict['all'] = avgAllList
        return dataDict
    except IOError:
        print >> sys.stderr, "Failed to read datafile %s"%(infile)
        return {}
#}}}
def GetDataList(datapath, thlist, thfixed, thHigh, isEvodist, mode):#{{{
    dataList = []
    for thvary in thlist:
        if mode  == 0:
            th_TM2TM = thvary
            th_TM2GAP = thfixed
        else:
            th_TM2TM = thfixed
            th_TM2GAP = thvary

        if isEvodist:
            datafile = (datapath + os.sep +
                    "method1.TM2TM%s_TM2GAP%s"%(th_TM2TM, th_TM2GAP)
                    +"_type_all_dg100.0_gap0.0_ps0.0.evodist-helixcmpclass.table.txt")
        else:
            datafile = (datapath + os.sep +
                    "method1.TM2TM%s_TM2GAP%s"%(th_TM2TM, th_TM2GAP)
                    +"_type_all_dg100.0_gap0.0_ps0.0-helixcmpclass.table.txt")
        if not os.path.exists(datafile):
            print >> sys.stderr, "datafile %s does not exist. Ignore" %datafile
#         print datafile
        dataDict = ReadDataFile(datafile, thHigh)
        if dataDict != {}:
            dataDict['th_TM2TM'] = th_TM2TM
            dataDict['th_TM2GAP'] = th_TM2GAP
            dataDict['thHigh'] = thHigh
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
    mode = 0
    th_TM2GAPList = ["0.4", "0.5", "0.6", "0.7", "0.8" ]
    th_TM2TMList = ["0.1", "0.2", "0.3", "0.33", "0.4", "0.5", "0.6"]
    th_TM2TM = "0.33"
    th_TM2GAP = "0.5"

    thHigh = None
    thHigh_evodist = 1.5
    thHigh_seqidt = 30.0
    
    
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
            elif argv[i] in ["-mode", "--mode"]:
                mode = int(argv[i+1])
                i += 2
            elif argv[i] in ["-datapath", "--datapath"]:
                datapath = argv[i+1]
                i += 2
            elif argv[i] in ["-thlist", "--thlist"]:
                th_TM2GAPList = [argv[i+1].split()]
                i += 2
            elif argv[i] in ["-thhigh", "--thhigh"]:
                thHigh = float(argv[i+1])
                i += 2
            elif argv[i] in ["-th", "--th"]:
                th = argv[i+1]
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
            print >> sys.stderr, "Error! Wrong argument:", argv[i]
            return 1

    if datapath == "":
        print >> sys.stderr, "datapath not set. exit"
        return 1

    if mode == 0:
        thlist = th_TM2TMList
        th = th_TM2GAP
    elif mode == 1:
        thlist = th_TM2GAPList
        th = th_TM2TM

    if thHigh == None:
        if isEvodist:
            thHigh = thHigh_evodist
        else:
            thHigh = thHigh_seqidt

    if outpath == "":
        outpath = datapath
    if not os.path.exists(outpath):
        os.system("mkdir -p %s"%outpath)

    dataList = GetDataList(datapath, thlist, th, thHigh, isEvodist, mode)
    WriteData(dataList, outpath, th, thHigh, isEvodist, mode)


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

#!/usr/bin/env python
# Description:
import os
import sys
import myfunc
usage = """
usage:  %s -datapath DIR [-evodist] 

Description: merge cmpclass table for various reliability score thresholds

Options:
  -outpath DIR    Set output path, (default: $datapath)
  -datapath DIR   Set datapath
  -evodist        Use evodist
  -pslist  STR    Threshold list for reliability score, 
                  default: 0.0 10.0 20.0 30.0 40.0 50.0 60.0 70.0 75.0 80.0
                           85.0 90.0 95.0
  -tm2tm  STR     Threshold for TM2TM, (default: 0.33)
  -tm2gap STR     Threshold for TM2GAP, (default: 0.5)
  -thhigh  FLOAT  Threshold to divided sequence identity or evolutionary
                  distance into high or low, default
                  30 for sequence identity
                  1.5 for evolutionary distance
  -q              Quiet mode
  -h, --help      Print this help message and exit

Created 2012-09-18, updated 2012-09-20, Nanjiang Shu 
"""%(os.path.basename(sys.argv[0]))

def PrintHelp():
    print usage
def WriteData(dataList, outpath, thHigh, isEvodist):
    if len(dataList) < 1:
        print >> sys.stderr, "No data read in. Exit"
        return 1

    th_TM2TM = g_params['th_TM2TM']
    th_TM2GAP = g_params['th_TM2GAP']

    itemList = ["low", "high", "all"]
    itemStrList = [".lt%g"%thHigh, ".ge%g"%thHigh, ".All"]
    addname = "method1.TM2TM%s.TM2GAP%s"%(th_TM2TM, th_TM2GAP)

    for i in xrange(len(itemList)):
        item = itemList[i]
        itemStr = itemStrList[i]
        if isEvodist:
            outfile = (outpath + os.sep + addname +
                    ".varyps.evodist%s.cmpclass.txt"%(itemStr))
            xlabel = "Evolutionary distance"
        else:
            outfile = (outpath + os.sep + addname +
                    ".varyps.seqidt%s.cmpclass.txt"%(itemStr))
            xlabel = "Sequence identity"
        try:
            fpout = open(outfile, "w")
            cnt = 0
            fpout.write("#%3s %5s %6s %6s %6s %6s %6s %10s\n"%("Idx", "rlty", "IDT",
                 "INV", "TM2GAP", "TM2SEQ", "TM2GAP_AND_TM2SEQ", "Count"))
            for dt in dataList:
                if len(dt[item]) > 0:
#                     print item, dt[item]
                    fpout.write("%-4d %5s %6.3f %6.3f %6.3f %6.3f %6.3f %10.0f\n"%(
                        cnt, dt['th_ps'],  dt[item][0], dt[item][1], dt[item][2],
                        dt[item][3], dt[item][4], dt[item][6]))
                    cnt += 1
            fpout.close()
#plotfigure
            print "%s output"%outfile
            cmd = "%s/plotCmpClass_mp1_merged_varyps.sh %s -xlabel \"%s\" -outpath %s"%(rundir,
                    outfile, xlabel, outpath)
            os.system(cmd)
            cmd = "%s/plotCmpClass_mp1_merged_varyps.sh %s -xlabel \"%s\" -outpath %s -norm"%(rundir,
                    outfile, xlabel, outpath)
            os.system(cmd)
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
                if len(strs) >= 9:
                    [th1, th2] = [float(x) for x in strs[1].split("-")]
                    li1 = [float(strs[j]) for j in range(2,9)]
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
def GetDataList(datapath, pslist, thHigh, isEvodist):#{{{
    dataList = []
    th_TM2TM = g_params['th_TM2TM']
    th_TM2GAP = g_params['th_TM2GAP']
    for ps in pslist:
        if isEvodist:
            datafile = (datapath + os.sep +
                    "method1.TM2TM%s_TM2GAP%s"%(th_TM2TM, th_TM2GAP)
                    +"_type_all_dg100.0_gap0.0_ps%s.evodist-cmpclass.table.txt"%(ps))
        else:
            datafile = (datapath + os.sep +
                    "method1.TM2TM%s_TM2GAP%s"%(th_TM2TM, th_TM2GAP)
                    +"_type_all_dg100.0_gap0.0_ps%s-cmpclass.table.txt"%(ps))
        if not os.path.exists(datafile):
            print >> sys.stderr, "datafile %s does not exist. Ignore" %datafile
#         print datafile
        dataDict = ReadDataFile(datafile, thHigh)
        if dataDict != {}:
            dataDict['th_TM2TM'] = th_TM2TM
            dataDict['th_ps'] = ps
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
    pslist = ["0.0", "10.0", "20.0", "30.0", "40.0", "50.0", "60.0", "70.0",
            "75.0", "80.0", "85.0", "90.0", "95.0"]

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
            elif argv[i] in ["-datapath", "--datapath"]:
                datapath = argv[i+1]
                i += 2
            elif argv[i] in ["-pslist", "--pslist"]:
                pslist = [argv[i+1].split()]
                i += 2
            elif argv[i] in ["-thhigh", "--thhigh"]:
                thHigh = float(argv[i+1])
                i += 2
            elif argv[i] in ["-tm2tm", "--tm2tm"]:
                g_params['th_TM2TM'] = argv[i+1]
                i += 2
            elif argv[i] in ["-tm2gap", "--tm2gap"]:
                g_params['th_TM2GAP'] = argv[i+1]
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


    if thHigh == None:
        if isEvodist:
            thHigh = thHigh_evodist
        else:
            thHigh = thHigh_seqidt

    if outpath == "":
        outpath = datapath
    if not os.path.exists(outpath):
        os.system("mkdir -p %s"%outpath)

    dataList = GetDataList(datapath, pslist,  thHigh, isEvodist)
    WriteData(dataList, outpath, thHigh, isEvodist)


#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['th_TM2TM'] = "0.33"
    g_params['th_TM2GAP'] = "0.5"
    g_params['isEvodist'] = False
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    rundir = os.path.dirname(sys.argv[0])
    sys.exit(main(g_params))

#!/usr/bin/env python
# Description:
import os
import sys

progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))
usage = """
usage:  %s  FILE [-outpath DIR]
Description:

Options:
  -outpath DIR    Set ouput path
  -win INT 
  -format STR     auto, cmpclass, helixcmpclass
  -h, --help      Print this help message and exit

Created 2012-12-14, updated 2012-12-14, Nanjiang Shu 
"""%(progname)

rundir = os.path.dirname(sys.argv[0])
binpath = rundir

def PrintHelp():
    print usage

def ReadFileRltyCmpclass(infile):
    try:
        recordList = []
        fpin = open(infile,"r")
        lines = fpin.readlines()
        fpin.close()
        for line in lines:
            strs = line.split()
            try:
                ps = float(strs[0])
                cmpclass = strs[1]
                recordList.append((ps, cmpclass))
            except (IndexError, TypeError):
                pass
        return recordList
    except IOError:
        print >> sys.stderr, "Failed to read infile %s"%(infile)
        return []
def InitDataTable(dataTable, classList):#{{{
    numClass = len(classList)
    dataTable['lx'] = []
    dataTable['subsum'] = []
    dataTable['ly'] = []
    for i in xrange(numClass):
        dataTable['ly'].append([])
    #}}}

def CalRunningAverage(recordList, cmpclassList, dataTable, win, step, classList):
    numRecord = len(recordList)
    numClass = len(classList)
    beg = 0
    while beg < numRecord - win:
        count = [0] * len(classList)
        end = beg + win
        pivot = int((beg + end)/2.0+0.5)
        md = recordList[pivot][0]

        beg_tmp = beg + step
        while beg < numRecord-win:
            end_tmp = beg_tmp + win
            pivot = int((beg_tmp + end_tmp)/2.0+0.5)
            md_tmp = recordList[pivot][0]
            if md_tmp <= md:
                beg_tmp += step
            else:
                break
        end = beg_tmp + win
        subCmpclassList = cmpclassList[beg:end]
        subnum = len(subCmpclassList)
        dataTable['lx'].append(md)
        dataTable['subsum'].append(end - beg)
        for j in xrange(numClass):
            count[j] = subCmpclassList.count(classList[j])
            dataTable['ly'][j].append(count[j]/float(len(subCmpclassList)))
        beg = end + step - win

def WriteTable2D(dataTable, classList, outfile):
    try:
        fpout = open(outfile, "w")
        numClass = len(classList)
        fpout.write("%4s %7s"%("#Idx","RLTY"))
        for cls in classList:
            fpout.write(" %9s"%(cls))
        fpout.write(" %10s"%"Occurrence")
        fpout.write("\n")
        for i in xrange(len(dataTable['lx'])):
            fpout.write("%-4d %7.1f"%(i, dataTable['lx'][i]))
            for j in xrange(numClass):
                fpout.write(" %9.3f"%(dataTable['ly'][j][i]*100))
            fpout.write(" %10d"%dataTable['subsum'][i])
            fpout.write("\n")
        fpout.close()
        print outfile, "output."
    except IOError:
        print >> sys.stderr, "Failed to write to file %s"%(outfile)

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    classList_helixcmpclass = ["0", "1", "2"]
    classList_cmpclass = ["IDT", "INV", "TM2GAP", "TM2SEQ",
            "TM2GAP_AND_TM2SEQ"]

    outpath = "./"
    infile = ""
    win = 5000
    file_format = "auto"

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
            elif argv[i] in ["-outpath", "--outpath"]:
                outpath = argv[i+1]
                i += 2
            elif argv[i] in ["-win", "--win"]:
                win = int( argv[i+1])
                i += 2
            elif argv[i] in ["-format", "--format"]:
                file_format = argv[i+1]
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

    recordList = ReadFileRltyCmpclass(infile)
    recordList = sorted(recordList, key=lambda x:x[0], reverse=False)
    cmpclassList = [x[1] for x in recordList]

    rootname = os.path.basename(os.path.splitext(infile)[0]);
    if file_format == "auto":
        if infile.find("helixcmpclass") != -1:
            file_format = "helixcmpclass"
        elif infile.find("cmpclass") != -1:
            file_format = "cmpclass"
    
        else:
            print >> sys.stderr, "unknow file format. exit"
            return 1
    if file_format == "cmpclass":
        classList = classList_cmpclass
    elif file_format == "helixcmpclass":
        classList = classList_helixcmpclass

    dataTable = {}
    step = int(win*0.05)
#0-1000, 10-1010
    InitDataTable(dataTable, classList)
    CalRunningAverage(recordList, cmpclassList, dataTable, win, step, classList)
    outfile = outpath + os.sep + rootname + "win%d_step%d.table.txt"%(win, step)
    WriteTable2D(dataTable, classList, outfile)
    if file_format == "helixcmpclass":
        cmd = ("%s/plotHelixCmpClass_mp1_xy.sh %s"%(binpath, outfile))
    elif file_format == "cmpclass":
        cmd = ("%s/plotCmpClass_mp1_xy.sh %s"%(binpath, outfile))
    os.system(cmd)
    
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

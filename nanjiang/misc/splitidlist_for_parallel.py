#!/usr/bin/env python
# Description:
import os
import sys
import myfunc
from math import ceil
usage = """
usage:  splitidlist_for_parallel.py infile -ns INT [-outpath DIR]
Description:

Options:
  -ns      INT    Number of split
  -outpath DIR    Set ouput path, (default: ./)
  -m, -method INT Set method for splitting, (default: 0)
                  0: evenly split by the number of ID
                  1: evenly split according to the filesize, i.e. 
                     sum(filesize) of each split = sum(filesize)/N
                  2: evenly split according to the square of filesize, i.e.
                     sum(filesize^2) of each split = sum(filesize^2)/N
  -datapath DIR   Datapath for file
  -ext      STR   Extension for file, filename = $id$ext
  -q              Quiet mode
  -h, --help      Print this help message and exit

Created 2012-09-03, updated 2012-09-03, Nanjiang Shu 

Example:
    splitidlist_for_parallel.py -ns 4 idlist.txt -outpath result -m 1 -datapath seqdir/ -ext .fa
"""

def PrintHelp():
    print usage
def SplitIDList(infile, nsplit, method, datapath, ext, outpath):#{{{
    idList = myfunc.ReadIDList(infile)
    rootname = os.path.basename(os.path.splitext(infile)[0])
    numID = len(idList)
    if numID <= 0:
        print >> sys.stderr, "no ID in the idlist file %s"%(infile)
        return 1

    if method == 0:
        nfile_per_split = int(ceil(numID / float(nsplit)))
        i = 0
        cntfile = 0
        while i < numID:
            outfile = outpath + os.sep + rootname + ".split_%d"%(cntfile)
            fpout = myfunc.myopen(outfile, None, "w", True)
            cntID_per_split = 0
            for j in xrange(i, i + nfile_per_split):
                if j < numID:
                    fpout.write("%s\n"%(idList[j]))
                    cntID_per_split += 1
            myfunc.myclose(fpout)
            print "split to %s \t %4d IDs"%(outfile, cntID_per_split)
            cntfile += 1
            i += nfile_per_split
    elif method in  [1,2]:
        sumFileSize = 0.0
        fsizeList = []
        for idd in idList:
            fname = datapath + os.sep + idd + ext
            if os.path.exists(fname):
                fsize = os.path.getsize(fname)
                if fsize > 0:
                    if method == 1:
                        sumFileSize += float(fsize)
                        fsizeList.append((idd, float(fsize)))
                    elif method == 2:
                        sumFileSize += float(fsize)*float(fsize)
                        fsizeList.append((idd, float(fsize)*fsize))
        sumfilesize_per_split = ceil(sumFileSize / float(nsplit))
        fsizeList = sorted(fsizeList, key = lambda x:x[1], reverse=True)
        print "sumfilesize_per_split = %g"% (sumfilesize_per_split)
        i = 0
        cntfile = 0
        numID = len(fsizeList)
        while i < numID:
            outfile = outpath + os.sep + rootname + ".split_%d"%(cntfile)
            fpout = myfunc.myopen(outfile, None, "w", True)
            j = 0
            cntID_per_split = 0
            sumFileSize = 0.0
            while sumFileSize <= sumfilesize_per_split:
                idx = i+j
                if i+j > numID -1:
                    break
                idd = fsizeList[idx][0]
                fsize = fsizeList[idx][1]
                sumFileSize += fsize
                fpout.write("%s\n"%(idd))
                cntID_per_split += 1
                j += 1
            myfunc.myclose(fpout)
            print "split to %s \t %4d IDs sumsize = %g"%(outfile, cntID_per_split, sumFileSize)
            cntfile += 1
            i += j
        return 1
#}}}

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = "./"
    infile = ""
    method = 0
    datapath = ""
    ext = ""
    nsplit = 0

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
            elif argv[i] in ["-datapath", "--datapath"]:
                datapath = argv[i+1]
                i += 2
            elif argv[i] in ["-ext", "--ext"]:
                ext = argv[i+1]
                i += 2
            elif argv[i] in ["-m", "--m", "-method", "--method"] :
                method = int(argv[i+1])
                i += 2
            elif argv[i] in ["-ns", "--ns"] :
                nsplit = int(argv[i+1])
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
        print >> sys.stderr, "infile %s not set or does not exist. Exit"
        return 1
    if nsplit == 0:
        print >> sys.stderr, "nsplit not set. Exit"
        return 1
    if method in [0, 1,2]:
        if method in [1,2]:
            if datapath == "" or not os.path.exists(datapath):
                print >> sys.stderr, \
                        "method is set to %d, but datapath not set. Exit" %(
                        method)
                return 1
            if ext == "":
                print >> sys.stderr, \
                        "method is set to %d, but ext not set. Exit" %(
                        method)
                return 1
    else:
        print >> sys.stderr, "Unrecognized method %d" %( method)
        return 1

    if not os.path.exists(outpath):
        os.system("mkdir -p %s"%(outpath))

    return SplitIDList(infile, nsplit, method, datapath, ext, outpath)
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

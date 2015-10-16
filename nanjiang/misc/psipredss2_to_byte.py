#!/usr/bin/env python
# Description:
# convert the ASCII format psipred ss2 file to binary format
import os
import sys
import myfunc
from array import array
import numpy as np
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))
usage = """
usage:  %s [-l LISTFILE] [-outpath DIR] [-q]
        %s FILE [FILE ...]
Description: convert the ASCII format psipred ss2 file to binary format

Options:
  -outpath DIR    Set ouput path
  -l      FILE    Set the fileListFile
  -q              Quiet mode
  -h, --help      Print this help message and exit

Created 2013-02-25, updated 2013-02-25, Nanjiang Shu 
"""%(progname, wspace)

DEBUG = 0;
DEBUG_READBINARY = 0;

def PrintHelp():
    print usage

def PrintTextSS2(aaSeq, ssSeq, probArrayList):
    print "# PSIPRED VFORMAT (PSIPRED V3.2)\n"
    seqLength = len(aaSeq)
    for i in xrange(seqLength):
        msg = "%4d %s %s  %6.3f %6.3f %6.3f"
        print msg%(i+1, aaSeq[i], ssSeq[i],
                probArrayList[0][i]/1000.0,
                probArrayList[1][i]/1000.0,
                probArrayList[2][i]/1000.0)
def ReadPSIPREDSS2(infile):
    hdl = myfunc.ReadLineByBlock(infile)
    if hdl.failure:
        return (None, None, None)
    aaSeqList = []
    ssSeqList = []
    arrayList = []
    for i in range(3):
        arrayList.append(array('h'))

    lines = hdl.readlines()
    while lines != None:
        for line in lines:
            strs = line.split()
            if len(strs) == 6 and strs[0].isdigit():
                aaSeqList.append(strs[1])
                ssSeqList.append(strs[2])
                for i in range(3):
                    try:
                        value = int(float(strs[i+3])*1000)
                        arrayList[i].append(value)
                    except (ValueError, IndexError):
                        msg = "Bad record \"%s\" in file %s"
                        print >> sys.stderr, msg%(line, infile)
                        return (None, None, None)
        lines = hdl.readlines()
    hdl.close()
    aaSeq = "".join(aaSeqList)
    ssSeq = "".join(ssSeqList)
    return (aaSeq, ssSeq, arrayList)
#}}}

def WriteBinarySS2(aaSeq, ssSeq, probArrayList, outfile):
    try:
        fpout = open(outfile, "wb")
    except IOError:
        print >> sys.stderr, "Failed to write to file %s"%(outfile)
        return 1
    seqLength = len(aaSeq)
    vI = array('I')
    vI.append(seqLength)
    vI.tofile(fpout)
    fpout.write(aaSeq)
    fpout.write(ssSeq)
    for i in range(3):
        probArrayList[i].tofile(fpout)

    fpout.close()

def ReadBinarySS2(infile):
    try:
        aaSeq = ""
        ssSeq = ""
        arrayList = []
        fpin = open(infile, "rb")
        vI = np.fromfile(fpin, dtype=np.int32, count=1)

        numrecord = vI[0];
        if DEBUG_READBINARY:
            print "numrecord=",numrecord

        vI = np.fromfile(fpin, dtype=np.int32, count=1)
        sizeAnno = vI[0];
        if DEBUG_READBINARY:
            print "sizeAnno=",sizeAnno

        anno = fpin.read(sizeAnno);
        if DEBUG_READBINARY:
            print "anno=",anno

        vI = np.fromfile(fpin, dtype=np.int32, count=1)
        seqlen = vI[0];
        if DEBUG_READBINARY:
            print "seqlen=",seqlen

        aaSeq = fpin.read(seqlen)
        ssSeq = fpin.read(seqlen)
        if DEBUG_READBINARY:
            print "aaSeq=",aaSeq
            print "ssSeq=",ssSeq
#        sys.exit(1)
        for i in range(3):
            arrayList.append(np.fromfile(fpin, dtype=np.int16,
                count=seqlen))

        fpin.close()
        return (aaSeq, ssSeq, anno, arrayList)
    except IOError:
        print >> sys.stderr, "Failed to read file %s"%(infile)
        return (None, None, None, None)


def PSIPRED_ss2_to_binary(infile, g_outpath):
    (aaSeq, ssSeq, probArrayList) = ReadPSIPREDSS2(infile)
    if aaSeq == None:
        print >> sys.stderr, "Failed to read file %s"%(infile)
        return ""
    seqLength = len(aaSeq)
    if g_outpath != "":
        outpath = g_outpath
    else:
        outpath = os.path.dirname(infile)
        if outpath == "":
            outpath = "."
    basename = os.path.basename(infile)
    outfile = outpath + os.sep + basename + ".bin"
    WriteBinarySS2(aaSeq, ssSeq, probArrayList, outfile)
    return outfile

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = ""
    fileListFile = None
    fileList = []

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            fileList.append(argv[i])
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
            elif argv[i] in ["-l", "--l"] :
                fileListFile = argv[i+1]
                i += 2
            elif argv[i] in ["-q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            fileList.append(argv[i])
            i += 1

    if fileListFile != None:
        try:
            fp = open(fileListFile,"r")
            fileList += fp.read().split()
            fp.close()
        except IOError:
            print >> sys.stderr, "file %s does not exist." %fileListFile

#     for i in xrange(len(fileList)):
#         infile = fileList[i];
#         (aaSeq, ssSeq, seqAnno, probArrayList) = ReadBinarySS2(infile)
#         #print "seqAnno=",seqAnno
#         if aaSeq != None:
#             PrintTextSS2(aaSeq, ssSeq, probArrayList);
#     return 1
    for i in xrange(len(fileList)):
        ss2file = fileList[i]
        outfile = PSIPRED_ss2_to_binary(ss2file, outpath)
# debug, read binary file to see if it is correct
        if DEBUG and outfile != "":
            (aaSeq, ssSeq, probArrayList) = ReadBinarySS2(outfile)
            if aaSeq != None:
                print "# PSIPRED VFORMAT (PSIPRED V3.2)\n"
                seqLength = len(aaSeq)
                for i in xrange(seqLength):
                    msg = "%4d %s %s  %6.3f %6.3f %6.3f"
                    print msg%(i+1, aaSeq[i], ssSeq[i],
                            probArrayList[0][i]/1000.0,
                            probArrayList[1][i]/1000.0,
                            probArrayList[2][i]/1000.0)
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

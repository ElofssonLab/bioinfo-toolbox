#!/usr/bin/env python
# Description:
import os
import sys
import myfunc
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))
usage = """
usage:  %s topofile -sp signalp_file [-o OUTFILE]
Description: mask predicted TM helices also predicted as signal peptide
             i/o status is taken as the i/o status of the residue after
             the masked TM region, so that the status of the rest TM
             regions will not be changed
Options:
  -o  OUTFILE  output result to file
  -sp  FILE    signalp result file, short format
  -f   STR     format of the signal peptide prediction, 
               (default: signalp), can be phobius, singalp
  -q           Quiet mode
  -h, --help   Print this help message and exit

Created 2013-01-14, updated 2013-01-16, Nanjiang Shu 
"""%(progname)

DEBUG = 1

def PrintHelp():
    print usage

def MaskTopologyBySignalPeptide(idList, topoList, signalpDict):
    newTopoList = []
    for i in xrange(len(idList)):
        topo = topoList[i]
        if idList[i] in signalpDict:
            posTMList = myfunc.GetTMPosition(topo)
            try:
                posSigP = signalpDict[idList[i]]
                (b,e) = (posTMList[0][0],posTMList[0][1])
                cov = myfunc.coverage(0, posSigP, b, e)
                if float(cov)/(e-b) > 0.5:
#mask
                    masked_state = topo[e]
                    newTopo = ( "".join([masked_state]*(e)) +
                            topo[e:])
                    newTopoList.append(newTopo)
                    if DEBUG:
                        print
                        print "posTM", (b,e), "SignalPeptide", posSigP
                        print topo
                        print newTopo
                else:
                    newTopoList.append(topo)
            except (KeyError, IndexError):
                newTopoList.append(topo)
        else:
            newTopoList.append(topo)
    return newTopoList
def ReadSignalPeptide_signalp(infile):
    try:
        signalpDict = {}
        fpin = open(infile, "r")
        lines = fpin.readlines()
        fpin.close()
        for line in lines:
            if line[0] != "#":
                strs = line.split()
                try:
                    status = strs[9]
                    if status == "Y":
                        seqid = myfunc.GetSeqIDFromAnnotation(strs[0])
                        signalpDict[seqid] = int(strs[2])
                except IndexError:
                    pass
        return signalpDict
    except IOError:
        print >> sys.stderr, "Failed to read infile %s"%infile
        return {}

def ReadSignalPeptide(infile, format_sp):#{{{
    if format_sp == "signalp":
        return ReadSignalPeptide_signalp(infile)
    elif format_sp == "phobius":
        return ReadSignalPeptide_phobius(infile)
#}}}
def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile = ""
    infile = ""
    signalp_file = ""
    format_sp = "signalp"

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
            elif argv[i] in ["-o", "--o"]:
                outfile = argv[i+1]
                i += 2
            elif argv[i] in ["-sp", "--sp"] :
                signalp_file = argv[i+1]
                i += 2
            elif argv[i] in ["-f", "--f", "-format", "--format"] :
                format_sp = argv[i+1]
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
        print >> sys.stderr, "infile not set or does not exist"
        return 1
    if signalp_file == "" or not os.path.exists(signalp_file):
        print >> sys.stderr, "signalp file not set or does not exist"
        return 1
    if not format_sp in ["signalp", "phobius"]:
        print >> sys.stderr, "format_sp = %s is not supported. Exit." %(
                format_sp)

    
    signalpDict = ReadSignalPeptide(signalp_file, format_sp)
    (idList, annoList, topoList) = myfunc.ReadFasta(infile)
    
    newTopoList = MaskTopologyBySignalPeptide(idList, topoList, signalpDict)

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)

    for i in xrange(len(idList)):
        fpout.write(">%s\n"%(annoList[i]))
        fpout.write("%s\n"%(newTopoList[i]))

    myfunc.myclose(fpout)

#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

#!/usr/bin/env python
# Description:  
#    Run multiple sequence alignment on uppmax by calling the
#    script "runmsa_multiprog.sh"
import os
import sys
import tempfile
import subprocess
import myfunc

default_progList=[
"kalign",
"kalignP",
"muscle_maxiter2",
"mafft",
"mafft_auto",
"clustalo",
"clustalo_auto",
"clustalw",
"muscle",
"t_coffee",
"probcons"
        ];
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))
usage = """
usage:  %s seqfilelist-with-numseq  [-prog PROG1 -prog PROG2]

DESCRIPTION: 
    Sumit runmsa_multiprog jobs to uppmax

OPTIONS:
  -outpath  DIR   Set ouput path, (default: the same as sequence file)
  -prog     STR   Set the program name in command line.
  -proglist FILE  Set the program list file, one item per line
  -gzip  yes|no   Whether gzip the alignment file, (default: no)
  -num       INT  Number of average sequenec per split,(default: 1000)
  -l        FILE  Set the idListFile
  -q              Quiet mode
  -h, --help      Print this help message and exit

Created 2013-02-13, updated 2013-11-27, Nanjiang Shu

Examples:
    %s t1.seqfilelist.withnumseq -prog kalignP -prog clustalo10 -outpath msa -gzip no
"""%(progname, progname)
projnumber = "snic2013-11-2"
# it takes about one hour to run a query with 1000 sequences
nodename = os.uname()[1]
if nodename.find("uppmax.uu.se") != -1:
    print "OK, you are running from uppmax node"
else:
    msg = "Error. This script should be run from uppmax login node. Exit"
    print >> sys.stderr, msg
#    sys.exit(1)


def PrintHelp():
    print usage

def GetTimeThreshold(prog):
    # running msa for 1000 sequences take 1 hour
    if prog in ["kalign"]:
        threshold = g_params['num_per_split']*300*300/2*4.0
    elif prog in ["kalignP", "muscle_maxiter2", "mafft", "mafft_auto"]:
        threshold = g_params['num_per_split']*300*300/2*1.0
    elif prog in ["clustalo", "clustalo_auto"]:
        threshold = g_params['num_per_split']*300*300/2*0.6
    elif prog in ["clustalw", "muscle"]:
        threshold = g_params['num_per_split']*300*300/2*0.4
    elif prog in ["t_coffee", "probcons"]:
        threshold = g_params['num_per_split']*300*300/2*0.2
    elif prog.find("clustalo") != -1:
        intr = int(prog.replace("clustalo",""))
        threshold =  g_params['num_per_split']*300*300/2*0.6/intr
    return threshold

def ReadInputList(infile):#{{{
    try:
        inputList = []
        fpin = open(infile, "r")
        lines = fpin.readlines()
        fpin.close()
        for line in lines:
            if not line or line[0] == "#":
                continue
            strs = line.split()
            try:
                inputList.append((strs[0], int(strs[1])))
            except (IndexError, TypeError, ValueError):
                print >> sys.stderr, "bad record: %s"%(line)
                pass
        return inputList
    except IOError:
        return []
#}}}
def WriteIDListFile(outfile, idList):
    try:
        fpout = open(outfile, "w")
        for idd in idList:
            print >> fpout, idd
        fpout.close()
        return 0
    except IOError:
        print >> sys.stderr, "Failed to write to file %s"%(outfile)
        return 1


def WriteScriptFile(scriptfile, idlistfile, logfile, prog):
    try:
        if prog == "t_coffee":
            numparallel = 1
        else:
            numparallel = 8

        numparallel = 1
        fpout = open(scriptfile, "w")
        if g_params['isGzip'] == True:
            str_gzip = "yes"
        else:
            str_gzip = "no"
        msg = """\
#!/bin/bash
$DATADIR3/wk/MPTopo/src/runmsa_multiprog.sh -prog {} -x {}  -l {} > {} 2>&1
"""
        print >> fpout, msg.format(prog, str_gzip, idlistfile, logfile)
        fpout.close()
        return 0
    except IOError:
        print >> sys.stderr, "Failed to write to file %s"%(outfile)
        return 1

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    infile = ""
    progList = []
    progListFile = ""
    outpath = ""

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
            elif argv[i] in ["-prog", "--prog"]:
                tmpstr, i = myfunc.my_getopt_str(argv, i)
                progList.append(tmpstr)
            elif argv[i] in ["-gzip", "--gzip"]:
                tmpstr, i = myfunc.my_getopt_str(argv, i)
                if tmpstr.upper()[0] == "-":
                    print >> sys.stderr, "Bad argument, -gzip should be"\
                            " followed by yes or no"
                    return 1
                elif tmpstr.upper()[0] == "Y":
                    g_params['isGzip'] = True
                else:
                    g_params['isGzip'] = False
            elif argv[i] in ["-num", "--num"]:
                g_params['num_per_split'], i = myfunc.my_getopt_int(argv, i)
            elif argv[i] in ["-proglist", "--proglist"]:
                progListFile, i = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-outpath", "--outpath"]:
                outpath, i = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q"]:
                g_params['isQuiet'] = True; i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            infile = argv[i]
            i += 1

    if myfunc.checkfile(infile,"infile") != 0:
        return 1

    inputList = ReadInputList(infile) # [(filename, numseq)]
    inputList = sorted(inputList, key=lambda x:x[1], reverse=False)
    rtname_infile = os.path.basename(os.path.splitext(infile)[0])

# get progList
    if len(progList) == 0 and progListFile == "":
        progList = default_progList
    else:
        if progListFile != "":
            tmp_list = myfunc.ReadIDList(progListFile)
            if len(tmp_list) == 0:
                print >> sys.stderr, "progListFile %s does not exist or empty"%(
                        progListFile)
                return 1
            else:
                progList += tmp_list
        if len(progList) == 0:
            print >> sys.stderr, "progList is empty. exit"
            return 1

    if outpath != "" and not os.path.exists(outpath):
        try:
            subprocess.check_output(["mkdir", "-p", outpath])
        except subprocess.CalledProcessError, e:
            print e
            return 1

    numInput = len(inputList)
    cwd = os.getcwd()
    workdir = tempfile.mkdtemp(prefix="%s/workdir.%s."%(cwd, rtname_infile))
    print "workdir=", workdir
    for prog in progList:
        print "\n============== PROG: %s =============\n"%(prog)
        threshold = GetTimeThreshold(prog)
        i = 0
        cntf = 0
        while i < numInput:
            splitList = []
            cntscore = 0
            j = 0
            while i+j < numInput:
                pair = inputList[i+j]
                j += 1
                seqfilename = pair[0]
                numseq = pair[1]
                if numseq >= 2: # There should be at least two sequences to run msa
                    rootname = os.path.basename(os.path.splitext(seqfilename)[0])
                    if outpath != "":
                        cur_outpath = outpath
                    else:
                        cur_outpath = myfunc.my_dirname(seqfilename)
                    gzfile = cur_outpath + os.sep + "%s.%s.mfa.gz"%(rootname, prog)
                    msafile = cur_outpath + os.sep + "%s.%s.mfa"%(rootname, prog)
                    if ((g_params['isGzip'] == True and not os.path.exists(gzfile)) or 
                            (g_params['isGzip'] == False and not os.path.exists(msafile))):
                        cntscore += pair[1]**2
                        splitList.append(pair[0])
                        if cntscore > threshold:
                            break
                    else:
                        print "%s already exist. Ignore"%(gzfile);

            if len(splitList) > 0:
                splitidlistfile = workdir + os.sep \
                        + "splitidlist.%s.%d.idlist"%(prog, cntf)
                splitscriptfile = workdir + os.sep \
                        + "splitscript.%s.%d.sh"%(prog, cntf)
                logfile = workdir + os.sep + "log.%s.%d.txt"%(prog, cntf)
                jobname = "msa.%s.%d"%(prog, cntf)
                WriteIDListFile(splitidlistfile, splitList)
                WriteScriptFile(splitscriptfile, splitidlistfile, logfile, prog)

                sbatch_cmd_str = "sbatch -A {} -p core -n 1 -t 1-24:00:00 -J {} {}"
                sbatch_cmd =  sbatch_cmd_str.format(projnumber, jobname, splitscriptfile)
                print sbatch_cmd
                os.system(sbatch_cmd)
                cntf += 1
            i += j
    print "workdir =", workdir

#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['isGzip'] = False
    g_params['num_per_split'] = 1000
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

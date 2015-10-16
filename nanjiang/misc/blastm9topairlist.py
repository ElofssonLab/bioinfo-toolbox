#!/usr/bin/env python
# get pairlist from blast m9 output
import sys,re,os
import myfunc
import tempfile

progname =  os.path.basename(sys.argv[0])

usage="""
Usage: %s [Options] FILE [FILE ...] [-l filelist]

Description:
    get the matched pairlist form blast m9 output

OPTIONS:
  -i         FILE      input file
  -o         FILE      outputfile
  -evalue    float     set maximum evalue threshold, (default: 1e-10)
  -seqidt    FLOAT     set the minimum sequence identity threshold, (default: 95)
  -round     int       to use the hits of which iteration, default is using the last iteration
  -h,--help            print this help message and exit

Created 2013-12-11, updated 2013-12-11, Nanjiang Shu 
"""%(progname)

def PrintHelp():
    print usage

def ReadBlastOutput(infile, iteration=99999, fmt=9):
    """
    Read output of blastpgp
    return recordlist in dictionary
    """
    try:
        linenumber_iteration=[]
        fpin = open(infile, "r")
        lines = fpin.read().split('\n')
        fpin.close()
    except IOError:
        print >> sys.stderr, "Failed to read file %s"%(infile)
        return []

    recordList = []
    cntLine=0
    for line in lines:
        if not line:
            continue
        if line[0] == "#" and line.find("Iteration") >=0:
            linenumber_iteration.append(cntLine)
        cntLine +=1

    if iteration > len(linenumber_iteration):
        iteration = len(linenumber_iteration)

    cntLine=0
    for line in lines:
        if not line:
            continue
        if cntLine > linenumber_iteration[iteration-1]:
            if line.find("BLASTP") >=0 :
                break
            elif line[0] != "#":
                strs=line.split('\t') 
                queryID = strs[0]
                hitID=strs[1]
                seqidt = float(strs[2])
                evalue=float(strs[10])
                if (evalue <= g_params['evalue_th'] and
                        seqidt >= g_params['seqidt_th']):
                    hit = {}
                    seqidt =float(strs[2])
                    alnLength=int(strs[3])
                    misMatch=int(strs[4])
                    gapOpens=int(strs[5])
                    qstart=int(strs[6])
                    qend=int(strs[7])
                    sstart=int(strs[8])
                    send=int(strs[9])
                    evalue=float(strs[10])
                    bitscore=float(strs[11])
                    hit['queryID'] = queryID
                    hit['hitID'] = hitID
                    hit['seqidt'] = seqidt
                    hit['alnLength']=alnLength
                    hit['misMatch']=misMatch
                    hit['gapOpens']=gapOpens
                    hit['qstart']=qstart
                    hit['qend']=qend
                    hit['sstart']=sstart
                    hit['send']=send
                    hit['evalue']=evalue
                    hit['bitscore']=bitscore
                    recordList.append(hit)
        cntLine+=1
#first get line number of the Iteration record
    return recordList
#}}}
def BlastM9toPairlist(infile, fpout):#{{{
    recordList = ReadBlastOutput(infile, iteration=g_params['iteration'], fmt=9)
    #create a new dict to sort the idlist by evalue in  ascending order
    for i in range(min(1, len(recordList))):
        rd = recordList[i]
        print >> fpout, rd['queryID'], myfunc.GetSeqIDFromAnnotation(rd['hitID'])
#}}}


def main(g_params):#{{{
    # Check argv
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile = ""
    fileList = []
    fileListFile = ""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            fileList.append(sys.argv[i])
            isNonOptionArg=False
            i += 1
        elif sys.argv[i] == "--":
            isNonOptionArg=True
            i += 1
        elif sys.argv[i][0] == "-":
            if sys.argv[i] in [ "-h", "--help"]:
                PrintHelp()
                return 1
            elif sys.argv[i] in  ["-l", "--l"]:
                fileListFile, i = myfunc.my_getopt_str(sys.argv,i)
            elif sys.argv[i] in [ "-o" ,"--o", "-outfile"]:
                outfile, i = myfunc.my_getopt_str(sys.argv,i)
            elif sys.argv[i] in ["-evalue" ,  "--evalue"]:
                g_params['evalue_th'], i  = myfunc.my_getopt_float(sys.argv, i)
            elif sys.argv[i] in ["-seqidt" ,  "--seqidt"]:
                g_params['seqidt_th'], i  = myfunc.my_getopt_float(sys.argv, i)
            elif sys.argv[i] in ["-round",  "--round"]:
                g_params['iteration'] = myfunc.my_getopt_int(sys.argv,i)
            else:
                print >> sys.stderr,("Error! Wrong argument: '%s'" % sys.argv[i])
                return 1
        else:
            fileList.append(sys.argv[i])
            i += 1


    if fileListFile != "":
        fileList += myfunc.ReadIDList(fileListFile, delim="\n")
    if len(fileList) < 1:
        print >> sys.stderr, "No input set. exit"
        return 1

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)


    for infile in fileList:
        BlastM9toPairlist(infile, fpout)

    myfunc.myclose(fpout)

#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['evalue_th'] = 1e-10
    g_params['seqidt_th'] = 95.0
    g_params['iteration'] = 99999 # using the last round by default
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

#!/usr/bin/env python
# Description:
import os
import sys
import myfunc
import tempfile
from Bio import SeqIO

progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))
usage = """
usage:  %s fastaseq [-o OUTFILE]
Description: Retrieve sequence by refseq id or uniprot id, target or query
sequence are output as it is.

Options:
  -o  OUTFILE    Output the result to file
  -l  LISTFILE   Set list file, two str per line, delimited by tab
                 Format: infile TAB outfile
  -seqdb STR     Supply sequence database, formated by indexfasta.py
  -q             Quiet mode
  -h, --help     Print this help message and exit

Created 2013-02-04, updated 2013-02-13, Nanjiang Shu 
"""%(progname)


def PrintHelp():
    print usage

def GetFileType(infile):#{{{
    try:
        fpin = open(infile, "r")
        line = fpin.readline()
        filetype = ""
        while line:
            line.rstrip("\n")
            if line:
                if line[0] == ">":
                    filetype = "fasta"
                    break
                else:
                    filetype = "idlist"
                    break
            line = fpin.readline()
        fpin.close()
        return filetype
    except IOError:
        print >> sys. stderr, "Failed to read infile %s"%(infile)
        return ""
#}}}
def GetDatabaseIDList(annoList):#{{{
    idList = [] 
    for anno in annoList:
        if anno == "" or anno[0] == "#":
            continue
        firstword = myfunc.GetFirstWord(anno)
        lengthword = len(firstword)
        p1 = firstword.find('(')
        if p1 == -1: 
            p1 = lengthword
        p2 = firstword.find('/')
        if p2 == -1: 
            p2 = lengthword

        firstword = firstword[:min(p1,p2)]
        if firstword.find("target") != -1:
            pass
        else:
            seqid = myfunc.GetSeqIDFromAnnotation(firstword)
            idList.append(seqid)

    #print len(myfunc.uniquelist(idList))
    #print len(set(idList))

    idList = myfunc.uniquelist(idList)
    return idList
#}}}

def GetFullSeq(infile, hdl_seqdb, fpout):#{{{
    hdl = myfunc.ReadLineByBlock(infile)
    if hdl.failure:
        return (1, 0, 0)

    cntRetrieved = 0

    idList = []
    lines = hdl.readlines()
    while lines != None:
        idList += GetDatabaseIDList(lines)
        lines = hdl.readlines()
    hdl.close()

    idList = myfunc.uniquelist(idList)
    numID = len(idList)
    for seqid in idList:
        record = hdl_seqdb.GetRecord(seqid)
        if record:
            fpout.write(record)
            cntRetrieved += 1
        else:
            msg = "Failed to retrieve record for ID %s"
            print >> sys.stderr, msg%(seqid)
    return (0, numID, cntRetrieved)
#}}}
def GetSeqFromMSA(infile, outfile, hdl_seqdb):#{{{
    if not os.path.exists(infile):
        print >> sys.stderr, "infile %s does not exist."%(infile)
        return 1

    outdir = os.path.dirname(outfile)
    if outdir != "" and not os.path.exists(outdir):
        os.system("mkdir -p %s"%(outdir))

    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)

    filetype = GetFileType(infile)
    inputfile = ""
    if filetype == "fasta":
        tmpf = tempfile.mktemp()
        os.system("grep '^>' %s | sed 's/>//' > %s"%(infile, tmpf))
        inputfile = tmpf
        try:
            first_record = SeqIO.parse(open(infile, "rU"), "fasta").next()
            if first_record.id == "target":
                seq = (first_record.seq._data).replace("-","")
                fpout.write(">%s\n"%(first_record.description))
                fpout.write("%s\n"%(seq))
        except (IOError,ValueError,KeyError):
            pass
        
    elif filetype == "idlist":
        inputfile = infile
    else:
        print >> sys.stderr, "Unrecognized infile type"
        return 1

    (status, numseq, numRetrieved) = GetFullSeq(inputfile, hdl_seqdb, fpout)
    if status == 1:
        msg = "%s retrieved %d out of %d sequences. Failed to read."
    else:
        msg = "%s retrieved %d out of %d sequences. Succeeded."
    print msg%(infile, numRetrieved, numseq)

    if filetype == "fasta":
        os.system("rm -f %s"%(tmpf))

    myfunc.myclose(fpout)

    return 0
#}}}

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outfile = ""
    infile = ""
    seqdb = ""
    listfile = ""

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
            elif argv[i] in ["-o", "--o", "-outfile", "--outfile"]:
                outfile = argv[i+1]
                i += 2
            elif argv[i] in ["-i", "--i"] :
                infile = argv[i+1]
                i += 2
            elif argv[i] in ["-l", "--l"] :
                listfile = argv[i+1]
                i += 2
            elif argv[i] in ["-seqdb", "--seqdb"] :
                seqdb = argv[i+1]
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

    runList = [] # runList is a list of tuples (infile, outfile)


    if infile != "":
        runList.append((infile, outfile))
    if listfile != "":
        if os.path.exists(listfile):
            try:
                fpin = open(listfile, "rU")
                lines = fpin.readlines()
                fpin.close()
                for line in lines:
                    if not line or line[0]== "#":
                        continue
                    strs = line.split("\t")
                    infile = strs[0].strip()
                    try:
                        outfile = strs[1].strip()
                    except IndexError:
                        outfile = ""
                    runList.append((infile, outfile))
            except IOError:
                print >> sys.stderr, "Failed to read file %s"%(listfile)
        else:
            print >> sys.stderr, "listfile %s does not exist"%(listfile)

    numInput = len(runList)
    if numInput < 1:
        print >> sys.stderr, "No input set. Exit"
        return 1

    if seqdb == "":
        print >> sys.stderr, "seqdb not set."
        return 1
    elif not os.path.exists(seqdb+"0.db"):
        print >> sys.stderr, "seqdb %s does not exist."%(seqdb)
        return 1

    hdl_seqdb = myfunc.MyDB(seqdb)
    if hdl_seqdb.failure:
        print >> sys.stderr, "Failed to open seqdb %s"%(seqdb)
        return 1

    for (infile, outfile) in runList:
        GetSeqFromMSA(infile, outfile, hdl_seqdb)

    hdl_seqdb.close()
    
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


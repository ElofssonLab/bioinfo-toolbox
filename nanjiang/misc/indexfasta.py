#!/usr/bin/env python
# ChangeLog 2013-08-22 
#   solved the absolute path in the soft link
import os
import sys
import myfunc
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))
usage="""
Usage:  %s FASTA-FILE

Description: Index the fasta file, so that it can be read by my_extract.py.
             Sequences with redundant seqIDs will be removed.
Options:
  -dbname   STR  Set output dbname, the result will be output as
                 (default: basename(FILE))
                 $dbname.index
                 $dbname0.db
  -idtype   INT  Set how seqid is obtained from annotation line
                 0: using GetSeqIDFromAnnotation
                 1: first word delimited by " \\t"
  -q             Quiet mode
  -h, --help     Print this help message and exit

Created 2011-10-20, updated 2013-08-22, Nanjiang Shu 
"""%(progname)
MAX_DBFILE_SIZE = 50*1024*1024*1024
BLOCK_SIZE = 100000

def PrintHelp():
    print usage

def WriteIndexFasta(seqWithAnno, fpdb, dbname, fpindex, cntdbfile, #{{{
        record_offset, idSet, idtype):
    """Write sequence to indexed fasta file, sequences with redundant IDs are
    ignored"""
    if idtype == 0:
        seqid = myfunc.GetSeqIDFromAnnotation(seqWithAnno)
    elif idtype == 1:
        seqid = myfunc.GetFirstWord(seqWithAnno.lstrip(">"))
    if seqid in idSet:
        return (fpdb, record_offset)
    else:
        seqWithAnno+="\n"
        if fpdb == None:
            dbfile=dbname+"%d.db"%(cntdbfile)
            fpdb=open(dbfile, "wb")
            print "dbfile %s is created."%dbfile
        fpindex.write("%s %d %d %d\n"%(seqid, cntdbfile, record_offset,
            len(seqWithAnno)))
        fpdb.write("%s"%seqWithAnno)
        record_offset += len(seqWithAnno)
        idSet.add(seqid)
        return (fpdb,record_offset)
#}}}
def IndexFastaFile(infile, dbname, idtype): #{{{
    path_of_dbname = os.path.dirname(dbname)
    if path_of_dbname != "" and not os.path.exists(path_of_dbname):
        os.system("mkdir -p %s"%path_of_dbname)
    fpin = None
    try:
        fpin = open(infile,"rb")
    except IOError:
        print >> sys.stderr, "Failed to open file %s for read"%(infile)
        return 1

    cntdbfile = 0
    record_offset = 0
    
    dbfile = dbname+"%d.db"%(cntdbfile)
    indexfile = dbname+".index"
    try:
        fpindex = open(indexfile,"wb")
    except IOError:
        msg = "Failed to open indexfile {} for write"
        print >> sys.stderr, msg.format(indexfile)
        return 1
    fpdb = None

    fpindex.write("DEF_DBNAME %s\n"%dbname)
    idSet = set([])
    isFirstSeq = True
    totalLength = 0
    buff = fpin.read(BLOCK_SIZE)
    brokenSeqWithAnnoLine = ""; ##for the annotation line broken by BLOCK read
    while buff:
        beg=0
        end=0
        while 1:
            if brokenSeqWithAnnoLine:
                if brokenSeqWithAnnoLine[len(brokenSeqWithAnnoLine)-1] == "\n":
                    end = buff.find(">")
                else:
                    end = buff.find("\n>")
                if end >= 0:
                    seqWithAnno = brokenSeqWithAnnoLine + buff[0:end]
                    (fpdb, record_offset) = WriteIndexFasta(seqWithAnno,
                            fpdb,dbname, fpindex, cntdbfile, record_offset,
                            idSet, idtype)
                    brokenSeqWithAnnoLine = ""
                    beg=end
                else:
                    brokenSeqWithAnnoLine += buff
                    break

            beg = buff.find(">",beg)
            end = buff.find("\n>",beg+1)
            if beg >= 0:
                if end >=0:
                    seqWithAnno=buff[beg:end]
                    (fpdb, record_offset) = WriteIndexFasta(seqWithAnno, fpdb,
                            dbname, fpindex, cntdbfile, record_offset, idSet,
                            idtype)
                    beg=end
                else:
                    brokenSeqWithAnnoLine=buff[beg:]
                    break
            else:
                break

        if record_offset > MAX_DBFILE_SIZE:
            fpdb.close()
            fpdb = None
            cntdbfile +=1
            record_offset=0
        buff = fpin.read(BLOCK_SIZE)
    
    if brokenSeqWithAnnoLine:
        seqWithAnno = brokenSeqWithAnnoLine
        (fpdb, record_offset) = WriteIndexFasta(seqWithAnno, fpdb, dbname,
                fpindex, cntdbfile, record_offset, idSet, idtype)
    fpin.close()   
    fpindex.close()
    if fpdb != None:
        fpdb.close()

# post processing
    numIndexedSeq = len(idSet)
    numSeq = myfunc.CountFastaSeq(infile)

    cmd = "%s/my_indexformatconvert.py -f %s.index"%(g_params['binpath'], dbname)
    os.system(cmd)

    if numIndexedSeq == numSeq and cntdbfile == 0:
        print "%d sequences indexed successfully"

        dbfile = "%s0.db"%(dbname)
        dbfile_base = os.path.basename(dbfile)
        dbfile_path = os.path.dirname(dbfile)
        infile_path = os.path.dirname(infile)
        if dbfile_path == "": 
            dbfile_path = "."
        if infile_path == "": 
            infile_path = "."
        relpath = os.path.relpath(dbfile_path, infile_path)

        cmd = "rm -f %s"%infile
        print cmd
        os.system(cmd)

        cmd = "ln -s %s%s%s %s"%(relpath, os.sep, dbfile_base, infile)
        print cmd
        os.system(cmd)
    else:
        print >> sys.stderr, "numIndexedSeq (%d) conflicts with numSeq (%d)" %(
                numIndexedSeq, numSeq)
    return 0
#}}}
def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = False
    return g_params
#}}}
def main(g_params):#{{{
    numArgv = len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    dbname = ""
    infile = ""
    idtype = 0

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            infile = sys.argv[i]
            isNonOptionArg=False
            i += 1
        elif sys.argv[i] == "--":
            isNonOptionArg=True
            i += 1
        elif sys.argv[i][0] == "-":
            if sys.argv[i] ==  "-h" or  sys.argv[i] == "--help":
                PrintHelp()
                sys.exit()
            elif (sys.argv[i] in ["-outname", "--outname" , "-dbname",
                "--dbname"]):
                dbname = sys.argv[i+1]
                i += 2
            elif (sys.argv[i] in ["-idtype", "--idtype"]):
                idtype = int(sys.argv[i+1])
                i += 2
            elif sys.argv[i] == "-q":
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", sys.argv[i]
                return 1
        else:
            infile = sys.argv[i]
            i += 1

    if infile == "" or not os.path.exists(infile):
        print >> sys.stderr, "infile not set"
        return 1
    if dbname == "":
        dbname = infile

    if not idtype in [0,1]:
        print >> sys.stderr, "unrecognized idtype: %d"%(idtype)
        return 1

    rundir = os.path.dirname(sys.argv[0]) 
    g_params['binpath'] = rundir

    return IndexFastaFile(infile, dbname, idtype)
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter() 
    sys.exit(main(g_params))

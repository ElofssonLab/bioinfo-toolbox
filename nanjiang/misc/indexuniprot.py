#!/usr/bin/env python
# 
import os
import sys
usage="""
usage:  indexuniprot.py FILE

index the uniprot dat file, so that it can be read by my_extract.py
  
Options:
  -q             Quiet mode
  -h, --help     Print this help message and exit

Created 2012-06-05, updated 2012-06-05, Nanjiang Shu 
"""
MAX_DBFILE_SIZE=2*1024*1024*1024
BLOCK_SIZE=100000

def PrintHelp():
    print usage

def GetACList(recordContent):
    lines = recordContent.split("\n")
    str_accession = ""     # AC
    for line in lines:
        if len(line) > 2:
            tag = line[0:2]
            if tag == "AC":
                str_accession += line[5:]
            elif tag == "SQ":
                break
    if str_accession != "":
        acList = str_accession.split(";")
        acList = filter(None, acList)
        acList = [x.strip() for x in acList]
        return acList
    else:
        return []


    
def WriteIndex(fpindex, record_offset, recordCotent):#{{{
    "Write indexing file, "
    acList = GetACList(recordCotent)
    for ac in acList:
        fpindex.write("%s %d %d %d\n"%(ac, 0, record_offset, len(recordCotent)))
    record_offset += len(recordCotent)
    return record_offset
    
#}}}
def Read_UniprotData_from_buffer(buff, recordList, isEOFreached):#{{{
    if not buff:
        return ""
    unprocessedBuffer = ""
    beg = 0
    end = 0
    while 1:
        beg=buff.find("ID ",beg)
        if beg >= 0:
            end = buff.find("\n//\n",beg+1)
            if end >= 0:
                recordContent = buff[beg:end+4]
                recordList.append(recordContent)
                beg = end
            else:
                unprocessedBuffer = buff[beg:]
                break
        else:
            unprocessedBuffer = buff[end:]
            break
    if isEOFreached and unprocessedBuffer:
        recordContent = unprocessedBuffer
        recordList.append(recordContent)
        unprocessedBuffer = ""
    return unprocessedBuffer
#}}}
def IndexUniprotFile(infile, dbname): #{{{
    path_of_dbname = os.path.dirname(dbname)
    if path_of_dbname != "" and not os.path.exists(path_of_dbname):
        os.system("mkdir -p %s"%path_of_dbname)

    fpin = None
    try:
        fpin=open(infile,"rb")
    except IOError:
        print >> sys.stderr, "Failed to open file %s for read"%(infile)
        raise

    record_offset = 0
    
    indexfile = dbname + ".index"
    try:
        fpindex = open(indexfile,"wb")
    except IOError:
        print >> sys.stderr, "Failed to open indexfile %s for write"%(indexfile)
        raise

    fpindex.write("DEF_DBNAME %s\n"%dbname)
    cntRecord = 0

    unprocessedBuffer = ""; ##for the annotation line broken by BLOCK read
    isEOFreached = False

    while 1:
        buff = fpin.read(BLOCK_SIZE)
        if len(buff) < BLOCK_SIZE:
            isEOFreached = True
        buff = unprocessedBuffer + buff
        recordList = []
        unprocessedBuffer = Read_UniprotData_from_buffer(
                buff, recordList, isEOFreached)
        if len(recordList) > 0: 
            cntRecord += len(recordList)
            for rd in recordList:
                record_offset = WriteIndex(fpindex, record_offset, rd)
        if isEOFreached == True:
            break
    fpin.close()

    print "%d record indexed" % (cntRecord)
    print "index file %s output" % (indexfile)

    fpindex.close()
    dbfile = infile + "0.db"
    os.system("ln -sf %s %s" %(infile, dbfile))
    
    return 0

#}}}

def main(g_params):
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    isQuiet=False
    dbname=""
    infile=""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            infile=sys.argv[i]
            isNonOptionArg=False
            i += 1
        elif sys.argv[i] == "--":
            isNonOptionArg=True
            i += 1
        elif sys.argv[i][0] == "-":
            if sys.argv[i] ==  "-h" or  sys.argv[i] == "--help":
                PrintHelp()
                return 1
            elif sys.argv[i] == "-outname" or sys.argv[i] == "--outname" or sys.argv[i] == "-dbname" or sys.argv[i] == "--dbname":
                dbname=sys.argv[i+1]
                i += 2
            elif sys.argv[i] == "-q":
                isQuiet=True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", sys.argv[i]
                return 1
        else:
            infile=sys.argv[i]
            i += 1

    if infile == "" or not os.path.exists(infile):
        print >> sys.stderr, "infile not set"
        return 1
    if dbname == "":
        dbname = os.path.basename(infile)

    IndexUniprotFile(infile,dbname)

if __name__ == '__main__' :
    g_params = {}
    sys.exit(main(g_params))


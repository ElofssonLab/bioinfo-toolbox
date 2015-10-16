#!/usr/bin/env python
# 
import os
import sys
import myfunc
usage="""
usage: index_topconssingle_result.py  FILE

Index the result of topocons_single all.info
  
Options:
  -dbname   STR  Set output dbname, the result will be output as
                 (default: basename(FILE))
                 $dbname.index
                 $dbname0.db
  -q             Quiet mode
  -h, --help     Print this help message and exit

Created 2011-12-09, updated 2011-12-09, Nanjiang Shu 
"""
MAX_DBFILE_SIZE=2*1024*1024*1024;
BLOCK_SIZE=100000;

def PrintHelp():
    print usage;

def WriteIndex(fpdb, record_offset) =  WriteIndex(recordList, fpdb, dbname, #{{{
                    fpindex, cntdbfile, record_offset, processedTopoIDSet);
    "Write sequence to indexed fasta file, sequences with redundant IDs are ignored"
    seqid = myfunc.GetSeqIDFromAnnotation(seqWithAnno);
    if seqid in idSet:
        return (fpdb, record_offset);
    else:
        seqWithAnno+="\n";
        if fpdb == None:
            dbfile=dbname+"%d.db"%(cntdbfile);
            fpdb=open(dbfile, "wb");
            print "dbfile %s is created."%dbfile;
        fpindex.write("%s %d %d %d\n"%(seqid, cntdbfile, record_offset, len(seqWithAnno)));
        fpdb.write("%s"%seqWithAnno);
        record_offset += len(seqWithAnno);
        idSet.add(seqid);
        return (fpdb,record_offset);
    
#}}}
def Index_topconssingle_result(infile, dbname): #{{{
    path_of_dbname = os.path.dirname(dbname);
    if path_of_dbname != "":
        os.system("mkdir -p %s"%path_of_dbname);

    fpin = None;
    try:
        fpin=open(infile,"rb");
    except IOError:
        print >> sys.stderr, "Failed to open file %s for read"%(infile);
        raise;

    cntdbfile=0;
    record_offset=0;

    dbfile=dbname+"%d.db"%(cntdbfile);
    indexfile=dbname+".index";
    try:
        fpindex=open(indexfile,"wb");
    except IOError:
        print >> sys.stderr, "Failed to open indexfile %s for write"%(indexfile);
        raise;
    fpdb = None;

    fpindex.write("DEF_DBNAME %s\n"%dbname);

    unprocessedBuffer="";
    isEOFreached = False;
    processedTopoIDSet = set([]);
    while 1:
        buff = fpin.read(BLOCK_SIZE);
        if len(buff) < BLOCK_SIZE:
            isEOFreached=True;
        buff = unprocessedBuffer + buff;
        recordList = [];
        unprocessedBuffer = Read_topconssingle_result_from_buffer(
                buff, recordList, isEOFreached);
        if len(recordList) > 0: 
            (fpdb, record_offset) =  WriteIndex(recordList, fpdb, dbname,
                    fpindex, cntdbfile, record_offset, processedTopoIDSet);
        if isEOFreached == True:
            break;
    fpin.close();

    fpindex.close();
    if fpdb != None:
        fpdb.close();
    
    return 0;

#}}}
def main(g_params):#{{{
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1
    
    infile="";

    i = 1;
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            infile=sys.argv[i];
            isNonOptionArg=False;
            i += 1;
        elif sys.argv[i] == "--":
            isNonOptionArg=True;
            i += 1;
        elif sys.argv[i][0] == "-":
            if sys.argv[i] ==  "-h" or  sys.argv[i] == "--help":
                PrintHelp();
                return 1
            elif sys.argv[i] in ["-outname", "--outname", "-dbname", "--dbname"]:
                g_params['dbname'] = sys.argv[i+1];
                i += 2;
            elif sys.argv[i] == "-q":
                g_params['isQuiet'] = True;
                i += 1;
            else:
                print >> sys.stderr, "Error! Wrong argument:", sys.argv[i];
                return 1;
        else:
            infile = sys.argv[i];
            i += 1

    if infile == "" or not os.path.exists(infile):
        print >> sys.stderr, "infile not set";
        return 1;
    if g_params['dbname'] == "":
        dbname = os.path.basename(infile);
    else:
        dbname = g_params['dbname'];

    Index_topconssingle_result(infile, dbname);
#}}}

def InitGlobalPrameter():#{{{
    g_params = {};
    g_params['isQuiet'] = False;
    g_params['dbname'] = "";
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalPrameter();
    sys.exit(main(g_params));

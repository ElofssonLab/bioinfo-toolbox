#!/usr/bin/env python
# extract data by tag ID from formatted db by my_formatdb.py
import os
import sys
import myfunc
from array import array
from mydb_common import *
#import numpy as np
#import cProfile

# ChangeLog 2011-11-15 #{{{
#    added the following options
#    -sn 
#    -splitall 
#    removed the global declaration, use g_params instead. 
# ChangeLog 2013-02-05 
#   1. for text format index file, read in to array, memory efficient
#   2. for binary format index file, can read offset of both array('I') and
#   array('L')
#}}}

progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

usage="""
Usage: %s [ID [ID ... ]] or  [-l IDLISTFILE ]
       %s [-splitall -outpath DIR]
       %s [-split -outpath DIR [-dataext STR] [-dataprefix STR]]
       %s [-mode STR]
       %s -dbname DBNAME
Description:
    Extract individual data from the formatted database if -split or -splitall
    is enabled, the individual file will be output to
    $outpath/$dataprefix$ID$dataext

Options:
  -l  LISTFILE  Set the idListFile
  -dbname  STR  Set the input database name
  -o      FILE  Set output file
  -outpath DIR  Set output path, (default: "./")
  -de, -dataext STR  
                Set the extension of datafile, including '.', e.g. '.fa'
                (default: read from database or "")
  -dp, -dataprefix STR  
                Set the prefix of datafile (default: read from database or "")
  -format  STR  Set the format of the index file, binary or text (default: text)
  -split        Split the output to individual files for every supplied ID
  -splitall     Split the output to individual files for all IDs in the
                database
  -sn, -shownumrecord
                Show the number of records in the database
  -wall         Print Warning information
  -h, --help    Print this help message and exit

Additional options
  -maxinput INT    Max input for using linear string search (default, 10)

Created 2011-09-18, updated 2013-02-06, Nanjiang Shu

Examples:

# extract the entry ID1 from the database mydb
    my_extractdb.py -dbname mydb ID1

# split entries in the idlist.txt from mydb to the folder outdir
    my_extractdb.py -dbname mydb -l idlist.txt -split -outpath outdir

# split all entries in the database mydb to folder outdir
    my_extractdb.py -dbname mydb -splitall -outpath outdir
"""%(progname, wspace, wspace, wspace, wspace)

def PrintHelp():
    print usage

def ReadNumRecord_binary(indexfile):#{{{
# return numRecord
    try:
        fpin=open(indexfile, "rb");
        vI=array('I');
        vI.fromfile(fpin,1);
        dumpedtext=fpin.read(vI[0]);
        #read in other information
        vI=array('I');
        vI.fromfile(fpin,1);
        dumpedidlist=fpin.read(vI[0]);
        vI=array('I');
        vI.fromfile(fpin,1);
        numRecord=vI[0];
        fpin.close();
        return numRecord;
    except IOError:
        msg = "{}: Failed to read indexfile {}"
        print >> sys.stderr, msg.format(sys.argv[0],indexfile);
        raise;
#}}}
def ReadRecordNumber(indexfile):#{{{
    numRecord=0;
    if g_params['formatindex']==FORMAT_TEXT:
        (headerinfo, dbfileindexList, index) = ReadIndex_text(indexfile);
        numRecord = len(index);
    else:
        numRecord = ReadNumRecord_binary(indexfile);
    return numRecord;
#}}}
def GetOutfileName(idd,origext,origprefix):#{{{
    prefix="";
    ext="";
    if g_params['dataprefix'] != "":
        prefix=g_params['dataprefix'];
    elif origprefix != "":
        prefix=origprefix;

    if  g_params['dataext'] != "":
        ext=g_params['dataext'];
    elif origext != "":
        ext=origext;
    outfile=g_params['outpath']+os.sep+prefix+idd+ext;
    return outfile;
#}}}
def ReadIndex(indexfile):#{{{
# return (indexList, headerinfo, dbfileindexList)
    if g_params['formatindex'] == FORMAT_TEXT:
        return ReadIndex_text(indexfile, g_params['isPrintWarning'])
    else:
        return ReadIndex_binary(indexfile, g_params['isPrintWarning'])
#}}}
def ExtractDBForEachIDWithIndexDict(idList, origext, origprefix, #{{{
        indexList, indexDict, fpdbList, indexfile, fpout):
# @params
# idd               ID
# indexDict         A dictionary storing record location for each ID
# fpdbList          A list of FILE handles for dbfiles
# indexedKeysSet    Unique set of IDs in the index file
# indexfile         Index file
    outfile="";

    indexedIDList = indexList[0]
    v1 = indexList[1]
    v2 = indexList[2]
    v3 = indexList[3]
    isSplit = g_params['isSplit']
    isSplitAll = g_params['isSplitAll']
    isQuiet = g_params['isQuiet']
    for idd in idList:
        try: 
            idxItem = indexDict[idd]
            if isSplit or isSplitAll:
                outfile=GetOutfileName(idd,origext,origprefix);
                fpout = open (outfile,"wb");

            fpdb = fpdbList[v1[idxItem]];
            fpdb.seek(v2[idxItem]);
            data = fpdb.read(v3[idxItem]);

            fpout.write(data) ;
            if isSplit or isSplitAll:
                fpout.close()
                if not isQuiet:
                    print >> sys.stdout, "%s output"%(outfile);
        except (KeyError, IndexError):
            msg = "{}: ID {} not found in the index file {}"
            print >> sys.stderr, msg.format(sys.argv[0], idd, indexfile );
        except IOError:
            msg = "{}: Failed to read or write record {}"
            print >> sys.stderr, msg.format(sys.argv[0], idd);
#}}}
def ExtractDBForEachIDWithIndexList(idList, origext, origprefix, #{{{
        indexList, fpdbList, indexfile, fpout):
# @params
# idd               ID
# indexList         A list of 4-tuples (ID, dbfileindex, offset, blocksize)
# fpdbList          A list of FILE handles for dbfiles
# indexedKeysSet    Unique set of IDs in the index file
# indexfile         Index file
    outfile = ""
    indexedIDList = indexList[0]
    v1 = indexList[1]
    v2 = indexList[2]
    v3 = indexList[3]
    isSplit = g_params['isSplit']
    isSplitAll = g_params['isSplitAll']
    isQuiet = g_params['isQuiet']
    for idd in idList:
        try:
            idxItem = indexedIDList.index(idd)
            if isSplit or isSplitAll:
                outfile = GetOutfileName(idd,origext,origprefix)
                fpout = open (outfile,"wb")
            fpdb = fpdbList[v1[idxItem]];
            fpdb.seek(v2[idxItem]);
            data = fpdb.read(v3[idxItem]);

            fpout.write(data) ;
            if isSplit or isSplitAll:
                fpout.close()
                if not isQuiet:
                    print >> sys.stdout, "%s output"%(outfile);
        except ValueError:
            msg = "{}: ID {} not found in the index file {}"
            print >> sys.stderr, msg.format(sys.argv[0], idd, indexfile );
        except IOError:
            msg = "{}: Failed to read or write record {}"
            print >> sys.stderr, msg.format(sys.argv[0], idd);
#}}}
def ExtractDB(idList, g_params, fpout):#{{{
    dbname = g_params['dbname'];
    formatindex = g_params['formatindex']
    (indexfile, formatindex) = GetIndexFile(dbname, formatindex);
    g_params['formatindex'] = formatindex
    if indexfile == "":
        msg = "{}: can not find indexfile for database {}. Exit" 
        print >> sys.stderr, msg.format(sys.argv[0],dbname);
        return 1;

    (indexList, headerinfo, dbfileindexList) = ReadIndex(indexfile)
    (origdbname, origversion, origext, origprefix) = headerinfo
    numRecord = len(indexList[0])
    if numRecord <= 0:
        msg = "{}: Read index file {} failed. Exit."
        print >> sys.stderr, msg.format(sys.argv[0],indexfile);
        return 1;
    if g_params['isSplitAll']:
        idList = indexList[0]

    fpdbList = []
    for i in dbfileindexList:
        dbfile = dbname+"%d.db"%(i);
        try:
            fpdbList.append(open(dbfile,"rb"));
        except IOError:
            msg = "{}: Failed to read dbfile {}"
            print >> sys.stderr, msg.format(sys.argv[0],dbfile);
            raise

    if g_params['typeindex'] == TYPE_LIST:
        ExtractDBForEachIDWithIndexList(idList, origext, origprefix, indexList,
                fpdbList, indexfile, fpout)
    else:
        indexDict = {}
        #t = np.arange(0, numRecord, 1, dtype=np.int32)
        for i in xrange(numRecord):
            indexDict[indexList[0][i]] = i
        ExtractDBForEachIDWithIndexDict(idList, origext, origprefix, indexList,
                indexDict, fpdbList, indexfile, fpout);

    for fp in fpdbList:
        fp.close();
    return 0;

#}}}
def main(g_params):#{{{
    numArgv=len(sys.argv);
    if numArgv < 2:
        PrintHelp();
        return 1;

    idListFile = "" ;
    idList=[];

    i = 1;
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            idList.append(sys.argv[i]);
            isNonOptionArg=False;
            i += 1;
        elif sys.argv[i] == "--":
            isNonOptionArg=True;
            i += 1;
        elif sys.argv[i][0] == "-":
            if sys.argv[i] in ["-h", "--help"]:
                PrintHelp();
                return 1
            elif (sys.argv[i] in ["-db","--db","-dbname","--dbname"]):
                g_params['dbname']=sys.argv[i+1];
                i += 2;
            elif (sys.argv[i] in ["-outpath" ,"--outpath"]):
                g_params['outpath']=sys.argv[i+1];
                i += 2;
            elif (sys.argv[i] in ["-de", "--de", "-dataext" , "--dataext"]):
                g_params['dataext']=sys.argv[i+1];
                i += 2;
            elif (sys.argv[i] in [ "-dp", "--dp", "-dataprefix" ,
                "--dataprefix"]):
                g_params['dataprefix']=sys.argv[i+1];
                i += 2;
            elif (sys.argv[i] in ["-maxinput" ,  "--maxinput"]):
                g_params['MAX_SMALL_NUMBER_INPUT']=int(sys.argv[i+1]);
                i += 2;
            elif (sys.argv[i] in ["-format" , "--format"]):
                if sys.argv[i+1].lower()[0]== "b":
                    g_params['formatindex']=FORMAT_BINARY;
                else:
                    g_params['formatindex']=FORMAT_TEXT;
                i += 2;
            elif (sys.argv[i] in ["-l" , "--l"]) :
                idListFile=sys.argv[i+1];
                i += 2;
            elif (sys.argv[i] in [ "-o" , "--o" ]):
                g_params['outfile']=sys.argv[i+1];
                i += 2;
            elif (sys.argv[i] in ["-split" , "--split"]):
                g_params['isSplit']=True;
                i += 1;
            elif (sys.argv[i] in ["-splitall" , "--splitall"]):
                g_params['isSplitAll']=True;
                i += 1;
            elif (sys.argv[i] in ["-sn" ,"--sn", "-shownumrecord",
                "--shownumrecord"]):
                g_params['isShowNumRecord']=True;
                i += 1;
            elif (sys.argv[i] in ["-q","--q","-quiet","--quiet"]):
                g_params['isQuiet']=True;
                i += 1;
            else:
                print >> sys.stderr, "Error! Wrong argument:", sys.argv[i];
                return 1;
        else:
            idList.append(sys.argv[i])
            i += 1

    if g_params['dbname']=="":
        print >> sys.stderr, "dbname is not set. Exit.";
        return 1;
    elif g_params['isShowNumRecord'] == True:
        (indexfile, g_params['formatindex']) = GetIndexFile(g_params['dbname'],
                g_params['formatindex']);
        if indexfile == "":
            msg = "{}: can not find indexfile for database {}. Exit"
            print >> sys.stderr, msg.format(sys.argv[0],g_params['dbname'])
            return 1
        else:
            numRecord = ReadRecordNumber(indexfile);
            if numRecord < 0:
                return 1;
            else :
                print numRecord;
                return 0;
    else:
        if idListFile != "":
            try:
                fp=open(idListFile,"r");
                idList+=fp.read().split();
                fp.close();
            except IOError:
                print >> sys.stderr, "Failed to read idlistfile %s."
        if (not g_params['isSplitAll']) and len(idList) <= 0:
            print >> sys.stderr,"no ID has been set, Exit." ;
            return 1;
        if g_params['isSplit'] or g_params['isSplitAll']:
            os.system("mkdir -p %s"%(g_params['outpath']));

        fpout = sys.stdout
        if (not (g_params['isSplit'] or g_params['isSplitAll'])):
            fpout = myfunc.myopen(g_params['outfile'], sys.stdout, "wb", False)

        idList = myfunc.uniquelist(idList)

        if len(idList) <= g_params['MAX_SMALL_NUMBER_INPUT']:
            g_params['typeindex']=TYPE_LIST;

        if ExtractDB(idList, g_params, fpout) != 0:
            print >> sys.stderr, "ExtractDB failed.";
        myfunc.myclose(fpout)
        return 0
#}}}

def InitGlobalParameter():#{{{
    g_params = {};   # define global parameters
    g_params['isQuiet'] = False;
    g_params['dbname']="";
    g_params['outpath']="./";
    g_params['dataext']="";
    g_params['dataprefix']="";
    g_params['isSplit']=False;
    g_params['isSplitAll']=False;
    g_params['isShowNumRecord']=False;
    g_params['formatindex']=FORMAT_BINARY;
    g_params['outfile']="";
    g_params['typeindex']=TYPE_DICT;

    # if the number of input ID is < 10, use list instead of dictionary
    g_params['MAX_SMALL_NUMBER_INPUT'] = 10; 

    g_params['isPrintWarning'] = False
    return g_params
#}}}
if __name__ == "__main__":
    g_params = InitGlobalParameter()
#    cProfile.run("main()");
    sys.exit(main(g_params));

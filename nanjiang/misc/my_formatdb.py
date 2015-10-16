#!/usr/bin/env python
# format a number of separated files with unique IDs to dumped file with
# formatted index

# Log:#{{{
# Log 2011-09-21 
# Version 1.0:
#   Initial version
# Version 1.1:  2011-09-21 
#   Add binary format for index file
# Version 1.2:  2011-09-21
#   in header, add prefix and ext, so that when doing extraction, they do not
#   need to be supplied
# Version 1.3: 2011-10-10
#   dataext should also contain the '.', the file name will therefore be
#   $dataprefix$ID$dataext, do not add .
# ChangeLog 2011-10-20
#   the prefix of dbname will be created is not exist
#   dbname can be supplied with prefix path
# ChangeLog 2011-11-09 
#     dataext can be supplied as ""
# ChangeLog 2011-11-19
#     The default lastDBFileIndex is set to 0 (before was -1)
#     when numRecord == 0, do not read further
# ChangeLog 2013-02-06
#     Faster reading of index file and more memory efficient. 
#     Also for binary format, array('L') is used instead of array('I') when
#     maximal offset exceeds 1.5GB
# ChangeLog 2014-12-17
#     MAXDBFILESIZE is changed to 8Gb, so the allowed database size will be 
#     255*8 = 2TB
#     The max number of splitted db files is 255, that is = array('B').itemsize
#}}}

import os
import sys
from array import array
from mydb_common import *
import myfunc

progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))
usage="""
Usage: %s [ID [ID ... ]] or  [-l IDLISTFILE ]
       %s [-mode STR]  [-dataprefix STR] [-dataext STR]
       %s -datapath DIR -dbname DBNAME

Description: Format a number of separated files with unique IDs to dumped file
             with formatted index. Outfile will be $DBNAME0.db and
             $DBNAME.index, if the size of db file exceeds 1GB, new files will
             be created and named as $DBNAME1.db $DBNAME2.db ...
             datafile=$datapath/$dataprefix$ID$dataext

             Note that the largest size of individual file is 1.5Gb 
             The maximum size of the database is 2TB
Options:
  -l   IDLISTFILE   Set the idListFile
  -datapath   DIR   Set the path containing data
  -dataext    STR   Set the extension that determines the datafile
                    including '.', e.g. '.fa', (default: "")
  -dataprefix STR   Set the prefix of datafile, (default: "")
  -dbname  DBNAME   Set the output database name
  -mode       STR   Set the mode, new or append (default: new)
  -format     STR   Set the format of the index file, binary or text 
                    (default: binary)
  -wall             Print warning message, (default: no)
  -h, --help        Print this help message and exit

Created 2011-09-18, updated 2014-12-17, Nanjiang Shu

Examples:
    my_formatdb.py -l idlist.txt -datapath test -dbname out1/dumpdb
"""%(progname, wspace, wspace)

MAXDBFILESIZE=1024*1024*1024*8; # 8GB # changed 2014-12-17
MODE_NEW=0
MODE_APPEND=1

def PrintHelp():
    print usage
def WriteDBForEachID(idd, origIdListSet, datapath, dataext, dataprefix, #{{{
        cntdbfile, size_dbfile, dbname, fpdb, indexList):
    isQuiet = g_params['isQuiet']
    datafile = datapath+os.sep+dataprefix+idd+dataext
    dbfile = dbname+"%d.db"%(cntdbfile)
    if idd in origIdListSet:
        msg = "ID {} already exists in the database. Ignore."
        print >> sys.stderr, msg.format(idd)
    elif not os.path.exists(datafile):
        print >> sys.stderr,"datafile %s does not exist. Ignore." % (datafile)
    else:
        try:
            if fpdb == None:
                fpdb=open(dbfile, "wb")
                if not isQuiet:
                    print "dbfile %s is created."%dbfile
            fpin=open(datafile,"rb")
            data=fpin.read()
            fpin.close()
            sizeblock=len(data)
            offset=fpdb.tell()

            indexList[0].append(idd)
            indexList[1].append(cntdbfile)
            indexList[2].append(offset)
            indexList[3].append(sizeblock)
#            indexList.append((idd,cntdbfile, offset, sizeblock))

            print >> fpdb, data
            size_dbfile += sizeblock
            if size_dbfile >= MAXDBFILESIZE:
                fpdb.close()
                cntdbfile += 1
                fpdb = None
                size_dbfile = 0
        except IOError:
            msg = "{}: Failed to write to dbfile {}"
            print >> sys.stderr, msg.format(sys.argv[0],datafile)
            raise
    return (fpdb, cntdbfile, size_dbfile)
#}}}
def WriteDB(idList, origIdListSet, datapath, dataext, dataprefix, #{{{
        cntdbfile, size_dbfile, dbname, fpdb, indexList):
    indexList.append([])
    indexList.append(array('B'))
    indexList.append(array('L'))
    indexList.append(array('I'))
    for idd in idList:
        (fpdb, cntdbfile, size_dbfile) = WriteDBForEachID(idd,origIdListSet,
                datapath, dataext, dataprefix, cntdbfile, size_dbfile, dbname,
                fpdb, indexList)
    if fpdb != None:
        fpdb.close()
    return 0
#}}}

def CreateNewFormattedDB(idList,datapath,dataext, dataprefix, #{{{
        indexfile, dbname):
    isQuiet = g_params['isQuiet']
    size_dbfile=0
    cntdbfile=0
    dbfile=dbname+"%d"%cntdbfile+".db"
    fpindex=open(indexfile,"wb")
    fpdb=open(dbfile,"wb")
    indexList = []
    if not isQuiet:
        print "dbfile %s is created."%dbfile

    indexFileHeaderText=["DEF_VERSION %s"%version, "DEF_DBNAME %s"%dbname,
            "DEF_EXTENSION %s"%dataext,"DEF_PREFIX %s"%dataprefix ]
    formatindex = g_params['formatindex']
    WriteIndexHeader(indexFileHeaderText, formatindex, fpindex)
    WriteDB(idList, set([]), datapath, dataext, dataprefix, cntdbfile,
            size_dbfile, dbname, fpdb, indexList)
#Write indexList
    WriteIndexContent(indexList, formatindex, fpindex)
    fpindex.close()
    if not isQuiet:
        msg = "{} records have been added to the database {}."
        print msg.format(len(indexList[0]), dbname)
        print "indexfile %s output."%(indexfile)
#}}}
def AppendFormattedDB(idList,datapath,dataext, dataprefix, #{{{
        indexfile,dbname):
    isQuiet = g_params['isQuiet']
    formatindex = g_params['formatindex']
    origIDListSet=set([])
    lastDBFileIndex = 0
    indexList=[]
    indexFileHeaderText=[]
    if formatindex==FORMAT_TEXT:
        (indexList, headerinfo, dbfileindexList) = ReadIndex_text(indexfile,
                g_params['isPrintWarning'])
        origIdListSet = set(indexList[0])
        lastDBFileIndex = indexList[1][len(indexList[0])-1]
    else:
        (indexList, headerinfo, dbfileindexList) = ReadIndex_binary(indexfile,
                g_params['isPrintWarning'])
        origIdListSet = set(indexList[0])
        indexFileHeaderText = GetIndexFileHeaderText(headerinfo)
        lastDBFileIndex = indexList[1][len(indexList[0])-1]
    if lastDBFileIndex < 0:
        msg = "Fatal: Read index file {} failed in function {}. Exit."
        print >> sys.stderr, msg.format(indexfile, sys._getframe().f_code.co_name)
        return 1

    numOrigRecord = len(indexList[0])
    cntdbfile=lastDBFileIndex
    dbfile=dbname+"%d"%cntdbfile+".db"
    fpdb=open(dbfile,"ab+")
    fpdb.seek(0,os.SEEK_END)
    size_dbfile=fpdb.tell()
    WriteDB(idList,origIdListSet, datapath, dataext, dataprefix, cntdbfile,
            size_dbfile, dbname, fpdb, indexList)
    if len(indexList[0]) > numOrigRecord:
        fpindex=None
        if formatindex==FORMAT_TEXT:
            fpindex=open(indexfile,"ab+")
        else:
            fpindex=open(indexfile,"wb")
        if formatindex == FORMAT_BINARY:
            WriteIndexHeader(indexFileHeaderText, formatindex, fpindex)
        WriteIndexContent(indexList, formatindex, fpindex)
        fpindex.close()
        if not isQuiet:
            msg = "%d records have been appended to the database %s"
            print msg.format(len(indexList[0])-numOrigRecord, dbname)
#}}}

def FormatDB(idList, datapath, dataext, dataprefix, dbname):#{{{
    mode = g_params['mode']
    formatindex = g_params['formatindex']
    path_of_dbname = os.path.dirname(dbname)
    if path_of_dbname != "" and not os.path.exists(path_of_dbname):
        os.system("mkdir -p %s"%path_of_dbname)
    indexfile=""
    if formatindex == FORMAT_BINARY:
        indexfile=dbname+".indexbin"
    else:
        indexfile=dbname+".index"
    if mode == MODE_NEW or not os.path.exists(indexfile):
        CreateNewFormattedDB(idList,datapath, dataext,dataprefix, indexfile,
                dbname)
    elif mode == MODE_APPEND:
        AppendFormattedDB(idList,datapath,dataext,dataprefix, indexfile,
                dbname)
    else:
        print >> sys.stderr,("%s: Unrecognized mode (%d). Exit." %
                (sys.argv[0],mode))

#}}}

def main(g_params):#{{{
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    dbname=""
    idListFile=""
    datapath=""
    dataext=""
    dataprefix=""
    idList=[]

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            idList.append(sys.argv[i])
            isNonOptionArg=False
            i += 1
        elif sys.argv[i] == "--":
            isNonOptionArg=True
            i += 1
        elif sys.argv[i][0] == "-":
            if sys.argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif sys.argv[i] in ["-dbname", "--dbname"]:
                dbname=sys.argv[i+1]
                i += 2
            elif sys.argv[i] in ["-datapath",  "--datapath"]:
                datapath = sys.argv[i+1]
                i += 2
            elif sys.argv[i] in ["-dataext", "--dataext"]:
                dataext=sys.argv[i+1]
                i += 2
            elif sys.argv[i] in ["-dataprefix", "--dataprefix"]:
                dataprefix=sys.argv[i+1]
                i += 2
            elif sys.argv[i] in ["-mode", "--mode"]:
                if sys.argv[i+1].lower()[0]== "n":
                    g_params['mode'] = MODE_NEW
                else:
                    g_params['mode'] = MODE_APPEND
                i += 2
            elif sys.argv[i] in ["-format", "--format"]:
                if sys.argv[i+1].lower()[0]== "b":
                    g_params['formatindex'] = FORMAT_BINARY
                else:
                    g_params['formatindex'] = FORMAT_TEXT
                i += 2
            elif sys.argv[i] in ["-l", "--l", '-list', '--list']:
                idListFile=sys.argv[i+1]
                i += 2
            elif sys.argv[i] in ["-q", '--q', '-quiet', '--quiet'] :
                g_params['isQuiet'] = True
                i += 1
            elif sys.argv[i] in ["-wall", '--wall'] :
                g_params['isPrintWarning'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", sys.argv[i]
                return 1
        else:
            idList.append(sys.argv[i])
            i += 1

    if dbname=="":
        print >> sys.stderr, "dbname is not set. Exit."
        return -1
#     if dataext=="":
#         print >> sys.stderr, "dataext is not set. Exit."
#         sys.exit()
    if datapath == "" or not os.path.exists(datapath):
        print >> sys.stderr, "datapath is not set or not existing. Exit."
        return -1

    if idListFile != "":
        try:
            fp=open(idListFile,"r")
            idList+=fp.read().split()
            fp.close()
        except IOError:
            print >> sys.stderr, "Failed to read idListFile %s"%idListFile
            pass
    if len(idList) <= 0:
        print >> sys.stderr,"no ID has been set, Exit." 
        return -1
    idList = myfunc.uniquelist(idList)
    FormatDB(idList, datapath, dataext, dataprefix, dbname)
#}}}

def InitGlobalParameter():#{{{
    g_params = {}    # define global parameters, in this way, one avoid of
                     # using global statement
    g_params['mode'] = MODE_NEW
    g_params['formatindex'] = FORMAT_BINARY; # format of the index file
    g_params['isQuiet'] = False
    g_params['isPrintWarning'] = False
    return g_params
#}}}

if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

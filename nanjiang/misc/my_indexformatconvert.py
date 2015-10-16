#!/usr/bin/env python
# convert index file format 
import os,sys
from array import array
import myfunc
from mydb_common import *
#ChangeLog 2013-02-05 
#   array('I') or array('L'), solve the problem of which the db file size is
#   greater than 2GB
#   after version 1.4,
#   vIarray changed to [array('B'), array('L'),array('I')]
#                       dbfileindex   offset    filesize
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))
usage="""
Usage: %s FILE [FILE ...]

Description: Convert the format of index file created by my_formatdb.py
    -f, -force       Force index file conversion
    -wall            Print Warning information

Created 2011-09-21, updated 2013-02-05, Nanjiang Shu  
"""%(progname)
def PrintHelp():
    print usage
def ConvertTextToBinary(infile):#{{{
    (indexList, headerinfo, dbfileindexList) = ReadIndex_text(infile,
            g_params['isPrintWarning'])
    indexFileHeaderText = GetIndexFileHeaderText(headerinfo)
    outfile = os.path.splitext(infile)[0]+".indexbin"
    if os.path.exists(outfile):
        if (not g_params['isforcewritten'] and myfunc.confirm(
            prompt ='outfile %s already exists, force overwrite'%(outfile),
            resp=False) == False):
            return
    fpindex = open(outfile,"wb")
    WriteIndexHeader(indexFileHeaderText, FORMAT_BINARY, fpindex)
    WriteIndexContent(indexList, FORMAT_BINARY, fpindex)
    fpindex.close()
    if not g_params['isQuiet']:
        print "%s ==> %s"%(infile,outfile)
#}}}
def ConvertBinaryToText(infile):#{{{
    (indexList, headerinfo, dbfileindexList) = ReadIndex_binary(infile,
            g_params['isPrintWarning'])
    indexFileHeaderText = GetIndexFileHeaderText(headerinfo)
    outfile = os.path.splitext(infile)[0]+".index"
    if os.path.exists(outfile):
        if (not g_params['isforcewritten'] and myfunc.confirm(
            prompt='outfile %s already exists, force overwrite'%(outfile),
            resp=False) == False):
            return
    fpindex=open(outfile,"wb")
    WriteIndexHeader(indexFileHeaderText, FORMAT_TEXT, fpindex)
    WriteIndexContent(indexList, FORMAT_TEXT, fpindex)
    fpindex.close()
    if not g_params['isQuiet']:
        print "%s ==> %s"%(infile,outfile)
#}}}
def ConvertFormat(infile):#{{{
    (rootnamewithpath, ext) = os.path.splitext(infile)
    dbfile = "%s0.db"%(rootnamewithpath)
    if not os.path.exists(dbfile):
        msg = "dbfile for indexfile {} does not exist. Ignore."
        print >> sys.stderr, msg.format(infile)
        return 1
    dbfilesize = os.path.getsize(dbfile)
    if dbfilesize > LargeFileThresholdSize:
        g_params['isLargeFile'] = True
    if ext == ".index":
        ConvertTextToBinary(infile)
    elif ext == ".indexbin":
        ConvertBinaryToText(infile)
    else:
        msg = "{}: Unrecognized file type for file {}. Ignore."
        print >> sys.stderr, msg.format%(sys.argv[0], infile)
#}}}
def main(g_params):#{{{
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1


    fileList = []
    i = 1
    isNonOptionArg = False
    while i < numArgv:
        if isNonOptionArg == True:
            fileList.append(sys.argv[i])
            isNonOptionArg=False
            i += 1
        elif sys.argv[i] == "--":
            isNonOptionArg=True
            i += 1
        elif sys.argv[i][0] == "-":
            if sys.argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif sys.argv[i] == "-q":
                g_params['isQuiet'] = True
                i += 1
            elif sys.argv[i] in ["-wall", "--wall"]:
                g_params['isPrintWarning'] = True
                i += 1
            elif sys.argv[i] in ["-f", "-force", "--f", "--force" ]:
                g_params['isforcewritten'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", sys.argv[i]
                return 1
        else:
            fileList.append(sys.argv[i])
            i += 1

    if len(fileList) <= 0:
        print >> sys.stderr,"no FILE has been set, Exit." 
        return 1

    for f in fileList:
        ConvertFormat(f)

    return 0
#}}}
def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isforcewritten'] = False
    g_params['isQuiet'] = False
    g_params['isPrintWarning'] = False
    g_params['formatindex'] = FORMAT_BINARY
    g_params['isLargeFile'] = False
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

#!/usr/bin/python
# Description:
import os
import sys
import myfunc
import subprocess
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

usage_short="""
Usage: %s -workdir DIR [-idlist IDLISTFILE -ext STR] [-filelist LISTFILE] 
"""%(progname)

usage_ext="""
Description:
    Split files to sub-folders

OPTIONS:
  -workdir         DIR   Output the result to OUTFILE
  -idlist   IDLISTFILE   Set the idlist file
  -filelist   LISTFILE   Set the filelistfile
  -ext             STR   Set extension if idlist is given
  -method          INT   Method to split folder (default: 0)
  -max             INT   maximum files per folder, (default: 2000)
  -q            Quiet mode
  -h, --help    Print this help message and exit

Created 2013-07-02, updated 2013-07-02, Nanjiang Shu
"""
usage_exp="""
Examples:
# split $id.fa $id.png under \"result/\" to sub-folders for id in
# idlistfile.txt, sub-folders are split_$N
    %s -workdir result -idlist idlistfile.txt -ext .fa -ext .png

# split files in filelist.txt to subfolders
    %s -workdir result -filelist filelist.txt

"""%(progname, progname)


def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}
def GetVacantSplitIndex(workdir, split_idx, maxfile_per_folder):#{{{
    while 1:
        subfolder = workdir + os.sep + "split_%d"%(split_idx)
        if not os.path.exists(subfolder):
            cmd = ["mkdir", "-p", subfolder]
            subprocess.call(cmd)
            return (split_idx, maxfile_per_folder)
        else:
            numfile_exist_subfolder = len(os.listdir(subfolder))
            if numfile_exist_subfolder < maxfile_per_folder:
                return (split_idx, maxfile_per_folder-numfile_exist_subfolder)
            else:
                split_idx += 1
#}}}
def MoveFileToFolder(filelist, folder, MAX_PER_MOVE=100):#{{{
    numfile = len(filelist)
    if numfile < 1:
        return 1
    else:
        li = range(0, numfile, MAX_PER_MOVE) + [numfile]
        for i in xrange(len(li)-1):
            subfilelist = filelist[li[i]:li[i+1]]
            cmd = ["/bin/mv", "-f"] + subfilelist + [folder]
            subprocess.call(cmd)
        return 0
#}}}
def SplitToFolder_idlist(idList, workdir, extList, maxfile_per_folder):#{{{
    split_idx = 0
    numExt = len(extList)
    numID = len(idList)
    i = 0
    while i < numID:
        (split_idx, numfile_vacant) = GetVacantSplitIndex(workdir, split_idx,
                maxfile_per_folder)
        numID_for_subfolder = max(1, numfile_vacant/numExt)
        subfolder = workdir + os.sep + "split_%d"%(split_idx)
        filelist_to_move = []
        for j in xrange(i, min(numID, i+numID_for_subfolder)):
            for ext in extList:
                f = workdir + os.sep + idList[j] + ext
                filelist_to_move.append(f)
        MoveFileToFolder(filelist_to_move, subfolder)
        i += numID_for_subfolder
    return 0
#}}}
def SplitToFolder_filelist(fileList, workdir, maxfile_per_folder):#{{{
    split_idx = 0
    numFile = len(fileList)
    i = 0
    while i < numFile:
        (split_idx, numfile_vacant) = GetVacantSplitIndex(workdir, split_idx,
                maxfile_per_folder)
        subfolder = workdir + os.sep + "split_%d"%(split_idx)
        filelist_to_move = fileList[i:i+numfile_vacant]
        MoveFileToFolder(filelist_to_move, subfolder)
        i += numfile_vacant
    return 0
#}}}
def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    workdir = ""
    fileListFile = ""
    idListFile = ""
    extList = []
    maxfile_per_folder = 2000
    method = 0

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            print >> sys.stderr, "Error! Wrong argument:", argv[i]
            return 1
            isNonOptionArg = False
            i += 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-idlist", "--idlist"]:
                (idListFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-filelist", "--filelist"]:
                (fileListFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-workdir", "--workdir"]:
                (workdir, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-max", "--max"]:
                (maxfile_per_folder, i) = myfunc.my_getopt_int(argv, i)
            elif argv[i] in ["-method", "--method"]:
                (method, i) = myfunc.my_getopt_int(argv, i)
            elif argv[i] in ["-ext", "--ext"]:
                (tmpstr, i) = myfunc.my_getopt_str(argv, i)
                extList.append(tmpstr)
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            print >> sys.stderr, "Error! Wrong argument:", argv[i]
            return 1

    if myfunc.checkfile(workdir) != 0:
        return 1

    if idListFile == "" and fileListFile == "":
        print >> sys.stderr, "At least one of idListFile and fileListFile need to be set"
        return 1

    if idListFile != "":
        if os.path.exists(idListFile):
            idList = myfunc.ReadIDList(idListFile)
            if len(idList) <= 0:
                print >> sys.stderr, "No ID in idListFile %s"%(idListFile)
            elif len(extList) <= 0:
                print >> sys.stderr, "No extension set when idList is used."

            else:
                SplitToFolder_idlist(idList, workdir, extList, maxfile_per_folder)
        else:
            print >> sys.stderr, "idListFile %s does not exist"%(idListFile)

    if fileListFile != "":
        if os.path.exists(fileListFile):
            fileList = open(fileListFile, "r").read().split("\n")
            fileList = filter(None, fileList)
            if len(fileList) <= 0:
                print >> sys.stderr, "No file in fileListFile %s"%(fileListFile)
            else:
                SplitToFolder_filelist(fileList, workdir, maxfile_per_folder)
        else:
            print >> sys.stderr, "fileListFile %s does not exist"%(fileListFile)

#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

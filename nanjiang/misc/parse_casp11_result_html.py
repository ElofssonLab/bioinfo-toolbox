#!/usr/bin/python
# Description:
# Parse the html result file downloaded from the CASP website

import os
import sys
import myfunc
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))

usage_short="""
Usage: %s FILE [-o outpath]
"""%(progname)

usage_ext="""
Description:

OPTIONS:
  -l LISTFILE   Set the listfile
  -q            Quiet mode
  -h, --help    Print this help message and exit

Created 2014-11-25, updated 2014-11-25, Nanjiang Shu
"""
usage_exp="""
Examples:
"""

def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}
def ParseHTMLCASP(infile, outpath):
    try:
        fpin = open(infile, "r")
        buff = fpin.read()
        fpin.close()

        p1 = buff.find("<textarea")
        if p1 == -1:
            print >> sys.stderr, "no <textarea> found in %s. Ignore"%(infile)
            return 1
        p2 = buff[p1:].find(">")
        p3 = buff[p1+p2:].find("</textarea>")

        if not (p1>=0 and p2>=0 and p3>=0):
            print >> sys.stderr, "not data found in %s"%(infile)
            return 1

        data = buff[p1+p2+1:p1+p2+p3]

        lines = data.split("\n")
        pred_type = ""
        target_name = ""
        model_nr = ""
        for line in lines:
            if line.find("PFRMAT") != -1:
                pred_type = line.strip().split()[1]
            elif line.find("TARGET") != -1:
                target_name =  line.strip().split()[1]
            elif line.find("MODEL") != -1:
                model_nr =  line.strip().split()[1]

            if pred_type != "" and target_name != "" and model_nr != "":
                break

#         print pred_type
#         print target_name
#         print model_nr
        if not (pred_type != "" and target_name != "" and model_nr != ""):
            print >> sys.stderr,"Bad data, pred_type=%s, target_name=%s, model_nr=%s. Ingore %s"%(pred_type, target_name, model_nr, infile) 
            return 1

        outfile = "%s/%s_%s_%s.txt"%(outpath, target_name, pred_type, model_nr)
        fpout = open(outfile, "w")
        fpout.write("%s"%(data))
        fpout.close()
        print "%s output"%(outfile)
    except IOError:
        print >> sys.stderr, "Failed to read file %s"%(infile)
        return 1

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = "./"
    fileListFile = ""
    fileList = []

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            fileList.append(argv[i])
            isNonOptionArg = False
            i += 1
        elif argv[i] == "--":
            isNonOptionArg = True
            i += 1
        elif argv[i][0] == "-":
            if argv[i] in ["-h", "--help"]:
                PrintHelp()
                return 1
            elif argv[i] in ["-outpath", "--outpath"]:
                (outpath, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-l", "--l"] :
                (fileListFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            fileList.append(argv[i])
            i += 1

    if not os.path.exists(outpath):
        os.system("mkdir -p %s"%(outpath))

    if fileListFile != "":
        fileList += myfunc.ReadIDList(fileListFile)

    for i in xrange(len(fileList)):
        ParseHTMLCASP(fileList[i], outpath)
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

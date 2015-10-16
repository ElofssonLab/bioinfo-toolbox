#!/usr/bin/env python
# Get Tree List Order (sequential order)

import os
import sys

itol_path = os.environ['HOME'] + os.sep + ".local/lib/python2.7/itol"
if not os.path.exists(itol_path):
    itol_path = "/afs/pdc.kth.se/home/n/nanjiang" + os.sep +  "Public/lib/python2.7/itol" 

sys.path.append(itol_path)

import Itol
import ItolExport
import myfunc

progname =  os.path.basename(sys.argv[0])


usage="""
Usage: %s TreeFile [TreeFile ...]
    Get the order of items in the phylo tree by given the Newick tree file
Options:
  -outpath  DIR   Set outpath
  -q              Quiet mode
  -h|--help       Print this help message and exit

Created 2012-03-19, updated 2013-09-19, Nanjiang Shu  
"""%(progname)

def PrintHelp():
    print usage
def GetTreeListOrder(treefile, outpath):
#Create the Itol class
    itl = Itol.Itol()
#Set the tree file
    tree = treefile

    dirname = os.path.dirname(treefile)
    if dirname == "":
        dirname = "."
    if outpath == "":
        outpath = dirname
    elif not os.path.exists(outpath):
        os.system("mkdir -p %s"%(outpath))

    rootname = os.path.basename(os.path.splitext(treefile)[0])

#===================================
    itl.add_variable('treeFile',tree)
    itl.add_variable('treeName',rootname)
    itl.add_variable('treeFormat','newick')
#===================================
# Check parameters
# itl.print_variables()
#Submit the tree
    print ''
    good_upload = itl.upload()
    if good_upload == False:
        print 'There was an error:'+itl.comm.upload_output
        return 1

#Export to pdf
    print 'Exporting to pdf'
    tree_id = itl.comm.tree_id

    itol_exporter = itl.get_itol_export()
    itol_exporter.set_export_param_value('format',"eps")
    itol_exporter.set_export_param_value('displayMode',"normal")
    itol_exporter.set_export_param_value('lineWidth',"1")
    epsfile = outpath + os.sep + rootname + '.itolnormal.eps'
    pdffile = outpath + os.sep + rootname + '.itolnormal.pdf'
    itol_exporter.export(epsfile)
    orderfile = outpath + os.sep + rootname + '.listorder.txt'

    os.system("epstopdf %s" % epsfile )
    os.system("pdftotext -nopgbrk %s %s" %(pdffile, orderfile))
    os.system("rm -f %s"% pdffile)
    os.system("rm -f %s"% epsfile)
    print 'list order output to ',orderfile    

def main(g_params):#{{{
    argv = sys.argv
    datapath = ''
    ext = ""
    outpath = ''
    fileList = []
    fileListFile = ''

    i = 1
    numArgv=len(argv)
    if numArgv < 2:
        PrintHelp()
        return()
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            isNonOptionArg=False
            fileList.append(argv[i])
            i = i + 1
        elif argv[i] == "--":
            isNonOptionArg=True
            i = i + 1
        elif argv[i][0] == "-":
            if argv[i] ==  "-h" or  argv[i] == "--help":
                PrintHelp()
                return(0)
            elif argv[i] == "-outpath" or argv[i] == "--outpath":
                (outpath, i) = myfunc.my_getopt_str(argv, i)
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % argv[i])
                return(1)
        else:
            fileList.append(argv[i])
            i+=1

    if fileListFile != "":
        fileList += myfunc.ReadIDList(fileListFile)

    if len(fileList) <= 0:
        print >> sys.stderr, "No input set. Exit"
        return(1)
    else:
        cnt = 0
        for treefile in fileList:
            print "================== ", cnt , treefile, " ===================="
            GetTreeListOrder(treefile, outpath)
            cnt += 1
def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))


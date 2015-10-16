#!/usr/bin/env python

import os
import sys
from ete2 import Tree
import myfunc

itol_path = os.environ['HOME'] + os.sep + ".local/lib/python2.7/itol"
if not os.path.exists(itol_path):
    itol_path = "/afs/pdc.kth.se/home/n/nanjiang" + os.sep +  "Public/lib/python2.7/itol" 
#itol_path = "/data3/downloads/itol/"
sys.path.append(itol_path)

import Itol
import ItolExport

progname = os.path.basename(sys.argv[0])

usage="""
Usage: %s  -l pfamidlist [ID [ID ...]]
    
Description:
    Generate ITOL tree with Sofia's data

OPTIONS:

  -datapath DIR   Set datapath, (default: ./)
  -outpath  DIR   Set outpath
  -q              Quiet mode
  -h|--help       Print this help message and exit

Created 2013-12-12, updated 2013-12-12, Nanjiang Shu 
"""%(progname)

def PrintHelp():
    print usage
def GetFontSize(numLeave):
    fontsize = int(numLeave/50.0)
    fontsize = max(16, fontsize)
    fontsize = min(200, fontsize)

    return fontsize

def ReadIDList(infile):
    try:
        fpin = open(infile,"r")
        li = fpin.read().split()
        fpin.close()
        return li
    except IOError:
        print "Failed to read listfile ", infile


def Itol_Tree(pfamid, datapath, outpath):
#Create the Itol class
    itl = Itol.Itol()
#Set the tree file
    treefile = "%s%s%s%s"%(datapath ,  os.sep , pfamid ,
            ".TMpro.clustalo10.fasttree")
    taxofile = "%s%s%s%s"%(datapath , os.sep , pfamid ,
            "-taxonomy-file.txt")
    TMdeffile = "%s%s%s%s"%(datapath , os.sep , pfamid ,
            "tms.txt")

    if not os.path.exists(treefile):
        print >> sys.stderr, "tree file %s does not exist. Ignore" %(treefile)
        return 1
    if myfunc.checkfile(taxofile, "taxofile") != 0:
        return 1

    
    numLeave = len(open(taxofile, "r").readlines())

#     t = Tree(treefile)
#     leaves = t.get_leaves()


    fontsize = GetFontSize(numLeave)

    colordeffile = taxofile
    branchlabelfile = ""
    dataset1 = ""
    dataset2 = ""
    dataset3 = ""
    dataset4 = ""

    dataset1 = TMdeffile
#===================================
    itl.add_variable('treeFile',treefile)
    itl.add_variable('treeName', pfamid)
    itl.add_variable('treeFormat','newick')
    if os.path.exists(colordeffile):
        itl.add_variable('colorDefinitionFile', colordeffile)
    if os.path.exists(branchlabelfile):
        itl.add_variable('branchLabelsFile', branchlabelfile)

    if os.path.exists(dataset1):
        itl.add_variable('dataset1File',dataset1)
        itl.add_variable('dataset1Label','numTM_and_repeat')
        itl.add_variable('dataset1Separator','comma')
        itl.add_variable('dataset1Type','domains')
        itl.add_variable('dataset1PreventOverlap','1')
        itl.add_variable('dataset1Color','#FF0000')
#        itl.add_variable('dataset1BarSizeMax','300')

#===================================
#     itl.add_variable('dataset1File',dataset1)
#     itl.add_variable('dataset1Label','numTM')
#     itl.add_variable('dataset1Separator','comma')
#     itl.add_variable('dataset1Type','simplebar')
#     itl.add_variable('dataset1Color','#FF0000')

#===================================
    if os.path.exists(dataset2):
        itl.add_variable('dataset2File', dataset2)
        itl.add_variable('dataset2Label', 'taxonomy')
        itl.add_variable('dataset2Separator','comma')
        itl.add_variable('dataset2Type','colorstrip')
        itl.add_variable('dataset2StripWidth','300')
        itl.add_variable('dataset2PreventOverlap','1')
        itl.add_variable('dataset2ColoringType','both')
        itl.add_variable('dataset2CirclesSpacing','100')


#===================================
    if os.path.exists(dataset3):
        itl.add_variable('dataset3File', dataset3)
        itl.add_variable('dataset3Label', 'pfam')
        itl.add_variable('dataset3Separator','tab')
        itl.add_variable('dataset3Type','colorstrip')
#        itl.add_variable('dataset3Type','ColorDefinitionFile')
#        itl.add_variable('dataset3StripWidth','300')
#        itl.add_variable('dataset3PreventOverlap','1')
#        itl.add_variable('dataset3ColoringType','both')
#        itl.add_variable('dataset3CirclesSpacing','100')

#===================================
    if os.path.exists(dataset4):
        itl.add_variable('dataset4File', dataset4)
        itl.add_variable('dataset4Label', 'cluster')
        itl.add_variable('dataset4Separator','comma')
        itl.add_variable('dataset4Type','colorstrip')
        itl.add_variable('dataset4StripWidth','200')
        itl.add_variable('dataset4PreventOverlap','1')
        itl.add_variable('dataset4ColoringType','both')
#itl.add_variable('dataset1BarSizeMax','1')

#===================================
# Check parameters
# itl.print_variables()


#Submit the tree
    print ''
    print 'Uploading the tree.  This may take some time depending on how large the tree is and how much load there is on the itol server'
    good_upload = itl.upload()
    if good_upload == False:
        print 'There was an error:'+itl.comm.upload_output
        sys.exit(1)

#Read the tree ID
    print 'Tree ID: '+str(itl.comm.tree_id)

#Read the iTOL API return statement
    print 'iTOL output: '+str(itl.comm.upload_output)

#Website to be redirected to iTOL tree
    print 'Tree Web Page URL: '+itl.get_webpage()

# Warnings associated with the upload
    print 'Warnings: '+str(itl.comm.warnings)

#Export to pdf
    print 'Exporting to pdf'
    itol_exporter = itl.get_itol_export()
#itol_exporter = itolexport.ItolExport()
#itol_exporter.set_export_param_value('tree','18793532031912684633930')
    itol_exporter.set_export_param_value('format', 'eps')
    itol_exporter.set_export_param_value('displayMode',"circular")
    itol_exporter.set_export_param_value('showBS',"0")
    itol_exporter.set_export_param_value('fontSize',fontsize)
    itol_exporter.set_export_param_value('alignLabels',"1")
    itol_exporter.set_export_param_value('datasetList','dataset1,dataset2,dataset3,dataset4')
    epsfile = outpath + os.sep + pfamid + '-itol.eps'
    pdffile = outpath + os.sep + pfamid + '-itol.pdf'
    jpgfile = outpath + os.sep + pfamid + '-itol.jpg'
    svgfile = outpath + os.sep + pfamid + '-itol.svg'
    itol_exporter.export(epsfile)
    itol_exporter.export(svgfile)
    os.system("epstopdf %s" % epsfile )
    os.system("convert  %s %s" % (epsfile, jpgfile) )
    print 'exported tree to ', pdffile

if __name__ == '__main__' :

#    datapath = '/data3/wk/MPTopo/pfamAna/pairwise/withinPfam/test'
    datapath = "."
    outpath = './'
    idList = []
    idListFile = ''
    
    i = 1;
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        sys.exit()
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            isNonOptionArg=False;
            idList.append(sys.argv[i])
            i = i + 1;
        elif sys.argv[i] == "--":
            isNonOptionArg=True;
            i = i + 1;
        elif sys.argv[i][0] == "-":
            if sys.argv[i] ==  "-h" or  sys.argv[i] == "--help":
                PrintHelp();
                sys.exit(0);
            elif sys.argv[i] in [ "-datapath", "--datapath"]:
                datapath = sys.argv[i+1]
                i = i + 2;
            elif sys.argv[i] in [ "-l", "--l"]:
                idListFile = sys.argv[i+1]
                i = i + 2;
            elif sys.argv[i] == "-outpath" or sys.argv[i] == "--outpath":
                outpath=sys.argv[i+1];
                i = i + 2;
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i]);
                sys.exit(1);
        else:
            idList.append(sys.argv[i]);
            i+=1;

    if idListFile != "":
        idList += ReadIDList(idListFile)
    numID = len(idList)
    if numID < 1:
        print >> sys.stderr, "No ID set. exit"
        sys.exit(1)

    if datapath == "" or not os.path.exists(datapath):
        print >> sys.stderr, "datapath (%s) does not exist "%(datapath)
        sys.exit(1)

    os.system("mkdir -p %s"%outpath)
    cnt = 0
    for pfamid in idList:
        print "================== ", cnt , pfamid, " ===================="
        Itol_Tree(pfamid, datapath, outpath)
        cnt += 1


#!/usr/bin/env python

import os
import sys
from ete2 import Tree

itol_path = os.environ['HOME'] + os.sep + ".local/lib/python2.7/itol"
if not os.path.exists(itol_path):
    itol_path = "/afs/pdc.kth.se/home/n/nanjiang" + os.sep +  "Public/lib/python2.7/itol" 
#itol_path = "/data3/downloads/itol/"
sys.path.append(itol_path)

import Itol
import ItolExport

usage="""
Usage:  itol_pfamtree.py -l pfamidlist [ID [ID ...]]
    Display pfam trees
Options:

  -datapath DIR   Set datapath, (default: ./)
  -outpath  DIR   Set outpath
  -q              Quiet mode
  -h|--help       Print this help message and exit

Created 2012-03-13, updated 2012-10-22, Nanjiang Shu 
"""

def PrintHelp():
    print usage
def GetFontSize(numLeave):
    fontsize = int(numLeave/100.0*80)
    fontsize = max(30, fontsize)
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
    tree = datapath + os.sep + pfamid + '.kalignp.fasttree'
    if not os.path.exists(tree):
        print >> sys.stderr, "tree file %s does not exist. Ignore" %(tree)
        return 1
    t = Tree(tree)
    leaves = t.get_leaves()
    numLeave = len(leaves)

    fontsize = GetFontSize(numLeave)

    colordeffile = datapath + os.sep + pfamid + '.pfam.colordef.txt'
    branchlabelfile = datapath + os.sep + pfamid + '.branchlabel.txt'
    dataset1 = ""
    dataset2 = ""
    dataset3 = ""
    dataset4 = ""

    dataset1 = datapath + os.sep + pfamid + '.numTM_and_io.txt'
    dataset2 = datapath + os.sep + pfamid + '.taxo.colordef.txt'
    #dataset3 = datapath + os.sep + pfamid + '.ntermstate.colordef.txt'
    dataset4 = datapath + os.sep + pfamid + 'cluster.colordef.txt'

#===================================
    itl.add_variable('treeFile',tree)
    itl.add_variable('treeName', pfamid)
    itl.add_variable('treeFormat','newick')
    if os.path.exists(colordeffile):
        itl.add_variable('colorDefinitionFile', colordeffile)
    if os.path.exists(branchlabelfile):
        itl.add_variable('branchLabelsFile', branchlabelfile)

    if os.path.exists(dataset1):
        itl.add_variable('dataset1File',dataset1)
        itl.add_variable('dataset1Label','numTM_and_io')
        itl.add_variable('dataset1Separator','comma')
        itl.add_variable('dataset1Type','multibar')
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
    thumbfile = outpath + os.sep + "thumb." + pfamid + '-itol.jpg'
    itol_exporter.export(epsfile)
    os.system("epstopdf %s" % epsfile )
    os.system("convert  %s %s" % (epsfile, jpgfile) )
    os.system("convert -thumbnail 200 %s %s" % (jpgfile, thumbfile))
    print 'exported tree to ',pdffile

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
    if len(idList) > 0:
        os.system("mkdir -p %s"%outpath)
        cnt = 0
        for pfamid in idList:
            print "================== ", cnt , pfamid, " ===================="
            Itol_Tree(pfamid, datapath, outpath)
            cnt += 1


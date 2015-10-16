#!/usr/bin/env python

import os
import sys
from ete2 import Tree
import math
import myfunc

itol_path = os.environ['HOME'] + os.sep + ".local/lib/python2.7/itol"
#itol_path = "/data3/downloads/itol/"
if not os.path.exists(itol_path):
    itol_path = "/afs/pdc.kth.se/home/n/nanjiang" + os.sep +  "Public/lib/python2.7/itol" 
sys.path.append(itol_path)

import Itol
import ItolExport

usage="""
USAGE:  itol_pfamtree.py [-datapath DIR] -l pfamidlist [ID [ID ...]]
    Visualize phylogenetic tree of Pfam family, highlighting important features
    of membrane proteins
OPTIONS:
  -m, -method INT Method for visualization, (default: 0)
  -datapath DIR   Set datapath, (default: ./)
  -outpath  DIR   Set outpath, (default: ./)
  -q              Quiet mode
  -h, --help      Print this help message and exit

Created 2012-03-13, updated 2014-10-06, Nanjiang Shu 
"""

def PrintHelp():#{{{
    print usage
    #}}}
def GetFontSize(numLeave):#{{{
    fontsize = 500/math.sqrt(numLeave)
    fontsize = max(30, fontsize)
    fontsize = min(200, fontsize)

    return fontsize
#}}}

def Itol_Tree_m0(pfamid, datapath, outpath):#{{{
#Create the Itol class
    itl = Itol.Itol()
#Set the tree file
    tree = datapath + os.sep + pfamid + '.kalignp.fasttree'
    (dataset1, dataset2, dataset3, dataset4) = ("", "", "", "")
    if not os.path.exists(tree):
        print >> sys.stderr, "tree file %s does not exist. Ignore" %(tree)
        return 1
    t = Tree(tree)
    leaves = t.get_leaves()
    numLeave = len(leaves)

    fontsize = GetFontSize(numLeave)

    dataset1 = datapath + os.sep + pfamid + '.numTM_and_io.txt'
    dataset2 = datapath + os.sep + pfamid + '.cmpclass.colordef.txt'
#    dataset3 = datapath + os.sep + pfamid + '.ntermstate.colordef.txt'
    dataset4 = datapath + os.sep + pfamid + '.cluster.colordef.txt'

    colordeffile = datapath + os.sep + pfamid + '.pfam.colordef.txt'
    branchlabelfile = datapath + os.sep + pfamid + '.branchlabel.txt'

#===================================
    itl.add_variable('treeFile',tree)
    itl.add_variable('treeName','PF00854')
    itl.add_variable('treeFormat','newick')
    if os.path.exists(colordeffile):
        itl.add_variable('colorDefinitionFile', colordeffile)
    if os.path.exists(branchlabelfile):
        itl.add_variable('branchLabelsFile', branchlabelfile)

#===================================
    if os.path.exists(dataset1):
        itl.add_variable('dataset1File',dataset1)
        itl.add_variable('dataset1Label','numTM_and_io')
        itl.add_variable('dataset1Separator','comma')
        itl.add_variable('dataset1Type','multibar')
        itl.add_variable('dataset1PreventOverlap','1')
        itl.add_variable('dataset1Color','#FF0000')

#===================================
    if os.path.exists(dataset2):
        itl.add_variable('dataset2File', dataset2)
        itl.add_variable('dataset2Label', 'cmpclass')
        itl.add_variable('dataset2Separator','comma')
        itl.add_variable('dataset2Type','colorstrip')
        itl.add_variable('dataset2StripWidth','200')
        itl.add_variable('dataset2PreventOverlap','1')
        itl.add_variable('dataset2ColoringType','both')

#===================================
    if os.path.exists(dataset3):
        itl.add_variable('dataset3File', dataset3)
        itl.add_variable('dataset3Label', 'NtermState')
        itl.add_variable('dataset3Separator','comma')
        itl.add_variable('dataset3Type','colorstrip')
        itl.add_variable('dataset3StripWidth','200')
        itl.add_variable('dataset3PreventOverlap','1')
        itl.add_variable('dataset3ColoringType','both')

#===================================
    if os.path.exists(dataset4):
        itl.add_variable('dataset4File', dataset4)
        itl.add_variable('dataset4Label', 'cluster')
        itl.add_variable('dataset4Separator','comma')
        itl.add_variable('dataset4Type','colorstrip')
        itl.add_variable('dataset4StripWidth','200')
        itl.add_variable('dataset4PreventOverlap','1')
        itl.add_variable('dataset4ColoringType','both')
        itl.add_variable('dataset4BranchColoringType','both')
        
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
    os.system("convert %s %s" % (epsfile, jpgfile) )
    os.system("convert -thumbnail 200 %s %s" % (jpgfile, thumbfile))
    print 'exported tree to ',pdffile
#}}}
def Itol_Tree_m1(pfamid, datapath, outpath):#{{{
# TM helices are treated as domains
#Create the Itol class
    itl = Itol.Itol()
#Set the tree file
    tree = datapath + os.sep + pfamid + '.kalignp.fasttree'
    (dataset1, dataset2, dataset3, dataset4) = ("", "", "", "")
    if not os.path.exists(tree):
        print >> sys.stderr, "tree file %s does not exist. Ignore" %(tree)
        return 1
    t = Tree(tree)
    leaves = t.get_leaves()
    numLeave = len(leaves)

    fontsize = GetFontSize(numLeave)

    dataset1 = datapath + os.sep + pfamid + '.numTM_and_io.txt'
    dataset2 = datapath + os.sep + pfamid + '.cmpclass.colordef.txt'
#    dataset3 = datapath + os.sep + pfamid + '.ntermstate.colordef.txt'
    dataset4 = datapath + os.sep + pfamid + '.cluster.colordef.txt'

    colordeffile = datapath + os.sep + pfamid + '.pfam.colordef.txt'
    branchlabelfile = datapath + os.sep + pfamid + '.branchlabel.txt'

#===================================
    itl.add_variable('treeFile',tree)
    itl.add_variable('treeName','PF00854')
    itl.add_variable('treeFormat','newick')
    if os.path.exists(colordeffile):
        itl.add_variable('colorDefinitionFile', colordeffile)
    if os.path.exists(branchlabelfile):
        itl.add_variable('branchLabelsFile', branchlabelfile)

#===================================
    if os.path.exists(dataset1):
        itl.add_variable('dataset1File',dataset1)
        itl.add_variable('dataset1Label','numTM_and_io')
        itl.add_variable('dataset1Separator','comma')
        itl.add_variable('dataset1Type','multibar')
        itl.add_variable('dataset1PreventOverlap','1')
        itl.add_variable('dataset1Color','#FF0000')

#===================================
    if os.path.exists(dataset2):
        itl.add_variable('dataset2File', dataset2)
        itl.add_variable('dataset2Label', 'cmpclass')
        itl.add_variable('dataset2Separator','comma')
        itl.add_variable('dataset2Type','colorstrip')
        itl.add_variable('dataset2StripWidth','200')
        itl.add_variable('dataset2PreventOverlap','1')
        itl.add_variable('dataset2ColoringType','both')

#===================================
    if os.path.exists(dataset3):
        itl.add_variable('dataset3File', dataset3)
        itl.add_variable('dataset3Label', 'NtermState')
        itl.add_variable('dataset3Separator','comma')
        itl.add_variable('dataset3Type','colorstrip')
        itl.add_variable('dataset3StripWidth','200')
        itl.add_variable('dataset3PreventOverlap','1')
        itl.add_variable('dataset3ColoringType','both')

#===================================
    if os.path.exists(dataset4):
        itl.add_variable('dataset4File', dataset4)
        itl.add_variable('dataset4Label', 'cluster')
        itl.add_variable('dataset4Separator','comma')
        itl.add_variable('dataset4Type','colorstrip')
        itl.add_variable('dataset4StripWidth','200')
        itl.add_variable('dataset4PreventOverlap','1')
        itl.add_variable('dataset4ColoringType','both')
        itl.add_variable('dataset4BranchColoringType','both')
        
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
    os.system("convert %s %s" % (epsfile, jpgfile) )
    os.system("convert -thumbnail 200 %s %s" % (jpgfile, thumbfile))
    print 'exported tree to ',pdffile
#}}}

def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    datapath = "."
    outpath = './'
    idList = []
    idListFile = ''

    i = 1;
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
            if sys.argv[i] in [ "-h", "--help"]:
                PrintHelp();
                return 1
            elif sys.argv[i] in [ "-datapath", "--datapath"]:
                datapath = sys.argv[i+1]
                i += 2;
            elif argv[i] in [ "-method", "--method"]:
                g_params['method'], i = myfunc.my_getopt_int(sys.argv, i)
            elif sys.argv[i] in [ "-l", "--l"]:
                idListFile = sys.argv[i+1]
                i = i + 2;
            elif sys.argv[i] in ["-outpath", "--outpath"]:
                outpath = sys.argv[i+1];
                i = i + 2;
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i]);
                return 1
        else:
            idList.append(sys.argv[i]);
            i+=1;

    if idListFile != "":
        idList += myfunc.ReadIDList(idListFile)
    if len(idList) > 0:
        os.system("mkdir -p %s"%outpath)
        cnt = 0
        for pfamid in idList:
            print "================== ", cnt , pfamid, " ===================="
            if g_params['method'] == 0:
                Itol_Tree_m0(pfamid, datapath, outpath)
            elif g_params['method'] == 1:
                Itol_Tree_m1(pfamid, datapath, outpath)
            cnt += 1
#}}}
def InitGlobalParameter():#{{{
    g_params = {}
    g_params['method'] = 0
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))



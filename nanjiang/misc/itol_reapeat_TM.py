#!/usr/bin/env python
# Filename: itol_reapeat_TM.py
import os
import sys
import myfunc
#from ete2 import Tree

itol_path = os.environ['HOME'] + os.sep + ".local/lib/python2.7/itol"
if not os.path.exists(itol_path):
    itol_path = "/afs/pdc.kth.se/home/n/nanjiang" + os.sep +  "Public/lib/python2.7/itol" 
sys.path.append(itol_path)

import Itol
import ItolExport

usage="""
Usage:  itol_pfamtree.py treefile
    Display phylogenetic trees
    All datafiles needed should be under the same folder as treefile
OPTIONS:
  -q              Quiet mode
  -m INT          set methods, (default: 0)
  -h|--help       Print this help message and exit

Created 2013-09-19, updated 2013-09-19, Nanjiang Shu 
"""

def PrintHelp():
    print usage
def GetFontSize(numLeave):#{{{
    fontsize = int(numLeave/100.0*80)
    fontsize = max(30, fontsize)
    fontsize = min(200, fontsize)

    return fontsize
#}}}

def Itol_Tree_m0(treefile, outpath, method): #{{{
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

    treefile_without_ext = os.path.splitext(treefile)[0]
    rootname = os.path.basename(os.path.splitext(treefile)[0])


    if not os.path.exists(tree):
        print >> sys.stderr, "tree file %s does not exist. Ignore" %(tree)
        return 1

    #t = Tree(tree)
    #leaves = t.get_leaves()
    #numLeave = len(leaves)


    colordeffile = treefile_without_ext + '.pfam.colordef.txt'
    branchlabelfile = treefile_without_ext + '.branchlabel.txt'
    dataset1 = ""
    dataset2 = ""
    dataset3 = ""
    dataset4 = ""

    dataset1 = treefile_without_ext + '.numTM_and_io.txt'
    dataset2 = treefile_without_ext + '.domain.txt'
    dataset3 = treefile_without_ext + '.taxo.colordef.txt'
    dataset4 = treefile_without_ext + '.cluster.colordef.txt'

    numLeave = len(open(dataset1, "r").readlines()) - 2
    fontsize = GetFontSize(numLeave)

#===================================
    itl.add_variable('treeFile',tree)
    itl.add_variable('treeName', rootname)
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
        itl.add_variable('dataset1MultiBarAlign','0')
        itl.add_variable('dataset1Color','#FF0000')
        itl.add_variable('dataset1BarSizeMax','1000')
#===================================
    if os.path.exists(dataset2):
        itl.add_variable('dataset2File', dataset2)
        itl.add_variable('dataset2Label', 'Domain')
        itl.add_variable('dataset2Separator','comma')
        itl.add_variable('dataset2Type','domains')
        itl.add_variable('dataset2ProtSizeMax','1000')
        itl.add_variable('dataset2PreventOverlap','1')
        itl.add_variable('dataset2CirclesSpacing','100')
#===================================
    if os.path.exists(dataset3):
        itl.add_variable('dataset3File', dataset3)
        itl.add_variable('dataset3Label', 'taxonomy')
        itl.add_variable('dataset3Separator','comma')
        itl.add_variable('dataset3Type','colorstrip')
        itl.add_variable('dataset3StripWidth','300')
        itl.add_variable('dataset3PreventOverlap','1')
        itl.add_variable('dataset3ColoringType','both')
        itl.add_variable('dataset3CirclesSpacing','100')
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
        return 1

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
    epsfile = treefile_without_ext + '.m%d.itol.eps'%(method)
    pdffile = treefile_without_ext + '.m%d.itol.pdf'%(method)
    jpgfile = treefile_without_ext + '.m%d.itol.jpg'%(method)
    svgfile = treefile_without_ext + '.m%d.itol.svg'%(method)
    thumbfile = treefile_without_ext + '.m%d.itol.thumb.jpg'%(method)
    itol_exporter.export(epsfile)
    os.system("epstopdf %s" % epsfile )
    os.system("convert  %s %s" % (epsfile, jpgfile) )
    os.system("convert -thumbnail 200 %s %s" % (jpgfile, thumbfile))
    print 'exported tree to ',pdffile

    itol_exporter.set_export_param_value('format', 'svg')
    itol_exporter.export(svgfile)
#}}}
def Itol_Tree_m1(treefile, outpath, method): #{{{
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

    treefile_without_ext = os.path.splitext(treefile)[0]
    rootname = os.path.basename(os.path.splitext(treefile)[0])


    if not os.path.exists(tree):
        print >> sys.stderr, "tree file %s does not exist. Ignore" %(tree)
        return 1

    #t = Tree(tree)
    #leaves = t.get_leaves()
    #numLeave = len(leaves)


    #colordeffile = treefile_without_ext + '.pfam.colordef.txt'
    colordeffile = ""
    branchlabelfile = treefile_without_ext + '.branchlabel.txt'
    dataset1 = ""
    dataset2 = ""
    dataset3 = ""
    dataset4 = ""

    dataset1 = treefile_without_ext + '.numTM.txt'
    dataset2 = treefile_without_ext + '.domain.txt'
    dataset3 = treefile_without_ext + '.taxo.colordef.txt'
    dataset4 = treefile_without_ext + '.cluster.colordef.txt'

    numLeave = len(open(dataset1, "r").readlines()) - 2
    fontsize = GetFontSize(numLeave)

#===================================
    itl.add_variable('treeFile',tree)
    itl.add_variable('treeName', rootname)
    itl.add_variable('treeFormat','newick')
    if os.path.exists(colordeffile):
        itl.add_variable('colorDefinitionFile', colordeffile)
    if os.path.exists(branchlabelfile):
        itl.add_variable('branchLabelsFile', branchlabelfile)



#===================================
    if os.path.exists(dataset1):
        itl.add_variable('dataset1File',dataset1)
        itl.add_variable('dataset1Label','numTM')
        itl.add_variable('dataset1Separator','comma')
        itl.add_variable('dataset1Type','multibar')
        itl.add_variable('dataset1PreventOverlap','1')
        itl.add_variable('dataset1MultiBarAlign','1')
        itl.add_variable('dataset1Color','#FF0000')
        itl.add_variable('dataset1BarSizeMax','1500')
#===================================
    if os.path.exists(dataset2):
        itl.add_variable('dataset2File', dataset2)
        itl.add_variable('dataset2Label', 'Domain')
        itl.add_variable('dataset2Separator','comma')
        itl.add_variable('dataset2Type','domains')
        itl.add_variable('dataset2PreventOverlap','1')
        itl.add_variable('dataset2CirclesSpacing','100')
        itl.add_variable('dataset2ProtSizeMax','1000')
#===================================
    if os.path.exists(dataset3):
        itl.add_variable('dataset3File', dataset3)
        itl.add_variable('dataset3Label', 'taxonomy')
        itl.add_variable('dataset3Separator','comma')
        itl.add_variable('dataset3Type','colorstrip')
        itl.add_variable('dataset3StripWidth','300')
        itl.add_variable('dataset3PreventOverlap','1')
        itl.add_variable('dataset3ColoringType','both')
        itl.add_variable('dataset3CirclesSpacing','100')
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
        return 1

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

    epsfile = treefile_without_ext + '.m%d.itol.eps'%(method)
    pdffile = treefile_without_ext + '.m%d.itol.pdf'%(method)
    jpgfile = treefile_without_ext + '.m%d.itol.jpg'%(method)
    svgfile = treefile_without_ext + '.m%d.itol.svg'%(method)
    thumbfile = treefile_without_ext + '.m%d.itol.thumb.jpg'%(method)

    itol_exporter.export(epsfile)
    os.system("epstopdf %s" % epsfile )
    os.system("convert  %s %s" % (epsfile, jpgfile) )
    os.system("convert -thumbnail 200 %s %s" % (jpgfile, thumbfile))
    print 'exported tree to ',pdffile

    itol_exporter.set_export_param_value('format', 'svg')
    itol_exporter.export(svgfile)
#}}}

def main(g_params): #{{{
    argv = sys.argv
    outpath = ''


    numArgv=len(argv)

    if numArgv < 2:
        PrintHelp()
        return 1

    fileList = []
    method = 0

    i = 1
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
                sys.exit(0)
            elif argv[i] in [ "-outpath", "--outpath"]:
                (outpath, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in [ "-m", "--m"]:
                (method, i) = myfunc.my_getopt_int(argv, i)
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % argv[i])
                return 1
        else:
            fileList.append(argv[i])
            i+=1

    if len(fileList) < 1:
        print >> sys.stderr, "No input set"
        return 1
    cnt = 0
    for treefile in fileList:
        print "================== ", cnt , treefile, " ===================="
        if method == 0:
            Itol_Tree_m0(treefile, outpath, method)
        elif method == 1:
            Itol_Tree_m1(treefile, outpath, method)
        cnt += 1
#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))


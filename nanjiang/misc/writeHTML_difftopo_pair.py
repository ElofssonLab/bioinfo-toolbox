#!/usr/bin/python
# Description:
import os
import sys
sys.path.append("%s/wk/MPTopo/src"%(os.environ['DATADIR3']))
import myfunc
progname =  os.path.basename(sys.argv[0])
wspace = ''.join([" "]*len(progname))
import subprocess

usage_short="""
Usage: %s -l pairlistfile -datapath DIR -outpath DIR
"""%(progname)

usage_ext="""
Description:

OPTIONS:
  -l LISTFILE   Set the listfile
  -q            Quiet mode
  -h, --help    Print this help message and exit

Created 2009-06-08, updated 2013-03-14, Nanjiang Shu
"""
usage_exp="""
Examples:
"""

def PrintHelp(fpout=sys.stdout):#{{{
    print >> fpout, usage_short
    print >> fpout, usage_ext
    print >> fpout, usage_exp#}}}

def WriteHTMLHeader(title, fpout):#{{{
    print >> fpout, "<HTML>"
    print >> fpout, "<head>"
    print >> fpout, "<title>%s</title>"%(title)
    print >> fpout, "<script src=\"sorttable.js\"></script>"
    print >> fpout, "<style type=\"text/css\" media=\"all\"> @import \"layout.css\";"
    print >> fpout, "<!--"
    #print >> fpout, "td {font-family: \"SansSerif\", \"SansSerif\", mono; font-size: 10px; text-align: center; }"

    print >> fpout, "table.sortable thead {"
    print >> fpout, "    background-color:#eee;"
    print >> fpout, "    color:#666666;"
    print >> fpout, "    font-weight: bold;"
    print >> fpout, "    cursor: default;"
    print >> fpout, "}"

    print >> fpout, "table.maininfo"
    print >> fpout, "{"
    print >> fpout, "   border: 1px solid black;"
    print >> fpout, "   margin: 10px 10px 10px 10px;"
    print >> fpout, "   padding: 2px 10px 2px 10px;"
    print >> fpout, "}"
    print >> fpout, "table.maininfo td"
    print >> fpout, "{"
    print >> fpout, "   border: 0px solid black;"
    print >> fpout, "   padding: 2px 20px 2px 20px;"
    print >> fpout, "}"
    print >> fpout, "table.content"
    print >> fpout, "{"
    print >> fpout, "   border: 1px solid black;"
    print >> fpout, "   margin: 5px 5px 5px 5px;"
    print >> fpout, "   padding: 5px 5px 5px 5px;"
    print >> fpout, "}"
    print >> fpout, "table.content td"
    print >> fpout, "{"
    print >> fpout, "   font-family:Arial; font-size: 14px;"
    print >> fpout, "   vertical-align: top;"
    print >> fpout, "   text-align: center;"
    print >> fpout, "   border: 3px solid black;"
    print >> fpout, "   padding: 15px 5px 15px 5px;"
    print >> fpout, "}"
    print >> fpout, "table.record"
    print >> fpout, "{"
    print >> fpout, "   border: 0px solid black;"
    print >> fpout, "   padding: 2px 10px 2px 10px;"
    print >> fpout, "}"
    print >> fpout, "table.record td, th"
    print >> fpout, "{"
    print >> fpout, "   font-family:Arial; font-size: 14px;"
    print >> fpout, "   vertical-align: top;"
    print >> fpout, "   text-align: center;"
    print >> fpout, "   border: 0px solid black;"
    print >> fpout, "   padding: 5px 2px 5px 2px;"
    print >> fpout, "}"
    print >> fpout, "table.entry"
    print >> fpout, "{"
    print >> fpout, "   border: 1px solid black;"
    print >> fpout, "   border-collapse: collapse;"
    print >> fpout, "   padding: 2px 10px 2px 10px;"
    print >> fpout, "}"
    print >> fpout, "table.entry td, th"
    print >> fpout, "{"
    print >> fpout, "   font-family:Arial; font-size: 14px;"
    print >> fpout, "   vertical-align: middle;"
    print >> fpout, "   text-align: center;"
    print >> fpout, "   max-width: 200px;"
    print >> fpout, "   border: 1px solid black;"
    print >> fpout, "   padding: 2px 2px 2px 2px;"
    print >> fpout, "}"
    print >> fpout, "table.hits"
    print >> fpout, "{"
    print >> fpout, "   border: 1px solid black;"
    print >> fpout, "   border-collapse: collapse;"
    print >> fpout, "   padding: 2px 10px 2px 10px;"
    print >> fpout, "}"
    print >> fpout, "table.hits td, th"
    print >> fpout, "{"
    print >> fpout, "   font-family:Arial; font-size: 14px;"
    print >> fpout, "   vertical-align: top;"
    print >> fpout, "   text-align: center;"
    print >> fpout, "   max-width: 150px;"
    print >> fpout, "   border: 1px solid black;"
    print >> fpout, "   padding: 2px 2px 2px 2px;"
    print >> fpout, "}"
    print >> fpout, "table.subtable"
    print >> fpout, "{"
    print >> fpout, "   border: 0px solid black;"
    print >> fpout, "   border-collapse: collapse;"
    print >> fpout, "   padding: 2px 10px 2px 10px;"
    print >> fpout, "}"
    print >> fpout, "table.subtable td, th"
    print >> fpout, "{"
    print >> fpout, "   font-family:Arial; font-size: 11px;"
    print >> fpout, "   vertical-align: top;"
    print >> fpout, "   text-align: center;"
    print >> fpout, "   border: 0px solid black;"
    print >> fpout, "   padding: 2px 2px 2px 2px;"
    print >> fpout, "}"
    print >> fpout, "-->"
    print >> fpout, "</style>"
    print >> fpout, "</head>"
    print >> fpout, "<BODY>"
#}}}
def WriteHTMLTail(fpout):#{{{
    print >> fpout, "</BODY>"
    print >> fpout, "</HTML>"
#}}}

def WriteHTMLTable(tablename, tabletitle, pairList, datapath, htmlname,
        outpath, fpout):
    numInputID = len (pairList)
    print >> fpout, "<a name=\"%s\"></a><h4>%s</h4>"%(tablename,tabletitle)
    print >> fpout, "<table class=\"sortable\" border=1>"

    targetpath = outpath + os.sep + "data"
    if not os.path.exists(targetpath):
        os.system("mkdir -p %s"%(targetpath))

    cntOutputID = 0

    headerItemList=[]
    headerItemList.append("No.")
    headerItemList.append("OPM ID")
    headerItemList.append("Uniprot ID")
    headerItemList.append("TopoAlign")
    headerItemList.append("Figure")

    print >> fpout, "<tr>"
    for item in headerItemList:
        print >> fpout, "<th>"
        print >> fpout, item
        print >> fpout, "</th>"
    print >> fpout, "</tr>"

    for i in xrange(numInputID):
        if len(pairList[i]) < 2:
            continue
        cntOutputID += 1
        print >> fpout, "<tr>"
#---------------------------
        print >> fpout, '<td>'
        print >> fpout, '%d'%(cntOutputID)
        print >> fpout, '</td>'
#---------------------------
        print >> fpout, '<td>'
        print >> fpout, '%s'%(pairList[i][0])
        print >> fpout, '</td>'
#---------------------------
        print >> fpout, '<td>'
        print >> fpout, '%s'%(pairList[i][1])
        print >> fpout, '</td>'
#---------------------------
        sourcefile = (datapath + os.sep + "%s_%s"%(pairList[i][0],
            pairList[i][1]) + ".topoaln.krbias.html")
        cmd = ["/bin/cp", "-f", sourcefile, targetpath]
        try:
            subprocess.check_call(cmd)
        except subprocess.CalledProcessError, e:
            print e
        print >> fpout, '<td>'
        print >> fpout, ("<a href=\"%s\" target=\"_blank\">"
                % ("data" + os.sep +
                    os.path.basename(sourcefile)))
        print >> fpout, os.path.basename(sourcefile)
        print >> fpout, "</a>"
        print >> fpout, '</td>'

#---------------------------
        sourcefile = (datapath + os.sep + "%s_%s"%(pairList[i][0],
            pairList[i][1]) + ".topoaln.krbias-crop.png")
        cmd = ["/bin/cp", "-f", sourcefile, targetpath]
        subprocess.check_call(cmd)
        print >> fpout, '<td>'
        print >> fpout, ("<a href=\"%s\" target=\"_blank\">"
                % ( "data" + os.sep +
                    os.path.basename(sourcefile)))
        print >> fpout, ("<img src=\"%s\" target=\"_blank\" height=\"200\" width=\"400\">"
                % ( "data" + os.sep +
                    os.path.basename(sourcefile)))
        print >> fpout, "</a>"
        print >> fpout, '</td>'

        print >> fpout, "</tr>"
    print >> fpout, "</table>"

def WriteHTML(pairList, datapath, outpath):#{{{
    """
    Write HTML
    """
    htmlname = g_params['htmlname']
# The name of the main html file is called index.html if no name is give
    htmlfilename = outpath + os.sep + htmlname + ".html"
    figuredir = outpath + os.sep + htmlname
    local_javalibpath = os.environ['HOME'] + os.sep + 'libjavascript'
    jsfileList = []
    jsfileList.append(local_javalibpath + os.sep + 'sorttable.js')
# copy jsfile to outpath 
    for f in jsfileList:
        if os.path.exists(f):
            os.system("cp -f %s %s"%(f, outpath))

    try:
        os.system("mkdir -p %s"%(figuredir)) 
        fpout = open(htmlfilename,"w")
        title="List of topologies predicted differently from OPM"
        WriteHTMLHeader(title, fpout)

        print >> fpout, "<dir id=\"Content\">"
        tablename = 'table1'
        tabletitle = ""
        WriteHTMLTable(tablename, tabletitle, pairList, datapath, htmlname,
                outpath, fpout)
        print >> fpout, "</dir>"
        WriteHTMLTail(fpout)
        fpout.close()
        print  "Result has been output to %s"%(htmlfilename)
    except IOError:
        raise

#}}}
def main(g_params):#{{{
    argv = sys.argv
    numArgv = len(argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outpath = ""
    datapath = ""
    pairListFile = ""
    pairList = []

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            pairList.append(argv[i].split())
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
            elif argv[i] in ["-datapath", "--datapath"]:
                (datapath, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-l", "--l"] :
                (pairListFile, i) = myfunc.my_getopt_str(argv, i)
            elif argv[i] in ["-q", "--q"]:
                g_params['isQuiet'] = True
                i += 1
            else:
                print >> sys.stderr, "Error! Wrong argument:", argv[i]
                return 1
        else:
            pairList.append(argv[i].split())
            i += 1


    if pairListFile != "":
        pairList += [x.split() for x in myfunc.ReadIDList(pairListFile,
            delim="\n")]

    numpair = len(pairListFile)
    if numpair < 1:
        print >> sys.stderr, "no pair set. exit"
        return 1
    if datapath == "":
        print >> sys.stderr, "datapath not set"
        return 1
    elif not os.path.exists(datapath):
        print >> sys.stderr, "datapath %s does not exist"%(datapath)
        return 1

    if outpath == "":
        print >> sys.stderr, "outpath not set"
        return 1
    elif not os.path.exists(outpath):
        cmd = ["mkdir", "-p", outpath]
        subprocess.check_call(cmd)


    WriteHTML(pairList, datapath, outpath)


#}}}

def InitGlobalParameter():#{{{
    g_params = {}
    g_params['isQuiet'] = True
    g_params['htmlname'] = "index"
    return g_params
#}}}
if __name__ == '__main__' :
    g_params = InitGlobalParameter()
    sys.exit(main(g_params))

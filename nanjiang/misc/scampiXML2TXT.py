#!/usr/bin/env python
# convert the modhmm xml file to text file, output probabilities as well

import sys,re,os
import os.path
import tempfile
import myfunc

#Note: 2010-08-26 
# still use xslt to read xml files. it is more reliable

usage="""
Usage:  compare_topos.py [Options] [-i] xmlfile
Options:
  -i         <file> : input file
  -aapath    <dir>  : set the path of aaSeqs for obtaining annotations and sequence information
  -o         <file> : output the result to file
  -h|--help         : print this help message and exit
  -m|--method 1|2   : method to convert xml to txt, method 1 uses xslt, method 2 uses python,default=1
Created 2010-08-19, updated 2010-08-25, Nanjiang
"""
def PrintHelp():
    print usage


def ReadSingleFasta_modhmm(inFile):
    seqID=""
    aaSeq=""
    annotation=""
    try:
        fpin = open(inFile, "r")
        line = fpin.readline()
        while line:
            line = line.rstrip('\n').strip()
            if line:
                if line[0] == ">":
                    seqID = myfunc.GetSeqIDFromAnnotation(line)
                    annotation = line.lstrip(">").strip()
                elif line[0] != "/": #neglecting membrane topology labelling
                    aaSeq=aaSeq+line
            line = fpin.readline()
        fpin.close()
    except:
        print >> sys.stderr, "Except for the input file ", inFile, "in the function ReadSingleFasta"
    return (seqID, annotation, aaSeq)

def SCAMPI_formattxtfile(inFile, aaPath, fpout):
#structure of predTopoSeq
#predTopoSeq={'seqID':"", 'annotation':"",'aaSeq':"",'topology':"",'length':0,'topoalphabet':"",'prob':[[]]}

    fpin = open(inFile, "r")
    line = fpin.readline()
    while line:
        line = line.rstrip('\n').strip()
        if line:
            strs=line.split(':')
            if len(strs) > 1 and strs[0] == "Seq ID":
                seqID = strs[1]
                annotation=""
                length=0
                normalizedLogLikelihood=""
                logodds=""
                reversi=""
                isTMProtein=""
                aaSeq=""
                topology =""
                topoalphabet =""
                prob=[[]]

                seqID=strs[1].strip()
                while 1:#{{{
                    line = fpin.readline().rstrip().strip()
                    if not line: 
                        break
                    #print line
                    strs=line.split(':')
                    if strs[0] == "Seq length":
                        length=int(strs[1].strip())
                    elif strs[0] == "NormalizedLogLikelihood":
                        normalizedLogLikelihood=strs[1]
                    elif strs[0]  == "Logodds":
                        logodds = strs[1]
                    elif strs[0]  == "Reversi":
                        reversi = strs[1]
                    elif strs[0]  == "Labeling":
                        topology=strs[1].strip()
                    elif strs[0]  == "Is TM protein":
                        isTMProtein="yes"
                    elif strs[0]  == "No TM protein":
                        isTMProtein="no"
                    elif strs[0][0] == "O":
                        topoalphabet=strs[0].replace(" ", "")
                        #init prob and allocate space
                        prob= [[0.0 for col in range(len(topoalphabet))] for row in range(length)]
                        for i in range(0,length):
                            line = fpin.readline()
                            strs=line.split()
                            for j in range(0,len(strs)):
                                prob[i][j]=float(strs[j])
                    #}}}
                if aaPath != "":
                    aaSeqFile=("%s/%s")%(aaPath, seqID)
                    if os.path.exists(aaSeqFile):
                        (seqID, annotation, aaSeq)=ReadSingleFasta_modhmm(aaSeqFile)
                else:
                    seqID=seqID
                    annotation=""
                    aaSeq="-"*length

                #print out#{{{
                print >> fpout, "SeqID:%s" % seqID
                print >> fpout, "SeqLength:%s" % length
                if normalizedLogLikelihood:
                    print >> fpout, "NormalizedLogLikelihood:%s" % normalizedLogLikelihood
                if logodds:
                    print >> fpout, "Logodds:%s" % logodds
                if reversi:
                    print >> fpout, "Reversi:%s" % reversi
                if isTMProtein:
                    print >> fpout, "IsTMProtein:%s" % isTMProtein
                if annotation:
                    print >> fpout, ">%s" %annotation
                if aaSeq:
                    print >> fpout,"Sequence:%s" % aaSeq
                if topology:
                    print >> fpout,"Topology:%s" % topology
                if topoalphabet:
                    fpout.write("%2s %2s" %("AA", "TP"))
                    for i in range(0,len(topoalphabet)):
                        fpout.write(" %5c" % topoalphabet[i])
                    fpout.write("\n")

                    for i in range (0,length):
                        fpout.write("%-2c %2c" %(aaSeq[i], topology[i]))
                        for j in range (0, len(topoalphabet)):
                            fpout.write(" %5.3f" %(prob[i][j]))
                        fpout.write("\n")
                fpout.write("\n")
                #}}}
        line = fpin.readline()
    fpin.close()

if __name__ == '__main__' :
    # Check argv
    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        sys.exit(1)

    aaPath = ""
    outFile=""
    inFile=""
    method=1

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            isNonOptionArg=False
            i = i + 1
        elif sys.argv[i] == "--":
            isNonOptionArg=True
            i = i + 1
        elif sys.argv[i][0] == "-":
            if sys.argv[i] ==  "-h" or  sys.argv[i] == "--help":
                PrintHelp()
                sys.exit(0)
            elif sys.argv[i] == "-i" or sys.argv[i] == "--infile":
                inFile=sys.argv[i+1]
                i = i + 2
            elif sys.argv[i] == "-o" or sys.argv[i] == "--outfile":
                outFile=sys.argv[i+1]
                i = i + 2
            elif sys.argv[i] == "-aapath" or sys.argv[i] == "--aapath" :
                aaPath=sys.argv[i+1]
                i = i + 2
            elif sys.argv[i] == "-m" or sys.argv[i] == "--method" :
                method=int(sys.argv[i+1])
                if method < 1 or method >2:
                    print >> sys.stderr,("Error! method should be 1 or 2")
                    sys.exit(1)
                i = i + 2
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i])
                sys.exit(1)
        else:
            inFile=sys.argv[i]
            i +=1

    if inFile == "":
        print >> sys.stderr,"Error! Input file not set."

    fpout = sys.stdout
    if outFile != "":
        fpout = open(outFile,"w")

    programpath=os.path.dirname(os.path.abspath(sys.argv[0]))
    try :
#first convert the xml file to txt using xlst script
        tmpFile=tempfile.mktemp()
        if method==1:
            os.system(("%s/my_modhmmxml2txt < %s > %s" % (programpath, inFile, tmpFile)))
        else:
            os.system(("%s/my_modhmmxml2txt.py  %s -o %s" % (programpath, inFile, tmpFile)))
        SCAMPI_formattxtfile(tmpFile, aaPath, fpout)
        os.remove(tmpFile)
        if fpout != sys.stdout:
            fpout.close()
    except :
        print >>sys.stderr, "except for the input file: %s" % inFile
        raise 

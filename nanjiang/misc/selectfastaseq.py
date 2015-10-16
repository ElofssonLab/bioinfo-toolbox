#!/usr/bin/env python
# Select sequences by give seqIDs
import sys
import os
import myfunc 

usage="""
usage:  selectfastaseq.py  -f fastaseqfile 
                           [ ID [ID ... ] ] [-l FILE] 
                           [ [-mine FLOAT ] [-maxe FLOAT] ]

Description: select fasta sequences by seqID from the given sequence file either 
             by the supplied seqID or by evalue threshold

  -l FILE       Set the ID list file
  -o FILE       Output the result to file, (default: stdout)
  -mine FLOAT   Set the minimal evalue threshold, (default: 0)
  -maxe FLOAT   Set the maximal evalue threshold, (default: 1e9)
  -h,--help     Print this help message and exit

Created 2011-11-16, updated 2012-05-30, Nanjiang Shu  

Examples:
    selectfastaseq.py -f seq.fa -l idlist.txt
    selectfastaseq.py -f seq.fa -maxe 1e-3
"""


def PrintHelp():
    print usage

def main(g_params):#{{{
    numArgv = len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        return 1

    outFile=""
    idList=[]
    idListFile=""
    fastaFile=""

    i = 1
    isNonOptionArg=False
    while i < numArgv:
        if isNonOptionArg == True:
            idList.append(sys.argv[i])
            isNonOptionArg=False
            i = i + 1
        elif sys.argv[i] == "--":
            isNonOptionArg=True
            i = i + 1
        elif sys.argv[i][0] == "-":
            if (sys.argv[i] in [ "-h", "--help"]):
                PrintHelp()
                return 1
            elif (sys.argv[i] in [ "-l", "--l", "-list", "--list"]):
                idListFile=sys.argv[i+1]
                i = i + 2
            elif (sys.argv[i] in [ "-f", "--f", "-fasta", "--fasta"]):
                fastaFile=sys.argv[i+1]
                i = i + 2
            elif (sys.argv[i] in [ "-o", "--o", "-outfile", "--outfile"]):
                outFile=sys.argv[i+1]
                i = i + 2
            elif (sys.argv[i] in [ "-mine", "--mine"]):
                g_params['min_evalue']=float(sys.argv[i+1])
                g_params['isEvalueSet'] = True
                i = i + 2
            elif (sys.argv[i] in [ "-maxe", "--maxe"]):
                g_params['max_evalue']=float(sys.argv[i+1])
                g_params['isEvalueSet'] = True
                i = i + 2
            else:
                print >> sys.stderr,("Error! Wrong argument:%s" % sys.argv[i])
                return 1
        else:
            idList.append(sys.argv[i])
            i+=1

    if fastaFile == "":
        print >> sys.stderr,"Fatal!  fasta file not set. Exit."
        return 1
    elif not os.path.exists(fastaFile):
        print >> sys.stderr,"Fatal! fasta file %s does not exist. Exit."%(fastaFile)
        return 1

    if os.path.exists(idListFile):
        idList += myfunc.ReadIDList(idListFile)
    
    if len(idList) > 0:
        isIDSet = True
    else:
        isIDSet = False
        
    if not g_params['isEvalueSet'] and not isIDSet:
        print >> sys.stderr, "Error! no ID nor evalue threshold is set. Eixt" 
        return 1

    idListSet = set(idList)
    fpout = myfunc.myopen(filename= outFile, default_fp = sys.stdout, mode="w", isRaise=False); 

    fpin = open (fastaFile, "rb")
    if not fpin:
        print >> sys.stderr, "Failed to open fastafile %s"%(fastaFile)
        return -1
    unprocessedBuffer=""
    isEOFreached = False
    BLOCK_SIZE = g_params['BLOCK_SIZE']
    isEvalueSet = g_params['isEvalueSet']
    min_evalue = g_params['min_evalue']
    max_evalue = g_params['max_evalue']
    while 1:
        buff = fpin.read(BLOCK_SIZE)
        if len(buff) < BLOCK_SIZE:
            isEOFreached=True
        buff = unprocessedBuffer + buff
        recordList = []
        unprocessedBuffer = myfunc.ReadFastaFromBuffer(buff,recordList, isEOFreached)
        if len(recordList) > 0: 
            for r in recordList:
                if ((not isIDSet) or (r[0] in idListSet)):
                    if (not isEvalueSet or r[1].lower().find('evalue') < 0):
                        fpout.write(">%s\n"%r[1])
                        fpout.write("%s\n"%r[2])
                    else:
                        evalue=myfunc.GetEvalueFromAnnotation(r[1])
                        if (evalue == None or (evalue >= min_evalue and
                            evalue <= max_evalue)):
                            fpout.write(">%s\n"%r[1])
                            fpout.write("%s\n"%r[2])

        if isEOFreached == True:
            break
    fpin.close()
    myfunc.myclose(fpout)

#}}}
if __name__ == '__main__' :
    # Check argv
    g_params = {}
    g_params['BLOCK_SIZE'] = 100000
    g_params['min_evalue'] = 0.0
    g_params['max_evalue'] = 1e9
    g_params['isEvalueSet'] = False;#whether evalue threshold is set as a criteria
    sys.exit(main(g_params))

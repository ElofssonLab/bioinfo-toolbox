#!/usr/bin/env python

# one arg: arg[1] is a top file that we count helices in
# two args: arg[1] is an alitop-file and we compare the two tops, arg[2] ignored
# OR just use the functions :)

# Modified by Nanjiang 2010-08-11 
# changes:
# 1. For ReadAliTopo function, unnecessary gaps existing in both sequences are
# removed. (this is caused by the multiple sequence alignment and the pairwise
# alignment are directly copied from the MSA)
# 2. The former CompareToposLocally is renamed as CompareToposGaplessly. see
# the description below
# 3. The CompareToposLocally is modified, so that the gaps in the sequences are
# not removed but only the gapping at both ends are removed. 

import sys,re
import libtopologycmp as lcmp

def DetectFileFormat(inFile):#{{{
#file format:
# 0. Fasta format
# 1. one record per line, annotation: seq
# 2. pairwise alignments
    fileFormat = 0
    try:
        fpin = open(inFile, "r")
        line = fpin.readline()
        while line:
            line = line.rstrip('\n')
            if line:
                if line.find("#Number of alignments") == 0 or line.find("#Topology alignment") == 0:
                    fileFormat=2
                elif line[0] == '>': 
                    fileFormat = 0
                else:
                    fileFormat = 1
                break
            line = fpin.readline()
        fpin.close()
    except:
        print >> sys.stderr, "Except for the input file ", inFile, "in the function", sys._getframe().f_code.co_name
    return fileFormat
#}}}
def counttopo(strTop) :#{{{
    intNumMem = len(re.findall('[io]M+[io]', strTop))
    return (intNumMem, strTop[0])
#}}}
def counttopo2(strTop):#{{{
    strTop = re.sub('M+','M', strTop)
    n = strTop.count('M')
    return n
#}}}
def filterTopo(strTop):#{{{
#modify some unreasonble topology status
#     strTop = re.sub('^M+i','oMi', strTop)
#     strTop = re.sub('^M+o','iMo', strTop)
#     strTop = re.sub('oM+$','oMi', strTop)
#     strTop = re.sub('iM+$','iMo', strTop)
    strTop = re.sub('^M+i', lambda x:'o'+'M'*(len(x.group(0))-2)+'i', strTop)
    strTop = re.sub('^M+o', lambda x:'i'+'M'*(len(x.group(0))-2)+'o', strTop)
    strTop = re.sub('oM+$', lambda x:'o'+'M'*(len(x.group(0))-2)+'i', strTop)
    strTop = re.sub('iM+$', lambda x:'i'+'M'*(len(x.group(0))-2)+'o', strTop)
    strTop = re.sub('MoM','ooo',strTop)
    strTop = re.sub('MiM','iii',strTop)
    return strTop
#}}}
def trimTopo_old(strTop) :#{{{
    strTop = re.sub('[ORCcg]','o', strTop)
    strTop = re.sub('[IrnT]','i', strTop)
    strTop = re.sub('[Wh]','M', strTop)
    #  i..... => iiiiii, o.....=>ooooo, '-', '.', 'X' all means gaps
    strTop = re.sub('([io])[-.X]+',lambda x:x.group(1)*len(x.group()) ,strTop)
    #  ......i => iiiiii, .....o=>ooooo, '-', '.', 'X' all means gaps
    strTop = re.sub('[-.X]+([io])',lambda x:x.group(1)*len(x.group()) ,strTop)
    #  MMMMM.....MMMMMMM => MMMMMMMMMMMMMMMMM
    strTop = re.sub('M+([-.])+M+',lambda x:'M'*len(x.group()) ,strTop)
    # ....MMMMiiii => ooooMMMMiiii
    strTop = re.sub('([-.X]+)(M+[i]+)', lambda x:'o'*len(x.group(1))+x.group(2), strTop)
    # ....MMMMoooo => iiiiMMMMoooo
    strTop = re.sub('([-.X]+)(M+[o]+)', lambda x:'i'*len(x.group(1))+x.group(2), strTop)
    # iiiiMMMM.... => iiiiMMMMoooo
    strTop = re.sub('([i]+M+)([-.X]+)', lambda x:x.group(1)+'o'*len(x.group(2)), strTop)
    # ooooMMMM.... => ooooMMMMiiii
    strTop = re.sub('([o]+M+)([-.X]+)', lambda x:x.group(1)+'i'*len(x.group(2)), strTop)
    return strTop
#}}}
def trimTopo(strTop) :#{{{
#bug fixed 2011-08-15
    strTop = re.sub('[ORCcg]','o', strTop)
    strTop = re.sub('[IrnT]','i', strTop)
    strTop = re.sub('[Wh]','M', strTop)
    #  i..... => iiiiii, o.....=>ooooo, '-', '.', 'X' all means gaps
    while True:
        newstr = re.sub('([io])[-.X]+',lambda x:x.group(1)*len(x.group()) ,strTop)
        if newstr == strTop:
            break
        else:
            strTop=newstr

    #  ......i => iiiiii, .....o=>ooooo, '-', '.', 'X' all means gaps
    while True:
        newstr = re.sub('[-.X]+([io])',lambda x:x.group(1)*len(x.group()) ,strTop)
        if newstr == strTop:
            break
        else:
            strTop=newstr

    #  MMMMM.....MMMMMMM => MMMMMMMMMMMMMMMMM
    while True:
        newstr = re.sub('M+([-.])+M+',lambda x:'M'*len(x.group()) ,strTop)
        if newstr == strTop:
            break
        else:
            strTop=newstr
    # ....MMMMiiii => ooooMMMMiiii
    while True:
        newstr = re.sub('([-.X]+)(M+[i]+)', lambda x:'o'*len(x.group(1))+x.group(2), strTop)
        if newstr == strTop:
            break
        else:
            strTop=newstr
    # ....MMMMoooo => iiiiMMMMoooo
    while True:
        newstr = re.sub('([-.X]+)(M+[o]+)', lambda x:'i'*len(x.group(1))+x.group(2), strTop)
        if newstr == strTop:
            break
        else:
            strTop=newstr
    # iiiiMMMM.... => iiiiMMMMoooo
    while True:
        newstr = re.sub('([i]+M+)([-.X]+)', lambda x:x.group(1)+'o'*len(x.group(2)), strTop)
        if newstr == strTop:
            break
        else:
            strTop=newstr
    # ooooMMMM.... => ooooMMMMiiii
    while True:
        newstr = re.sub('([o]+M+)([-.X]+)', lambda x:x.group(1)+'i'*len(x.group(2)), strTop)
        if newstr == strTop:
            break
        else:
            strTop=newstr
    return strTop
#}}}
def top2short(strTop) :#{{{
    strShortTop = strTop
    strShortTop = re.sub('M+','M', strShortTop)
    strShortTop = re.sub('i{1,20}','i', strShortTop)
    strShortTop = re.sub('o{1,20}','o', strShortTop)
    return strShortTop
#}}}
def shiftcheck(strTop1,strTop2):#{{{
    for i in range(1,len(strTop1)):
        if (strTop1[i] == 'i' and strTop2[i] == 'o') or  (strTop1[i] == 'o' and strTop2[i] == 'i'):
            return False
    return True
#}}}
def readTopo(strFile) :#{{{
    flh1 = file(strFile)
    strTitle1 = flh1.readline().lstrip('>').strip()
    strNextLine = flh1.readline()
    strTop1 = ""
    while ( not strNextLine.startswith('>') ) and ( len(strNextLine) > 0 ) :
        strTop1 += strNextLine
        strNextLine = flh1.readline()
    flh1.close()
    return strTop1
#}}}
def GetSeqIDFromAnnotation(line):#{{{
# get the seqID from the annotation line of the fasta  file
    seqID = ""
    line = line.lstrip('>').split()[0]; #get the first word after '>'
    if line and line.find('|') >= 0:# if the annotation line has |, e.g. >sp|P0AE31|ARTM_ECOL6 Arginine ABC transporter permease
        strs = line.split("|")
        if (strs[0] == "sp" or strs[0] == "lcl" or strs[0] == "tr") : seqID = strs[1]
        else                 : seqID = strs[0]
    else:
        seqID=line
    return seqID
#}}}
def GetPIDFromAnnotation(line):#{{{
# get the sequence identity from the annotation line of the fasta  file
    m = re.search("pid=[^, ]*", line)
    if m: pid = m.group(0).split('=')[1].rstrip('%')
    else: pid = ""
    return pid
#}}}
def GetEvalueFromAnnotation(line):#{{{
    m=re.search("evalue=[^, ]*", line)
    if m: evalue = m.group(0).split('=')[1]
    else: evalue = ""
    return evalue
#}}}
def ReadAliMSA(inFile):#{{{
#read in the multiple sequence alignment in the format
#text1:KALSGJSFJA
#text2:LAS--SFJQO
    idList=[]
    topoSeqList=[]
    aaSeq=""
    seqID = ""
    cntSeq = 0
    fpin = open(inFile, "r")
    line = fpin.readline()
    while line:
        line = line.rstrip('\n')
        if line:
            strs=line.split(":")
            if (len(strs)== 2):
                aaSeq=strs[1] 
                seqID = GetSeqIDFromAnnotation(strs[0])
                if seqID == "": seqID = ("seq_%s" % cntSeq)
            idList.append(seqID)
            topoSeqList.append(aaSeq)
            cntSeq = cntSeq +1
        line = fpin.readline()
    fpin.close()
    return (idList, topoSeqList)
#}}}
def ReadFasta2(inFile):#{{{
#try to extract sequence identity and evalues from the annotation as well
    idList=[]
    topoSeqList=[]
    pidList=[]
    evalueList=[]
    aaSeq=""
    seqID = ""
    pid =""
    evalue=""
    cntSeq = 0

    fpin = open(inFile, "r")
    line = fpin.readline()
    while line:
        line = line.rstrip('\n')
        if line:
            if line[0] == ">":
                cntSeq = cntSeq +1
                if (cntSeq > 1): # if not the first occurrence of the leading '>', it means a new sequence has been read in
                    idList.append(seqID)
                    topoSeqList.append(aaSeq)
                    pidList.append(pid)
                    evalueList.append(evalue)
                    seqID = ""
                    pid = ""
                    evalue = ""
                    aaSeq = ""


                seqID = GetSeqIDFromAnnotation(line)
                pid = GetPIDFromAnnotation(line)
                evalue = GetEvalueFromAnnotation(line)
                if seqID == ""     : seqID = ("seq_%s" % cntSeq)
                if pid == ""    : pid="-"
                if evalue == "" : evalue = "-"
            else:
                aaSeq=aaSeq+line
        line = fpin.readline()

    fpin.close()

    if (cntSeq > 0): # when reaching the end of file, a new sequence should also be added
        idList.append(seqID)
        topoSeqList.append(aaSeq)
        pidList.append(pid)
        evalueList.append(evalue)

    return (idList, topoSeqList, pidList, evalueList)
#}}}
def ReadAliMSA2(inFile):#{{{
#try to extract sequence identity and evalues from the annotation as well
#read in the multiple sequence alignment in the format
#text1:KALSGJSFJA
#text2:LAS--SFJQO
    idList=[]
    topoSeqList=[]
    pidList=[]
    evalueList=[]
    aaSeq=""
    seqID = ""
    pid=""
    evalue=""
    cntSeq = 0

    fpin = open(inFile, "r")
    line = fpin.readline()
    while line:
        line = line.rstrip('\n')
        if line:
            strs=line.split(":")
            if (len(strs)== 2):
                aaSeq=strs[1] 
                seqID = GetSeqIDFromAnnotation(strs[0])
                pid = GetPIDFromAnnotation(strs[0])
                evalue = GetEValueFromAnnotation(strs[0])
                if seqID == "": seqID = ("seq_%s" % cntSeq)
                if pid == ""    : pid="-"
                if evalue == "" : evalue = "-"
            idList.append(seqID)
            topoSeqList.append(aaSeq)
            pidList.append(pid)
            evalueList.append(evalue)
            cntSeq = cntSeq +1
        line = fpin.readline()
    fpin.close()
    return (idList, topoSeqList, pidList, evalueList)
        #}}}
def ReadAliTopoNew(inFile, fileFormat):#{{{
    idList=[];     # list of seqIDs
    topoSeqList=[];# list of topology sequences
    pidList=[]; # list of sequence identities
    evalueList=[]; # list of evalues 
    if fileFormat == 0:
        (idList, topoSeqList, pidList, evalueList) = ReadFasta2(inFile)
    else:
        (idList, topoSeqList, pidList, evalueList) = ReadAliMSA2(inFile)
    return (idList, topoSeqList, pidList, evalueList)
#}}}
def readAliTopo(strFile) : # pairwise alignment#{{{
    flh1 = file(strFile)
    strTitle1 = flh1.readline().lstrip('>').strip()
    strNextLine = flh1.readline()
    strTop1 = ""
    while ( not strNextLine.startswith('>') ) and ( len(strNextLine) > 0 ) :
        strNextLine = strNextLine.rstrip('\n').strip();#debug
        strTop1 += strNextLine
        strNextLine = flh1.readline()
    strTitle2 = strNextLine.lstrip('>').strip()
    strNextLine = flh1.readline()
    strTop2 = ""
    while ( not strNextLine.startswith('>') ) and ( len(strNextLine) > 0 ) :
        strNextLine = strNextLine.rstrip('\n').strip();#debug
        strTop2 += strNextLine
        strNextLine = flh1.readline()
    flh1.close()
    return strTop1, strTop2
#}}}
def CompareToposGloballyNew(strTop1, strTop2, strProtein1, strProtein2, fpLog):#{{{
    [strTop1, strTop2] = lcmp.RemoveUnnecessaryGap([strTop1, strTop2])
    strTop1 = trimTopo(strTop1)
    strTop2 = trimTopo(strTop2)
    strTop1 = filterTopo(strTop1)
    strTop2 = filterTopo(strTop2)
    if fpLog != 0:
        print >> fpLog, "Global"
        print >> fpLog, "%-20s:%s"%(strProtein1, strTop1)
        print >> fpLog, "%-20s:%s"%(strProtein2, strTop2)
        print >> fpLog

    if len(strTop1) <= 0 and len(strTop2) <=0:
        return ("DIFF", 0,0 )
    elif len(strTop1)*len(strTop2) == 0 and len(strTop1)+len(strTop2) > 0:
        print >> sys.stderr,"%s %s global length does not match" % (strProtein1, strProtein2)
        sys.exit(1)

    (intNumMem1,Nterm1)=counttopo(strTop1)
    (intNumMem2,Nterm2)=counttopo(strTop2)
    return compareTopos(intNumMem1,intNumMem2,strTop1,strTop2, Nterm1, Nterm2);#}}}
def CompareToposGaplesslyNew(strTop1, strTop2, strProtein1, strProtein2, fpLog):#{{{
# -----iiiMMMooo ___\ iMMMooo
# iiiiii--MMMMoo    / iMMMMoo
    [strTop1, strTop2] = lcmp.RemoveUnnecessaryGap([strTop1, strTop2])
    if fpLog != 0:
        print >> fpLog, "Unnecessary gaps removed"
        print >> fpLog, "%-20s:%s"%(strProtein1, strTop1)
        print >> fpLog, "%-20s:%s"%(strProtein2, strTop2)
        print >> fpLog
    strNewTop1 = ''
    strNewTop2 = ''
    for i in range(len(strTop1)) :
        if not ( strTop1[i] == '-' or strTop2[i] == '-' ) :
            strNewTop1 += strTop1[i]
            strNewTop2 += strTop2[i]
    strNewTop1 = filterTopo(strNewTop1)
    strNewTop2 = filterTopo(strNewTop2)
    if fpLog != 0:
        print >> fpLog, "Gapless"
        print >> fpLog, "%-20s:%s"%(strProtein1, strNewTop1)
        print >> fpLog, "%-20s:%s"%(strProtein2, strNewTop2)
        print >> fpLog

    if len(strNewTop1) <= 0 and len(strNewTop2) <=0:
        return ("DIFF", 0,0 )
    elif len(strNewTop1)*len(strNewTop2) == 0 and len(strNewTop1)+len(strNewTop2) > 0:
        print >> sys.stderr,"%s %s gapless length does not match" % (strProtein1, strProtein2)
        sys.exit(1)
    (intMems1,N1) = counttopo(strNewTop1)
    (intMems2,N2) = counttopo(strNewTop2)
    return compareTopos(intMems1,intMems2,strNewTop1,strNewTop2,N1,N2)
        #}}}
def CompareToposLocallyNew(strTop1, strTop2, strProtein1, strProtein2, fpLog):#{{{
# -----iiiMMMooo-- ___\ iiiMMMooo
# iiiiii--MMMMoooo    / i--MMMMoo
    [strTop1, strTop2] = lcmp.RemoveUnnecessaryGap([strTop1, strTop2])
# 1. treat the beginning
    nbegin=0
    if (strTop1[0] == '-' or strTop2[0] == '-'):
        i=0
        while(i < len(strTop1) and (strTop1[i] == '-' or strTop2[i] == '-')):
                i = i +1
        nbegin = i
# 2. treat the ending                                  
    nend = len(strTop1)
    if (strTop1[len(strTop1)-1] == '-' or strTop2[len(strTop1)-1] == '-'):
        i=len(strTop1)-1
        while(i >= 0 and (strTop1[i] == '-' or strTop2[i] == '-')):
                i = i - 1
        nend = i+1

    tmpStrTop1 = strTop1[nbegin:nend]
    tmpStrTop2 = strTop2[nbegin:nend]
# 3. remove unnecessary gaps
    strNewTop1 = ''
    strNewTop2 = ''
    for i in range(len(tmpStrTop1)) :
        if not ( tmpStrTop1[i] == '-' and tmpStrTop2[i] == '-' ) :
            strNewTop1 += tmpStrTop1[i]
            strNewTop2 += tmpStrTop2[i]
    strNewTop1 = trimTopo(strNewTop1); #after local treatment, gaps may still exist in the alignment, use the function trimTopo to remove these gaps
    strNewTop2 = trimTopo(strNewTop2)

    strNewTop1 = filterTopo(strNewTop1)
    strNewTop2 = filterTopo(strNewTop2)

    if fpLog != 0:
        print >> fpLog, "Locally"
        print >> fpLog, "%-20s:%s"%(strProtein1, strNewTop1)
        print >> fpLog, "%-20s:%s"%(strProtein2, strNewTop2)
        print >> fpLog

    if len(strNewTop1) <= 0 and len(strNewTop2) <=0:
        return ("DIFF", 0,0 )
    elif len(strNewTop1)*len(strNewTop2) == 0 and len(strNewTop1)+len(strNewTop2) > 0:
        print >> sys.stderr,"%s %s local length does not match" % (strProtein1, strProtein2)
        sys.exit(1)

    (intMems1,N1) = counttopo(strNewTop1)
    (intMems2,N2) = counttopo(strNewTop2)
    return compareTopos(intMems1,intMems2,strNewTop1,strNewTop2,N1,N2);#}}}
def compareTopos(intNumMem1,intNumMem2,strTop1,strTop2, Nterm1, Nterm2) :#{{{
    #OK:4, INV:3, SHIFT:2, DIFF:1
    if intNumMem1 == intNumMem2:
        if Nterm1 == Nterm2:
            if shiftcheck(strTop1,strTop2):
                return ("OK",intNumMem1,intNumMem2)
            else:
                return ("SHIFT",intNumMem1,intNumMem2)
        else:
            return  ("INV",intNumMem1,intNumMem2)
    return ("DIFF",intNumMem1,intNumMem2);#}}}
def compareTopologiesCategories(strTop1,strTop2) :#{{{
    #OK:4, INV:3, SHIFT:2, DIFF:1
    strTop1 = filterTopo(strTop1)
    strTop2 = filterTopo(strTop2)
    intNumMem1, Nterm1 = counttopo(strTop1)
    intNumMem2, Nterm2 = counttopo(strTop2)
    if intNumMem1 == intNumMem2:
        if Nterm1 == Nterm2:
            if shiftcheck(strTop1,strTop2):
                return "OK"
            else:
                return "SHIFT"
        else:
            return "INV"
    return "DIFF";#}}}
def compareToposGlobally(strAlitopFile) :#{{{
    print "########\nCompareToposGlobally"; #debug
    (strTop1, strTop2) = readAliTopo(strAlitopFile)
    [strTop1, strTop2] = lcmp.RemoveUnnecessaryGap([strTop1, strTop2])
    print "strTop1:%s" % (strTop1); #debug
    print "strTop2:%s" % (strTop2); #debug
    strTop1 = trimTopo(strTop1)
    strTop2 = trimTopo(strTop2)
    print "After trimming";#debug
    print "strTop1:%s" % (strTop1); #debug
    print "strTop2:%s" % (strTop2); #debug
    strTop1 = filterTopo(strTop1)
    strTop2 = filterTopo(strTop2)
    (intNumMem1,Nterm1)=counttopo(strTop1)
    (intNumMem2,Nterm2)=counttopo(strTop2)
    return compareTopos(intNumMem1,intNumMem2,strTop1,strTop2, Nterm1, Nterm2);#}}}
def compareToposGaplessly(strAlitopFile) :#{{{
# 1st version:
# Cut off unaligned ends: 
# By Nanjiang 2010-08-11: this actually cut off all unaligned regions, not only
# the two endings. The code really deal with the function described below in
# written in another def, and compareToposLocally is renamed as
# compareToposGaplessly
# that is 
# -----iiiMMMooo ___\ iMMMooo
# iiiiii--MMMMoo    / iMMMMoo
    print "##########\nCompareToposGaplessly\n"; #debug
    (strTop1, strTop2) = readAliTopo(strAlitopFile)
    [strTop1, strTop2] = lcmp.RemoveUnnecessaryGap([strTop1, strTop2])

    print "strTop1:%s" % (strTop1); #debug
    print "strTop2:%s" % (strTop2); #debug

    strNewTop1 = ''
    strNewTop2 = ''
    for i in range(len(strTop1)) :
        if not ( strTop1[i] == '-' or strTop2[i] == '-' ) :
            strNewTop1 += strTop1[i]
            strNewTop2 += strTop2[i]

    print "After Gapless treatment\n";#debug
    print "strTop1:%s" % (strNewTop1); #debug
    print "strTop2:%s" % (strNewTop2); #debug

    strNewTop1 = filterTopo(strNewTop1)
    strNewTop2 = filterTopo(strNewTop2)
    (intMems1,N1) = counttopo(strNewTop1)
    (intMems2,N2) = counttopo(strNewTop2)
    return compareTopos(intMems1,intMems2,strNewTop1,strNewTop2,N1,N2);#}}}
def compareToposLocally(strAlitopFile) :#{{{
# Cut off unaligned ends: 
# Modified by Nanjiang from the original code 2010-08-11
# compareToposLocally, that is 
# -----iiiMMMooo-- ___\ iiiMMMooo
# iiiiii--MMMMoooo    / i--MMMMoo
    print "###########\nCompareToposLocally\n"; #debug
    (strTop1, strTop2) = readAliTopo(strAlitopFile)

    [strTop1, strTop2] = lcmp.RemoveUnnecessaryGap([strTop1, strTop2])
    print "strTop1:%s" % (strTop1); #debug
    print "strTop2:%s" % (strTop2); #debug

# 1. treat the beginning
    nbegin=0
    if (strTop1[0] == '-' or strTop2[0] == '-'):
        i=0
        while(i < len(strTop1) and (strTop1[i] == '-' or strTop2[i] == '-')):
                i = i +1
        nbegin = i
# 2. treat the ending                                  
    nend = len(strTop1)
    if (strTop1[len(strTop1)-1] == '-' or strTop2[len(strTop1)-1] == '-'):
        i=len(strTop1)-1
        while(i >= 0 and (strTop1[i] == '-' or strTop2[i] == '-')):
                i = i - 1
        nend = i+1

    tmpStrTop1 = strTop1[nbegin:nend]
    tmpStrTop2 = strTop2[nbegin:nend]
# 3. remove unnecessary gaps
    strNewTop1 = ''
    strNewTop2 = ''
    for i in range(len(tmpStrTop1)) :
        if not ( tmpStrTop1[i] == '-' and tmpStrTop2[i] == '-' ) :
            strNewTop1 += tmpStrTop1[i]
            strNewTop2 += tmpStrTop2[i]

    print "After local treatment\n";#debug
    print "strTop1:%s" % (strNewTop1); #debug
    print "strTop2:%s" % (strNewTop2); #debug

    strNewTop1 = trimTopo(strNewTop1); #after local treatment, gaps may still exist in the alignment, use the function trimTopo to remove these gaps
    strNewTop2 = trimTopo(strNewTop2)

    strNewTop1 = filterTopo(strNewTop1)
    strNewTop2 = filterTopo(strNewTop2)
    (intMems1,N1) = counttopo(strNewTop1)
    (intMems2,N2) = counttopo(strNewTop2)
    return compareTopos(intMems1,intMems2,strNewTop1,strNewTop2,N1,N2)
#}}}

if __name__ == '__main__' :
    if len(sys.argv) == 2:
        strTop1 = trimTopo(readTopo(sys.argv[1]))
        strTop1 = filterTopo(strTop1)
        print "%s\t%s" % counttopo(strTop1)
    else:
        print compareToposGlobally(sys.argv[1])

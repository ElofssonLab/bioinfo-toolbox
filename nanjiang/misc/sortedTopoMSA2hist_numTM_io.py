#!/usr/bin/env python

import sys
import os
import myfunc

infile = sys.argv[1]
red="#FF0000"
green="#00FF00"
blue="#0000FF"
GAP = '-'

# read in taxonomy def
if not os.path.exists(infile):
    print >> sys.stderr, "Error! file infile (%s) does not exist." %infile;
    sys.exit(1)

def GetNtermState(topo):#{{{
    if topo[0] != GAP:
        return topo[0];
    else:
        topo=topo.lstrip(GAP);
        if topo != "":
            return topo[0];
        else:
            return None;
#}}}

def IsAlternate(liStatus, numAltSerie):
    num = len(liStatus)
    for i in range(0, num-numAltSerie+1):
        isAlt = True
        for j in range(i, i+numAltSerie-1):
            if liStatus[j] * liStatus[j+1] != -1:
                isAlt = False
                break
        if isAlt:
            return True
    return False

(idList, seqList) = myfunc.ReadFasta_without_annotation(infile)

# write out taxdef
fpout = sys.stdout
numSeq = len(idList)


countDict = {}
maxNumTM = 0
for i in xrange(numSeq):
    gid = idList[i]
    if gid != 'Consensus':
        NtermState = GetNtermState(seqList[i])
        numTM = myfunc.CountTM(seqList[i])
        maxNumTM = max(maxNumTM, numTM)
        if not numTM in countDict:
            countDict[numTM] = [0,0]
        if NtermState == 'i':
            countDict[numTM][0] += 1
        else:
            countDict[numTM][1] += 1


# check whether it is alternating
liStatus = []
th = 0.75 # only when 75% of the nterm status is either i or o are recognized
minCount = 10
# i : 1
# o : -1
# non-determined: 0

for i in xrange(1, maxNumTM+1):
    try:
        (n1, n2) = (countDict[i][0], countDict[i][1])
    except (KeyError, IndexError):
        (n1, n2) = (0, 0)
    if myfunc.FloatDivision(n1, n1+n2) >= th and n1 >= minCount:
        liStatus.append(1)
    elif myfunc.FloatDivision(n2, n1+n2) >= th and n2 >= minCount:
        liStatus.append(-1)
    else:
        liStatus.append(0)

numAltSerie = 3
isAlternate = IsAlternate(liStatus, numAltSerie)


#write the result
fpout.write("#isAlternate: %d\n"%(isAlternate))
fpout.write("%6s %8s %8s\n"%("#numTM", "i(Nterm)",  "o(Nterm)"))
for i in xrange(1, maxNumTM+1):
    try:
        (n1, n2) = (countDict[i][0], countDict[i][1])
    except (KeyError, IndexError):
        (n1, n2) = (0, 0)
    fpout.write("%6d %8d %8d\n"%(i, n1, n2))

if fpout != sys.stdout:
    fpout.close()

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

(idList, seqList) = myfunc.ReadFasta_without_annotation(infile)

# write out taxdef
fpout = sys.stdout
numSeq = len(idList)
fpout.write("LABELS,mylabel_1,mylabel_2\n")
fpout.write("COLORS,%s,%s\n"%(red, blue))

for i in xrange(numSeq):
    gid = idList[i]
    if gid != 'Consensus':
        n_i = 0
        n_o = 0
        NtermState = GetNtermState(seqList[i])
        numTM = myfunc.CountTM(seqList[i])
        if NtermState == 'o':
            n_i = 0
            n_o = numTM
        else:
            n_i = numTM
            n_o = 0
        fpout.write("%s,%d,%d\n"%(gid, n_i,n_o ))
fpout.write("\n")

if fpout != sys.stdout:
    fpout.close()

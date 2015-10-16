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
numSeq = len(idList)
for i in xrange(numSeq):
    gid = idList[i]
    if gid != 'Consensus':
        color = red
        NtermState = GetNtermState(seqList[i])
        if NtermState == 'o':
            color = blue
        sys.stdout.write("%s,%s\n"%(gid, color ))
sys.stdout.write("\n")

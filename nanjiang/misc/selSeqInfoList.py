#!/usr/bin/env python

import sys
import os
import myfunc
def ReadSeqDefInfo(infile):#{{{
    try:
        seqInfoDict = {}
        fpin = open(infile,"r")
        line = fpin.readline()
        line = fpin.readline()
        while line:
            strs = line.split('|')
            if len(strs) == 4:
                gid = strs[0].strip()
                refseqid = strs[1].strip()
                pfamid = strs[2].strip()
                seqdef = strs[3].strip()
                seqInfoDict[gid] = {}
                seqInfoDict[gid]['pfamid'] = pfamid 
                seqInfoDict[gid]['refseqid'] = refseqid 
                seqInfoDict[gid]['seqdef'] = seqdef 
            line = fpin.readline()
        fpin.close();
        return seqInfoDict
    except IOError:
        print >> sys.stderr, "Error! file seqDefFile (%s) does not exist." %infile;
        return {};
#}}}

def main(g_params):
    seqinfofile = "/data3/wk/MPTopo/pfamAna/pfam2-giid-refseqid-pfamid-description.txt"
    seqidlistfile = "/data3/wk/MPTopo/pfamAna/pairwise/all/pfamfullseq.selTM_uniq.seqidlist"
    seqInfoDict = {}
    seqInfoDict = ReadSeqDefInfo(seqinfofile)
    idList = myfunc.ReadIDList(seqidlistfile)

    print "#gi_id | refseq_id | pfamid | sequence_description"
    for idd in idList:
        print "%s | %s | %s | %s" % ( idd,
                seqInfoDict[idd]['refseqid'],seqInfoDict[idd]['pfamid'],
                seqInfoDict[idd]['seqdef'])
    return 0
        

if __name__ == '__main__' :
    g_params = {};
    sys.exit(main(g_params));


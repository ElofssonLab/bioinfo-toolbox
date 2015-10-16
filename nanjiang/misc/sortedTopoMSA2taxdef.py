#!/usr/bin/env python

import sys
import os

seqDefFile = '/data3/wk/MPTopo/pfamAna/pfam2-selTM-giid-refseqid-pfamid-description.txt'
infile = sys.argv[1]

# read in taxonomy def
if not os.path.exists(infile):
    print >> sys.stderr, "Error! file infile (%s) does not exist." %infile;
    sys.exit(1)


if not os.path.exists(seqDefFile):
    print >> sys.stderr, "Error! file seqDefFile (%s) does not exist." %seqDefFile;
    sys.exit(1)
fpin = open(seqDefFile,"r")
seqInfoDict = {}
line = fpin.readline()
line = fpin.readline()
while line:
    strs = line.split('|')
    if len(strs) == 4:
        gid = strs[0].strip()
        refseqid = strs[1].strip()
        pfamid = strs[2].strip()
        seqdef = strs[3].strip()
#         print seqdef
        taxdef = ""
        if seqdef.find('[') !=  -1:
            taxdef = seqdef.split('[')[1].split(']')[0]
        seqInfoDict[gid] = {}
        seqInfoDict[gid]['pfamid'] = pfamid 
        seqInfoDict[gid]['refseqid'] = refseqid 
        seqInfoDict[gid]['seqdef'] = seqdef 
        seqInfoDict[gid]['taxdef'] = taxdef
    line = fpin.readline()
fpin.close();

# write out taxdef
fpin = open(infile,"r")
line = fpin.readline()
while line:
    if line[0] == '>':
        gid = line.split()[0].lstrip('>')
        if gid != 'Consensus':
            taxdef = ''
            if gid in seqInfoDict:
                taxdef = seqInfoDict[gid]['taxdef'] 
            sys.stdout.write("%s,%s\n"%(gid,taxdef ))
    line = fpin.readline()
fpin.close()

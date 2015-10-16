#!/usr/bin/env python
import os
import sys

file_pairlistwithclanid = "Pfam-A-full.perTM75_nseq20.nr100.filtered.withinclan.pairlistwithclanid"

file_sel_pairlist = "Pfam-A-full.perTM75_nseq20.nr100.filtered.withinclan.max30000.pairlist"

outfile_sel_pairlistwithclanid = "Pfam-A-full.perTM75_nseq20.nr100.filtered.withinclan.max30000.pairlistwithclanid"

fpin = open(file_pairlistwithclanid, "r")
lines = fpin.readlines()
fpin.close()
dict1 = {}
for line in lines:
    if line and line[0] != "#":
        strs = line.split()
        if len(strs) >= 3:
            id1 = strs[0]
            id2 = strs[1]
            famid = strs[2]
            dict1[(id1,id2)] = famid

fpin = open(file_sel_pairlist, "r")
lines = fpin.readlines()
fpin.close()
fpout = open(outfile_sel_pairlistwithclanid, "w")
for line in lines:
    if line and line[0] != "#":
        strs = line.split()
        if len(strs) >= 2:
            id1 = strs[0]
            id2 = strs[1]
            famid = ""
            if (id1, id2) in dict1:
                famid = dict1[(id1,id2)]
            elif (id2,id1) in dict1:
                famid = dict1[(id2,id1)]
            if famid != "":
                fpout.write("%s %s %s\n"%(id1, id2, famid))
            else:
                print >> sys.stderr, "%s - %s not found" % (id1, id2)

fpout.close()
            





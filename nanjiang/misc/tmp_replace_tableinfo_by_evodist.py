#!/usr/bin/env python

import os, sys, myfunc

def ReadEvoDist(infile):
    fpin = open(infile, "r")
    lines = fpin.readlines()
    fpin.close()
    dict1 = {}
    for line in lines:
        if line and line[0] != "#":
            strs = line.split()
            if len(strs) >= 3:
                id1 = strs[0]
                id2 = strs[1]
                evodist = float(strs[2])
                dict1[(id1,id2)] = evodist
    return dict1


in_tableinfo_file= "/data3/wk/MPTopo/pfamAna_refpro/cellular_filter_all/pairwise/withinClan/Pfam-A-full.perTM75_nseq20.nr100.filtered.withinclan.max30000.kalignP.tableinfo"
in_evodist_file = "/data3/wk/MPTopo/pfamAna_refpro/cellular_filter_all/pairwise/withinClan/evodist/allpair.evodist.txt"

outfile_evodist_seqidt1 = "evodist-vs-seqidt1.txt"
outfile_tableinfo = "Pfam-A-full.perTM75_nseq20.nr100.filtered.withinclan.max30000.evodist.tableinfo"
# evodist *= 10

evodistDict = ReadEvoDist(in_evodist_file)

fpout_table = open(outfile_tableinfo, "w")
fpout_dist = open(outfile_evodist_seqidt1, "w")


fpin = open (in_tableinfo_file, "r")
lines = fpin.readlines()
fpin.close()

for line in lines:
    if line and line[0] != "#":
        strs = line.split()
        if len(strs) >= 13:
            id1 = strs[0]
            id2 = strs[1]
            seqidt0 = float(strs[2])
            alnLength = int(strs[4])
            seqLength1 = int(strs[5])
            seqLength2 = int(strs[6])
            numIDT = int ( strs[8])
            numGap = int ( strs[10])
            seqidt1 = float(strs[11])
            seqidt2 = float(strs[12])

            evodist = -100
            if (id1,id2) in evodistDict:
                evodist = evodistDict[(id1,id2)]
            elif (id2,id1) in evodistDict:
                evodist = evodistDict[(id2,id1)]

            if evodist != -100:
                fpout_table.write("%-16s %-15s %6.1f %6.1f %9d %6d %6d %9.1f %6d %6d %6d %6.1f %6.1f\n"% (
                    id1, id2, seqidt0, -1.0,
                    alnLength,
                    seqLength1, seqLength2,
                    -1.0,
                    numIDT, -1, numGap, evodist, seqidt2))
                fpout_dist.write("%f %f\n"%(seqidt1, evodist))


fpout_table.close()
fpout_dist.close()

print "%s output "%outfile_tableinfo
print "%s output "%outfile_evodist_seqidt1


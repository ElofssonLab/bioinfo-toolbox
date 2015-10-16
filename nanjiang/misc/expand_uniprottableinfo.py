#!/usr/bin/env python

import os
import sys
import myfunc
usage="""
USAGE: %s tableinfoFile expanded_tableinfoFile

Created 2014-11-06, updated 2014-11-06, Nanjiang Shu
"""%(sys.argv[0])

try:
    infile=sys.argv[1]
except IndexError:
    print usage
    sys.exit(1)

try:
    outfile=sys.argv[2]
except IndexError:
    print usage
    sys.exit(1)

fpout = myfunc.myopen(outfile, sys.stdout,"w", False)

hdl = myfunc.ReadLineByBlock(infile)
dt = {}
if hdl.failure:
    sys.exit(1)
lines = hdl.readlines()
while lines != None:
    for line in lines:
        if not line:
            continue
        if line[0] == "#":
            print >> fpout, line
        else:
            strs = line.split("\t")
            strs2 = strs[0].split(";")
            content = "\t".join(strs[1:])
            for j in xrange(len(strs2)):
                idd = strs2[j].strip()
                if idd != "":
                    print >> fpout, "%s\t%s"%(idd,content)
    lines = hdl.readlines()
hdl.close()

myfunc.myclose(fpout)

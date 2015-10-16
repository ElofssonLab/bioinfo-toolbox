#!/usr/bin/env python
import os
import sys
import re
import myfunc
DEBUG = 1

usage = """
usage: %s infile outfile
"""%(sys.argv[0])
def InitXY(N):#{{{
    lx = range(N+1)
    ly = [0]*(N+1)
    return (lx, ly)
    #}}}

def GetSegPos(string, keyC):#{{{
    """
    Get segment of a continue keyC state
    e.g. given a string "0001100022000111100"
    and keyC = '1'
    return [(3,5), (13,17)]
    """
    posList = []
    ex = "(%s+)"%(keyC)
    m = re.finditer(ex,string)
    for i in m:
        posList.append((i.start(0), i.end(0)))
    return posList
#}}}
def FilterSegPos(posList, string, neighbour_char):
### return only list of "0110"
    newList = []
    N = len(string)
    for (b,e) in posList:
        if b>0 and string[b-1] != neighbour_char:
            continue
        if e < N-1 and string[e] != neighbour_char:
            continue
        newList.append((b,e))
    return newList

def Stat1(lines, outfile):
    MAX_NUMTM = 14
    dt = {}
    for item in [0.0, 0.5,1.0]:
        dt[item] = InitXY(MAX_NUMTM)
    cnt = 0
    for line in lines:
        if line:
            strs = line.split()
            if strs[0] == "TMMap":
                cnt += 1
                numTM = int(strs[4].rstrip(":"))
                if numTM <= 1 or numTM >= MAX_NUMTM :
                    continue
                mapArray = [int(x) for x in strs[5:]]
                st = 1
                str_maparray_list = ["%d"%x for x in mapArray]
                str_maparray = "".join(str_maparray_list)
                posContList = GetSegPos(str_maparray, "%d"%st)
                neighbour_char = "0"
                posContList = FilterSegPos(posContList, str_maparray, neighbour_char)
                if len(posContList) >= 1:
                    for (b,e) in posContList:
                        if b == 0:
                            pp = 0.0
                        elif e == numTM:
                            pp = 1.0
                        else:
                            pp = 0.5
                            if DEBUG and (e-b) == 1:
                                print "Maparray", mapArray, posContList
                        dt[pp][1][e-b] += 1
    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    fpout.write("#%2s %8s %3s %8s %3s %8s\n"%("Nx", "Ny", "Ix", "Iy", "Cx", "Cy"))
    for i in xrange(1, MAX_NUMTM+1):
        for item in [0.0, 0.5, 1.0]:
            (lx, ly) =  dt[item]
            fpout.write("%3d %8d "%(lx[i], ly[i]))
        fpout.write("\n")
    myfunc.myclose(fpout)

def Stat2(lines, outfile):
    MAX_NUMTM = 14
    dt = {}
    for item in [1, 12, 2]:
        dt[item] = InitXY(MAX_NUMTM)
    cnt = 0
    for line in lines:
        if line:
            strs = line.split()
            if strs[0] == "TMMap":
                cnt += 1
                numTM = int(strs[4].rstrip(":"))
                if numTM <= 1 or numTM >= MAX_NUMTM :
                    continue
                mapArray = [int(x) for x in strs[5:]]
                st = 1
                str_maparray_list = ["%d"%x for x in mapArray]
                str_maparray = "".join(str_maparray_list)
                num_1 = str_maparray.count("1")
                num_2 = str_maparray.count("2")
                num_12 = num_1 + num_2
                dt[1][1][num_1] += 1
                dt[2][1][num_2] += 1
                dt[12][1][num_12] += 1
    fpout = myfunc.myopen(outfile, sys.stdout, "w", False)
    fpout.write("#%2s %8s %3s %8s %3s %8s\n"%("1x", "1y", "2x", "2y", "12x", "12y"))
    for i in xrange(1, MAX_NUMTM+1):
        for item in [1, 2, 12]:
            (lx, ly) =  dt[item]
            fpout.write("%3d %8d "%(lx[i], ly[i]))
        fpout.write("\n")
    myfunc.myclose(fpout)

try: 
    infile = sys.argv[1]
    outfile = sys.argv[2]
    fpin = open (infile, "r") 
    lines = fpin.readlines()
    fpin.close()

    Stat1(lines, outfile)
#    Stat2(lines, outfile)

except (IOError, IndexError):
    print >> sys.stderr, usage



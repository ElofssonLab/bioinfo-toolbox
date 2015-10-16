#!/usr/bin/python
import os
import sys

try:
    infile = sys.argv[1]
except IndexError:
    sys.exit(1)

fpin = open(infile, "r")
lines = fpin.readlines()
fpin.close()

li = []
i = 0
for line in lines:
    strs = line.split()
    for j in range(len(strs)):
        li.append((float(strs[j]), i, j))
    i += 1

li2 = sorted(li, key=lambda x:x[0], reverse=True)

for item in li2:
    if item[0] > 0 and item[1] != item[2]:
        sys.stdout.write("%2d %2d %6.3f\n"%(item[1],item[2], item[0]))

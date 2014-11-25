#!/usr/bin/env python

#A sslight modification of script by Maricn
from localconfig import *
import os, sys, subprocess

affixes = []
for method in ['hh', 'jh']:
        for cutoff in ['E0', 'E4', 'E10', 'E40']:
                affixes.append( method+cutoff )

bin = PconsC3 + '/bin/'
l = sys.argv[1]

# Trying to create inpu to scripts/constructsetG.py

q = [bin+'/constructsetG.py']
for a in affixes:
        q.append( l + "." + a + ".plmdca20")
        q.append( l + "." + a + ".gdca")


q.append(l + ".phycmap")
for a in affixes:
        q.append( l + "." + a + ".stats")
        q.append( l + "." + a + ".trimmed")

q.append(l + ".rsa")
q.append(l + ".ss2")
#q.append(l + ".contactCB")
q.append(l + ".PconsC3.training")

count = 0
subprocess.call(q)
if not os.path.exists(q[-1]):
        subprocess.call(q)
        count += 1


try:
        subprocess.call(bin + '/pconsc26-p0-e.py ' + l + '.PconsC3.training', shell=True)
        subprocess.call(bin + '/pconsc26-p1-e.py 1 ' + l + '.PconsC3.training', shell=True)
#        subprocess.call(bin + '/pconsc26-p1-e.py 2 ' + l + '.PconsC3.training', shell=True)
#        subprocess.call(bin + '/pconsc26-p1-e.py 3 ' + l + '.PconsC3.training', shell=True)
#        subprocess.call(bin + '/pconsc26-p1-e.py 4 ' + l + '.PconsC3.training', shell=True)
#        subprocess.call(bin + '/pconsc26-p1-e.py 5 ' + l + '.PconsC3.training', shell=True)
#        subprocess.call(bin + '/pconsc26-p1-e.py 6 ' + l + '.PconsC3.training', shell=True)
except:
        pass



sys.exit(0)


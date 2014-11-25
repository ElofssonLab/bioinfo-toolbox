#!/usr/bin/env python

import sys, subprocess, os

alignment = sys.argv[1]

if not os.path.exists(alignment):
	sys.exit(0)

stem = alignment[:alignment.rfind('.trimmed')]

if os.path.exists(stem + '.gneff') and os.path.getsize(stem + '.gneff') > 3:
	sys.exit(0)

a = subprocess.check_output('julia /scratch/arne/PconsC3/bin/rungDCA.jl ' + alignment + ' ' + stem + '.gdca', shell=True)

f = open(stem + '.gneff', 'w')
f.write(a)
f.close()

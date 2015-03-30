#!/usr/bin/env python
from localconfig import *
import sys, subprocess, os
import string as s



def check_output(command):
	return subprocess.Popen(command, stdout=subprocess.PIPE).communicate()[0]
	
if not os.path.exists(matlab):
        print "Matlab not found"
	sys.exit(1)



cpu = 1
if '-c' in sys.argv:
    idx = sys.argv.index('-c')
    try:
         cpu = int(sys.argv[idx+1])
    except:
        print 'Number of cores -c must be an integer, %r is not. Default is %s.' % (sys.argv[idx+1], cpu)
        sys.exit(1)
    del sys.argv[idx]
    del sys.argv[idx]

if len(sys.argv) < 2:
	print sys.argv[0], '<target>'
	sys.exit(0)

infile = sys.argv[1]


os.chdir(os.path.abspath(infile)[:os.path.abspath(infile).rfind('/')]) 
infilestem = infile.split('/')[-1]
infilestem = infilestem[:infilestem.rfind('.')]

print "Running plmDCA"
        
print ([matlab, '-nodesktop', '-r', "path(path, '"+ plmdcapath +"'); path(path, '"+ plmdcapath +"/functions'); path(path, " + plmdcapath + "/3rd_party_code/minFunc/'); plmDCA_asymmetric ( '" + infilestem + ".trimmed', '" + infilestem + ".plmdca20', 0.1, 4); exit;"])
t = check_output([matlab, '-nodesktop', '-nosplash' , '-r', "path(path, '"+ plmdcapath +"'); path(path, '"+ plmdcapath +"/functions'); path(path, ' " + plmdcapath + "/3rd_party_code/minFunc/'); disp  'Starting plmDCA_asymmetric' ; plmDCA_asymmetric ( '" + infilestem + ".trimmed', '" + infilestem + ".plmdca20', 0.1, "+ str(cpu) +"); exit;"])
print t

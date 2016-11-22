#!/usr/bin/env python
import sys

sys.path.append("/home/x_arnel/git/bioinfo-toolbox/")
from parsing import parse_fasta
from parsing import parse_contacts

sfile = sys.argv[1]
cfile = sys.argv[2]
target = sys.argv[3]
server = sys.argv[4]
ofilepath = sys.argv[5]
minsep= sys.argv[6]
minscore = sys.argv[7]

seq = parse_fasta.read_fasta(open(sfile)).items()[0][1][0]

contacts = parse_contacts.parse(open(cfile), min_dist=0)

print len(contacts)
print contacts[0]
print seq

ofile = open(ofilepath, 'w')

if server == "Pcons-net":
    ofile.write("PFRMAT RR\nTARGET %s\nAUTHOR 5450-4562-0389\nMETHOD Pcons-net\nREMARK PconsC3\nMETHOD Improved contact predictions on\nMETHOD small protein families.\nMODEL  1\n" % target)
elif server == "PconsC2":
    ofile.write("PFRMAT RR\nTARGET %s\nAUTHOR 4146-6019-9011\nMETHOD PconsC2\nREMARK PconsC2\nMETHOD Improved contact predictions using the\nMETHOD recognition of protein like contact\nMETHOD patterns.\nMODEL  1\n" % target)
elif server == "PconsC31":
    ofile.write("PFRMAT RR\nTARGET %s\nAUTHOR 2066-8709-0532\nMETHOD PconsC31\nREMARK PconsC31\nMETHOD Improved contact predictions on\nMETHOD small protein families.\nMODEL  1\n" % target)
else:
    ofile.write("PFRMAT RR\nTARGET %s\nAUTHOR XXXX-XXXX-XXXX\nMETHOD %s\nREMARK PconsC31\nMODEL  1\n" % (target, server))



tmp_i = 0
for aa in seq:
    ofile.write(aa)
    tmp_i += 1
    if tmp_i % 50 == 0:
        ofile.write("\n")
print tmp_i
if tmp_i % 50 != 0:
    ofile.write("\n")

n_lines = 0
for c in contacts:
    if ( int(c[2]-c[1]) >= int(minsep) and float(c[0]) > float(minscore)):
        ofile.write("%d %d 0 8 %.5f\n" % (c[1], c[2], c[0]))
        n_lines += 1
#    else:
#        print "skipping: ",c[2]-c[1],c[1], c[2], c[0]
    if n_lines > 30000:
        break


ofile.write("END")

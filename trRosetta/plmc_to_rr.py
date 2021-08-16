#!/usr/bin/env python3

#import numpy as np
import os
import sys
from Bio import SeqIO
import pandas as pd
#from collections import namedtuple


cwd = os.getcwd()

if len(sys.argv) != 3:
   print("Usage: {}   seq.fasta contacts.plmc ".format(os.path.basename(__file__)))
   sys.exit()


fasta = os.path.join(cwd,sys.argv[1])
contacts = os.path.join(cwd,sys.argv[2])

with open(fasta, "r") as handle:
   for record in SeqIO.parse(handle, "fasta"):
      seq=record



fasta_seq=str(seq.seq)


#target = '.'.join(input_file.split('/')[-1].split('.')[:-1])


print ("PFRMAT RR")
print ("TARGET {}".format(fasta))
print ("AUTHOR Markslab")
print ("METHOD plmc")
print ("REMARK contacts created by plmc")
print ("MODEL 1")
n=50
for i in range(0, len(fasta_seq), n):
   print (fasta_seq[i:i+n])
factor=2.5
numcontacts=factor*len(fasta_seq)
contacts_handle = open(contacts, 'r')
# Strips the newline character
data=pd.read_csv(contacts,sep=" ",names=("i","foo","j","bar","num","score"))

x=0

for i,data in data.sort_values(by=["score"]).iterrows():
   print (data["i"],data["j"],0.0,8.0,data["score"])
   x+=1
   if (x>numcontacts): break

print ("END")


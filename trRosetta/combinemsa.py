#!/usr/bin/python3

import sys, getopt,re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

#from Bio.SeqFeature import SeqFeature, FeatureLocation


file1=sys.argv[1]
file2=sys.argv[2]
outfile=sys.argv[3]

type="fasta"
handle1 = open(file1, 'r')
handle2 = open(file2, 'r')


for record in SeqIO.parse(handle1, type ) :
    #print (record.name)
    break
seq1=record


for record in SeqIO.parse(handle2, type ) :
    #print (record.name)
    break
seq2=record


seq=seq1+seq2
seq.id=seq1.id+"-"+seq2.id
seq.name="Merged MSA sequence"
seq.description="Merged MSA sequence"
#print (seq.format("fasta"))
extension1=""
extension2=""
for i in range(len(seq1)):
    extension1+="-"
for i in range(len(seq2)):
    extension2+="-"
#align=MultipleSeqAlignment([seq])
align=[]
align+=[seq]
#print (len(seq),len(seq1),len(seq2),len(extension1),len(extension2))
for record in SeqIO.parse(handle1, type ) :
    record.seq+=extension2
    #print (record.format("fasta"))
    align+=[record]
    #print (len(record))
for record in SeqIO.parse(handle2, type ) :
    tempseq=SeqRecord(seq=extension1+record.seq,id=record.id,name=record.name,description=record.description)
    #print (tempseq.format("fasta"))
    align+=[tempseq]
    #print (len(tempseq))
#print (seq.format("fasta"))
#print (align)
outobject=MultipleSeqAlignment(align)
AlignIO.write(outobject, outfile , "fasta")

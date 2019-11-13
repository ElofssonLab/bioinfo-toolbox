#!/usr/bin/env python

import fileinput
import re

inchain=False
for line in fileinput.input():
    #print (line)
    if re.search("^File:",line):
        temp=re.split("\s+",line)
        filename=temp[1]
        inchain=False
        i=0
        score=[]
    elif re.search("^Chain",line):
        inchain=True
    elif (inchain and re.search("^[A-Z ]\s",line)):
        temp=re.split("\s+",line)
        #print (i,temp,temp[4])
        score+=[temp[4]]
        i+=1
    elif len(line)==1:
        inchain=False
        print (filename,score)
        

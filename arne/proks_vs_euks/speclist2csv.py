#!/usr/bin/env python

import re
import sys
name=''
string=''
with open('../data/speclist.txt') as fp:
    for line in fp:
        #print (line)
        p = re.compile(r'^\W\W\W\W\W\s[ABEXV]\s')
        c = re.compile(r'^\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\sC=')
        c = re.compile(r'^\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\s\sN=')
        if p.match(line):
            elsif
        p = re.compile(r'\W+')
        s=p.split(line)
        #print (s[0],name)
        if s[0] != name:
            print (string)
            string=s[0]+" "+s[1]
            name=s[0]
        else:
            string+=" "+s[1]
print string


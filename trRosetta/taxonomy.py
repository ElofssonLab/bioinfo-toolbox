#!/usr/bin/env python3

import argparse


arg_parser = argparse.ArgumentParser(description="summarize taxonomy of a MSA")
arg_parser.add_argument("-i","--input", type=str)
#arg_parser.add_argument("-d","--database", type=str,default="~/git/bioinfo-toolbox/trRosetta/taxonomy_subset.tab")
#arg_parser.add_argument("-o","--out", type=str)
#arg_parser.add_argument("-O","--output", type=str)


args = arg_parser.parse_args()
std = 1

f=open (args.input,"r")


tiny=1.e-20
cutoff=0.75
maxvalue=0
sumvalue=tiny
maxkey=''
maxnum=0
i=0
for l in f.readlines():
    i+=1
    key,v=l.split(":")
    value=int(v)
    sumvalue+=value
    if value>maxvalue:
        maxkey=key
        maxnum=i
        maxvalue=value
fraction=maxvalue/sumvalue
if (maxvalue/sumvalue > cutoff):
    print (maxnum,"\t",maxkey,"\t",fraction)
else:
    print (0,"\t","Mixed ",maxkey,"\t",fraction)


#!/usr/bin/env python3


import argparse
import sys
import os
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description = '''format graph-QA for CASP''')
parser.add_argument('-loc', nargs=1, type= str, default=sys.stdin, help = 'Path to outdir.')
parser.add_argument('-glob', nargs=1, type= str, default=sys.stdin, help = 'Path to outdir.')
parser.add_argument('-out', nargs=1, type= str, default=sys.stdin, help = 'Path to outdir.')
args = parser.parse_args()
localfile = args.loc[0]
globalfile = args.glob[0]
outdir=args.out[0]


localdata = pd.read_csv(localfile)
globaldata = pd.read_csv(globalfile)

targets=globaldata['target_id'].drop_duplicates()
header="PFRMAT QA"
method="METHOD Graph-QA\nMODEL 2\nQMODE 2"
remark="REMARK Graph-QA"
author="AUTHOR 5229-7541-3942"


def convert(x):
    d0=5
    score=d0*np.sqrt(1/x-1)
    if score < 15.0:
        return (score)
    else:
        return(15.0)


for target in targets:
    tempglob=globaldata.loc[globaldata['target_id']==target]
    temploc=localdata.loc[localdata['target_id']==target]
    models=tempglob['decoy_id'].drop_duplicates()
    f = open(outdir+"/"+target+".QA", 'w')
    f.write(header+"\n")
    f.write("TARGET "+str(target)+"\n")
    f.write(author+"\n"+remark+"\n"+method+"\n")
    count=0
    for model in models:
        globscore=tempglob.loc[tempglob['decoy_id']==model]['gdtts'].max()
        temp=temploc.loc[temploc['decoy_id']==model]['lddt']
        locscore=temp.apply(lambda x:convert(x)).to_list()
        f.write(str(target)+str(model)+" ")
        f.write(str(round(globscore,3))+" ")
        
        for i in range(0,len(locscore)):
            f.write(str(round(locscore[i],3))+" ")
            count+=1
            if count>20:
                f.write("\n")
                count=0
        f.write("\nEND\n")
    f.close()


    # Convert local-CA-CA using
#			  rmsd=d0*sqrt(1/Sstr[i][j]-1);
#			  if(rmsd>15)
#			    {
#			      rmsd=15;
#			    }



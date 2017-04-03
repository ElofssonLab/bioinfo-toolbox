#!/usr/bin/env python
from __future__ import print_function

import re
import sys
import os
import csv
  

# we have to guess the method from various type of names. inpout is a .tar.gz name
def parse_method(fname):
    temp=file
    p=re.compile('([jh]h[Ee]\d+)')
    m=p.search(fname)
    if m:
        ali=m.group(1)
    else:
        ali='UNKNOWN'
        #
    p=re.compile('_(\d\d)_(\d\d)')
    m=p.search(fname)
    if m:
        mindist=m.group(1)
        maxdist=m.group(2)
    else:
        mindist=0
        maxdist=0
        #
    p=re.compile('\_(\d\.\d)\_')
    m=p.search(fname)
    if m:
        num=m.group(1)
    else:
        num=9999
        #
    return(ali,num,mindist,maxdist)

def parse_pcons(id,fname):
    with open(fname) as f:
        score = {}
        for l in f:
            l_arr = l.strip().split()
            if not l_arr:
                continue
            target = l_arr[0]
            pdb = l_arr[2]
            model = l_arr[3]
            pcons = l_arr[4]
            score[model]=pcons
    f.close
    return score,pdb

def parse_cns(id,fname):
    with open(fname) as f:
        score_total = {}
        score_noe = {}
        for l in f:
            l_arr = l.strip().split()
            if not l_arr:
                continue
            target = l_arr[0]
            model = l_arr[3]
            total = l_arr[5]
            bond = l_arr[6]
            angle = l_arr[7]
            imp = l_arr[8]
            vdw = l_arr[9]
            noe = l_arr[10]
            score_total[model]=total
            score_noe[model]=noe
    f.close
    return (score_total,score_noe)

def parse_TM(id,fname):
    with open(fname) as f:
        score = {}
        for l in f:
            l_arr = l.strip().split()
            if not l_arr:
                continue
            p=re.compile('(fa\_[0-9]+).*\.pdb')
            m=p.search(l_arr[0])
            if m:
                model=m.group(1)
            else:
                model="UNKNOWN"
            TM=l_arr[3]
            score[model]=TM
    f.close
    return score

def parse_proq(fname):
    with open(fname) as f:
        for l in f:
            if not re.match('^\d',l):
                continue
            l_arr = l.strip().split()
            if not l_arr:
                continue
            proq2=l_arr[0]
            proqC=l_arr[1]
            proqF=l_arr[2]
            proq3=l_arr[3]
    f.close
    return proq2,proqC,proqF,proq3

def length_PDB(fname):
    length=0
    #print fname
    with open(fname) as f:
        for l in f:
            p=re.compile('ATOM.*(CA ).*')
            m=p.search(l)
            if m:
                length=length+1
    f.close
    return length





if __name__=="__main__":
    ProQ2D={}
    ProQ3D={}
    tarfname = os.path.basename(sys.argv[1])
    p=re.compile('.tar.gz$')
    fname=p.sub("",tarfname)
    dname=os.path.dirname(os.path.realpath(sys.argv[1]))
    p=re.compile(".*/PF")
    target=p.sub("PF",dname)
    #print fname,target
    (ali,num,mindist,maxdist)=parse_method(fname)
    #print ali,num,mindist,maxdist
    (pcons,pdb)=parse_pcons(target,dname+"/"+fname+"_pcons.out")
    # print pcons,pdb
    TM=parse_TM(target,dname+"/"+fname+"_TM.out")
    #print tm
    (cns,noe)=parse_cns(target,dname+"/"+fname+"_cns.out")
    #print cns
    #print noe
    for model in pcons:
        proq=parse_proq(dname+"/"+fname+"_proq3/"+target+"."+pdb+"."+model+".pdb.proq3.global")
        ProQ2D[model]=proq[0]
        ProQ3D[model]=proq[3]
        #    print ProQ3D
    length=length_PDB(dname+"/"+fname+"_proq3/"+target+"."+pdb+".fa_1.pdb")
    #    print length
    #    now we need to 
    print ('target , ','ali , ','num , ','mindist , ','maxdist , ','length , ','model , ','TM , ','Pcons , ','cns , ','noe , ','ProQ2D , ','ProQ3D')

                
    
    for model in pcons:
        try:
            print(target," , ",ali," , ",num," , ",mindist," , ",maxdist," , ",length," , ",model," , ",TM[model]," , ",pcons[model]," , ",cns[model]," , ",noe[model]," , ",ProQ2D[model]," , ",ProQ3D[model])
        except:
            print('Error printng output\n', file=sys.stderr)

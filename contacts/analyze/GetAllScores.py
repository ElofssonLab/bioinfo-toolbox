#!/usr/bin/env python
from __future__ import print_function

import re
import sys
import os
import csv
import glob, os

def filename(dir,pattern):
    os.chdir(dir)
    file='NO-HIT"'
    for f in glob.glob(pattern):
        file=f
    return file
        

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
        num=99
        #
    return(ali,num,mindist,maxdist)

def parse_pcons(id,fname):
    IGNORE_LST = ['PFRMAT','TARGET','AUTHOR','REMARK','METHOD','MODEL','QMODE','END']
    score = {}
    with open(fname) as f:
        for l in f:
            l_arr = l.strip().split()
            first = l_arr[0]
            if first in IGNORE_LST or not '.pdb' in first:
                continue
            pdb=l_arr[0]
            model = id+"/stage1/"+pdb
            pcons = l_arr[1]
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
            pdb = l_arr[0]
            model = id+"/stage1/"+pdb
            total = l_arr[1]
            bond = l_arr[2]
            angle = l_arr[3]
            imp = l_arr[4]
            vdw = l_arr[5]
            noe = l_arr[6]
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
            model=l_arr[0]
            TM=l_arr[3]
            score[model]=TM
    f.close
    return score

def parse_gneff(id,dname,pattern):
    file=filename(dname,pattern)
    fname=dname+"/"+file
    with open(fname) as f:
        gneff = {}
        for l in f:
            l_arr = l.strip().split()
            if not l_arr[0] == "M":
                continue
            score=l_arr[8]
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
    (pcons,pdb)=parse_pcons(fname,dname+"/"+fname+".raw")
    # print pcons,pdb
    TM=parse_TM(target,dname+"/"+fname+"_TM.out")
    # print pcons,pdb
    gneff=parse_gneff(target,dname,"*"+ali+"*.gneff")
#    print ("GNEFF: ",gneff)
    #print tm
    (cns,noe)=parse_cns(fname,dname+"/"+fname+"_cns.out")
    #print cns
    #print noe
    for model in noe:
        proq=parse_proq(dname+"/"+fname+"_proq3/"+pdb+".proq3.global")
        ProQ2D[model]=proq[0]
        ProQ3D[model]=proq[3]
        #    print ProQ3D
    length=length_PDB(dname+"/"+fname+"_proq3/"+pdb)
    #    print length
    #    now we need to 
    print ('target , ','model , ','ali , ','num , ','mindist , ','maxdist , ','length , ','TM , ','Pcons , ','cns , ','noe , ','ProQ2D , ','ProQ3D ,','Meff')

                
    
    for model in pcons:
        try:
            print(target," , ",model," , ",ali," , ",num," , ",mindist," , ",maxdist," , ",length," , ",TM[model]," , ",pcons[model]," , ",cns[model]," , ",noe[model]," , ",ProQ2D[model]," , ",ProQ3D[model]," , ",gneff)
        except:
            print('Error printng output\n', file=sys.stderr)

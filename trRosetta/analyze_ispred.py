#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse
from argparse import RawTextHelpFormatter
from Bio.PDB import *
        
##args = docopt.docopt(__doc__)
#out_dir = args['--output_folder']
 

p = argparse.ArgumentParser(description = '- plotting trRosetta maps-',
                            formatter_class=RawTextHelpFormatter)
p.add_argument('-data1','--input1','-i', required= True, help='Input ispred file 1')
p.add_argument('-data2','--input2','-j', required= True, help='Input ispred file 2')
p.add_argument('-pdb','--pdb','-p', required= False, help='Pdb file for analysis od distances')
ns = p.parse_args()

input_file1 = np.load(ns.input1)
input_file2 = np.load(ns.input2)

three2one = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D',
             'CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
             'HIS':'H','ILE':'I','LEU':'L','LYS':'K',
             'MET':'M','PHE':'F','PRO':'P','SER':'S',
             'THR':'T','TRP':'W','TYR':'Y','VAL':'V',
             'MSE':'M'}

def pdb_scan(pdb):
    prv = ''
    seq = ''
    for line in pdb:
        if line.startswith('ATOM'):
            if line[22:27].strip() != prv:
                seq += three2one[line[17:20]]
                prv = line[22:27].strip()
    return seq




def find_shortest_distance(res1, res2):

    min_dist = 0
    for atom1 in res1:
        for atom2 in res2:
            seta=atom1.coord
            setb=atom2.coord
            dab = math.sqrt((seta[0]-setb[0])**2+(seta[1]-setb[1])**2+(seta[2]-setb[2])**2)
            if min_dist==0 or min_dist>dab: min_dist=dab

    return min_dist



borders=[]
if ns.sequence:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    #import A3MIO
    # We can read it as fasta as we only care about the first sequence withouth gaps
    import re
    with open(ns.sequence, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq=record
            #print (record)
            break
    borders+=[len(seq)]
    
# If we have two inputs put one at each diagonal but             
if ns.inputB:
    input_fileB = np.load(ns.inputB)
    distA=dist
    distB = input_fileB["dist"]
    p_lenB = distB.shape[0]
    #resB = np.zeros((p_len, p_len))
    #resB.fill(20)
    #np.fill_diagonal(resB, 4)
    if (p_len != p_lenB):
        print ("NPZ files of different lengths")
        sys.exit(1)
    if (len(borders)!=1):
        print ("Not two chains",borders)
        sys.exit(1)
    shiftA=borders[0]
    shiftB=p_len-1-borders[0]-seplen
    for i in range(shiftA+seplen,p_len-1):
        for j in range(i+1,p_len-1):
            dist[i, j, 0:]=distB[i-shiftA+1-seplen, j-shiftA+1-seplen, 0:]
    for i in range(shiftB):
        for j in range(i+1,shiftB):
            dist[i+shiftA+seplen, j+shiftA+seplen, 0:]=distB[i, j, 0:]
    for i in range(shiftB+seplen,p_len-1):
        for j in range(shiftB):
            dist[i, j, 0:]=distB[j, i, 0:]
             
            
#print (ns.domain)
for i in range(p_len-1):
    #for j in range(i+1):
    for j in range(p_len-1):
        prob = dist[i, j, 0]
        if prob > 0.5:
            continue
        d_slice = dist[i, j, 1:]
        mean_dist = np.sum(np.multiply(bins, d_slice/np.sum(d_slice)))
        res[i, j] = mean_dist
        #res[j, i] = mean_dist

maxdist=20
if ns.pdb:
    p = PDBParser()
    structure = p.get_structure('', ns.pdb)

    #chains=[]
    pdblen=[]
    for chain in structure[0]:
        #print (chain,len(chain))
        #chains+=[chain]
        pdblen+=[len(chain)]

    dimerlen=pdblen[0]+pdblen[1]
    if (dimerlen+seplen != p_len):
        print ("PDB file is of different lengths")
        print (dimerlen, p_len,pdblen,seplen)
        sys.exit(1)
        #seplen=0
        #borders+=[pdblen[0]]
    
    pdbdist = np.zeros((dimerlen+seplen, dimerlen+seplen))

    #print (chains,pdblen)
    i=-seplen
    for chain1 in structure[0]:
        i+=seplen
        for residue1 in chain1:
        #print (residue1.get_resname())
            if (residue1.get_resname()=="GLY"):
                c1=residue1["CA"]
            else:
                try:
                    c1=residue1["CB"]
                except:
                    break
            j=-seplen    
            for chain2 in structure[0]:
                j+=seplen
                for residue2 in chain2:
                #print (residue2.get_resname())
                    if (residue2.get_resname()=="GLY"):
                        c2=residue2["CA"]
                    else:
                        try:
                            c2=residue2["CB"]
                        except:
                            break
                    pdbdist[i,j]=c1-c2
                    if i<j: res[i,j]=min(maxdist,pdbdist[i,j])
                    #print(i,j,c1,c2,c1-c2)
                    j+=1
            i+=1


#sys.exit()    
        
        
borders+=[p_len-1]        
startx=0
average=[]
mindist=[]
numdist=[]
mindist=[]
averagedist=[]
nummedcontacts=[]
fractionmedcontacts=[]
numshortcontacts=[]
fractionshortcontacts=[]
numlongcontacts=[]
fractionlongcontacts=[]
numprob=[]
fractionprob=[]
x=0
y=0
skip=5
short=5
med=8
long=12
probcut=0.5
# We only do this for two domains at the moment

def get_area(i,j,cut):

    if i<j:
        if i<cut and j<cut : x= 0
        elif i<cut and j>cut : x= 2
        else: x= 1
    else:
        if i<cut and j<cut : x= 3
        elif i>cut and j>cut : x= 4
        else: x= 5
    #print (i,j,cut,x)
    return x
extradist=2
z=[0,0,0,0,0,0]
mindist=[9999,9999,9999,9999,9999,9999]
average=[0,0,0,0,0,0]
numdist=[0,0,0,0,0,0]
numprob=[0,0,0,0,0,0]
numshortcontacts=[0,0,0,0,0,0]
numlongcontacts=[0,0,0,0,0,0]
nummedcontacts=[0,0,0,0,0,0]
averagedist=[0,0,0,0,0,0]
shortTP=[0,0,0,0,0,0]
shortFP=[0,0,0,0,0,0]
shortPPV=[0,0,0,0,0,0]
medTP=[0,0,0,0,0,0]
medFP=[0,0,0,0,0,0]
medPPV=[0,0,0,0,0,0]
longTP=[0,0,0,0,0,0]
longFP=[0,0,0,0,0,0]
longPPV=[0,0,0,0,0,0]
if (ns.sequence):
    for m in borders:
        starty=0
        for n in borders:
            #print (x,mindist,average,startx,starty,m,n)
            for i in range(startx,m):
                for j in range(starty,n):
                    # Avoid sequences separated by less than 5 resideus
                    if np.abs(i-j)<skip: continue
                    # We have six groups, let's call them 0-5
                    x=get_area(i,j,borders[0])
                    #print (i,j)
                    prob = dist[i, j, 0]
                    average[x]+=1-prob
                    z[x]+=1
                    if prob > probcut:  # Should we include all or only those with prob>probcut?
                        numprob[x]+=1
                    numdist[x]+=1
                    #d_slice = dist[i, j, 1:]
                    #mean_dist = np.sum(np.multiply(bins, d_slice/np.sum(d_slice)))
                    mean_dist=res[i,j]
                    mindist[x]=min(mean_dist,mindist[x])
                    if (mean_dist<short):
                        numshortcontacts[x]+=1
                    if (mean_dist<med):
                        nummedcontacts[x]+=1
                    if (mean_dist<long):
                        numlongcontacts[x]+=1
                    averagedist[x]+=mean_dist
                    # Now we need to check agreement for all predicted long contacts
                    if i>j:
                        if res[i,j]<short:
                            if res[j,i]<short+extradist:
                                shortTP[x]+=1
                            else:
                                shortFP[x]+=1
                        if res[i,j]<med:
                            if res[j,i]<med+extradist:
                                medTP[x]+=1
                            else:
                                medFP[x]+=1
                        if res[i,j]<long:
                            if res[j,i]<long+extradist:
                                longTP[x]+=1
                            else:
                                longFP[x]+=1
            starty=n+len(sepseq)
        startx=m+len(sepseq)

for x in range(0,6):
    print (x,z[x])
    average[x]=average[x]/z[x]
    averagedist[x]=average[x]/z[x]
    fractionprob+=[numprob[x]/z[x]]
    Z=np.sqrt(z[x])
    fractionshortcontacts+=[numshortcontacts[x]/Z]
    fractionlongcontacts+=[numlongcontacts[x]/Z]
    fractionmedcontacts+=[nummedcontacts[x]/Z]
    shortPPV[x]=shortTP[x]/(shortTP[x]+shortFP[x]+1.e-20)        
    medPPV[x]=medTP[x]/(medTP[x]+medFP[x]+1.e-20)        
    longPPV[x]=longTP[x]/(longTP[x]+longFP[x]+1.e-20)        
# o        
#sys.exit()
        
fig = plt.figure()
ax = fig.add_subplot(111)
#fig, (ax1, ax2) = plt.subplots(ncols=2)
#ax2=ax.twin()
cax = ax.matshow(res, cmap="hot")
#print (res)
if type(ns.domain) is list:
    for cut in ns.domain:
        #x=[0,p_len-1,cut,cut]
        #y=[cut,cut,0,p_len-1]
        x=[0,p_len-1]
        y=[float(cut),float(cut)]
        #print (x,y)
        ax.plot(x,y,lw=3,c="b",alpha=0.2)
        ax.plot(y,x,lw=3,c="b",alpha=0.2)
    #ax.set(xlim=[0,500],ylim=[0,500])
line=ns.input+" NumLongContacts: " + str(numlongcontacts[5])  + " LongPPV: "+str(np.round(longPPV[5],3))
ax.set(title=line)
fig.colorbar(cax)
if ns.output:
    fig.savefig(ns.output)
#plt.show()
print ("AverageProb",ns.input,np.round(average,3))
print ("Mindist",ns.input,np.round(mindist,3))
print ("Numdist",ns.input,np.round(numdist,3))
print ("AverageDistance",ns.input,np.round(averagedist,3))
print ("NumShortContacts",ns.input,np.round(numshortcontacts,3))
print ("FractionShortContacts",ns.input,np.round(fractionshortcontacts,3))
print ("NumMediumContacts",ns.input,np.round(nummedcontacts,3))
print ("FractionMediumContacts",ns.input,np.round(fractionmedcontacts,3))
print ("NumLongContacts",ns.input,np.round(numlongcontacts,3))
print ("FractionLongContacts",ns.input,np.round(fractionlongcontacts,3))
print ("NumProb",ns.input,np.round(numprob,3))
print ("FractionProb",ns.input,np.round(fractionprob,3))
print ("ShortTP",shortTP)
print ("ShortFP",shortFP)
print ("ShortPPV",np.round(shortPPV,3))
print ("MedTP",medTP)
print ("MedFP",medFP)
print ("MedPPV",np.round(medPPV,3))
print ("LongTP",longTP)
print ("LongFP",longFP)
print ("LongPPV",np.round(longPPV,3))

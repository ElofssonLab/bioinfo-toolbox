#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import argparse
from argparse import RawTextHelpFormatter
        
##args = docopt.docopt(__doc__)
#out_dir = args['--output_folder']
 

p = argparse.ArgumentParser(description = '- Merging three NPZ files (each file + interaction area)',
                            formatter_class=RawTextHelpFormatter)
p.add_argument('-data','--input','-i', required= True, help='Input trRossetta NPZ file')
#p.add_argument('-seq','--sequence','-s', required= True, help='sequence file to identify domain baorders')
#p.add_argument('-out','--output','-o', required= True, help='output NPX file')
p.add_argument('-dom','--domain','-d', required= False, help='positions of domain borders', nargs='+')
p.add_argument('-seq','--sequence','-s', required= True, help='sequence file to identify domain baorders')
#parser.add_argument('--nargs', nargs='+')
p.add_argument("--sepseq","-sep","-S",required=False, help='Separation sequence between protein in MSA' ,default="GGGGGGGGGGGGGGGGGGGG")

ns = p.parse_args()


bin_step = 0.5
bins = np.array([2.25+bin_step*i for i in range(36)])

rst = np.load(ns.input)
dist = rst["dist"]
borders=[]
if ns.sequence:
    sepseq=ns.sepseq
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    #import A3MIO
    with open(ns.sequence, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq=record
            #print (record)
            break
    if (len(seq)==dist.shape[0])
    for m in re.finditer(sepseq,str(seq.seq)):
        borders+=[m.start()]
        #print(m.start(), m.group())
        for i in range(m.start(),m.start()+len(sepseq)):
            ns.domain+=[i]

# We have three areas to calculate scores in

skip=5
# Scores to calculate


mindist=[]
averagedist=[]
numcontacts=[]
fractioncontacts=[]
numshortcontacts=[]
fractionshortcontacts=[]
numlongcontacts=[]
fractionlongcontacts=[]

def getdistdata(dist,i,j,mindist,averagedist,numcontacts,numshortcontacts,numlongcontacts):
    d_slice = dist[i, j, 1:]
    mean_dist = np.sum(np.multiply(bins, d_slice/np.sum(d_slice)))
    short=5
    contact=8
    long=12
    mindist=min(mindist,mean_dist)
    averagedist+=mean_dist
    if (mean_dist<short):
        numshortcontacts+=1
    elif (mean_dist<contact):
        numcontacts+=1
    if (mean_dist<long):
        numlongcontacts+=1
    return mindist,averagedist,numcontacts,numshortcontacts,numlongcontacts
# Seq 1
N=0
pos=0
mindist+=[9999]
averagedist+=[0]
numcontacts+=[0]
numshortcontacts+=[0]
numlongcontacts+=[0]

for i in range(len(seq)):
    for j in range(i+skip,len(seq)):
        mindist[pos],averagedist[pos],numcontacts[pos],numshortcontacts[pos],numlongcontacts[pos]=getdistdata(dist,i,j,mindist[pos],averagedist[pos],numcontacts[pos],numshortcontacts[pos],numlongcontacts[pos])
        N+=1
averagedist[pos]=averagedist[pos]/N
N=np.sqrt(2*N)
fractionshortcontacts+=[2*numshortcontacts[pos]/N]
fractioncontacts+=[2*numcontacts[pos]/N]
fractionlongcontacts+=[2*numlongcontacts[pos]/N]


# Seq 2
N=0
pos=1
mindist+=[9999]
averagedist+=[0]
numcontacts+=[0]
numshortcontacts+=[0]
numlongcontacts+=[0]
for i in range(len(seq)+1,dist.shape[0]):
    for j in range(i+skip,dist.shape[0]):
        mindist[pos],averagedist[pos],numcontacts[pos],numshortcontacts[pos],numlongcontacts[pos]=getdistdata(dist,i,j,mindist[pos],averagedist[pos],numcontacts[pos],numshortcontacts[pos],numlongcontacts[pos])
        N+=1
averagedist[pos]=averagedist[pos]/N
N=np.sqrt(2*N)
fractionshortcontacts+=[2*numshortcontacts[pos]/N]
fractioncontacts+=[2*numcontacts[pos]/N]
fractionlongcontacts+=[2*numlongcontacts[pos]/N]

# Inter-contacts
N=0
pos=2
mindist+=[9999]
averagedist+=[0]
numcontacts+=[0]
numshortcontacts+=[0]
numlongcontacts+=[0]
for i in range(len(seq)):
    for j in range(len(seq)+1,dist.shape[0]):
        mindist[pos],averagedist[pos],numcontacts[pos],numshortcontacts[pos],numlongcontacts[pos]=getdistdata(dist,i,j,mindist[pos],averagedist[pos],numcontacts[pos],numshortcontacts[pos],numlongcontacts[pos])
        N+=1
averagedist[pos]=averagedist[pos]/N
N=np.sqrt(N)
fractionshortcontacts+=[numshortcontacts[pos]/N]
fractioncontacts+=[numcontacts[pos]/N]
fractionlongcontacts+=[numlongcontacts[pos]/N]


for i in range(3):
    print (i,mindist[i],averagedist[i],numcontacts[i],fractioncontacts[i],
               numshortcontacts[i],fractionshortcontacts[i],numlongcontacts[i],
               fractionlongcontacts[i])
        



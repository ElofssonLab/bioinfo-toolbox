#!/usr/bin/env python3
import matplotlib.pyplot as plt
#import matplotlib
#matplotlib.use('Agg')
import numpy as np
import argparse
from argparse import RawTextHelpFormatter
from Bio import pairwise2  # Should perhaps replace with newer pairwise aligner - but speed is minimal.
from Bio import SeqIO
from Bio.SubsMat.MatrixInfo import blosum62
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

##args = docopt.docopt(__doc__)
#out_dir = args['--output_folder']
 

p = argparse.ArgumentParser(description = '- plotting trRosetta maps-',
                            formatter_class=RawTextHelpFormatter)
p.add_argument('-data','--input','-i', required= True, help='Input trRossetta NPZ file')
p.add_argument('-seq','--sequence','-s', required= True, help='sequence file to identify domain borders')
#p.add_argument('-ali','--alignment','-a', required= True, help='Alignment of first and second sequence')
p.add_argument("--sepseq","-sep","-S",required=False, help='Separation sequence between protein in MSA' ,default="GGGGGGGGGGGGGGGGGGGG")
p.add_argument('-out','--output','-o', required= True, help='output NPZ file')
p.add_argument('-png','--heatmap','-p', required= False, help='also save PNG files.', action='store_true')
#parser.add_argument('--nargs', nargs='+')
ns = p.parse_args()

rst = np.load(ns.input)
bin_step = 0.5
numbin=36
numomega=24
numtheta=24
numphi=12
bins = np.array([2.25+bin_step*i for i in range(numbin)])
#max_dist=2.25+bin_step*numbin+1
max_dist=20
dist = rst["dist"]
p_len = dist.shape[0]
res = np.zeros((p_len, p_len))
res.fill(max_dist)
np.fill_diagonal(res, 4)
#new_res = np.zeros((p_len, p_len))
#new_res.fill(max_dist)
#np.fill_diagonal(new_res, 4)

#print(rst)


borders=[]
if ns.sequence:
    sepseq=ns.sepseq
    #import A3MIO
    # We can read it as fasta as we only care about the first sequence withouth gaps
    import re
    with open(ns.sequence, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq=record
            #print (record)
            break
    ns.domain=[]
    for m in re.finditer(sepseq,str(seq.seq)):
        borders+=[m.start()]
        #print(m.start(), m.group())
        for i in range(m.start(),m.start()+len(sepseq)):
            ns.domain+=[i]
seplen=len(sepseq)
seqA=seq[0:borders[0]]
seqB=seq[borders[0]+seplen:]

#print (seqA.seq,seqB.seq)
alignments = pairwise2.align.localds(seqA.seq, seqB.seq, blosum62, -10, -0.5)
#print(alignments[0])
#print(alignments[0][0])
#print(alignments[0][1])
#print(alignments[0][2])
#print(alignments[0][3])
#print(alignments[0][4])


xmap=np.zeros(len(seqA.seq)+1)
ymap=np.zeros(len(seqB.seq)+1)
x=0
y=0
#print (len(seqA.seq))
#print (len(seqB.seq))
for i in (range(len(alignments[0][0]))):
    #print (i,x,y,alignments[0][0][i],alignments[0][1][i])
    X=x
    Y=y
    if alignments[0][0][i] == "-":
        xmap[x]=y
        ymap[y]=-1
        y+=1
        #x+=1
    elif alignments[0][1][i] == "-":
        xmap[x]=-1
        ymap[y]=x
        #y+=1
        x+=1
    else:
        xmap[x]=y
        ymap[y]=x
        y+=1
        x+=1
    #print (X,Y,xmap[X],ymap[Y])

        
#print(xmap,ymap)
#for i in range(len(alignments[0])

#print(pairwise2.format_alignment(*alignments[0])) # doctest:+ELLIPSIS

#sys.exit()


            
shift=0
new_rst = {'dist' : [], 'omega' : [], 'theta' : [], 'phi' : [] } # , 'rep' : []}
#new_rst=rst

#print (borders,rst["dist"].shape[0])
length=dist.shape[0]

#print("dist",rst["dist"][0,0,0:])
#print("omega",rst["omega"][0,0,0:])
#print("theta",rst["theta"][0,0,0:])
#print("phi",rst["phi"][0,0,0:])

zero_rst={'dist' : [], 'omega' : [], 'theta' : [], 'phi' : [] }
zero_rst['dist']=np.zeros(numbin+1, dtype=np.float32)
zero_rst['omega']=np.zeros(numomega+1, dtype=np.float32)
zero_rst['theta']=np.zeros(numtheta+1, dtype=np.float32)
zero_rst['phi']=np.zeros(numphi+1, dtype=np.float32)
# We do not set it to 0 as this cause problems with log,
zero_rst["dist"].fill(.1/numbin)
zero_rst["omega"].fill(.1/numomega)
zero_rst["theta"].fill(.1/numtheta)
zero_rst["phi"].fill(.1/numphi)
zero_rst['dist'][0]=0.9
zero_rst['omega'][0]=0.9
zero_rst['theta'][0]=0.9
zero_rst['phi'][0]=0.9

max_intradist=16
intradist_rst={'dist' : [], 'omega' : [], 'theta' : [], 'phi' : [] }
intradist_rst['dist']=np.zeros(numbin+1, dtype=np.float32)
intradist_rst['omega']=np.zeros(numomega+1, dtype=np.float32)
intradist_rst['theta']=np.zeros(numtheta+1, dtype=np.float32)
intradist_rst['phi']=np.zeros(numphi+1, dtype=np.float32)
# We do not set it to 0 as this cause problems with log,
intradist_rst["dist"]=np.array([0.499,
                                0.,0.,0.,0.,0.,0.,
                                0.,0.,0.,0.,0.,0.02,
                                0.02,0.02,0.02,0.02,0.02,0.02,
                                0.02,0.02,0.02,0.02,0.02,0.02,
                                0.02,0.02,0.02,0.02,0.02,0.02,
                                0.02,0.02,0.02,0.02,0.02,0.02
],dtype=np.float32)

print (intradist_rst["dist"],zero_rst["dist"])
d_slice = intradist_rst["dist"][ 1:]
print (intradist_rst["dist"].sum(),zero_rst["dist"].sum(),np.sum(np.multiply(bins, d_slice/np.sum(d_slice))))

sys.exit()
#intradist_rst["dist"].fill(.1/numbin)
intradist_rst["omega"].fill(.1/numomega)
intradist_rst["theta"].fill(.1/numtheta)
intradist_rst["phi"].fill(.1/numphi)
#intradist_rst['dist'][0]=0.6
intradist_rst['omega'][0]=0.9
intradist_rst['theta'][0]=0.9
intradist_rst['phi'][0]=0.9

new_rst["dist"]=np.copy(rst["dist"])
new_rst["omega"]=np.copy(rst["omega"])
new_rst["theta"]=np.copy(rst["theta"])
new_rst["phi"]=np.copy(rst["phi"])

#print("dist",zero_rst["dist"][0:])
#print("omega",zero_rst["omega"][0:])
#print("theta",zero_rst["theta"][0:])
#print("phi",zero_rst["phi"][0:])


for i in range(p_len-1):
    #for j in range(i+1):
    for j in range(p_len-1):
        prob = dist[i, j, 0]
        if prob > 0.5:
            #res[i, j] = max_dist
            continue
        d_slice = dist[i, j, 1:]
        mean_dist = np.sum(np.multiply(bins, d_slice/np.sum(d_slice)))
        res[i, j] = mean_dist
        #res[j, i] = mean_dist
new_res=np.copy(res)

#i=403
#j=403
#print ("TEST",i,j,res[i,j],dist[i,j])
m=borders[0]
cutoff=19
for x in range(0,m-1):
    for y in range(m+seplen,length-1):
    # intra-chain contacts - keep as is
        x2=x
        y2=int(ymap[y-(m+seplen)])
        x3=int(xmap[x])+m+seplen
        y3=y
        #print ("Test",m,seplen,x,y,x2,y2,x3,y3)
        #print ("Dist",res[x,y],res[x2,y2],res[x3,y3])
        if (  res[x,y]<cutoff):
            if (res[x2,y2]<cutoff and res[x3,y3]<cutoff): # and  (res[x,y2]<cutoff and res[x2,y]<cutoff):
                #print ("Found",x,y,x2,y2,x3,y3,res[x,y],res[x2,y2],res[x3,y3])
                new_res[x,y]=max_dist
                new_res[y,x]=max_dist
                for d in rst.files:
                    new_rst[d][x,y]=zero_rst[d]
                    new_rst[d][y,x]=zero_rst[d]
                    #new_rst[d][x,y]=new_rst[d][x,y]
                    #new_rst[d][y,x]=new_rst[d][y,x]
            elif(res[x2,y2]<cutoff or res[x3,y3]<cutoff): # and  (res[x,y2]<cutoff and res[x2,y]<cutoff):
                #print ("Found",x,y,x2,y2,x3,y3,res[x,y],res[x2,y2],res[x3,y3])
                if (res[x,y]<maxintradist):  # what should we put this to?
                    new_res[x,y]=max_intradist
                    new_res[y,x]=max_intradist
                    print ("homocontact",x,y,res[x,y],new_rst["dist"][x,y],intradist_rst["dist"])
                    for d in rst.files:
                        new_rst[d][x,y]=intradist_rst[d]
                        new_rst[d][y,x]=intradist_rst[d]
                        #new_rst[d][x,y]=new_rst[d][x,y]
                        #new_rst[d][y,x]=new_rst[d][y,x]
            else:
                print ("contact",res[x,y],x,y,new_rst["dist"][x,y],intradist_rst["dist"])
            #    #print ("Found",x,y,x2,y2,x3,y3,res[x,y],res[x2,y2],res[x3,y3])
            #    if (res[x,y]<maxintradist):
            #        new_res[x,y]=max_intradist
            #        new_res[y,x]=max_intradist
            #        for d in rst.files:
            #            new_rst[d][x,y]=intradist_rst[d]
            #            new_rst[d][y,x]=intradist_rst[d]
            #        #new_rst[d][x,y]=new_rst[d][x,y]
            #        #new_rst[d][y,x]=new_rst[d][y,x]
        #else:
        #    print (x,y,x2,y2,res[x,y],res[x2,y2],res[x,y2],res[x2,y])
        #    for d in rst.files:
        #        new_rst[d][x,y]=rst[d][x,y,0:]
        #        new_rst[d][y,x]=rst[d][y,x]
        #    new_res[x,y]=res[x,y]
        #    new_res[y,x]=res[y,x]

#print (new_rst)

# Now we delete the sepseq
#shift=0
#for m in borders:
#    print ("Deleting :",m,m+len(sepseq))
#    for x in rst.files:
#        #print (x)
#        new_rst[x]=np.delete(new_rst[x],slice(m+shift,m+shift+len(sepseq)),1)
#        new_rst[x]=np.delete(new_rst[x],slice(m+shift,m+shift+len(sepseq)),0)
#    shift+=len(sepseq)


# Save
np.savez_compressed(ns.output, dist=new_rst['dist'], omega=new_rst['omega'], theta=new_rst['theta'], phi=new_rst['phi'])# , rep=new_rst['rep'])

#sys.exit()

if (ns.heatmap):
    outfig1=ns.output+"-org.png"
    outfig2=ns.output+"-new.png"
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cax = ax.matshow(res, cmap="hot")
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
    line="Original plot"
    ax.set(title=line)
    fig.colorbar(cax)
    fig.savefig(outfig1)
    
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cax = ax.matshow(new_res, cmap="hot")
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
    line="updated plot"
    ax.set(title=line)
    fig.colorbar(cax)
    fig.savefig(outfig2)

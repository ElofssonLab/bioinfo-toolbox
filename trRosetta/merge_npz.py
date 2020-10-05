#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import argparse
from argparse import RawTextHelpFormatter
        
##args = docopt.docopt(__doc__)
#out_dir = args['--output_folder']
 

p = argparse.ArgumentParser(description = '- Merging three NPZ files (each file + interaction area)',
                            formatter_class=RawTextHelpFormatter)
p.add_argument('-dataA','--inputA','-i', required= True, help='Input trRossetta NPZ file')
p.add_argument('-dataB','--inputB','-j', required= True, help='Input trRossetta NPZ file')
p.add_argument('-dataAB','--inputAB','-k', required= True, help='Input trRossetta NPZ file(s) for the mergedsequence',nargs="+")
#p.add_argument('-seq','--sequence','-s', required= True, help='sequence file to identify domain baorders')
p.add_argument('-out','--output','-o', required= True, help='output NPZ file')
p.add_argument('-fast','--fast','-f', required= False, default="False",help='Do not merge NPZ files position by position',action='store_true')
#parser.add_argument('--nargs', nargs='+')
ns = p.parse_args()


bin_step = 0.5
bins = np.array([2.25+bin_step*i for i in range(36)])


rstA = np.load(ns.inputA)
rstB = np.load(ns.inputB)
#print ("reading ",ns.inputA)
#print ("reading ",ns.inputB)
rstAB=[]
for i in ns.inputAB:
    #print ("reading ",i)
    rstAB += [np.load(i)]

#print (rstAB)


distA = rstA["dist"].shape[0]
distB = rstB["dist"].shape[0]
distAB = rstAB[0]["dist"].shape[0]


if distA+distB != distAB:
    
    print ("Not correct sized of distnance matrices",distA,distB,distAB)
    exit(-1)

for r in rstAB:
    if r["dist"].shape[0] !=distAB:
        print ("Not correct sized of distance matrices",d,rstAB[r]["dist"].shape[0],distAB)
        exit(-1)

new_rst = {'dist' : [], 'omega' : [], 'theta' : [], 'phi' : [] } # , 'rep' : []}    

# first we merge the full length ones.
for f in rstAB[0].files:
    new_rst[f]=rstAB[0][f]


            
#bins = np.array([2.25+bin_step*i for i in range(36)])    
# Take the one with the lowsest  probability to not have a distance
# To speed up things we only do this for intrachan
if (ns.fast):
    orgprob=np.mean(new_rst["dist"][distA:distAB,0:distA,0])
    orgcontacts=(new_rst["dist"][distA:distAB,0:distA, 0]  < 0.5 ).sum()

    for d in range(1,len(rstAB)):
        newprob=np.mean(rstAB[d]["dist"][distA:distAB,0:distA,0])
        newcontacts=(rstAB[d]["dist"][distA:distAB,0:distA, 0]  < 0.5 ).sum()
        #print ("Test",d,orgprob,newprob,orgcontacts,newcontacts,
        #       np.mean(rstAB[d]["dist"][0:distA,0:distA,0]),
        #       np.mean(rstAB[d]["dist"][distA+1:distAB,distA+1:distAB,0]))
        #if newprob < orgprob:
        if newcontacts > orgcontacts:
            for f in rstAB[0].files:
                new_rst[f][distA:distAB,0:distA]=rstAB[d][f][distA:distAB,0:distA]
                new_rst[f][0:distB,distB:distAB]=rstAB[d][f][0:distB,distB:distAB]
            orgprob=newprob
            orgcontacts=newcontacts
else:
    for d in range(1,len(rstAB)):
        #print ("Using ",d,distA,distB)
        for i in range(distA):
            for j in range(distA+1,distA+distB):
                orgprob = new_rst["dist"][i, j, 0]
                newprob = rstAB[d]["dist"][i, j, 0]
                if newprob < orgprob:
                    for f in rstAB[0].files:
                        new_rst[f][i,j]=rstAB[d][f][i,j]
                        new_rst[f][j,i]=rstAB[d][f][j,i]

                    
# Then we add the distance for the individual proteins   
for i in rstAB[0].files:
    A=rstA[i]
    B=rstB[i]
    AB=new_rst[i]
    AB[0:distA,0:distA]=A    
    AB[distA:,distA:]=B    
    new_rst[i]=AB



np.savez_compressed(ns.output, dist=new_rst['dist'], omega=new_rst['omega'], theta=new_rst['theta'], phi=new_rst['phi']) #, rep=new_rst['rep'])

#!/usr/bin/env python3

import os
import re
import numpy as np
import math
import scipy
import pandas as pd
import matplotlib.pyplot as plt
#import seaborn as sns

def moving_average(a, n=50) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

dir="predictions/"
fig_dir="figures/"


dic_colors = {"E":"#004D40", "B":"#D81B60","A":"#1E88E5","All":"black",
              "E1":"green", "B1":"red", "A1":"blue","All1":"grey",
              "E2":"lightgreen","B2":"pink", "A2":"lightblue","All2":"lightgrey"}
dic_markers = {"E":"o", "B":"x", "A":"o"}

factor=10

flist = open('formatted_list','r')

# Make sure the keys are not duplicated..

GCfile={"outpred_noGC-cross-pro.csv":"outpred_GC-cross-pro.csv",
        "outpred_GC-cross-rna.csv":"outpred_noGC-cross-rna.csv",
        "outpred_GCgenomic-pro.csv":"outpred_GC-cross-pro.csv",
        "outpred_GC-cross-pro.csv":"outpred_GC-cross-rna.csv",
        "outpred_noGC-cross-rna.csv":"outpred_noGC-cross-pro.csv"}

for f in GCfile:
   df=pd.read_csv(dir+f, sep=',',header=0)
   df2=pd.read_csv(dir+GCfile[f], sep=',',header=0)
   #print (df)
   print (f,GCfile[f])
   file1=re.sub(r'.csv','',f)   
   file2=re.sub(r'.csv','',GCfile[f])   
   fig, ax = plt.subplots(figsize=(6,6))
   for kingdom in ["B","E","A"]:
       if kingdom=="All":
           d=df.sort_values('gc').dropna()
           d2=df2.sort_values('gc').dropna()
       else:
           d=df.loc[(df.kingdom == kingdom) ].sort_values('gc').dropna()
           d2=df2.loc[(df2.kingdom == kingdom) ].sort_values('gc').dropna()
       #print(d)
       numave=int(factor*math.sqrt(len(d)))
       x=(moving_average(d.gc.to_list(),numave))
       x2=(moving_average(d2.gc.to_list(),numave))
       pred=(moving_average(d.Pred.to_list(),numave))
       pred2=(moving_average(d2.Pred.to_list(),numave))
       iupred=(moving_average(d.IUPRED.to_list(),numave))
       diso=(moving_average(d.Diso.to_list(),numave))
       diso2=(moving_average(d2.Diso.to_list(),numave))
       err=pred-diso
       err2=pred2-diso2
       iuerr=iupred-diso
       ppv=(moving_average(d.PPV.to_list(),numave))
       mcc=(moving_average(d.MCC.to_list(),numave))
       f1=(moving_average(d.F1.to_list(),numave))
       #plt.plot(x,diso,label="Diso-"+kingdom,color=dic_colors[kingdom],lw=2)
       #plt.plot(x,pred,'--',label="Pred-"+kingdom,color=dic_colors[kingdom],lw=0.5)

       plt.plot(x,err,label=file1+"-"+kingdom,color=dic_colors[kingdom],lw=0.5)
       plt.plot(x2,err2,'--',label=file2+"-"+kingdom,color=dic_colors[kingdom],lw=2)
       #plt.plot(x,iuerr,'.',label="Error-IUpred-"+kingdom,color=dic_colors[kingdom],lw=2)
       
       #plt.plot(d.gc,d.Diso,'.',label="Diso-"+kingdom,color=dic_colors[kingdom])
       #plt.plot(d.gc,d.Pred,'.',label="Pred-"+kingdom,color=dic_colors[kingdom])
       
       #plt.plot(x,mcc,'.',label="MCC-"+kingdom,color=dic_colors[kingdom],lw=1)
       #plt.plot(x,f1,'-',label="F1-"+kingdom,color=dic_colors[kingdom],lw=1)
       #plt.plot(x,ppv,'--',label="PPV-"+kingdom,color=dic_colors[kingdom],lw=1)
   ax.set_title="Error in fraction disorder per protein"
   ax.set_xlabel("GC")
   ax.set_ylabel("Error in Disorder")
   ax.legend(fontsize=6)
   
   fig.savefig(fig_dir+file1+"-"+file2+"-comparison.eps",rasterized=True)

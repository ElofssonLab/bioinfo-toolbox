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




for f in os.listdir(dir):
    if f.endswith(".csv"):
        file=re.sub(r'.csv','',f)
        df=pd.read_csv(dir+file+".csv", sep=',',header=0)
        #print (df)

        fig, ax = plt.subplots(figsize=(6,6))
    
        for kingdom in ["All","B","E","A"]:
            if kingdom=="All":
                d=df.sort_values('gc').dropna()
            else:
                d=df.loc[(df.kingdom == kingdom) ].sort_values('gc').dropna()
            #print(d)
            numave=int(factor*math.sqrt(len(d)))
            x=(moving_average(d.gc.to_list(),numave))
            pred=(moving_average(d.Pred.to_list(),numave))
            iupred=(moving_average(d.IUPRED.to_list(),numave))
            diso=(moving_average(d.Diso.to_list(),numave))
            err=pred-diso
            iuerr=iupred-diso
            ppv=(moving_average(d.PPV.to_list(),numave))
            mcc=(moving_average(d.MCC.to_list(),numave))
            f1=(moving_average(d.F1.to_list(),numave))
            #plt.plot(x,diso,label="Diso-"+kingdom,color=dic_colors[kingdom],lw=2)
            #plt.plot(x,pred,'--',label="Pred-"+kingdom,color=dic_colors[kingdom],lw=0.5)

            plt.plot(x,err,label="Error-Diso-"+kingdom,color=dic_colors[kingdom],lw=0.5)
            plt.plot(x,iuerr,label="Error-IUpred-"+kingdom,color=dic_colors[kingdom],lw=2)
            
            #plt.plot(d.gc,d.Diso,'.',label="Diso-"+kingdom,color=dic_colors[kingdom])
            #plt.plot(d.gc,d.Pred,'.',label="Pred-"+kingdom,color=dic_colors[kingdom])
            
            #plt.plot(x,mcc,'.',label="MCC-"+kingdom,color=dic_colors[kingdom],lw=1)
            #plt.plot(x,f1,'-',label="F1-"+kingdom,color=dic_colors[kingdom],lw=1)
            #plt.plot(x,ppv,'--',label="PPV-"+kingdom,color=dic_colors[kingdom],lw=1)
        ax.set_title="Error in fraction disorder per protein"
        ax.set_xlabel("GC")
        ax.set_ylabel("Error in Disorder")
        ax.legend()
        
        fig.savefig(fig_dir+file+"-running-average.eps",rasterized=True)

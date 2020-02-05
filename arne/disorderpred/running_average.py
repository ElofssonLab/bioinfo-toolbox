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


dic_colors = {"E":"#004D40", "B":"#D81B60","A":"#1E88E5",
              "E1":"green", "B1":"red", "A1":"blue",
              "E2":"lightgreen","B2":"pink", "A2":"lightblue"}
dic_markers = {"E":"o", "B":"x", "A":"o"}

factor=30
for f in os.listdir(dir):
    if f.endswith(".csv"):
        file=re.sub(r'.csv','',f)
        df=pd.read_csv(dir+file+".csv", sep=',',header=0)
        #print (df)

        fig, ax = plt.subplots(figsize=(6,6))
    
        for kingdom in ["B","E","A"]:
            d=df.loc[(df.kingdom == kingdom) ].sort_values('gc').dropna()
            #print(d)
            numave=int(factor*math.sqrt(len(d)))
            x=(moving_average(d.gc.to_list(),numave))
            pred=(moving_average(d.Pred.to_list(),numave))
            diso=(moving_average(d.Diso.to_list(),numave))
            ppv=(moving_average(d.PPV.to_list(),numave))
            mcc=(moving_average(d.MCC.to_list(),numave))
            f1=(moving_average(d.F1.to_list(),numave))
            plt.plot(x,diso,label="Diso-"+kingdom,color=dic_colors[kingdom],lw=1)
            plt.plot(x,pred,label="Pred-"+kingdom,color=dic_colors[kingdom],lw=2)
            #plt.plot(x,mcc,'.',label="MCC-"+kingdom,color=dic_colors[kingdom],lw=1)
            #plt.plot(x,f1,'-',label="F1-"+kingdom,color=dic_colors[kingdom],lw=1)
            plt.plot(x,ppv,'--',label="PPV-"+kingdom,color=dic_colors[kingdom],lw=1)
            ax.legend()
        ax.set_title=file
        ax.set_xlabel("GC")
        ax.set_ylabel("Disorder")
        fig.savefig(fig_dir+file+"-running-average.eps",rasterized=True)

plt.close()
for f in os.listdir(dir):
    if f.endswith(".csv"):
        if re.match('.*noGC.*',f):
            file1=re.sub(r'.csv','',f)
            set1='formatted_list' in f
            type1='rna' in f 
            for g in os.listdir(dir):
                if g.endswith(".csv"):
                    if re.match(r'.*\_GC\_.*',g):
                        file2=re.sub(r'.csv','',g)
                        set2='formatted_list' in g
                        type2='rna' in g
                        #print ("Test",file1,file2,set1,set2,type1,type2)
                        if (set1==set2 and type1==type2):
                            df1=pd.read_csv(dir+file1+".csv", sep=',',header=0)
                            df2=pd.read_csv(dir+file2+".csv", sep=',',header=0)
                            fig, ax = plt.subplots(figsize=(6,6))
                            for kingdom in ["B","E","A"]:
                                d1=df1.loc[(df1.kingdom == kingdom) ].sort_values('gc').dropna()
                                d2=df2.loc[(df2.kingdom == kingdom) ].sort_values('gc').dropna()
                                #print(d)
                                numave=int(factor*math.sqrt(len(d1)))
                                x1=(moving_average(d1.gc.to_list(),numave))
                                x2=(moving_average(d2.gc.to_list(),numave))
                                pred1=(moving_average(d1.Pred.to_list(),numave))
                                pred2=(moving_average(d2.Pred.to_list(),numave))
                                diso1=(moving_average(d1.Diso.to_list(),numave))
                                ppv1=(moving_average(d1.PPV.to_list(),numave))
                                mcc1=(moving_average(d1.MCC.to_list(),numave))
                                f11=(moving_average(d1.F1.to_list(),numave))
                                ppv2=(moving_average(d2.PPV.to_list(),numave))
                                mcc2=(moving_average(d2.MCC.to_list(),numave))
                                f12=(moving_average(d2.F1.to_list(),numave))
                                plt.plot(x1,diso1,label="Diso-"+kingdom,color=dic_colors[kingdom],lw=1)
                                plt.plot(x1,pred1,label="Pred-noGC-"+kingdom,color=dic_colors[kingdom+"1"],lw=1)
                                plt.plot(x2,pred2,label="Pred-GC-"+kingdom,color=dic_colors[kingdom+"2"],lw=1)
                                #fig.plot(x1,mcc1,'.',label="MCC-"+kingdom,color=dic_colors[kingdom],lw=1)
                                #plt.plot(x1,f11,'-',label="F1-"+kingdom,color=dic_colors[kingdom],lw=1)
                                #plt.plot(x2,f12,'-',label="F1-"+kingdom,color=dic_colors[kingdom],lw=1)
                                #fig.plot(x2,ppv1,'--',label="PPV_noGC-"+kingdom,color=dic_colors[kingdom+"1"],lw=1)
                                #fig.plot(x2,ppv2,'--',label="PPV_GC-"+kingdom,color=dic_colors[kingdom+"2"],lw=2)
                            ax.legend()
                            ax.set_title=file1+file2
                            ax.set_xlabel("GC")
                            ax.set_ylabel("Disorder")
                            fig.savefig(fig_dir+file1+"-"+file2+"-comparison.eps",rasterized=True)
                            plt.close()
                            fig, ax = plt.subplots(figsize=(6,6))
                            for kingdom in ["B","E","A"]:
                                d1=df1.loc[(df1.kingdom == kingdom) ].sort_values('gc').dropna()
                                d2=df2.loc[(df2.kingdom == kingdom) ].sort_values('gc').dropna()
                                #print(d)
                                numave=int(factor*math.sqrt(len(d1)))
                                x1=(moving_average(d1.gc.to_list(),numave))
                                x2=(moving_average(d2.gc.to_list(),numave))
                                pred1=(moving_average(d1.Pred.to_list(),numave))
                                pred2=(moving_average(d2.Pred.to_list(),numave))
                                diso1=(moving_average(d1.Diso.to_list(),numave))
                                ppv1=(moving_average(d1.PPV.to_list(),numave))
                                mcc1=(moving_average(d1.MCC.to_list(),numave))
                                f11=(moving_average(d1.F1.to_list(),numave))
                                ppv2=(moving_average(d2.PPV.to_list(),numave))
                                mcc2=(moving_average(d2.MCC.to_list(),numave))
                                f12=(moving_average(d2.F1.to_list(),numave))
                                #plt.plot(x1,diso1,label="Diso-"+kingdom,color=dic_colors[kingdom],lw=1)
                                #plt.plot(x1,pred1,label="Pred-noGC-"+kingdom,color=dic_colors[kingdom+"1"],lw=1)
                                #plt.plot(x2,pred2,label="Pred-GC-"+kingdom,color=dic_colors[kingdom+"2"],lw=1)
                                #fig.plot(x1,mcc1,'.',label="MCC-"+kingdom,color=dic_colors[kingdom],lw=1)
                                plt.plot(x1,f11,'-',label="F1_noGC-"+kingdom,color=dic_colors[kingdom+"1"],lw=1)
                                plt.plot(x2,f12,'-',label="F1_GC-"+kingdom,color=dic_colors[kingdom+"2"],lw=1)
                                #plt.plot(x2,ppv1,'--',label="PPV_noGC-"+kingdom,color=dic_colors[kingdom+"1"],lw=2)
                                #plt.plot(x2,ppv2,'--',label="PPV_GC-"+kingdom,color=dic_colors[kingdom+"2"],lw=2)
                            ax.legend()
                            ax.set_title=file1+file2
                            ax.set_xlabel("GC")
                            ax.set_ylabel("Disorder")
                            fig.savefig(fig_dir+file1+"-"+file2+"-PPV.eps",rasterized=True)
                                

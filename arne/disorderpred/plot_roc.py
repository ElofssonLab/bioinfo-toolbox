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
    if f.endswith(".roc"):
        file=re.sub(r'.roc','',f)
        df=pd.read_csv(dir+file+".roc", sep=',',header=0)
        #print (df)

        fig, ax = plt.subplots(figsize=(6,6))
    
        plt.plot(df.FPR,df.TPR,label=file,lw=0.5)
        ax.set_title=file
        ax.set_xlabel("FPR")
        ax.set_ylabel("TPR")
        ax.legend()
        
        fig.savefig(fig_dir+file+"-ROC.eps",rasterized=True)

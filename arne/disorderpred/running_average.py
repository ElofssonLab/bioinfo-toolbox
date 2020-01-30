#!/usr/bin/env python3

import numpy as np
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
for file in ["outpred_model_200-1-0.001_1_pro_GC_42","outpred_model_200-1-0.001_1_rna_GC_42","outpred_model_200-10-0.001_1_pro_noGC_72","outpred_model_200-10-0.001_1_rna_noGC_65"]:

    df=pd.read_csv(dir+file+".csv", sep=',',header=0)
    #print (df)

    fig, ax = plt.subplots(figsize=(6,6))
    
    for kingdom in ["B","E"]: #,"A"]:
        d=df.loc[(df.kingdom == kingdom) ].sort_values('gc').dropna()
        #print(d)
        temp=[d.Pred.to_list()]
        x=(moving_average(d.gc.to_list(),500))
        pred=(moving_average(d.Pred.to_list(),500))
        diso=(moving_average(d.Diso.to_list(),500))
        plt.plot(x,diso,label="Diso-"+kingdom)
        plt.plot(x,pred,label="Pred-"+kingdom)
        ax.legend()
    ax.set_title=file
    ax.set_xlabel("GC")
    ax.set_ylabel("Disorder")
    fig.savefig(fig_dir+file+"running-average.eps",rasterized=True)


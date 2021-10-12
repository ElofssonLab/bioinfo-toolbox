#!/usr/bin/env python3

import pandas as pd
df_dockqreorder=pd.read_csv("DockQ.csv",sep=",",names=["id","DockQreorder"])
df_MMreorder=pd.read_csv("MM.csv",sep=",",names=["id","MMreorder"])
df_dockq=pd.read_csv("DockQ.csv",sep=",",names=["id","DockQ"])
df_MM=pd.read_csv("MM.csv",sep=",",names=["id","MM"])
df_MMall=pd.read_csv("MMall.csv",sep=",",names=["id","MMall"])

columns=["id","NumRes","IF_plDDT","plDDT","SumIF_pLDDT"]

df_pLDDT=pd.read_csv("pLDDT.csv",sep=",",names=columns)
#df_pLDDT["SumIF_pLDDT"]=df_pLDDT.pd_plDDT_av*df_pLDDT.NumRes

df_temp=pd.merge(df_dockq,df_MM,on=["id"],how="outer")
df_temp2=pd.merge(df_temp,df_MMall,on=["id"],how="outer")
df_temp=pd.merge(df_temp2,df_MMreorder,on=["id"],how="outer")
df_temp2=pd.merge(df_temp,df_dockqreorder,on=["id"],how="outer")
df_merged=pd.merge(df_temp2,df_pLDDT,on=["id"],how="outer")
df_merged.to_csv("evaluation.csv")

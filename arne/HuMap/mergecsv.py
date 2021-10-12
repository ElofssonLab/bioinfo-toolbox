#!/usr/bin/env python3

import pandas as pd
df_dockqreorder=pd.read_csv("DockQreorder.csv",sep=",",names=["Name","DockQreorder"])
df_MMreorder=pd.read_csv("MMreorder.csv",sep=",",names=["Name","MMreorder"])
df_dockq=pd.read_csv("DockQ.csv",sep=",",names=["Name","DockQ"])
df_dockqall=pd.read_csv("DockQall.csv",sep=",",names=["Name","DockQall"])
df_MM=pd.read_csv("MM.csv",sep=",",names=["Name","MM"])
df_MMall=pd.read_csv("MMall.csv",sep=",",names=["Name","MMall"])
df_petras=pd.read_csv("map_petras.csv",sep=",")

columns=["Name","NumRes","IF_plDDT","plDDT"]

df_pLDDT=pd.read_csv("pLDDT.csv",sep=",",names=columns)
#df_pLDDT["SumIF_pLDDT"]=df_pLDDT.pd_plDDT_av*df_pLDDT.NumRes

df_temp=pd.merge(df_dockq,df_MM,on=["Name"],how="outer")
df_temp2=pd.merge(df_temp,df_MMall,on=["Name"],how="outer")
df_temp=pd.merge(df_temp2,df_MMreorder,on=["Name"],how="outer")
df_temp2=pd.merge(df_temp,df_dockqreorder,on=["Name"],how="outer")
df_temp=pd.merge(df_temp2,df_dockqall,on=["Name"],how="outer")
df_temp2=pd.merge(df_temp,df_petras,on=["Name"],how="outer")
df_merged=pd.merge(df_temp2,df_pLDDT,on=["Name"],how="outer")
df_merged["SumIFplDDT"]=df_merged.NumRes*df_merged.IF_plDDT
df_merged.to_csv("evaluation.csv")

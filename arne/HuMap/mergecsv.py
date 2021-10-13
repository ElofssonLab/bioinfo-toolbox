#!/usr/bin/env python3

import pandas as pd
df_dockqreorder=pd.read_csv("DockQreorder.csv",sep=",",names=["Name","DockQreorder"])
df_MMreorder=pd.read_csv("MMreorder.csv",sep=",",names=["Name","MMreorder"])
df_dockq=pd.read_csv("DockQ.csv",sep=",",names=["Name","DockQ"])
df_dockqall=pd.read_csv("DockQall.csv",sep=",",names=["Name","DockQall"])
df_MM=pd.read_csv("MM.csv",sep=",",names=["Name","MM"])
df_MMall=pd.read_csv("MMall.csv",sep=",",names=["Name","MMall"])
df_petras=pd.read_csv("map_petras.csv",sep=",")
df_seqlen=pd.read_csv("seqlen.csv",sep=",",names=["Name","SeqLen1","SeqLen2"])
df_seqlen["SeqLen"]=df_seqlen.SeqLen1+df_seqlen.SeqLen2

columns=["Name","NumRes","IF_plDDT","plDDT"]

df_plDDT=pd.read_csv("pLDDT.csv",sep=",",names=columns)
df_plDDT["SumIF"]=df_plDDT.IF_plDDT*df_plDDT.NumRes

df_temp=pd.merge(df_dockq,df_MM,on=["Name"],how="outer")
df_temp2=pd.merge(df_temp,df_MMall,on=["Name"],how="outer")
df_temp=pd.merge(df_temp2,df_MMreorder,on=["Name"],how="outer")
df_temp2=pd.merge(df_temp,df_dockqreorder,on=["Name"],how="outer")
df_temp=pd.merge(df_temp2,df_dockqall,on=["Name"],how="outer")
df_temp2=pd.merge(df_temp,df_petras,on=["Name"],how="outer")
df_temp=pd.merge(df_temp2,df_seqlen,on=["Name"],how="outer")
df_merged=pd.merge(df_temp,df_plDDT,on=["Name"],how="outer")
df_merged["SumIFplDDT"]=df_merged.NumRes*df_merged.IF_plDDT
df_merged.to_csv("evaluation.csv")



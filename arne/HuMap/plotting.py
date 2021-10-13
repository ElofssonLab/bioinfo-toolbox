#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sklearn.metrics as metrics


# In[2]:


from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline


# In[102]:


df=pd.read_csv("evaluation.csv")
df_negpLDDT=pd.read_csv("negatome-pLDDT.csv")
df_negseqlen=pd.read_csv("negatome-seqlen.csv",sep=",",names=["Name","SeqLen1","SeqLen2"])
df_negseqlen["SeqLen"]=df_negseqlen.SeqLen1+df_negseqlen.SeqLen2


# In[109]:


df_negatome=pd.merge(df_negpLDDT,df_negseqlen,on=["Name"],how="inner")
df_negatome


# In[110]:


df["SumLogID_pLDDT"]=df.IF_plDDT*np.log(df.NumRes)
df["SumIF"]=df.NumRes*df.IF_plDDT/100
df_negatome["SumIF"]=df_negatome.NumRes*df_negatome.IF_plDDT/100


# In[111]:


minreses=[1,2,5,10,20,30,50]
minrescols=[]
minIFcols=[]
for i in minreses:
    df["minres-"+str(i)]= np.where(df.NumRes>i, 1, 0)
    df_negatome["minres-"+str(i)]= np.where(df_negatome.NumRes>i, 1, 0)
    minrescols+=["minres-"+str(i)]
    df["IFmin-"+str(i)]=df["IF_plDDT"]*df["minres-"+str(i)]
    df_negatome["IFmin-"+str(i)]=df_negatome["IF_plDDT"]*df_negatome["minres-"+str(i)]
    minIFcols+=["IFmin-"+str(i)]


# In[112]:


df["Struct"]= np.where(df.DockQ.notna(), True, False)
df_struct=df[df.DockQ.notna()]
df_nostruct=df[df.DockQ.isna()]
#df_struct["Good"]= np.where(df_struct.DockQ>0.23, True, False)
#df_struct["TMGood"]= np.where(df_struct.MMall>0.7, True, False)
df_struct["GoodAll"]= np.where(df_struct.DockQ>0.23, True, False)
df_struct["TMGoodAll"]= np.where(df_struct.MMall>0.7, True, False)
df_corr=df_struct[df_struct.GoodAll==True]
df_incorr=df_struct[df_struct.GoodAll==False]


# In[113]:


sns.distplot(df_corr.IF_plDDT,label="PDB-correct")
sns.distplot(df_incorr.IF_plDDT,label="PDB-incorrect")
sns.distplot(df_nostruct.IF_plDDT,label="NoPDB")
sns.distplot(df_negatome.IF_plDDT,label="Negatome")
plt.legend()

plt.savefig("pLDDT.png")


# In[114]:


ax=sns.distplot(df_corr.NumRes,label="PDB-correct")
sns.distplot(df_incorr.NumRes,label="PDB-incorrect")
sns.distplot(df_nostruct.NumRes,label="NoPDB")
sns.distplot(df_negatome.NumRes,label="Negatome")
ax.set_xlim([0,500])
plt.legend()
plt.show()
plt.savefig("NumRes.png")


# In[115]:


ax=sns.distplot(df_corr.SumIF,label="PDB-correct")
sns.distplot(df_incorr.SumIF,label="PDB-incorrect")
sns.distplot(df_nostruct.SumIF,label="NoPDB")
sns.distplot(df_negatome.SumIF,label="Negatome")
#ax.set_xscale('log')
ax.set_xlim(0,300)
plt.legend()
plt.show()
plt.savefig("SumIF.png")


# In[129]:


f, ax = plt.subplots(figsize=(6.5, 6.5))
reg = LinearRegression()
cutoff=0.23
tempdf=df_struct.dropna(subset=["IF_plDDT","SumIF","NumRes"])
X=tempdf[["IF_plDDT","SumIF","NumRes","SeqLen"]]
Y=tempdf.DockQall
R=reg.fit(X,Y)
pred=reg.predict(X)
plt.scatter(Y,pred)
err=mean_squared_error(pred,Y)
correct=Y>cutoff            
#fig2.show()
fpr, tpr, threshold = metrics.roc_curve(correct, pred)
roc_auc = metrics.auc(fpr, tpr)
corr=np.corrcoef(Y,pred)[0,0]
ax.set_title('pDockQ, Cx:'+str(round(corr,2))+" ROC: "+str(round(roc_auc,2)))
    
plt.legend(loc = 'upper left')
plt.plot([0, 1], [0, 1],'r--')
plt.xlim([0, 1])
plt.ylim([0, 1])
plt.ylabel('pDockQ')
plt.xlabel('DockQ')
plt.savefig("pDockQ.png")
#print (corr,roc_auc)


# In[130]:


f, ax = plt.subplots(figsize=(6.5, 6.5))
tempdf=df_struct.dropna(subset=["IF_plDDT","SumIF","NumRes"])
correct=tempdf["GoodAll"]
X=tempdf[["IF_plDDT","SumIF","NumRes","SeqLen"]]
tempdf["pred"]=reg.predict(X)
for d in ["IF_plDDT","SumIF","NumRes","IFmin-50","pred"]:
    values=tempdf[d]
    fpr, tpr, threshold = metrics.roc_curve(correct, values)
    roc_auc = metrics.auc(fpr, tpr)
    plt.plot(fpr, tpr, label = d+': AUC = %0.2f' % roc_auc)
    ax.set_title('ROC: Correct vs Wrong Structure')
    
    plt.legend(loc = 'lower right')
    plt.plot([0, 1], [0, 1],'r--')
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
plt.savefig("ROC-error.png")


# In[131]:


df_struct[df_struct.MMall>0.9].dropna(subset=["IF_plDDT","SumIF","NumRes"]).sort_values("DockQall")


# In[132]:


df_negatome["GoodAll"]=False
df_negatome["Struct"]=False
df_negatome
df_concat=pd.concat([df_negatome,df_struct])[["SeqLen","Name","IF_plDDT","plDDT","NumRes","SumIF","GoodAll","Struct"]+minIFcols]
df_concat2=pd.concat([df_negatome,df_struct[df_struct.DockQ>0.23]])[["SeqLen","Name","IF_plDDT","plDDT","NumRes","SumIF","GoodAll","Struct"]+minIFcols]

df_concat


# In[133]:


f, ax = plt.subplots(figsize=(6.5, 6.5))
tempdf=df_concat.dropna()
correct=tempdf["Struct"]
X=tempdf[["IF_plDDT","SumIF","NumRes","SeqLen"]]
tempdf["pred"]=reg.predict(X)

for d in ["IF_plDDT","SumIF","NumRes","pred","IFmin-50"]:
    values=tempdf[d]
    fpr, tpr, threshold = metrics.roc_curve(correct, values)
    roc_auc = metrics.auc(fpr, tpr)
    plt.plot(fpr, tpr, label = d+': AUC = %0.2f' % roc_auc)
    ax.set_title('ROC: Interaction vs Negatome')
    
    plt.legend(loc = 'lower right')
    plt.plot([0, 1], [0, 1],'r--')
    plt.xlim([0, .1])
    plt.ylim([0, 1])
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
plt.savefig("ROC-negatome.png")


# In[134]:


f, ax = plt.subplots(figsize=(6.5, 6.5))
tempdf=df_concat2.dropna()
correct=tempdf["Struct"]
X=tempdf[["IF_plDDT","SumIF","NumRes","SeqLen"]]
tempdf["pred"]=reg.predict(X)

for d in ["IF_plDDT","SumIF","NumRes",'IFmin-50',"pred"]:
    values=tempdf[d]
    fpr, tpr, threshold = metrics.roc_curve(correct, values)
    roc_auc = metrics.auc(fpr, tpr)
    plt.plot(fpr, tpr, label = d+': AUC = %0.2f' % roc_auc)
    ax.set_title('Receiver Operating Characteristic')
    
    plt.legend(loc = 'lower right')
    plt.plot([0, 1], [0, 1],'r--')
    plt.xlim([0, .1])
    plt.ylim([0, 1])
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')


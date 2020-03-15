#!/usr/bin/env python3


from __future__ import division
from __future__ import print_function

import pandas as pd  
import numpy as np  
import matplotlib.pyplot as plt  
import seaborn as seabornInstance 
from sklearn.model_selection import train_test_split 
from sklearn.linear_model import LinearRegression
from sklearn import metrics
import os
# %matplotlib inline
import docopt
from dateutil import parser

#args = docopt.docopt(__doc__)
#out = args['--output_folder']

out ="/Users/arnee/Desktop/covidify-output/"
def string_format(x):
    return '%s' % x

def float_format(x):
    return '%.0f' % x

def float_format2(x):
    return '%.3f' % x

def exp_format(x):
    return '%.3E' % x


def antest(df,k,l,columns):
    aagc = {}

    aagcnon = {}
    pr = {}
    for col in columns:
        lm = smf.ols( col + ' ~ kingdom +  GCcoding' , data=df).fit()
        #print(lm.summary())
        table = sm.stats.anova_lm(lm, typ=2)
        #print (table)
        aagc[dic_labels[col]] = [table["F"]["kingdom"]]
    Fvalue = pd.DataFrame(aagc,index =['F-test'])
    Ft=Fvalue.transpose()
    return Ft

dataset = pd.read_csv('daily-data.csv')

mindeath=10
maxdeath=225
bigset=dataset.loc[dataset['Death'] >mindeath].sort_values(by='Death')
bigset=bigset.loc[dataset['Death'] <maxdeath]

image_dir =  os.path.join(out,'reports', 'images')

if not os.path.exists(image_dir):
    print('Creating reports folder...')
    os.system('mkdir -p ' + image_dir)

inp=''
mindays=20
data_cols=[]
for j in range(mindays,0,-1):
    data_cols+= ["-"+str(j)]
#print (inp)
target_col=['Death']
#print (dataset.keys())
y = dataset[target_col]
X = dataset[data_cols]
model = LinearRegression()
model.fit(X, y)
r_sq = model.score(X, y)
#print (model.coef_,r_sq)
#print (X,y)
cof=model.coef_[0]



y = bigset[target_col]
X = bigset[data_cols]
model = LinearRegression()
model.fit(X, y)
r_sq = model.score(X, y)
#print (model.coef_,r_sq)
#print (X,y)
bigcof=model.coef_[0]


y = bigset[target_col]
X = bigset[data_cols]
model = LinearRegression()
model.fit(X, y)
r_sq = model.score(X, y)
#print (model.coef_,r_sq)
#print (X,y)
bigcof=model.coef_[0]


y = bigset[target_col]
x  = pd.DataFrame(data=bigset["-10"])
#print (x,y)
model = LinearRegression()
model.fit(x, y)
r_sq = model.score(x, y)
#print (model.coef_,r_sq)
#print (X,y)
#bigcof=model.coef_[0]
y2=model.predict(x)
fig, ax = plt.subplots(figsize=(20,10))
ax.plot(X["-10"],y2,label="8 days")
ax.scatter(X["-8"],y,label="8 days")
ax.scatter(X["-9"],y,label="9 days")
ax.scatter(X["-10"],y,label="10 days")
ax.scatter(X["-11"],y,label="11 days")
ax.set(xlabel="Cases X days before")
ax.set(ylabel="Deaths")
ax.set(Title="Correlation between confirmed and deaths")
ax.legend()
ax.set(xlim=(0, 5000), ylim=(0, 270))

fig = ax.get_figure()
fig.savefig(os.path.join(image_dir, 'Cases_vs_Death.png'))



cor=[]
bigcor=[]
for col in data_cols:
    x=dataset[col]
    target=dataset["Death"]
    #print (x,target)
    #print(col,np.corrcoef(x,target)[0,1] )
    cor+=[np.corrcoef(x,target)[0,1]]
    x=bigset[col]
    target=bigset["Death"]
    #print (x,target)
    #print(col,np.corrcoef(x,target)[0,1] )
    bigcor+=[np.corrcoef(x,target)[0,1]]

width=0.4
x=np.arange(-1*len(cof),0,1) # range does not work np.arange does
#print (x,cor,cof)

fig, ax = plt.subplots(figsize=(20,10))
plt.bar(x-width/2, cor, width=width, color='b', align='center',label="Correlation")
plt.bar(x+width/2, cof, width=width, color='r', align='center',label="LinRefCoef")
ax.set(xlabel="Days before")
ax.set(ylabel="R")
ax.set(Title="Correlation between confirmed and deaths")
ax.legend()
fig = ax.get_figure()
fig.savefig(os.path.join(image_dir, 'LinReg.png'))


fig, ax = plt.subplots(figsize=(20,10))
plt.bar(x-width/2, bigcor, width=width, color='b', align='center',label="Correlation")
plt.bar(x+width/2, bigcof, width=width, color='r', align='center',label="LinRefCoef")
ax.set(xlabel="Days before")
ax.set(ylabel="R")
ax.set(Title="Correlation between confirmed and deaths")
ax.legend()
fig = ax.get_figure()
fig.savefig(os.path.join(image_dir, 'LinRegBig.png'))

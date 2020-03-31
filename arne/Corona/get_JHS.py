#!/usr/bin/env python3
from __future__ import print_function
import pandas as pd
import numpy as np
import os
import sys, traceback
from pathlib import Path
import argparse
from argparse import RawTextHelpFormatter
import re
import glob
import math
import wget
import docopt
#import pickle
import os.path
import operator
from dominate import document
from dominate.tags import *

from dateutil.parser import parse
from datetime import datetime,date,time, timedelta
from dateutil import parser
#import pyarrow
#import matplotlib.pyplot as plt
#import matplotlib as mpl
# %matplotlib inline
from scipy.stats import linregress
from scipy.optimize import curve_fit

#mpl.rc('figure', max_open_warning = 0)

import preprocess as pp
import config as cf


##args = docopt.docopt(__doc__)
#out_dir = args['--output_folder']


#first we run data_pryp.py from covodify (slightly modifcied)


p = argparse.ArgumentParser(description =  '- get_ECDC.py - Extract data from ECDC and format it to merged.csv -',
            formatter_class=RawTextHelpFormatter)
p.add_argument('-f','--force', required= False, help='Force', action='store_true')
p.add_argument('-','--province', required= False, help='Force', action='store_true')
p.add_argument('-out','--output_folder', required= False, help='output folder')
ns = p.parse_args()

if ns.province:
    datafile="merged_province_JHS.csv"
else:
    datafile="merged_JHS.csv"

if ns.output_folder:
    out_dir = ns.output_folder
else:
    out_dir=home = str(Path.home())+"/Desktop/Corona/"



    
# Dynamic parameters
data_dir  = os.path.join(out_dir,'data') # os.path.join(out, 'data'  )+"/" # , str(datetime.date(datetime.now())))
#ECDC = "https://www.ecdc.europa.eu/sites/default/files/documents/"  # +2020-03-20+".xlsx
# import data
image_dir =  out_dir #os.path.join(out,'reports', 'images')
nation_dir =  os.path.join(out_dir,'nations')
reports_dir = os.path.join(out_dir,'data')

if not os.path.exists(image_dir):
    print('Creating reports folder...')
    os.system('mkdir -p ' + image_dir)
if not os.path.exists(nation_dir):
    print('Creating nation folder...')
    os.system('mkdir -p ' + nation_dir)
if not os.path.exists(data_dir):
    print('Creating reports folder...')
    os.system('mkdir -p ' + data_dir)
if not os.path.exists(reports_dir):
    print('Creating reports folder...')
    os.system('mkdir -p ' + reports_dir)

today=date.today()
yesterday=date.today() - timedelta(1)
aweekago=date.today() - timedelta(7)




# First we test if we already run the data for today:
try:
    merged_df=pd.read_csv(reports_dir+"/"+datafile, sep=',')
    date=merged_df['date'].max()
except:
    ns.force=True
    date=yesterday
if str(date)==str(today) and (not ns.force):
    print ("Exiting as todays plots are alredy run, use --force to rerun on yesterdays data")
    sys.exit(0)

# Get new data
os.system("python3 data_prep.py")
    
#excelfile="COVID-19-geographic-disbtribution-worldwide-"+str(today)+".xlsx"
#URL=ECDC+excelfile
infile=data_dir+"/"+str(today)+"/"+"agg_data_"+str(today)+".csv"

#print(infile)
if not os.path.isfile(infile):
    if (ns.force):
        infile=data_dir+"/"+str(yesterday)+"/"+"agg_data_"+str(yesterday)+".csv"
    else:
        print ("There is no new data yet")
        print(", use --force to rerun on yesterdays data")
        sys.exit(0)

# First we test if we already run the data for today:

    
print('Importing Data...')
print("Using: ",infile)
df=pd.read_csv(infile)
pp.fix_country_names(df)


merged_df = pd.DataFrame([])

agg_df=df.rename(columns={
    #"DateRep": "date",
    #"dateRep": "date",
    "date":"datestring",
    "Countries and territories":"country",
    "countriesAndTerritories":"country",
    #'Cases':"new_confirmed_cases",
    'cases':"confirmed_",
    #'confirmed':"new_confirmed_cases",
    #'recovered':"new_recovered_cases",
    'Deaths':"deaths",
    #'deaths':"deaths",
    })

# We need to turn date to right type (datetime.



#try:
agg_df['date']=agg_df.apply(lambda x:pp.FormatDateMerged(x.datestring), axis=1)
# We need to map american counties to states
pp.fix_states(agg_df)


del agg_df['datetime']
del agg_df['datestring']
del agg_df['file_date']

if ns.province:
    df["newcountry"] = df["country"] + df["province"]
    del agg_df['province']
    del agg_df['country']
    agg_df['country']=agg_df['newcountry']
    agg_df['province']=agg_df['newcountry']
    del agg_df['newcountry']
# We have to deal with missing data in only some provnces


dates=agg_df.groupby(['date'])['date'].first().dropna()

date=agg_df['date'].max()
if str(date.date())!=str(today) and (not ns.force):
    print ("Exiting as there is no new data for today, use --force to rerun on yesterdays data")
    sys.exit(0)

countries=agg_df['country'].drop_duplicates()
print ("Adding missing dates, province by province...")
for country in countries:
    provinces=agg_df.loc[agg_df['country']==country]["province"].drop_duplicates()
    for province in provinces:
        c=0
        r=0
        d=0
        print (country,province)
        for date in dates:
            if date=='' : continue
            if date>today:continue
            if agg_df.loc[ (agg_df['country']==country) & (agg_df['province']==province) & (agg_df['date']==date)].empty:
                agg_df=agg_df.append(
                {'date':date,
                    "country":country,
                    "province":province,
                    "confirmed":c,
                    "deaths":d,
                    "recovered":r},
                    ignore_index=True)
                #merged_df.append(data, ignore_index=True)
            else:
                c=int(agg_df.loc[ (agg_df['country']==country) & (agg_df['province']==province) & (agg_df['date']==pd.Timestamp(date))]['confirmed'].iloc[0])
                d=int(agg_df.loc[ (agg_df['country']==country) & (agg_df['province']==province) & (agg_df['date']==pd.Timestamp(date))]['deaths'].iloc[0])
                r=int(agg_df.loc[ (agg_df['country']==country) & (agg_df['province']==province) & (agg_df['date']==pd.Timestamp(date))]['recovered'].iloc[0])


    
del agg_df['province']
print (agg_df)


#sys.exit()

for country in countries:
    temp_df=agg_df.loc[agg_df['country'] == country].sort_values(by='date', ascending=True)
    country_df=pd.DataFrame([])
    date_df=pd.DataFrame([])
    country_df = temp_df.groupby(['date'])[['confirmed','deaths','recovered']].sum()
    date_df = temp_df.groupby(['date'])[['date']].first()
    country_df['new_deaths'] = country_df['deaths'].diff().fillna(0)
    country_df['new_confirmed_cases'] = country_df['confirmed'].diff().fillna(0)
    country_df['new_recovered_cases'] = country_df['recovered'].diff().fillna(0)
    country_df['country']=pd.Series(country,index=country_df.index)
    foo_df=pd.concat([date_df,country_df],axis=1)
    #print(foo_df)
    merged_df=pd.concat([merged_df, foo_df], ignore_index=True,sort=False)
# Remove a few dupliaed names
merged_df=merged_df.drop_duplicates(['country','date'], keep='last').dropna()
#print(merged_df)

# We shoudl complete the database with all missing dates.
first=merged_df["date"].to_list()[0]
firstdate=first
last=merged_df["date"].to_list()[-1]
lastdate=last



# We need to make two lists of countries
# Some parameters

countrylist={}
#sys.exit()
# This is all countries 
countries=merged_df['country'].drop_duplicates()
first={}
firstdate={}
startdate={}
firstdeaths={}
startdeaths={}
# Now we nedeathed to get the first date for each country (if <100 case last date)
for country in countries: # ["Afghanistan","Sweden","China"]: #countries:
    tempdf=merged_df.loc[merged_df['country'] == country]
    first[country]=tempdf["date"].to_list()[0]
    firstdate[country]=first[country]
    x=5
    try:
        start=tempdf[tempdf.confirmed > cf.minnum].iloc[0].date
    except:
        start=tempdf.date.tail(1).to_list()[0]
    try:
        deathsstart=tempdf[tempdf.deaths > cf.mindeaths].iloc[0].date
    except:
        deathsstart=tempdf.date.tail(1).to_list()[0]
    #if (type(start) is int):
    #    startdate[country]=parser.parse(start)
    #    startdeaths[country]=parser.parse(deathsstart)
    #else:
    startdate[country]=start
    startdeaths[country]=deathsstart

#print (startdeaths,startdate)
#tiny=0.000001

def Days(x,y):
    return (x-startdate[y]).days
def DeathsDays(x,y):
    return (x-startdeaths[y]).days

merged_df['Days']=merged_df.apply(lambda x:Days(x.date,x.country), axis=1)
merged_df['DeathsDays']=merged_df.apply(lambda x:DeathsDays(x.date,x.country), axis=1)

#print (merged_df)

dates=merged_df.groupby(['date'])['date'].first().dropna()

print ("Adding missing rows...")
for country in countries:
    c=0
    r=0
    d=0
    for date in dates:
        if date=='' : continue
        if date>today:continue
        #print ("TEST",date)
        if merged_df.loc[ (merged_df['country']==country) & (merged_df['date']==date)].empty:
            merged_df=merged_df.append(
            {'date':date,
                 "country":country,
                 #"new_confirmed_cases":0,
                 #"new_deaths":0,
                 #"new_recovered_cases":0,
                 "confirmed":c,
                 "deaths":d,
                 "recovered":r,
                      "Days":(date-startdate[country]).days,
                      "DeathsDays":(date-startdeaths[country]).days},
                                            ignore_index=True)
            #merged_df.append(data, ignore_index=True)
        else:
            c=int(merged_df.loc[ (merged_df['country']==country) & (merged_df['date']==pd.Timestamp(date))]['confirmed'])
            d=int(merged_df.loc[ (merged_df['country']==country) & (merged_df['date']==pd.Timestamp(date))]['deaths'])
            r=int(merged_df.loc[ (merged_df['country']==country) & (merged_df['date']==pd.Timestamp(date))]['recovered'])

#merged_df.to_csv(reports_dir+"/merged1.csv", sep=',')
merged_df['LogCases']=merged_df['confirmed'].apply(lambda x:(np.log2(max(x,cf.tiny))))
merged_df['LogDeaths']=merged_df['deaths'].apply(lambda x:(np.log2(max(x,cf.tiny))))
merged_df['Ratio'] = merged_df["deaths"]/merged_df["confirmed"]


linreg={}
countrylist={}
for country in countries:
    newdf=merged_df.loc[(merged_df['confirmed']>cf.minnum) & (merged_df['confirmed']<cf.maxnum) &(merged_df['country'] == country)]
    if (len(newdf)<4): # We need 6 points to fit a curve to the sigmoidal function.
        linreg[country]=linregress([0.0,1.0],[0.0,0.0])
        continue
    linreg[country]=linregress(newdf['Days'],newdf['LogCases'])
    countrylist[country]=linreg[country].slope

tmplist = sorted(countrylist.items() , reverse=True, key=lambda x: x[1])
sortedcountries=[]
for i in range(0,len(tmplist)):
    sortedcountries+=[tmplist[i][0]]

    # Sigmoidal (in log) funcion fit

#print (slopelist)

deathslist={}
deathsreg={}
for country in countries:
    newdf=merged_df.loc[(merged_df['deaths']>cf.mindeaths) & (merged_df['deaths']<cf.maxdeaths) & (merged_df['country'] == country)]
    if (len(newdf)<4):
        deathsreg[country]=linregress([0.0,1.0],[0.0,0.0])
        continue
    deathsreg[country]=linregress(newdf['DeathsDays'],newdf['LogDeaths'])
    #print(deathsreg[country])
    deathslist[country]=deathsreg[country].slope
tmplist = sorted(deathslist.items() , reverse=True, key=lambda x: x[1])
deathscountries=[]
for i in range(0,len(tmplist)):
    deathscountries+=[tmplist[i][0]]
    
    

def LinExp(x,y):
    return np.exp2(linreg[y].intercept+linreg[y].slope*x)

def DeathsExp(x,y):
    return np.exp2(deathsreg[y].intercept+deathsreg[y].slope*x)


merged_df['LinCases']=merged_df.apply(lambda x:LinExp(x.Days,x.country), axis=1)
merged_df['LinDeaths']=merged_df.apply(lambda x:DeathsExp(x.DeathsDays,x.country), axis=1)


#optimized2_ydata = my_func(xdata, ydata, est2_w, est2_k)

merged_df=merged_df.sort_values(by=['country', 'date'])
merged_df.to_csv(reports_dir+"/"+datafile, sep=',')

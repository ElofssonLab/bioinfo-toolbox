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
import matplotlib.pyplot as plt
import matplotlib as mpl
# %matplotlib inline
from scipy.stats import linregress
from scipy.optimize import curve_fit

mpl.rc('figure', max_open_warning = 0)


def fix_country_names(df):
    translations = {'United_States_of_America':'USA',
                    'United_Kingdom':'UK',
                    'Central_African_Republic':'CAR',
                    'United_Arab_Emirates':'UAE',
                    'United_Republic_of_Tanzania':'Tanzania',
                    'Democratic_Republic_of_the_Congo':'Congo',
                    'Others':'Cruise Ship',
                    'Iran (Islamic Republic of)':"Iran",
                    'Mainland China':'China',
                    'Korea, South':'South Korea',
                    'Republic of Korea':'South Korea',
                    'Cases_on_an_international_conveyance_Japan':'Cruise Ship'}
    df.replace(translations, inplace=True)



# Plot and save trendlinae graph
def nations_trend_line(tmp_df, name, cumconfirmed, cumdeath, ncases,ndeath,cdays,lincases,ddays,lindeaths):
    global curvefit,curvedeaths
    xdata=tmp_df[cdays].to_numpy()
    ydata=tmp_df[cumconfirmed].apply(lambda x:(np.log2(max(x,tiny)))).to_numpy()
    f, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]},figsize=(20,15))
    if name in curvefit.keys():
        yfit=sigmoidalfunction(xdata,ydata,curvefit[name][1],curvefit[name][2],curvefit[name][3],curvefit[name][4])
        ynum=np.exp2(yfit)
        ax1.plot(tmp_df['date'],ynum,label="sigmoidal curve fitting, max: "+str(round(np.exp2(curvefit[name][2]),0)))
    if name in curvedeath.keys():
        yfit=sigmoidalfunction0(xdata,ydata,curvedeath[name][1],curvedeath[name][2],curvedeath[name][3]) # ,curvedeath[name][4])
        ynum=np.exp2(yfit)
        ax1.plot(tmp_df['date'],ynum,label="sigmoidal death, max: "+str(round(np.exp2(curvedeath[name][2]),0)))
    tmp_df.groupby(['date'])[[cumconfirmed, cumdeath,lincases,lindeaths]].sum().plot(ax=ax1, marker='o')
    ax1.set_yscale('log')
    ax1.set(ylim=(0.5,maxcases))
    ax1.tick_params(axis='x', labelrotation=45 )
    days = tmp_df.groupby(['date'])[[cdays]].max()
    deathdays = tmp_df.groupby(['date'])[[ddays]].max()
    tmp = tmp_df.groupby(['date'])[[ncases]].sum()
    tmp2 = tmp_df.groupby(['date'])[[ndeath]].sum()
    tmp5 = tmp_df.groupby(['date'])[[cumconfirmed]].sum()
    tmp6 = tmp_df.groupby(['date'])[[cumdeath]].sum()
    ratio = tmp6[cumdeath]/tmp5[cumconfirmed]
    wratio =tmp6[cumdeath]/(tmp5[cumconfirmed].shift(7)+tiny)
    #print (ratio,wratio)
    if tmp[ncases].max()>0:
        ax1.bar(tmp.index,tmp[ncases], color="red",width=0.4, ls='dashed', lw=2)

        
    
    # We need to check if we have any deaths
    #print (tmp_df[ndeath])
    if tmp2[ndeath].max()>0:
        ax1.bar(tmp2.index,tmp2[ndeath], color="blue",width=0.8,alpha = 0.5, ls='dotted', lw=2  )
    #ax.set(xlabel="Days since > " + str(cutoff) + " cases")
    ax1.set(ylabel="Number of cases")
    ax1.set(Title="Covid-19 cases in " + name)
    #ax1.plot(tmp4.index, tmp4[lincases], 'red', label="Exponential curve fit cases")
    ax1.legend()
    
    fig = ax1.get_figure()
    x=tmp_df.groupby(['date'])['date'].first()
    ax2.bar(ratio.index,ratio,width=0.8, alpha=0.5,lw=2, color="blue")
    ax2.bar(wratio.index,wratio,width=0.4, color="red")
    ax2.tick_params(axis='x', labelrotation=45 )
    ax2.set(ylim=(0.,0.25))
    plt.xticks(rotation=45, ha='right')
    ax2.set(Title="Ratio of death (blue today red cmp with cases a week ago)")
    ax2.set(ylabel="Ratio of death")
    fig.savefig(os.path.join(nation_dir, name+'_trendline.png'.format(cumconfirmed)))
    plt.close(fig)
    plt.close('all')
    #sys.exit()

def create_trend_line(tmp_df, col, col2, col3):
    fig, ax = plt.subplots(figsize=(20,10))
    tmp_df.groupby(['date'])[[col, col2, col3]].sum().plot(ax=ax, marker='o')
    fig = ax.get_figure()
    fig.savefig(os.path.join(image_dir, '{}_trendline.png'.format(col)))

def create_bar(tmp_df, col, rgb):
    fig, ax = plt.subplots(figsize=(20,10))
    tmp = tmp_df.groupby(['date'])[[col]].sum()
    tmp.plot.bar(ax=ax, rot=45, color=rgb)
    fig = ax.get_figure()
    fig.savefig(os.path.join(image_dir, '{}_bar.png'.format(col)))
    
def create_stacked_bar(tmp_df, col1, col2, fig_title):
    tmp_df = tmp_df.set_index('date')
    fig, ax = plt.subplots(figsize=(20,10))
    tmp_df[[col2, col1]].plot.bar(ax=ax,
                                  rot=45,
                                  stacked=True,
                                  title=fig_title);
    fig = ax.get_figure()
    fig.savefig(os.path.join(image_dir, '{}_stacked_bar.png'.format(col2)))
    

    
# Set fonts etc.
font = {'weight' : 'bold',   'size'   : 22}
plt.rc('font', **font)

#set ggplot style
plt.style.use('ggplot')
markers = [ '.', ',', 'o', 'v', '^', '<', '>', '1', '2',
    '3', '4', '8', 's', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd', '|',
    '_', 'P', 'X' ,1,2,3,4,5,6,7,8,9]

#colours=['blue','green','red','cyan','magenta','yellow','black','grey','pink','brown','magenta']
colours=[]
for name, hex in mpl.colors.cnames.items():
    #print(name, hex)
    rc = hex[1:3]  # red      
    gc = hex[3:5]  # blue      
    bc = hex[5:7]  # green 
    if (int(rc,16)+int(gc,16)+int(bc,16)<512):
        colours+=[name]    
# Parameters for linreg
mindeaths=5
maxdeaths=1000
minnum=100
maxnum=5000

# Countries to select
minconfirmed=1000
mindeaths=10

# Cutoff to select startdate for expinential vurved
cutoff=500
# Plotting parameters
daysbefore=-10
daysafter=45
mincases=5
maxcases=250000
mindeathcases=1
maxdeathcases=25000
ddaysbefore=-5
ddaysafter=20
minslopedays=10
minslopeddays=-1

##args = docopt.docopt(__doc__)
#out_dir = args['--output_folder']


p = argparse.ArgumentParser(description = 
                                     '- AE-analysis.py - Extract data from date range and create models -',
            formatter_class=RawTextHelpFormatter)
p.add_argument('-f','--force', required= False, help='Force', action='store_true')
p.add_argument('-out','--output_folder', required= False, help='output folder')
ns = p.parse_args()

if ns.output_folder:
    out_dir = ns.output_folder
else:
    out_dir=home = str(Path.home())+"/Desktop/Corona/"

# Dynamic parameters
data_dir  = os.path.join(out_dir,'data') # os.path.join(out, 'data'  )+"/" # , str(datetime.date(datetime.now())))
ECDC = "https://www.ecdc.europa.eu/sites/default/files/documents/"  # +2020-03-20+".xlsx
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

#excelfile="COVID-19-geographic-disbtribution-worldwide-"+str(today)+".xlsx"
#URL=ECDC+excelfile
#infile=data_dir+"/"+excelfile
#
##print(infile)
#if not os.path.isfile(infile):
#    try:
#        #os.sys("wget -c " + URL + " -O " + infile )
#        excelfile=wget.download(URL, out=data_dir)
#    except:
#        if not ns.force:
#            print ("Exiting as there is no new data yet, use --force to rerun on yesterdays data")
#            sys.exit(0)
#        excelfile="COVID-19-geographic-disbtribution-worldwide-"+str(yesterday)+".xlsx"
#        infile=data_dir+"/"+excelfile
#        if not os.path.isfile(infile):
#            URL=ECDC+excelfile
#            excelfile=wget.download(URL, out=data_dir)
#else:
#    print ("Exiting as todays date is already run, use --force to rerun")
#    sys.exit(0)


# First we test if we already run the data for today:

merged_df=pd.read_csv(reports_dir+"/merged.csv", sep=',')
date=merged_df['date'].max()
if str(date)==str(today) and (not ns.force):
    print ("Exiting as todays plots are alredy run, use --force to rerun on yesterdays data")
    sys.exit(0)


infile=Path(data_dir+"/ECDC-data-"+str(today)+".csv")
if infile.is_file():
    print('Removing CVS file..')
    os.system('rm -f ' + str(infile))

URL="https://opendata.ecdc.europa.eu/covid19/casedistribution/csv"
csvfile=wget.download(URL, out=str(infile))
    
print('Importing Data...')
print("Using: ",infile)
df=pd.read_csv(csvfile,encoding='iso8859_16')
#df= pd.read_excel(infile)


merged_df = pd.DataFrame([])

agg_df=df.rename(columns={
    #"DateRep": "date",
    "dateRep": "DateRep",
    "Countries and territories":"country",
    "countriesAndTerritories":"country",
    'Cases':"new_confirmed_cases",
    'cases':"new_confirmed_cases",
    'Deaths':"new_deaths",
    'deaths':"new_deaths",
    })

# We need to turn date to right type (datetime.

def FormatDate(x):
    #return (parser.parse(x))
    return (datetime.strptime(x,"%d/%m/%Y"))

#try:
agg_df['date']=agg_df.apply(lambda x:FormatDate(x.DateRep), axis=1)


date=agg_df['date'].max()
if str(date.date())!=str(today) and (not ns.force):
    print ("Exiting as there is no new data for today, use --force to rerun on yesterdays data")
    sys.exit(0)

#sys.exit()
countries=agg_df['country'].drop_duplicates()
for country in countries:
    country_df=agg_df.loc[agg_df['country'] == country].sort_values(by=['year','month','day'], ascending=True)
    country_df['deaths'] = country_df['new_deaths'].cumsum()
    country_df['confirmed'] = country_df['new_confirmed_cases'].cumsum()
    merged_df=pd.concat([merged_df, country_df], ignore_index=True,sort=False)

# Remove a few dupliaed names
merged_df=merged_df.drop_duplicates(['country','date'], keep='last').dropna()


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
        start=tempdf[tempdf.confirmed > minnum].iloc[0].date
    except:
        start=tempdf.date.tail(1).to_list()[0]
    try:
        deathsstart=tempdf[tempdf.deaths > mindeaths].iloc[0].date
    except:
        deathsstart=tempdf.date.tail(1).to_list()[0]
    #if (type(start) is int):
    #    startdate[country]=parser.parse(start)
    #    startdeaths[country]=parser.parse(deathsstart)
    #else:
    startdate[country]=start
    startdeaths[country]=deathsstart

#print (startdeaths,startdate)
tiny=0.000001

def Days(x,y):
    return (x-startdate[y]).days
def DeathsDays(x,y):
    return (x-startdeaths[y]).days

merged_df['Days']=merged_df.apply(lambda x:Days(x.date,x.country), axis=1)
merged_df['DeathsDays']=merged_df.apply(lambda x:DeathsDays(x.date,x.country), axis=1)
merged_df['LogCases']=merged_df['confirmed'].apply(lambda x:(np.log2(max(x,tiny))))
merged_df['LogDeaths']=merged_df['deaths'].apply(lambda x:(np.log2(max(x,tiny))))
merged_df['Ratio'] = merged_df["deaths"]/merged_df["confirmed"]

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
                 "confirmed":c,
                 "deaths":d,
                 #"recovered":r,
                      "Days":(date-startdate[country]).days,
                      "DeathsDays":(date-startdeaths[country]).days,
                      #DeathsDate,
                      #StartDate,
                      #"FirstDate":first[country],
                      "LogCases":np.log2(c+tiny),
                      "LogDeaths":np.log2(d+tiny)},
                                            ignore_index=True)
            #merged_df.append(data, ignore_index=True)
        else:
            c=int(merged_df.loc[ (merged_df['country']==country) & (merged_df['date']==date)]['confirmed'])
            d=int(merged_df.loc[ (merged_df['country']==country) & (merged_df['date']==date)]['deaths'])
            #r=int(merged_df.loc[ (merged_df['country']==country) & (merged_df['date']==date)]['recovered'])

merged_df.to_csv(reports_dir+"/merged1.csv", sep=',')


linreg={}
countrylist={}
for country in countries:
    newdf=merged_df.loc[(merged_df['confirmed']>minnum) & (merged_df['confirmed']<maxnum) &(merged_df['country'] == country)]
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
def sigmoidalfunction(x,y,k,m,a,b):
    f=[]
    for i in x:
        f+=[(m-b) / (1 + np.exp(-k*(i-a))) + b ]
    return f
def sigmoidalfunction0(x,y,k,m,a):
    f=[]
    b=0
    for i in x:
        f+=[(m-b) / (1 + np.exp(-k*(i-a))) + b ]
    return f

print ("Plotting sigmoidal fit")

curvefit={}
totalcases={}
for country in countries:
    for test in [10,5,2,1,0]:
        newdf=merged_df.loc[(merged_df['country'] == country) & (merged_df['confirmed'] >test )]
        #print (country,mincases)
        xdata=newdf['Days'].to_numpy()
        ydata=newdf['LogCases'].to_numpy()
        try:
            curvefit[country], pcov = curve_fit(sigmoidalfunction,xdata,ydata) #,bounds=par_bounds)
            if curvefit[country][2]<22:
                totalcases[country]=np.exp2(curvefit[country][2])
                break
        except:
            continue
            #print ("error with ",country)
#sys.exit()
#print (totalcases)
#fig, ax = plt.subplots(figsize=(20,10))
sorted_td = sorted(totalcases.items(), key=operator.itemgetter(1),reverse=True)
x=[]
y=[]
#
curvedeath={}
totaldeaths={}
for country in countries:
    for test in [1,2,5,0,10,20,100]:
        newdf=merged_df.loc[(merged_df['country'] == country) & (merged_df['deaths'] >test )]
        #print (country,mincases)
        xdata=newdf['Days'].to_numpy()
        ydata=newdf['LogDeaths'].to_numpy()
        try:
            curvedeath[country], pcov = curve_fit(sigmoidalfunction0,xdata,ydata) #,bounds=par_bounds)
            if curvedeath[country][2]<18:
                totaldeaths[country]=np.exp2(curvedeath[country][2])
                break
        except:
            continue
            #print ("error with ",country)
#sys.exit()
#print (totaldeaths)
fig, ax = plt.subplots(figsize=(20,10))
sorted_d = sorted(totaldeaths.items(), key=operator.itemgetter(1),reverse=True)
x=[]
y=[]
z=[]
#print (sorted_td)
maxcountries=50
for i in range(0,min(maxcountries,len(sorted_td))):
    x+=[sorted_td[i][0]]
    z+=[sorted_td[i][1]]
    try:
        y+=[totaldeaths[sorted_d[i][0]]]
    except:
        y+=[0]
#print (x,y)
ax.bar(x,y,width=0.8,alpha=0.5,color="red")
ax.bar(x,z,width=0.4,color="blue")
ax.set_yscale('log')
ax.set(ylabel="Number of Cases/Deaths")
ax.set(title="Predicted total number of cases (blue) /deaths (red)")
plt.xticks(rotation=45, ha='right')
fig = ax.get_figure()
fig.savefig(os.path.join(image_dir, 'total_bar.png'))


# This is to get the change in slope
slopelist={}
newslopelist={}
for country in countries:
    ctoday=merged_df.loc[(merged_df['country']==country)]['Days'].max()
    if ctoday>7:
        slopelist[country]=[]
        newslopelist[country]=[]
    for days in range(7,ctoday):
        newdf=merged_df.loc[(merged_df['Days']>days-7) &(merged_df['Days']<=days) &(merged_df['country'] == country)]
        lr=linregress(newdf['Days'],newdf['LogCases'])
        slopelist[country]+=[lr.slope]
        dayone=merged_df.loc[(merged_df['Days']==(days-7)) &(merged_df['country'] == country)]['confirmed'].iloc[0]
        dayseven=merged_df.loc[(merged_df['Days']==(days)) &(merged_df['country'] == country)]['confirmed'].iloc[0]
        newslopelist[country]+=[np.power(float(dayseven)/float(dayone),(1/7))-1.0]

#print (slopelist,newslopelist)
#sys.exit()
# This is to get the change in slope
newdeathslopelist={}
deathslopelist={}
for country in countries:
    dtoday=merged_df.loc[(merged_df['country']==country)]['DeathsDays'].max()
    if dtoday>7:
        deathslopelist[country]=[]
        newdeathslopelist[country]=[]
    for days in range(7,dtoday):
        newdf=merged_df.loc[(merged_df['DeathsDays']>days-7) &(merged_df['DeathsDays']<=days) &(merged_df['country'] == country)]
        lr=linregress(newdf['DeathsDays'],newdf['LogCases'])
        dayone=merged_df.loc[(merged_df['DeathsDays']==(days-7)) &(merged_df['country'] == country)]['deaths'].iloc[0]
        dayseven=merged_df.loc[(merged_df['DeathsDays']==(days)) &(merged_df['country'] == country)]['deaths'].iloc[0]
        newdeathslopelist[country]+=[np.power(float(dayseven)/float(dayone),(1/7))-1.0]
        
        deathslopelist[country]+=[lr.slope]

#print (slopelist)

# This is to get the corrent (last week) linregression in cases
weeklist={}
weekreg={}
for country in countries:
    newdf=merged_df.loc[(merged_df['date']>pd.Timestamp(aweekago)) &(merged_df['country'] == country)]
    if (len(newdf)<4):
        weekreg[country]=linregress([0.0,1.0],[0.0,0.0])
        continue
    weekreg[country]=linregress(newdf['Days'],newdf['LogCases'])
    weeklist[country]=weekreg[country].slope
tmplist = sorted(weeklist.items() , reverse=True, key=lambda x: x[1])
#print (weekreg)

weekdeathslist={}
weekdeaths={}
for country in countries:
    newdf=merged_df.loc[(merged_df['date']>pd.Timestamp(pd.Timestamp(aweekago))) & (merged_df['country'] == country)]
    if (len(newdf)<4):
        weekdeaths[country]=linregress([0.0,1.0],[0.0,0.0])
        continue
    weekdeaths[country]=linregress(newdf['DeathsDays'],newdf['LogDeaths'])
    #print(deathsreg[country])
    weekdeathslist[country]=weekdeaths[country].slope
tmplist = sorted(weeklist.items() , reverse=True, key=lambda x: x[1])


deathslist={}
deathsreg={}
for country in countries:
    newdf=merged_df.loc[(merged_df['deaths']>mindeaths) & (merged_df['deaths']<maxdeaths) & (merged_df['country'] == country)]
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
    
    
linfit_df=merged_df.loc[(merged_df['confirmed']>minnum) & (merged_df['confirmed']<maxnum) & (merged_df['country'] != "China")]
deathsfit_df=merged_df.loc[(merged_df['deaths']>mindeaths) & (merged_df['deaths']<maxdeaths) & (merged_df['country'] != "China") ]    
newdf=linfit_df.groupby(['date']).sum()
linreg["RoW"]=linregress(newdf['Days'],newdf['LogCases'])
newdf=deathsfit_df.groupby(['date']).sum()
deathsreg["RoW"]=linregress(newdf['DeathsDays'],newdf['LogDeaths'])


def LinExp(x,y):
    return np.exp2(linreg[y].intercept+linreg[y].slope*x)

def DeathsExp(x,y):
    return np.exp2(deathsreg[y].intercept+deathsreg[y].slope*x)

merged_df['LinCases']=merged_df.apply(lambda x:LinExp(x.Days,x.country), axis=1)
merged_df['LinDeaths']=merged_df.apply(lambda x:DeathsExp(x.DeathsDays,x.country), axis=1)


#optimized2_ydata = my_func(xdata, ydata, est2_w, est2_k)

merged_df=merged_df.sort_values(by=['country', 'date'])
            
merged_df.to_csv(reports_dir+"/merged.csv", sep=',')

# ------------------------------------------------------------------------
##### Create Graphs #####
    
print('Creating graphs...')


print('... Time Series Trend Line')


fig, ax = plt.subplots(figsize=(20,10))
ax.set_yscale('log')
mark=0
col=0
minsize=10
for country in sortedcountries:
    tempdf=merged_df.loc[merged_df['country'] == country]
    first=tempdf["date"].to_list()[0]
    firstdate=first

    if tempdf["confirmed"].size<minsize: continue
    s=tempdf["confirmed"].max()
    try:
        start=tempdf[tempdf.confirmed > cutoff].iloc[0]
        startdate=start.date
    except:
        continue
    if s>cutoff:
        x=tempdf["date"].apply(lambda x:(x-startdate).days)
        y=tempdf['confirmed']
        ax.plot(x,y,label=country,marker=markers[mark])
        mark+=1
        if mark>=len(markers): mark=0
        col+=1
        if col>=len(colours): col=0
ax.set(xlim=(daysbefore, daysafter), ylim=(mincases, maxcases))
ax.set(xlabel="Days since " + str(mincases) + " cases")
ax.set(ylabel="Number of cases")
ax.set(Title="Covid-19 in countries")
ax.legend()
fig = ax.get_figure()
fig.savefig(os.path.join(image_dir, 'COUNTRIES_trendline.png'.format(col)))
#plt.show(block=False)
#time.sleep(5)
plt.close('all')



c=[]
r=[]
list={}
for country in countries:
    tempdf=merged_df.loc[merged_df['country'] == country]
    weekdf=tempdf.loc[(tempdf['date']<pd.Timestamp(pd.Timestamp(aweekago)))]

    y=tempdf[['confirmed','deaths']].max()
    q=weekdf[['confirmed','deaths']].max()
    if (y.deaths>1 or y.confirmed>100):
        c+=[country]
        r+=[y.deaths/y.confirmed]
        list[country]=[(tiny+y.deaths)/y.confirmed,y.deaths,y.confirmed,(tiny+y.deaths)/max(tiny,q.confirmed)]
fig, (ax1, ax2) = plt.subplots(2,1,gridspec_kw={'height_ratios': [1, 3]},figsize=(20,15))
tmplist=sorted(list.items() , reverse=True, key=lambda x: x[1])
x=[]
y=[]
z=[]
w=[]
v=[]
for i in range(0,len(tmplist)):
    x+=[tmplist[i][0]]
    y+=[float(tmplist[i][1][0])]
    z+=[float(tmplist[i][1][1])]
    w+=[float(tmplist[i][1][2])]
    v+=[float(tmplist[i][1][3])]
ax1.bar(np.arange(0,len(x))-0.2,z,width=0.4,color="green",label="Deaths")
ax1.bar(np.arange(0,len(x))+0.2,w,width=0.4,color="red",label="Cases")
plt.xticks(rotation=45, ha='right')
ax1.legend()
ax2.bar(x,y,color="blue",width=0.8,alpha=0.5,lw=2)
#ax2.bar(x,v,color="red",width=0.4,lw=2)  # This is actually confusing
#ax2.set(ylim=(0.,0.25))
plt.xticks(rotation=45, ha='right')
#ax.set(xlabel="Days since > " + str(cutoff) + " cases")
ax1.set(ylabel="Fraction of cases that are dead")
ax1.set(Title="Deaths Ratios in countries")
ax1.set_yscale('log')
fig = ax1.get_figure()
fig = ax2.get_figure()
fig.savefig(os.path.join(image_dir, 'ratio_bar.png'.format(col)))
plt.close('all')


# regression data for countr
#newdf=agg_df.loc[(agg_df['confirmed']>100) & (agg_df['confirmed']<10000)]
x=[]
y=[]
yerr=[]
z=[]
zerr=[]
mark=0
col=0
colorlist=[]

plt.close('all')
print('... Slope plots')
fig2, (ax2, ax3) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]},figsize=(20,15))
for country in sortedcountries:
    newdf=merged_df.loc[merged_df['country'] == country]
    if newdf['confirmed'].max()<minconfirmed: continue
    x+=[country]
    y+=[linreg[country].slope]
    yerr+=[linreg[country].stderr]
    z+=[weekreg[country].slope]
    zerr+=[weekreg[country].stderr]
    fig, ax = plt.subplots(figsize=(20,10))
    ax.set(ylabel="Log(Commulative cases)")
    ax.set(xlabel="Days from "+str(minnum)+" to "+str(maxnum) + " days")
    ax.set(Title=" Covid-19 log (cases) in different countries" )
    ax.scatter(newdf['Days'],newdf['confirmed'],label=country)
    #ax.plot(newdf['Days'], np.exp(linreg[country].intercept +
    #            linreg[country].slope*newdf['Days']), 'r',
    #            label=str(linreg[country]))
    ax.plot(newdf['Days'], newdf['LinCases'],label=str(linreg[country]))
    ax.set_yscale('log')
    ax.set(xlim=(daysbefore, daysafter), ylim=(mincases, maxcases))
    ax2.set(xlim=(daysbefore, daysafter), ylim=(mincases, maxcases))

    fig = ax.get_figure()
    ax.legend()

    ax2.plot(newdf['Days'],newdf['LinCases'] ,color=colours[col]) #, label='fitted line'+str(linreg[country])
    ax2.scatter(newdf['Days'],newdf['confirmed'],label=country,marker=markers[mark],color=colours[col])
    colorlist+=[colours[col]]
    mark+=1
    if mark>=len(markers): mark=0
    col+=1
    if col>=len(colours): col=0
        
    #plt.close('all')

ax2.set_yscale('log')
ax2.set(ylabel="Log(Commulative cases)")
ax2.set(Title=" Covid-19 log (cases) in different countries " )
fig2 = ax2.get_figure()
ax2.legend()
#plt.xticks(rotation=45, ha='right')
ax3.set(ylabel="Slope")
ax3.set(Title="Slope of Covid-19 log(cases) in different countries (last week in transparent)" )
#x=linreg.keys()
#y=linreg.slope
#yerr=linreg.stderr
ax3.bar(x,y,yerr=yerr,color=colorlist,width=0.4, ls='dashed', lw=2 )
ax3.bar(x,z,yerr=zerr,color=colorlist,width=0.8,alpha = 0.5, ls='dotted', lw=2 )
ax3.tick_params(axis='x', labelrotation=45 )
fig2 = ax3.get_figure()
plt.xticks(rotation=45, ha='right')
fig2.savefig(os.path.join(image_dir, 'slope.png'))
plt.close('all')


x=[]
y=[]
yerr=[]
z=[]
zerr=[]
mark=0
col=0
colorlist=[]
fig2, (ax2, ax3) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]},figsize=(20,15))
for country in deathscountries:
    newdf=merged_df.loc[merged_df['country'] == country]
    if newdf['deaths'].max()<mindeaths: continue
    x+=[country]
    y+=[deathsreg[country].slope]
    yerr+=[deathsreg[country].stderr]
    z+=[weekdeaths[country].slope]
    zerr+=[weekdeaths[country].stderr]
    ax2.plot(newdf['DeathsDays'],newdf['LinDeaths'],color=colours[col]) #, label='fitted line'+str(linreg[country])
    ax2.scatter(newdf['DeathsDays'],newdf['deaths'],label=country,marker=markers[mark],color=colours[col])
    colorlist+=[colours[col]]
    mark+=1
    if mark>=len(markers): mark=0
    col+=1
    if col>=len(colours): col=0
        
    #plt.close('all')
ax2.set(xlim=(ddaysbefore, ddaysafter),ylim=(mindeathcases, maxdeathcases))

ax2.set_yscale('log')
ax2.set(ylabel="Log(Commulative deaths)")
ax2.set(Title=" Covid-19 log (deaths) in different countries  (last week in transparent)" )
fig2 = ax2.get_figure()
ax2.legend()
#plt.xticks(rotation=45, ha='right')
ax3.set(ylabel="Slope")
ax3.set(Title="Slope of Covid-19 log(deaths) in different countries" )
ax3.bar(x,y,yerr=yerr,color=colorlist,width=0.4, ls='dashed', lw=2 )
ax3.bar(x,z,yerr=zerr,color=colorlist,width=0.8,alpha = 0.5, ls='dotted', lw=2 )
ax3.tick_params(axis='x', labelrotation=45 )
fig2 = ax3.get_figure()
plt.xticks(rotation=45, ha='right')
fig2.savefig(os.path.join(image_dir, 'deathslope.png'))
plt.close('all')

#sys.exit()

# Time Series Data Plots

print('... Slope figures')
colorlist=[]
markerlist=[]
col=0
mark=0
fig, ax = plt.subplots(figsize=(20,10))
for country in slopelist.keys():
    if (len(slopelist[country])>minslopedays):
        ax.plot(np.arange(len(slopelist[country])),slopelist[country],label=country,lw=2,marker=markers[mark],color=colours[col])
        colorlist+=[colours[col]]
        markerlist+=[markers[mark]]
        mark+=1
        if mark>=len(markers): mark=0
        col+=1
        if col>=len(colours): col=0
ax.legend() 
ax.set(title="Changes in slope from onset")
ax.set(ylabel="Slope (log2 base)")
ax.set(xlabel="Days since " + str(mincases) + " cases")
#ax.set(xlabel="Days")
ax.set(xlim=[0,daysafter])
#fig.show()
fig.savefig(os.path.join(image_dir, 'slope-evol.png'))

colorlist=[]
markerlist=[]
col=0
mark=0
fig, ax = plt.subplots(figsize=(20,10))
for country in newslopelist.keys():
    if (len(newslopelist[country])>minslopedays):
        ax.plot(np.arange(len(newslopelist[country])),newslopelist[country],label=country,lw=2,marker=markers[mark],color=colours[col])
        colorlist+=[colours[col]]
        markerlist+=[markers[mark]]
        mark+=1
        if mark>=len(markers): mark=0
        col+=1
        if col>=len(colours): col=0
ax.legend() 
ax.set(title="Changes in slope from onset")
ax.set(ylabel="Increase per day")
ax.set(xlabel="Days since " + str(mincases) + " cases")
#ax.set(xlabel="Days")
ax.set(xlim=[0,daysafter])
#fig.show()
fig.savefig(os.path.join(image_dir, 'newslope-evol.png'))


colorlist=[]
markerlist=[]
col=0
mark=0
fig, ax = plt.subplots(figsize=(20,10))
for country in deathslopelist.keys():
    if (len(deathslopelist[country])>minslopeddays):
        ax.plot(np.arange(len(deathslopelist[country])),deathslopelist[country],lw=3,label=country,marker=markers[mark],color=colours[col])
        colorlist+=[colours[col]]
        markerlist+=[markers[mark]]
        mark+=1
        if mark>=len(markers): mark=0
        col+=1
        if col>=len(colours): col=0
ax.legend() 
ax.set(title="Changes in death slope from onset")
ax.set(ylabel="Slope (log2 base)")
ax.set(xlabel="Days since " + str(mindeathcases) + " deaths")
#ax.set(xlabel="Days")
ax.set(xlim=[0,daysafter])
#fig.show()
fig.savefig(os.path.join(image_dir, 'deathslope-evol.png'))

colorlist=[]
markerlist=[]
col=0
mark=0
fig, ax = plt.subplots(figsize=(20,10))
for country in newdeathslopelist.keys():
    if (len(newdeathslopelist[country])>minslopeddays):
        ax.plot(np.arange(len(newdeathslopelist[country])),newdeathslopelist[country],lw=3,label=country,marker=markers[mark],color=colours[col])
        colorlist+=[colours[col]]
        markerlist+=[markers[mark]]
        mark+=1
        if mark>=len(markers): mark=0
        col+=1
        if col>=len(colours): col=0
ax.legend() 
ax.set(title="Changes in death slope from onset")
ax.set(ylabel="Increase per day")
ax.set(xlabel="Days since " + str(mindeathcases) + " deaths")
#ax.set(xlabel="Days")
ax.set(xlim=[0,daysafter])
#fig.show()
fig.savefig(os.path.join(image_dir, 'newdeathslope-evol.png'))


colorlist=[]
markerlist=[]
col=0
mark=0
fig, ax = plt.subplots(figsize=(20,10))
for country in countries:
    ctoday=merged_df.loc[(merged_df['country']==country)]['Days'].max()
    cases=merged_df.loc[(merged_df['country']==country)]['confirmed'].max()
    if cases<2000:
        continue
    X=[]
    Y=[]
    for day in range(7,ctoday):
        x=merged_df.loc[(merged_df['Days']==day) &(merged_df['country'] == country)]['confirmed'].iloc[0]
        y=x-merged_df.loc[(merged_df['Days']==day-7) &(merged_df['country'] == country)]['confirmed'].iloc[0]
        X+=[x]
        Y+=[y]
    ax.plot(X,Y,label=country,marker=markers[mark],color=colours[col])
    colorlist+=[colours[col]]
    markerlist+=[markers[mark]]
    mark+=1
    if mark>=len(markers): mark=0
    col+=1
    if col>=len(colours): col=0
x=[0,100000]
y=[0,50000]
ax.plot(x,y,label="Doubled each week",lw=4,color="black")
x=[0,100000]
y=[0,99000]
ax.plot(x,y,label="Doubled each day",lw=4,color="grey")
ax.legend() 

ax.set(title="Fraction of all infected last week from onset")
ax.set(ylabel="Increase from a week ago")
ax.set(xlabel="Total number of cases until today")
#fig.show()
fig.savefig(os.path.join(image_dir, 'weekly-increase.png'))
ax.set_xscale('log')
ax.set_yscale('log')
fig.savefig(os.path.join(image_dir, 'weekly-increase-log.png'))

colorlist=[]
markerlist=[]
col=0
mark=0
fig, ax = plt.subplots(figsize=(20,10))
for country in countries:
    ctoday=merged_df.loc[(merged_df['country']==country)]['Days'].max()
    cases=merged_df.loc[(merged_df['country']==country)]['deaths'].max()
    if cases<100:
        continue
    X=[]
    Y=[]
    for day in range(7,ctoday):
        x=merged_df.loc[(merged_df['Days']==day) &(merged_df['country'] == country)]['deaths'].iloc[0]
        y=x-merged_df.loc[(merged_df['Days']==day-7) &(merged_df['country'] == country)]['deaths'].iloc[0]
        X+=[x]
        Y+=[y]
    ax.plot(X,Y,label=country,marker=markers[mark],color=colours[col])
    colorlist+=[colours[col]]
    markerlist+=[markers[mark]]
    mark+=1
    if mark>=len(markers): mark=0
    col+=1
    if col>=len(colours): col=0
x=[0,25000]
y=[0,12500]
ax.plot(x,y,label="Doubled each week",lw=4,color="black")
x=[0,25000]
y=[0,24000]
ax.plot(x,y,label="Doubled each day",lw=4,color="grey")

ax.legend() 
ax.set(title="Fraction of all death occuring last week")
ax.set(ylabel="Increase from a week ago")
ax.set(xlabel="Total number of deaths until today")
#fig.show()
fig.savefig(os.path.join(image_dir, 'weekly-death.png'))
ax.set_xscale('log')
ax.set_yscale('log')
fig.savefig(os.path.join(image_dir, 'weekly-death-log.png'))



plt.close('All')
print('... Country Figures')
# Ratio plots
tempdf=merged_df.loc[merged_df['country'] != "China"]
nations_trend_line(tempdf, "RestOfWorld",  'confirmed', 'deaths', "new_confirmed_cases","new_deaths","Days","LinCases",'DeathsDays',"LinDeaths")

tempdf=merged_df.loc[merged_df['country'] == "China"]
nations_trend_line(tempdf, "China",  'confirmed', 'deaths', "new_confirmed_cases","new_deaths","Days","LinCases",'DeathsDays',"LinDeaths")

#sys.exit()

for country in countries:
    tempdf=merged_df.loc[merged_df['country'] == country]
    nations_trend_line(tempdf, country,  'confirmed', 'deaths', "new_confirmed_cases","new_deaths","Days","LinCases",'DeathsDays',"LinDeaths")

# Trend line for new cases
#create_trend_line(merged_df, 'new_confirmed_cases', 'new_deaths', 'new_recoveries')
    
    
#print('... Daily New Infections Differences')
#new_df = pd.DataFrame([])
#new_df['date'] = daily_df['date']
#new_df['confirmed_cases'] = merged_df.confirmed - daily_df.new_confirmed_cases
#new_df['new_confirmed_cases'] = daily_df.new_confirmed_cases
#create_stacked_bar(merged_df, 'new_confirmed_cases', 'confirmed_cases', "Stacked bar of confirmed and new cases by day")

#
print ('Makeing HTML file')

        
#topcountries=merged_df['country'].drop_duplicates()


selectedcountries="China|Italy|USA|UK|South_Korea|Germany|Sweden|Norway|Spain|Denmark|Finland"
overviewfiles="merged|slope|deaths|total"
allplots=[]
plots=[]
overview=[]
for f in os.listdir(out_dir):
    if (f.endswith(".png")):
        allplots +=[f]
    if (re.match(selectedcountries, f)):
        plots +=[f]
    if (re.match(overviewfiles, f)):
        overview +=[f]

#cvs = glob.glob(out_dir+'/*.csv')

with document(title='Arne Elofsson Corona') as doc:
    h1('Corona Data')
    for path in overview:
        div(img(src=path), _class='photo')

    h3('Selected countries')
    for path in plots:
        div(img(src=path), _class='photo')
            
    h3('All countries')
    for path in allplots:
        div(img(src=path), _class='photo')
    
with open(out_dir+'/corona.html', 'w') as f:
    f.write(doc.render())

print('Done!')

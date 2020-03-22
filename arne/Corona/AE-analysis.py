#!/usr/bin/env python3

"""
data_exploration.py - Extract data from date range and create models
Usage:
    data_exploration.py [options]
    data_exploration.py -h | --help

Options:
    -h --help             Show this message.
    --output_folder=OUT   Output folder for the data and reports to be saved
"""

from __future__ import print_function
import pandas as pd
import numpy as np
import os
import glob
import math
import wget
import docopt
#import pickle
import os.path
from dateutil.parser import parse
from datetime import datetime,date,time, timedelta
from dateutil import parser
#import pyarrow
import matplotlib.pyplot as plt
# %matplotlib inline
from scipy.stats import linregress

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
    f, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]},figsize=(20,15))
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
    if tmp[ncases].max()>0:
        ax1.bar(tmp.index,tmp[ncases], color="red",width=0.4)

        
    
    # We need to check if we have any deaths
    #print (tmp_df[ndeath])
    if tmp2[ndeath].max()>0:
        ax1.bar(tmp2.index,tmp2[ndeath], color="blue",width=0.9)
    #ax.set(xlabel="Days since > " + str(cutoff) + "cases")
    ax1.set(ylabel="Number of cases")
    ax1.set(Title="Covid-19 cases in " + name)
    #ax1.plot(tmp4.index, tmp4[lincases], 'red', label="Exponential curve fit cases")
    ax1.legend()
    
    fig = ax1.get_figure()
    x=tmp_df.groupby(['date'])['date'].first()
    ax2.bar(ratio.index,ratio,width=0.8, color="green")
    ax2.tick_params(axis='x', labelrotation=45 )
    ax2.set(ylim=(0.,0.1))
    plt.xticks(rotation=45, ha='right')
    ax2.set(ylabel="Ratio of death")
    fig.savefig(os.path.join(image_dir, name+'_trendline.png'.format(cumconfirmed)))
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

colours=['blue','green','red','cyan','magenta','yellow','black','grey','pink','brown']

# Parameters for linreg
mindeaths=5
maxdeaths=1000
minnum=50
maxnum=10000

# Countries to select
minconfirmed=1000

# Cutoff to select startdate for expinential vurved
cutoff=500
# Plotting parameters
daysbefore=-10
daysafter=45
mincases=5
maxcases=125000
mindeathcases=1
maxdeathcases=7500


args = docopt.docopt(__doc__)
out_dir = args['--output_folder']

# Dynamic parameters
data_dir  = out_dir # os.path.join(out, 'data'  )+"/" # , str(datetime.date(datetime.now())))
ECDC = "https://www.ecdc.europa.eu/sites/default/files/documents/"  # +2020-03-20+".xlsx
# import data
image_dir =  out_dir #os.path.join(out,'reports', 'images')
reports_dir = out_dir  #os.path.join(out,'reports')
if not os.path.exists(image_dir):
    print('Creating reports folder...')
    os.system('mkdir -p ' + image_dir)
if not os.path.exists(data_dir):
    print('Creating reports folder...')
    os.system('mkdir -p ' + data_dir)
if not os.path.exists(reports_dir):
    print('Creating reports folder...')
    os.system('mkdir -p ' + reports_dir)

today=date.today()
yesterday=date.today() - timedelta(1)
excelfile="COVID-19-geographic-disbtribution-worldwide-"+str(today)+".xlsx"
URL=ECDC+excelfile
infile=data_dir+"/"+excelfile

#print(infile)
if not os.path.isfile(infile):
    try:
        #os.sys("wget -c " + URL + " -O " + infile )
        excelfile=wget.download(URL, out=data_dir)
    except:
        excelfile="COVID-19-geographic-disbtribution-worldwide-"+str(yesterday)+".xlsx"
        infile=data_dir+"/"+excelfile
        if not os.path.isfile(infile):
            URL=ECDC+excelfile
            excelfile=wget.download(URL, out=data_dir)

print('Importing Data...')
print("Using: ",infile)
df= pd.read_excel(infile)

fix_country_names(df)
agg_df = pd.DataFrame([])

countries=df['Countries and territories'].drop_duplicates()
for country in countries:
    country_df=df.loc[df['Countries and territories'] == country].sort_values(by='DateRep')
    country_df['Cumulative deaths'] = country_df['Deaths'].cumsum()
    country_df['Cumulative cases'] = country_df['Cases'].cumsum()
    agg_df=pd.concat([agg_df, country_df], ignore_index=True)

merged_df=agg_df.rename(columns={
    "DateRep": "date",
    "Countries and territories":"country",
    'Cumulative cases':"confirmed",
    'Cases':"new_confirmed_cases",
    'Deaths':"new_deaths",
    'Cumulative deaths':"deaths"
    })
# Remove a few dupliaed names
merged_df=merged_df.drop_duplicates(['country','date'], keep='last')


# We shoudl complete the database with all missing dates.
first=merged_df["date"].to_list()[0]
firstdate=first
last=merged_df["date"].to_list()[-1]
lastdate=last



# We need to make two lists of countries
# Some parameters

countrylist={}
linreg={}
linregdeaths={}
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
    startdate[country]=start
    startdeaths[country]=deathsstart

#print (startdeaths,startdate)
tiny=0.001

def Days(x,y):
    return (x-startdate[y]).days
def DeathsDays(x,y):
    return (x-startdeaths[y]).days

merged_df['Days']=merged_df.apply(lambda x:Days(x.date,x.country), axis=1)
merged_df['DeathsDays']=merged_df.apply(lambda x:DeathsDays(x.date,x.country), axis=1)
merged_df['LogCases']=merged_df['confirmed'].apply(lambda x:(math.log(max(x,tiny))))
merged_df['LogDeaths']=merged_df['deaths'].apply(lambda x:(math.log(max(x,tiny))))
merged_df['Ratio'] = merged_df["deaths"]/merged_df["confirmed"]



dates=merged_df.groupby(['date'])['date'].first()
for country in countries:
    c=0
    r=0
    d=0
    for date in dates:
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
                      "LogCases":math.log(c+tiny),
                      "LogDeaths":math.log(d+tiny)},
                                            ignore_index=True)
            #merged_df.append(data, ignore_index=True)
        else:
            c=int(merged_df.loc[ (merged_df['country']==country) & (merged_df['date']==date)]['confirmed'])
            d=int(merged_df.loc[ (merged_df['country']==country) & (merged_df['date']==date)]['deaths'])
            #r=int(merged_df.loc[ (merged_df['country']==country) & (merged_df['date']==date)]['recovered'])

countrylist={}
for country in countries:
    newdf=merged_df.loc[(merged_df['confirmed']>minnum) & (merged_df['confirmed']<maxnum) &(merged_df['country'] == country)]
    if (len(newdf)<4):
        linreg[country]=linregress([0.0,1.0],[0.0,0.0])
        continue
    linreg[country]=linregress(newdf['Days'],newdf['LogCases'])
    countrylist[country]=linreg[country].slope
    
tmplist = sorted(countrylist.items() , reverse=True, key=lambda x: x[1])
sortedcountries=[]
for i in range(0,len(tmplist)):
    sortedcountries+=[tmplist[i][0]]


countrylist={}
deathsreg={}
for country in merged_df['country'].drop_duplicates():
    newdf=merged_df.loc[(merged_df['deaths']>mindeaths) & (merged_df['deaths']<maxdeaths) & (merged_df['country'] == country)]
    if (len(newdf)<4):
        deathsreg[country]=linregress([0.0,1.0],[0.0,0.0])
        continue
    deathsreg[country]=linregress(newdf['DeathsDays'],newdf['LogDeaths'])
    #print(deathsreg[country])
countrylist[country]=deathsreg[country].slope
tmplist = sorted(countrylist.items() , reverse=True, key=lambda x: x[1])
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
    return np.exp(linreg[y].intercept+linreg[y].slope*x)

def DeathsExp(x,y):
    return np.exp(deathsreg[y].intercept+deathsreg[y].slope*x)

merged_df['LinCases']=merged_df.apply(lambda x:LinExp(x.Days,x.country), axis=1)
merged_df['LinDeaths']=merged_df.apply(lambda x:DeathsExp(x.DeathsDays,x.country), axis=1)

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
ax.set(xlabel="Days since > " + str(cutoff) + "cases")
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
    y=tempdf[['confirmed','deaths']].max()
    if (y.deaths>1 or y.confirmed>100):
        c+=[country]
        r+=[y.deaths/y.confirmed]
        list[country]=[(0.0000000001+y.deaths)/y.confirmed,y.deaths,y.confirmed]


fig, (ax1, ax2) = plt.subplots(2,1,gridspec_kw={'height_ratios': [1, 3]},figsize=(20,15))
tmplist=sorted(list.items() , reverse=True, key=lambda x: x[1])
x=[]
y=[]
z=[]
w=[]
for i in range(0,len(tmplist)):
    x+=[tmplist[i][0]]
    y+=[float(tmplist[i][1][0])]
    z+=[float(tmplist[i][1][1])]
    w+=[float(tmplist[i][1][2])]
ax1.bar(np.arange(0,len(x))-0.2,z,width=0.4,color="green",label="Deaths")
ax1.bar(np.arange(0,len(x))+0.2,w,width=0.4,color="red",label="Cases")
plt.xticks(rotation=45, ha='right')
ax1.legend()
ax2.bar(x,y,color="green")
plt.xticks(rotation=45, ha='right')
#ax.set(xlabel="Days since > " + str(cutoff) + "cases")
ax1.set(ylabel="Fraction of cases that are dead (min 5 deaths)")
ax1.set(Title="Deaths Ration in countries")
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

mark=0
col=0
colorlist=[]

print('... Time Series Trend Line')
fig2, (ax2, ax3) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]},figsize=(20,15))
for country in sortedcountries:
    newdf=merged_df.loc[merged_df['country'] == country]
    if newdf['confirmed'].max()<minconfirmed: continue
    x+=[country]
    y+=[linreg[country].slope]
    yerr+=[linreg[country].stderr]
    fig, ax = plt.subplots(figsize=(20,10))
    ax.set(ylabel="Log(Commulative cases)")
    ax.set(xlabel="Days from "+str(minnum)+" to "+str(maxnum) + " days")
    ax.set(Title=" Covid-19 log (cases) in different countries" )
    ax.scatter(newdf['Days'],newdf['confirmed'],label=country)
    ax.plot(newdf['Days'], np.exp(linreg[country].intercept +
                linreg[country].slope*newdf['Days']), 'r',
                label=str(linreg[country]))
    ax.set_yscale('log')
    ax.set(xlim=(daysbefore, daysafter), ylim=(mincases, maxcases))
    ax2.set(xlim=(daysbefore, daysafter), ylim=(mincases, maxcases))

    fig = ax.get_figure()
    ax.legend()

    ax2.plot(newdf['Days'], np.exp(linreg[country].intercept +
                linreg[country].slope*newdf['Days']),color=colours[col]) #, label='fitted line'+str(linreg[country])
    ax2.scatter(newdf['Days'],newdf['confirmed'],label=country,marker=markers[mark],color=colours[col])
    colorlist+=[colours[col]]
    mark+=1
    if mark>=len(markers): mark=0
    col+=1
    if col>=len(colours): col=0
        
    #plt.close('all')

ax2.set_yscale('log')
ax2.set(ylabel="Log(Commulative cases)")
ax2.set(Title=" Covid-19 log (cases) in different countries" )
fig2 = ax2.get_figure()
ax2.legend()
#plt.xticks(rotation=45, ha='right')
ax3.set(ylabel="Slope")
ax3.set(Title="Slope of Covid-19 log(cases) in different countries" )
#x=linreg.keys()
#y=linreg.slope
#yerr=linreg.stderr
ax3.bar(x,y,yerr=yerr,color=colorlist)
ax3.tick_params(axis='x', labelrotation=45 )
fig2 = ax3.get_figure()
plt.xticks(rotation=45, ha='right')
fig2.savefig(os.path.join(image_dir, 'slope.png'))
plt.close('all')

x=[]
y=[]
yerr=[]
mark=0
col=0
colorlist=[]
fig2, (ax2, ax3) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]},figsize=(20,15))
mindeaths=5
for country in deathscountries:
    newdf=merged_df.loc[merged_df['country'] == country]
    if newdf['deaths'].max()<mindeaths: continue
    x+=[country]
    y+=[deathsreg[country].slope]
    yerr+=[deathsreg[country].stderr]
    ax2.plot(newdf['DeathsDays'],newdf['LinDeaths'],color=colours[col]) #, label='fitted line'+str(linreg[country])
    ax2.scatter(newdf['DeathsDays'],newdf['deaths'],label=country,marker=markers[mark],color=colours[col])
    colorlist+=[colours[col]]
    mark+=1
    if mark>=len(markers): mark=0
    col+=1
    if col>=len(colours): col=0
        
    #plt.close('all')

ax2.set_yscale('log')
ax2.set(ylabel="Log(Commulative deaths)")
ax2.set(Title=" Covid-19 log (deaths) in different countries" )
fig2 = ax2.get_figure()
ax2.legend()
#plt.xticks(rotation=45, ha='right')
ax3.set(ylabel="Slope")
ax3.set(Title="Slope of Covid-19 log(deaths) in different countries" )
ax3.bar(x,y,yerr=yerr,color=colorlist)
ax3.tick_params(axis='x', labelrotation=45 )
fig2 = ax3.get_figure()
plt.xticks(rotation=45, ha='right')
fig2.savefig(os.path.join(image_dir, 'deathslope.png'))
plt.close('all')

#sys.exit()

# Time Series Data Plots


print('... Daily Figures')

print('... Country Figures')
# Ratio plots
tempdf=merged_df.loc[merged_df['country'] != "China"]
nations_trend_line(tempdf, "RestOfWorld",  'confirmed', 'deaths', "new_confirmed_cases","new_deaths","Days","LinCases",'DeathsDays',"LinDeaths")

tempdf=merged_df.loc[merged_df['country'] == "China"]
nations_trend_line(tempdf, "CHINA",  'confirmed', 'deaths', "new_confirmed_cases","new_deaths","Days","LinCases",'DeathsDays',"LinDeaths")

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



print('Done!')

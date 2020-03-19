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
import docopt
import pickle
import os.path
from datetime import datetime
#from datetime import timedelta
from dateutil import parser
import pyarrow
import matplotlib.pyplot as plt
# %matplotlib inline
from scipy.stats import linregress

font = {'weight' : 'bold',
        'size'   : 22}
plt.rc('font', **font)


#set ggplot style
plt.style.use('ggplot')
 
args = docopt.docopt(__doc__)
out = args['--output_folder']

# Dynamic parameters
data_dir  = os.path.join(out, 'data', str(datetime.date(datetime.now())))
agg_file  = 'agg_data_{}.parquet.gzip'.format(datetime.date(datetime.now()))
trend_file  = 'trend_{}.csv'.format(datetime.date(datetime.now()))
report  = 'report_{}.xlsx'.format(datetime.date(datetime.now()))


# import data
print('Importing Data...')
agg_df = pd.read_parquet(os.path.join(data_dir, agg_file))
daily_df = pd.read_csv(os.path.join(data_dir, trend_file))


#Create place to save diagrams
image_dir =  os.path.join(out,'reports', 'images')
reports_dir =  os.path.join(out,'reports')

# Some remapping
#agg_df['newcountry']=agg_df['country'].apple x:re("United Kingdom"].#
agg_df.loc[agg_df['country'] == 'Others', 'country'] = "Cruise Ship"
agg_df.loc[agg_df['country'] == 'United Kingdom', 'country'] = "UK"
agg_df.loc[agg_df['country'] == 'Iran (Islamic Republic of)', 'country'] = "Iran"
agg_df.loc[agg_df['country'] == 'Mainland China', 'country'] = "China"
agg_df.loc[agg_df['country'] == 'Korea, South', 'country'] = "South Korea"
agg_df.loc[agg_df['country'] == 'Republic of Korea', 'country'] = "South Korea"



#print (agg_df)
# Merge to country..
# Convert types
for col in ['confirmed', 'deaths', 'recovered']:
    agg_df[col] = agg_df[col].replace('', 0).astype(int)

sum_df=agg_df.groupby(['date','country'])[['confirmed', 'deaths', 'recovered']].sum()
first_df=agg_df.groupby(['date','country'])['date','country'].first()


#print(sum_df,first_df)



merged=np.concatenate((first_df.to_numpy(),sum_df.to_numpy()),axis=1)
columns=['date','country','confirmed', 'deaths', 'recovered']
#merged_df = pd.DataFrame(data=merged,columns=columns),
merged_df = pd.DataFrame({'date': merged[:, 0], 'country': merged[:, 1], 'confirmed': merged[:, 2], 'deaths': merged[:, 3], 'recovered': merged[:, 4]})

# We shoudl complete the database with all missing dates.
first=merged_df["date"].to_list()[0]
firstdate=parser.parse(first)
last=merged_df["date"].to_list()[-1]
lastdate=parser.parse(last)

#print (firstdate,lastdate)
#



if not os.path.exists(image_dir):
    print('Creating reports folder...')
    os.system('mkdir -p ' + image_dir)

if not os.path.exists(reports_dir):
    print('Creating reports folder...')
    os.system('mkdir -p ' + reports_dir)

    # We need to make two lists of countries
markers = [ '.', ',', 'o', 'v', '^', '<', '>', '1', '2',
    '3', '4', '8', 's', 'p', '*', 'h', 'H', '+', 'x', 'D', 'd', '|',
    '_', 'P', 'X' ,1,2,3,4,5,6,7,8,9]

colours=['blue','green','red','cyan','magenta','yellow','black','grey','pink','brown']

# first we need to sort on slope
mindeaths=5
maxdeaths=1000
minnum=50
maxnum=10000
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
for country in countries:
    tempdf=merged_df.loc[merged_df['country'] == country]
    first[country]=tempdf["date"].to_list()[0]
    firstdate[country]=parser.parse(first[country])
    try:
        start=tempdf[tempdf.confirmed > minnum].iloc[0]
    except:
        start=tempdf.date.tail(1).to_list()[0]
    try:
        deathsstart=tempdf[tempdf.deaths > mindeaths].iloc[0]
    except:
        deathsstart=tempdf.date.tail(1).to_list()[0]
        #continue
    try:
        startdate[country]=parser.parse(start.date)
    except:
        startdate[country]=parser.parse(first[country])
    try:
        startdeaths[country]=parser.parse(deathsstart.date)
    except:
        startdeaths[country]=parser.parse(first[country])
    
tiny=0.0000000001
def Days(x,y):
    return (parser.parse(x)-startdate[y]).days
def DeathsDays(x,y):
    return (parser.parse(x)-startdeaths[y]).days

merged_df['Days']=merged_df.apply(lambda x:Days(x.date,x.country), axis=1)
merged_df['DeathsDays']=merged_df.apply(lambda x:DeathsDays(x.date,x.country), axis=1)
#merged_df['DeathsDate']=merged_df.apply(lambda x:startdeaths[x.country], axis=1)
#merged_df['StartDate']=merged_df.apply(lambda x:startdate[x.country], axis=1)
#merged_df['FirstDate']=merged_df.apply(lambda x:firstdate[x.country], axis=1)
merged_df['LogCases']=merged_df['confirmed'].apply(lambda x:(math.log(x+tiny)))
merged_df['LogDeaths']=merged_df['deaths'].apply(lambda x:(math.log(x+tiny)))



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
                 "recovered":r,
                      "Days":(parser.parse(date)-startdate[country]).days,
                      "DeathsDays":(parser.parse(date)-startdeaths[country]).days,
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
            r=int(merged_df.loc[ (merged_df['country']==country) & (merged_df['date']==date)]['recovered'])

countrylist={}
linreg={}
for country in countries:
    newdf=merged_df.loc[(merged_df['confirmed']>minnum) & (merged_df['confirmed']<maxnum) &(merged_df['country'] == country)]
    if (len(newdf)<4):
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
        continue
    deathsreg[country]=linregress(newdf['DeathsDays'],newdf['LogDeaths'])
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

merged_df=merged_df.sort_values(by=['country', 'date'])
            
merged_df.to_csv(reports_dir+"/merged.csv", sep=',')
#sys.exit()    



#
#sys.exit()


#merged_df['StartDate']=merged_df.loc[merged_df['country'] == country]["date"].apply(lambda x:startdate)
#merged_df['FirstDate']=merged_df.loc[merged_df['country'] == country]["date"].apply(lambda x:firstdate)
##### Define Graphs #####

# Plot and save trendlinae graph
def nations_trend_line(tmp_df, name, col, col2, col3,col4,col7,col5,slope,intercept,col6,deathsslope,deathsintercept):
    f, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]},figsize=(20,15))
    #print (tmp_df)
    #fig = plt.subplots()
    #ax=plt.subplot(2,1,1)
    tmp_df.groupby(['date'])[[col, col2]].sum().plot(ax=ax1, marker='o')
    ax1.set_yscale('log')
    ax1.set(ylim=(0.5,100000))
    ax1.tick_params(axis='x', labelrotation=45 )
    #ax1.set_xticklabels(labels=tmp_df.groupby(['date'])['date'], rotation=45 )
    #tmp_df.groupby(['date'])[[col4]].sum().plot.bar()
    days = tmp_df.groupby(['date'])[[col5]].max()
    deathdays = tmp_df.groupby(['date'])[[col6]].max()
    tmp = tmp_df.groupby(['date'])[[col4]].sum()
    tmp2 = tmp_df.groupby(['date'])[[col7]].sum()
    x=[]
    y=[]
    j=0
    for i in tmp[col4].keys():
        x+=[j]
        j+=1
        y+=[tmp[col4][i]]
    if tmp_df[col4].max()>0:
        ax1.bar(np.arange(0,len(x))-0.2,y,width=0.4, color="red")


    
    # We need to check if we have any deaths
    #print (tmp_df[col7])
    if tmp_df[col7].max()>0:
        x=[]
        y=[]
        j=0
        for i in tmp2[col7].keys():
            #print (i,j,tmp2[col7][i])
            x+=[j]
            j+=1
            y+=[tmp2[col7][i]]
        ax1.bar(np.arange(0,len(x))+0.2,y,width=0.4, color="blue")
        x=[]
        y=[]
        j=0
        if  (deathsslope>0):
            for i in deathdays[col6]:
                x+=[j]
                j+=1
                y+=[np.exp(deathsintercept+deathsslope*i)]
            #print (x,y,deathsintercept,deathsslope)
            ax1.plot(x, y, 'blue', label="Exponential curve fit death, Intercept:  ") # + str("%.5f" % deathsintercept) + " Slope: " + str("%.5" % deathsslope)   )


        
    #tiny=1
    
    
    #ax.set(xlabel="Days since > " + str(cutoff) + "cases")
    ax1.set(ylabel="Number of cases")
    ax1.set(Title="Covid-19 cases in " + name)
    x=[]
    y=[]
    j=0
    if (slope>0):
        for i in days[col5]:
            x+=[j]
            j+=1
            y+=[np.exp(intercept+slope*i)]
        ax1.plot(x, y, 'red', label="Exponential curve fit cases. Intercept: ") # + str("%.5f" % intercept) + " Slope: " + str("%.5f" % slope) )
        ax1.legend()
    
    fig = ax1.get_figure()

    #ax=plt.subplot(2,1,2)

    #ratio=[]
    #npa=tmp.to_numpy()
    #for i in range(1,len(npa-1)):
    #    if (npa[i][0]<tiny or npa[i-1][0]<tiny):
    #        ratio+=[1]
    #    else:
    #        ratio+=[max(npa[i][0],tiny)/max(npa[i-1][0],tiny)]
    #x=np.arange(1,len(npa))
    #ax2.set_yscale('log')
    #ax2.set(xlim=(-10, 25), ylim=(5, 750000))
    #ax2.set(ylim=(0., 0.1))
    x=tmp_df.groupby(['date'])['date'].first()
    #def Division(x,y):
    #    return (x/y)
    #y=tmp_df.apply(lambda x:Division(x.deaths,x.confirmed))
    j=0
    y=[]
    z=[]
    tmp = tmp_df.groupby(['date'])[[col]].sum()
    tmp2 = tmp_df.groupby(['date'])[[col2]].sum()
    for i in tmp[col].keys():
        y+=[tmp[col][i]]
    for i in tmp2[col2].keys():
        z+=[tmp2[col2][i]]
    for i in range(0,len(y)):
        y[i]=z[i]/(y[i]+tiny)
        
    ax2.bar(x,y,width=0.8, color="green")
    ax2.tick_params(axis='x', labelrotation=45 )
    ax2.set(ylim=(0.,0.1))
    plt.xticks(rotation=45, ha='right')
    #y=tmp_df.apply(lambda x:Division(x.deaths,x.confirmed))
    #ax2.bar(x,y, color="green",label="Ratio")
    ax2.set(ylabel="Ratio of death")
    fig.savefig(os.path.join(image_dir, name+'_trendline.png'.format(col)))
    plt.close()
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
    
    
##### Create Graphs #####
    
print('Creating graphs...')
print('... Time Series Trend Line')


cutoff=500
fig, ax = plt.subplots(figsize=(20,10))
ax.set_yscale('log')
mark=0
col=0
minsize=10
for country in sortedcountries:
    tempdf=merged_df.loc[merged_df['country'] == country]
    first=tempdf["date"].to_list()[0]
    firstdate=parser.parse(first)

    if tempdf["confirmed"].size<minsize: continue
    s=tempdf["confirmed"].max()
    try:
        start=tempdf[tempdf.confirmed > cutoff].iloc[0]
        startdate=parser.parse(start.date)
    except:
        continue
    if s>cutoff:
        x=tempdf["date"].apply(lambda x:(parser.parse(x)-startdate).days)
        y=tempdf['confirmed']
        #tempdf.groupby(['days'])[['confirmed']].sum().plot(ax=ax, marker='o')
        #ax.legend=([country])
        ax.plot(x,y,label=country,marker=markers[mark])
        mark+=1
        if mark>28: mark=0
        col+=1
        if col>10: col=0
ax.set(xlim=(-10, 25), ylim=(5, 75000))
ax.set(xlabel="Days since > " + str(cutoff) + "cases")
ax.set(ylabel="Number of cases")
ax.set(Title="Covid-19 in countries")
ax.legend()
fig = ax.get_figure()
fig.savefig(os.path.join(image_dir, 'COUNTRIES_trendline.png'.format(col)))
plt.close()


c=[]
r=[]
list={}
for country in countries:
    tempdf=merged_df.loc[merged_df['country'] == country]
    y=tempdf[['confirmed','deaths','recovered']].max()
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


# Data for predictions
mindays=20
#x=pd.DataFrame([])
file = open("daily-data.csv","w")  
string=''
for j in range(mindays,0,-1):
    string+= "-"+str(j)+","
string+="Deaths\n"
file.write(string)

for country in countries:
    tempdf=merged_df.loc[merged_df['country'] == country]
    if (len(tempdf) < mindays): continue
    tempdf["new_confirmed_cases"]=tempdf['confirmed'].diff()
    tempdf['new_deaths']=tempdf['deaths'].diff()
    tempdf['new_recoveries']=tempdf['recovered'].diff()
    tempdf["increased_confirmed_cases"]=tempdf['new_confirmed_cases'].diff()
    tempdf['increased_deaths']=tempdf['new_deaths'].diff()
    tempdf['increased_recoveries']=tempdf['new_recoveries'].diff()
    e=tempdf.increased_confirmed_cases.to_list()
    c=tempdf.new_confirmed_cases.to_list()
    d=tempdf.new_deaths.to_list()
    string=''
    for i in range(mindays+1,len(d)):
        #x=c[i-mindays:i-1],d[i])
        for j in range(i-mindays,i):
            string+= str(c[j])+","
        string+=str(d[i])+"\n"
    file.write(string)
    #x.append(y, ignore_index=True)
#x.to_csv("daily-data.csv")
file.close()
#sys.exit()

# regression data for countr
#newdf=agg_df.loc[(agg_df['confirmed']>100) & (agg_df['confirmed']<10000)]
x=[]
y=[]
yerr=[]

mark=0
col=0
colorlist=[]

    
minconfirmed=1000
fig2, (ax2, ax3) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]},figsize=(20,15))
for country in sortedcountries:
    newdf=merged_df.loc[merged_df['country'] == country]
    if newdf['confirmed']<minconfirmed: continue
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
    ax.set(xlim=(-10, 25), ylim=(5, 75000))
    ax2.set(xlim=(-10, 25), ylim=(5, 75000))
    #ax.set(ylim=(75, 80000))
    #ax2.set(ylim=(75, 80000))

    fig = ax.get_figure()
    ax.legend()
    #plt.xticks(rotation=45, ha='right')
    #fig.savefig(os.path.join(image_dir, country+'-slope.png'))

    ax2.plot(newdf['Days'], np.exp(linreg[country].intercept +
                linreg[country].slope*newdf['Days']),color=colours[col]) #, label='fitted line'+str(linreg[country])
    ax2.scatter(newdf['Days'],newdf['confirmed'],label=country,marker=markers[mark],color=colours[col])
    colorlist+=[colours[col]]
    mark+=1
    if mark>28: mark=0
    col+=1
    if col>9: col=0
        
    #plt.close()

ax2.set_yscale('log')
ax2.set(ylabel="Log(Commulative cases)")
ax2.set(Title=" Covid-19 log (cases) in different countries" )
fig2 = ax2.get_figure()
ax2.legend()
#plt.xticks(rotation=45, ha='right')
ax3.set(ylabel="Slope")
ax3.set(Title="Slope of Covid-19 log(cases) in differnt countries" )
#x=linreg.keys()
#y=linreg.slope
#yerr=linreg.stderr
ax3.bar(x,y,yerr=yerr,color=colorlist)
ax3.tick_params(axis='x', labelrotation=45 )
fig2 = ax3.get_figure()
plt.xticks(rotation=45, ha='right')
fig2.savefig(os.path.join(image_dir, 'slope.png'))

x=[]
y=[]
yerr=[]
mark=0
col=0
colorlist=[]
fig2, (ax2, ax3) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]},figsize=(20,15))
mindeaths=0
for country in deathscountries:
    newdf=merged_df.loc[merged_df['country'] == country]
    if newdf['deaths']<mindeaths: continue
    x+=[country]
    y+=[deathsreg[country].slope]
    yerr+=[deathsreg[country].stderr]
    fig, ax = plt.subplots(figsize=(20,10))
    ax.set(ylabel="Log(Commulative deaths)")
    ax.set(xlabel="Days from "+str(mindeaths)+" to "+str(maxdeaths) + " days")
    ax.set(Title=" Covid-19 log (deaths) in different countries" )
    ax.scatter(newdf['DeathsDays'],newdf['deaths'],label=country)
    ax.plot(newdf['DeathsDays'], np.exp(deathsreg[country].intercept +
                deathsreg[country].slope*newdf['DeathsDays']), 'r',
                label=str(deathsreg[country]))
    ax.set_yscale('log')
    ax.set(xlim=(-10, 25), ylim=(1, 5000))
    ax2.set(xlim=(-10, 25), ylim=(1, 5000))
    #ax.set(ylim=(75, 80000))
    #ax2.set(ylim=(75, 80000))

    fig = ax.get_figure()
    ax.legend()
    #plt.xticks(rotation=45, ha='right')
    #fig.savefig(os.path.join(image_dir, country+'-slope.png'))

    ax2.plot(newdf['DeathsDays'], np.exp(deathsreg[country].intercept +
                deathsreg[country].slope*newdf['DeathsDays']),color=colours[col]) #, label='fitted line'+str(linreg[country])
    ax2.scatter(newdf['DeathsDays'],newdf['deaths'],label=country,marker=markers[mark],color=colours[col])
    colorlist+=[colours[col]]
    mark+=1
    if mark>28: mark=0
    col+=1
    if col>9: col=0
        
    #plt.close()

ax2.set_yscale('log')
ax2.set(ylabel="Log(Commulative deaths)")
ax2.set(Title=" Covid-19 log (deaths) in different countries" )
fig2 = ax2.get_figure()
ax2.legend()
#plt.xticks(rotation=45, ha='right')
ax3.set(ylabel="Slope")
ax3.set(Title="Slope of Covid-19 log(deaths) in differnt countries" )
#x=linreg.keys()
#y=linreg.slope
#yerr=linreg.stderr
ax3.bar(x,y,yerr=yerr,color=colorlist)
ax3.tick_params(axis='x', labelrotation=45 )
fig2 = ax3.get_figure()
plt.xticks(rotation=45, ha='right')
fig2.savefig(os.path.join(image_dir, 'deathslope.png'))

#sys.exit()

# Time Series Data Plots


print('... Daily Figures')
# Daily Figures Data Plots
daily_figures_cols = ['new_confirmed_cases', 'new_deaths', 'new_recoveries', 'currently_infected']
for col, rgb in zip(daily_figures_cols, ['tomato', 'lightblue', 'mediumpurple', 'green']):
    create_bar(daily_df, col, rgb)    

print('... Country Figures')
# Ratio plots
tempdf=merged_df.loc[merged_df['country'] != "China"]
#print (country,y.diff())
tempdf["new_confirmed_cases"]=tempdf['confirmed'].diff()
tempdf['new_deaths']=tempdf['deaths'].diff()
tempdf['new_recoveries']=tempdf['recovered'].diff()
tempdf["increased_confirmed_cases"]=tempdf['new_confirmed_cases'].diff()
tempdf['increased_deaths']=tempdf['new_deaths'].diff()
tempdf['increased_recoveries']=tempdf['new_recoveries'].diff()
nations_trend_line(tempdf, "RestOfWorld",  'confirmed', 'deaths', 'recovered',"new_confirmed_cases","new_deaths","Days",linreg["RoW"].slope,linreg["RoW"].intercept,'DeathsDays',deathsreg["RoW"].slope,deathsreg["RoW"].intercept)

for country in countries:
    tempdf=merged_df.loc[merged_df['country'] == country]
    tempdf["new_confirmed_cases"]=tempdf['confirmed'].diff()
    tempdf['new_deaths']=tempdf['deaths'].diff()
    tempdf['new_recoveries']=tempdf['recovered'].diff()
    tempdf["increased_confirmed_cases"]=tempdf['new_confirmed_cases'].diff()
    tempdf['increased_deaths']=tempdf['new_deaths'].diff()
    tempdf['increased_recoveries']=tempdf['new_recoveries'].diff()
    try:
        print(country,linreg[country])
        nations_trend_line(tempdf, country,  'confirmed', 'deaths', 'recovered',"new_confirmed_cases","new_deaths","Days",linreg[country].slope,linreg[country].intercept,'DeathsDays',deathsreg[country].slope,deathsreg[country].intercept)
    except:
        try:
            nations_trend_line(tempdf, country,  'confirmed', 'deaths', 'recovered',"new_confirmed_cases","new_deaths","Days",linreg[country].slope,linreg[country].intercept,'DeathsDays',0.,0.)
        except:
            nations_trend_line(tempdf, country,  'confirmed', 'deaths', 'recovered',"new_confirmed_cases","new_deaths","Days",0.,0.,'DeathsDays',0.,0.)

        
# Trend line for new cases
create_trend_line(daily_df, 'new_confirmed_cases', 'new_deaths', 'new_recoveries')
    
    
#print('... Daily New Infections Differences')
new_df = pd.DataFrame([])
new_df['date'] = daily_df['date']
new_df['confirmed_cases'] = merged_df.confirmed - daily_df.new_confirmed_cases
new_df['new_confirmed_cases'] = daily_df.new_confirmed_cases
create_stacked_bar(new_df, 'new_confirmed_cases', 'confirmed_cases', "Stacked bar of confirmed and new cases by day")



#print('Creating excel spreadsheet report...')
workbook_writer = pd.ExcelWriter(os.path.join(reports_dir, report), engine='xlsxwriter')

# Add daily summary to spreadsheet
daily_df.to_excel(workbook_writer, sheet_name='daily figures')  


workbook = workbook_writer.book

def get_image_types(path):
    # get all the possible types of images in
    # the passed directory path
    types = []
    for fn in glob.glob(os.path.join(path, '*.png')):
        types.append(fn.split('_',)[-1].split('.')[0])
    
    return types

# Get all images for each type
def read_images(path, graph_type):
    image_list = []
    for fn in glob.glob(os.path.join(path, '*_{}.png'.format(graph_type))):
        image_list.append(fn)    
    images = {graph_type : image_list}
    return dict(images)

image_types = get_image_types(image_dir)

padding = 1 # Set padding for images in spreadsheet
for types in set(image_types):
    print('... reading images for:', types)
    type_dict = read_images(image_dir, types)
    
    # Add image to the worksheet
    worksheet = workbook.add_worksheet(name='{}_graphs'.format(types))
    for image in type_dict[types]:
        worksheet.insert_image('A' +str(padding), image) 
        padding += 50
    padding = 1
    
workbook.close()

print('Done!')

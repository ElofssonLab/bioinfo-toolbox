
# Parameters for linreg
mindeaths=5
maxdeaths=1000
minnum=100
maxnum=5000

# Countries to select (we try to include sweden)
minconfirmed=3000
mindeaths=50

# Selected countries
specialcountries=["China","USA","UK","Germany","America","Spain","Italy","Sweden","France","Denmark","Norway"]

# Cutoff to select startdate for exponential curves
cutoff=500

# Plotting parameters
daysbefore=-10
daysafter=75
mincases=5
maxcases=250000
mindeathcases=1
maxdeathcases=25000
ddaysbefore=-5
ddaysafter=20
minslopedays=10
minslopeddays=-1

maxcountries=20

# Misc
tiny=0.000001


import os

#
# CLI 
#
SCRIPT = '/pipeline.sh'
LIST_SCRIPT = '/pipeline.sh'


#
# DATA PREP
#
REPO = 'https://github.com/CSSEGISandData/COVID-19.git'
TMP_FOLDER = '/tmp/corona/'
TMP_GIT = os.path.join(TMP_FOLDER, REPO.split('/')[-1].split('.')[0])
DATA = os.path.join(TMP_GIT, 'csse_covid_19_data', 'csse_covid_19_daily_reports')

#Github cols
KEEP_COLS = ['country',
             'province', 
             'confirmed',
             'deaths',
             'datetime',
             'file_date',
             #'new_recovered_cases',
             #'new_confirmed_cases',
             #'new_deaths',
             'recovered',
             'date']

NUMERIC_COLS = ['confirmed', 
                'deaths', 
             #'new_recovered_cases',
             #'new_confirmed_cases',
             #'new_deaths'
                'recovered'
                    ]

import numpy as np
#from dateutil.parser import parse
from datetime import datetime # ,date,time, timedelta
#from dateutil import parser

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
                    'Cases_on_an_international_conveyance_Japan':'Cruise Ship',
                    'Mainland China':'China',
                    'Korea, South':'South Korea',
                    'Republic of Korea':'South Korea',
                    'Hong Kong SAR':'Hong Kong',
                    'Taipei and environs':'Taiwan',
                    'Taiwan*':'Taiwan',
                    'Macao SAR':'Macau',
                        'Iran (Islamic Republic of)':'Iran',
                        'Viet Nam':'Vietnam',
                        'UK':'United Kingdom',
                        ' Azerbaijan':'Azerbaijan',
                        'Bosnia and Herzegovina':'Bosnia',
                        'Czech Republic':'Czechia',
                        'Republic of Ireland':'Ireland',
                        'North Ireland':'Ireland',
                        'Republic of Moldova':'Moldova',
                        'Congo (Brazzaville:)':'Congo',
                        'Congo (Kinshasa)':'Congo',
                        'Republic of the Congo':'Congo',
                        'Gambia, The':'Gambia',
                        'The Gambia':'Gambia',
                        'USA':'Unites States',
                        'US':'United States',
                        'Bahamas, The':'The Bahamas',
                        'Bahamas':'The Bahamas',
                        'Cruise Ship':'Others'
                        }
    df.replace(translations, inplace=True)



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

def FormatDate(x):
    #return (parser.parse(x))
    return (datetime.strptime(x,"%d/%m/%Y"))

def FormatDateMerged(x):
    #return (parser.parse(x))
    return (datetime.strptime(x,"%Y-%m-%d"))

def replace_arg_space(country_str):
    return country_str.replace(' ', '_')

def replace_arg_score(country_str):
    return country_str.replace('_', ' ')

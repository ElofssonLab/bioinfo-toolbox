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
                        'Unites States':"USA",
                        'US':'USA',
                        "United Kindom":'UK',
                        'Bahamas, The':'The Bahamas',
                        'Bahamas':'The Bahamas',
                        'Cruise Ship':'Others'
                        }
    df.replace(translations, inplace=True)


def fix_states(df):
    df.replace(['.*AL'],['Alabama'],regex=True,inplace=True)
    df.replace(['.*AK'],['Alaska'],regex=True,inplace=True)
    df.replace(['.*AZ'],['Arizona'],regex=True,inplace=True)
    df.replace(['.*AR'],['Arkansas'],regex=True,inplace=True)
    df.replace(['.*CA'],['California'],regex=True,inplace=True)
    df.replace(['.*CO'],['Colorado'],regex=True,inplace=True)
    df.replace(['.*CT'],['Connecticut'],regex=True,inplace=True)
    df.replace(['.*DE'],['Delaware'],regex=True,inplace=True)
    df.replace(['.*DC'],['District of Columbia'],regex=True,inplace=True)
    df.replace(['.*FL'],['Florida'],regex=True,inplace=True)
    df.replace(['.*GA'],['Georgia'],regex=True,inplace=True)
    df.replace(['.*HI'],['Hawaii'],regex=True,inplace=True)
    df.replace(['.*ID'],['Idaho'],regex=True,inplace=True)
    df.replace(['.*IL'],['Illinois'],regex=True,inplace=True)
    df.replace(['.*IN'],['Indiana'],regex=True,inplace=True)
    df.replace(['.*IA'],['Iowa'],regex=True,inplace=True)
    df.replace(['.*KS'],['Kansas'],regex=True,inplace=True)
    df.replace(['.*KY'],['Kentucky'],regex=True,inplace=True)
    df.replace(['.*LA'],['Louisiana'],regex=True,inplace=True)
    df.replace(['.*ME'],['Maine'],regex=True,inplace=True)
    df.replace(['.*MD'],['Maryland'],regex=True,inplace=True)
    df.replace(['.*MA'],['Massachusetts'],regex=True,inplace=True)
    df.replace(['.*MI'],['Michigan'],regex=True,inplace=True)
    df.replace(['.*MN'],['Minnesota'],regex=True,inplace=True)
    df.replace(['.*MS'],['Mississippi'],regex=True,inplace=True)
    df.replace(['.*MO'],['Missouri'],regex=True,inplace=True)
    df.replace(['.*MT'],['Montana'],regex=True,inplace=True)
    df.replace(['.*NE'],['Nebraska'],regex=True,inplace=True)
    df.replace(['.*NV'],['Nevada'],regex=True,inplace=True)
    df.replace(['.*NH'],['New Hampshire'],regex=True,inplace=True)	
    df.replace(['.*NJ'],['New Jersey'],regex=True,inplace=True)	
    df.replace(['.*NM'],['New Mexico'],regex=True,inplace=True)	
    df.replace(['.*NY'],['New York'],regex=True,inplace=True)	
    df.replace(['.*NC'],['North Carolina'],regex=True,inplace=True)	
    df.replace(['.*ND'],['North Dakota'],regex=True,inplace=True)	
    df.replace(['.*OH'],['Ohio'],regex=True,inplace=True)
    df.replace(['.*OK'],['Oklahoma'],regex=True,inplace=True)
    df.replace(['.*OR'],['Oregon'],regex=True,inplace=True)
    df.replace(['.*PA'],['Pennsylvania'],regex=True,inplace=True)
    df.replace(['.*RI'],['Rhode Island'],regex=True,inplace=True)	
    df.replace(['.*SC'],['South Carolina'],regex=True,inplace=True)	
    df.replace(['.*SD'],['South Dakota'],regex=True,inplace=True)	
    df.replace(['.*TN'],['Tennessee'],regex=True,inplace=True)
    df.replace(['.*TX'],['Texas'],regex=True,inplace=True)
    df.replace(['.*UT'],['Utah'],regex=True,inplace=True)
    df.replace(['.*VT'],['Vermont'],regex=True,inplace=True)
    df.replace(['.*VA'],['Virginia'],regex=True,inplace=True)
    df.replace(['.*WA'],['Washington'],regex=True,inplace=True)
    df.replace(['.*WV'],['West Virginia'],regex=True,inplace=True)
    df.replace(['.*WI'],['Wisconsin'],regex=True,inplace=True)
    df.replace(['.*WY'],['Wyoming'],regex=True,inplace=True)
    df.replace(['.*WY'],['Wyoming'], regex=True, inplace=True)
    df.replace(['.*Diamond*'],['Cruise Ship'], regex=True, inplace=True)
    df.replace(['Virgin Islands, U.S.'],['Virgin Islands'], regex=True, inplace=True)
    
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

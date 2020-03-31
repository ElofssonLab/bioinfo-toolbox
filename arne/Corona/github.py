#!/usr/bin/env python3
from __future__ import print_function
import pandas as pd
import re
import os
import sys
import git
import numpy as np
from tqdm import tqdm
from time import strftime
from dateutil.parser import parse
from datetime import datetime, date, time 
from config import REPO, TMP_FOLDER, TMP_GIT, DATA, KEEP_COLS, NUMERIC_COLS

def clean_sheet_names(new_ranges):
    # Remove all sheets that dont have a numeric header
    return [x for x in new_ranges if re.search(r'\d', x)]

def clone_repo(TMP_FOLDER, REPO):
    print('Cloning Data Repo...')
    git.Git(TMP_FOLDER).clone(REPO)

def get_date(last_update):
    return parse(str(last_update).split(' ')[0]).strftime("%Y-%m-%d")

def get_csv_date(f):
    return get_date(f.split('.')[0] + ' ')


def fix_country_names(tmp_df):
    '''
    Cleaning up after JHU's bullshit data management
    '''
    # Asian Countries
    tmp_df['country'] = np.where((tmp_df['country']  == 'Mainland China'),'China', tmp_df['country'])
    tmp_df['country'] = np.where((tmp_df['country']  == 'Korea, South'),'South Korea', tmp_df['country'])
    tmp_df['country'] = np.where((tmp_df['country']  == 'Republic of Korea'),'South Korea', tmp_df['country'])
    tmp_df['country'] = np.where((tmp_df['country']  == 'Hong Kong SAR'),'Hong Kong', tmp_df['country'])
    tmp_df['country'] = np.where((tmp_df['country']  == 'Taipei and environs'),'Taiwan', tmp_df['country'])
    tmp_df['country'] = np.where((tmp_df['country']  == 'Taiwan*'),'Taiwan', tmp_df['country'])
    tmp_df['country'] = np.where((tmp_df['country']  == 'Macao SAR'),'Macau', tmp_df['country'])
    tmp_df['country'] = np.where((tmp_df['country']  == 'Iran (Islamic Republic of)'),'Iran', tmp_df['country'])
    tmp_df['country'] = np.where((tmp_df['country']  == 'Viet Nam'),'Vietnam', tmp_df['country'])

    #European Countries
    tmp_df['country'] = np.where((tmp_df['country']  == 'UK'),'United Kingdom', tmp_df['country'])
    tmp_df['country'] = np.where((tmp_df['country']  == ' Azerbaijan'),'Azerbaijan', tmp_df['country'])
    tmp_df['country'] = np.where((tmp_df['country']  == 'Bosnia and Herzegovina'),'Bosnia', tmp_df['country'])
    tmp_df['country'] = np.where((tmp_df['country']  == 'Czech Republic'),'Czechia', tmp_df['country'])
    tmp_df['country'] = np.where((tmp_df['country']  == 'Republic of Ireland'),'Ireland', tmp_df['country'])
    tmp_df['country'] = np.where((tmp_df['country']  == 'North Ireland'),'Ireland', tmp_df['country'])
    tmp_df['country'] = np.where((tmp_df['country']  == 'Republic of Moldova'),'Moldova', tmp_df['country'])

    #African Countries
    tmp_df['country'] = np.where((tmp_df['country']  == 'Congo (Brazzaville)'),'Congo', tmp_df['country'])
    tmp_df['country'] = np.where((tmp_df['country']  == 'Congo (Kinshasa)'),'Congo', tmp_df['country'])
    tmp_df['country'] = np.where((tmp_df['country']  == 'Republic of the Congo'),'Congo', tmp_df['country'])
    tmp_df['country'] = np.where((tmp_df['country']  == 'Gambia, The'),'Gambia', tmp_df['country'])
    tmp_df['country'] = np.where((tmp_df['country']  == 'The Gambia'),'Gambia', tmp_df['country'])

    # Western Countries
    tmp_df['country'] = np.where((tmp_df['country']  == 'USA'),'America', tmp_df['country'])
    tmp_df['country'] = np.where((tmp_df['country']  == 'US'),'America', tmp_df['country'])
    tmp_df['country'] = np.where((tmp_df['country']  == 'Bahamas, The'),'The Bahamas', tmp_df['country'])
    tmp_df['country'] = np.where((tmp_df['country']  == 'Bahamas'),'The Bahamas', tmp_df['country'])

    # Others
    tmp_df['country'] = np.where((tmp_df['country']  == 'Cruise Ship'),'Others', tmp_df['country'])

    return tmp_df

# Now that we have all the data we now need to clean it 
# - Fill null values
# - remore suspected values
# - change column names
def clean_data(df):
    tmp_df = df.copy()

    if 'Demised' in tmp_df.columns:
        tmp_df.rename(columns={'Demised':'deaths'}, inplace=True)

    if 'Country/Region' in tmp_df.columns:
        tmp_df.rename(columns={'Country/Region':'country'}, inplace=True)

    if 'Country_Region' in tmp_df.columns:
        tmp_df.rename(columns={'Country_Region':'country'}, inplace=True)
    
    if 'Province/State' in tmp_df.columns:
        tmp_df.rename(columns={'Province/State':'province'}, inplace=True)

    if 'Province_State' in tmp_df.columns:
        tmp_df.rename(columns={'Province_State':'province'}, inplace=True)

    if 'Last Update' in tmp_df.columns:
        tmp_df.rename(columns={'Last Update':'datetime'}, inplace=True)

    if 'Last_Update' in tmp_df.columns:
        tmp_df.rename(columns={'Last_Update':'datetime'}, inplace=True)

    #Lower case all col names
    tmp_df.columns = map(str.lower, tmp_df.columns) 

    for col in tmp_df[NUMERIC_COLS]:
        tmp_df[col] = tmp_df[col].fillna(0)
        tmp_df[col] = tmp_df[col].astype(int)

    return tmp_df

def get_data(cleaned_sheets):
    all_csv = []
    # Import all CSV's
    for f in tqdm(sorted(cleaned_sheets), desc='... loading data: '):
        if 'csv' in f:
            #print (f)
            try:
                tmp_df = pd.read_csv(os.path.join(DATA, f), index_col=None,header=0, parse_dates=['Last Update'])  
            except:
                # Temporary fix for JHU's bullshit data management
                tmp_df = pd.read_csv(os.path.join(DATA, f), index_col=None,header=0, parse_dates=['Last_Update'])  

            tmp_df = clean_data(tmp_df)
            tmp_df['date'] = tmp_df['datetime'].apply(get_date) # remove time to get date
            tmp_df['file_date'] = get_csv_date(f) #Get date of csv from file name
            tmp_df = tmp_df[KEEP_COLS]
            tmp_df['province'].fillna(tmp_df['country'], inplace=True) #If no region given, fill it with country
            all_csv.append(tmp_df)

    df_raw = pd.concat(all_csv, axis=0, ignore_index=True, sort=True)  # concatenate all csv's into one df
    df_raw = fix_country_names(df_raw)    # Fix mispelled country names
    df_raw = df_raw.sort_values(by=['datetime'])
    return df_raw


# use this function to fetch the data
def get():
        
    # Create Tmp Folder
    if not os.path.isdir(TMP_FOLDER):
        print('Creating folder...')
        print('...', TMP_FOLDER)
        os.mkdir(TMP_FOLDER)

    #Check if repo exists
    #git pull if it does
    if not os.path.isdir(TMP_GIT):
        clone_repo(TMP_FOLDER, REPO)
    else:
        try:
            print('git pull from', REPO)
            rep = git.Repo(TMP_GIT)
            rep.remotes.origin.pull()
        except:
            print('Could not pull from', REPO)
            sys.exit(1)

    sheets = os.listdir(DATA)
    
    # Clean the result to the sheet tabs we want
    print('Getting sheets...')
    cleaned_sheets = clean_sheet_names(sheets)

    # Aggregate all the data from sheets
    df = get_data(cleaned_sheets)
    #print("githubdata",df)
    #Clean the column names
    return df

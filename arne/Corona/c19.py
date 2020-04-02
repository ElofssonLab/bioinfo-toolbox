from lxml import html
#!/usr/bin/env python3

import os
from pathlib import Path
import re
import sys
import pandas as pd
from dateutil.parser import parse
from datetime import datetime,date,time, timedelta
from dateutil import parser

os.system('./c19.bash 2>/dev/null')
today=date.today()
yesterday=date.today() - timedelta(1)

inp_dir=home = str(Path.home())+"/Desktop/Corona/c19/data/"
out_file=home = str(Path.home())+"/Desktop/Corona/c19/data/c19-"+str(today)+".csv"
sys.stdout = open(out_file, "w")

columns=["country","province","DateRep", "confirmed", "new_confirmed_cases", "deaths", "IntensiveCare"]
country="Sweden"
alldata=[]
file_list=[]
print(f'{columns[0]:s},{columns[1]:s},{columns[2]:s},{columns[3]:s},{columns[4]:s},{columns[5]:s},{columns[6]:s}')

for f in os.listdir(inp_dir):
    if (f.find(".txt")!=-1):
        region=re.sub(r'.txt','',f)
        #print (region)
        row=0
        data={}
        with open(inp_dir+f) as file:
            lines =  [line.rstrip(']\n').lstrip('categorisdt: [') for line in file]
            for line in lines:
                #print (row,line)
                data[row]=line.split(',')
                row+=1
            for i in range(len(data[0])):
                #print(f'{f(x,n):d}, b={g(x,n):d}')
                d=re.sub(r" '|'","",data[0][i])+" 2020"
                date=datetime.strptime(d,"%b %d %Y").date()
                print(f'{country:s},{region:s},{date},{int(data[1][i]):d},{int(data[2][i]):d},{int(data[3][i]):d},{int(data[4][i]):d}')
                #printf  (country,region,i,data[0][i],data[1][i],data[2][i],data[3][i],data[4][i])
                #for j in range(len(data[0][0])):
                #alldata+=[[country],[region],[data[0][i]],[data[1][i]],[data[2][i]],[data[3][i]],[data[4][i]]]


            
        

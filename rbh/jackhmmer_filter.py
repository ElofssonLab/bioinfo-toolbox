import pandas as pd
import sys 

'''
Since jackhmmer does not have sequence identity, to choose the best hit we only use evalue as a parameter, which is: take the target has the best evalue as the hit of the query. And each for each query we only keep one hit.
'''
infile=sys.argv[1]

##overwrite the input file
df=pd.read_csv(infile,header=None,delimiter=' ')
df1=df.sort_values(2).groupby(1).first()
df1.reset_index(inplace=True)
df2=df1.sort_values(2).groupby(0).first()

df2.to_csv(infile,header=False,sep=' ')



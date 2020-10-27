import pandas as pd
import sys
import numpy as np
import argparse
import os


'''
given two jackhmmer rbh hit files(psiblast format or jackhmmer tblout format)get the target protein list based on species intersection.
psiblast format: 13 columns(last:species) 
jackhmmer format: 4 columns(query,target,evalue,species)

'''


parser = argparse.ArgumentParser(description = '''get target protein list based on species intersection.''')
parser.add_argument('--input', nargs=2, type= str, default=sys.stdin, help = 'hit file1 and hit file2')
parser.add_argument('--format',nargs=1, type= str, default=sys.stdin, help = 'psiblast hit format or jackhmmer format')
parser.add_argument('--out',nargs=1, type= str, default=sys.stdin, help = 'outputpath.')

args = parser.parse_args()

filename1=args.input[0].split('/')[-1]
filename2=args.input[1].split('/')[-1]

if args.format[0] == 'psiblast':
   df1=pd.read_csv(args.input[0],header=None,delimiter='\t')
   df2=pd.read_csv(args.input[1],header=None,delimiter='\t')
   mergedf=df1.merge(df2,on=12)

elif args.format[0] == 'hmmer':
   df1=pd.read_csv(args.input[0],header=None,delimiter=' ')
   df2=pd.read_csv(args.input[1],header=None,delimiter=' ')
   mergedf=df1.merge(df2,on=3)

else:
   print('wrong input format!')

if mergedf.empty:
   sys.exit("no proteins with common species!")

mergedf['1_x'].to_csv(args.out[0]+'/'+filename1,header=False,index=False)
mergedf['1_y'].to_csv(args.out[0]+'/'+filename2,header=False,index=False)





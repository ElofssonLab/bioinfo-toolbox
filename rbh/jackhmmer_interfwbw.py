import pandas as pd
import sys
import argparse
'''
get the intersection of forward hits and backward hits of jackhmmer ouput. 
'''

#Arguments for argparse module:
parser = argparse.ArgumentParser(description = "get the intersection of forward and backward hits of jackhmmer ouput and this is the last step of rbh using jackhmmer.")
parser.add_argument('-fwfile',nargs=1,required= True, type = str, default=sys.stdin, help = "forward hit file including path(after trimming).")
parser.add_argument('-bwfile',nargs=1,required= True,type = str, default=sys.stdin, help='backward hit file including path(after trimming).')
parser.add_argument('-o',nargs=1, required= True, type = str, default=sys.stdin, help="output file including path.")

##read file
args = parser.parse_args()
fw=pd.read_csv(args.fwfile[0],header=None,delimiter=' ')
bw=pd.read_csv(args.bwfile[0],header=None,delimiter=' ')
pd.merge(fw,bw,left_on=[0,1],right_on=[1,0])[['0_x','1_x','2_x']].to_csv(args.o[0],header=False,index=False,sep=' ')
##pname eg: UP123456





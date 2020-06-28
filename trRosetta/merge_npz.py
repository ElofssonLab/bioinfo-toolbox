#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import argparse
from argparse import RawTextHelpFormatter
        
##args = docopt.docopt(__doc__)
#out_dir = args['--output_folder']
 

p = argparse.ArgumentParser(description = '- Merging three NPZ files (each file + interaction area)',
                            formatter_class=RawTextHelpFormatter)
p.add_argument('-dataA','--inputA','-i', required= True, help='Input trRossetta NPZ file')
p.add_argument('-dataB','--inputB','-j', required= True, help='Input trRossetta NPZ file')
p.add_argument('-dataAB','--inputAB','-k', required= True, help='Input trRossetta NPZ file for combined sequence')
#p.add_argument('-seq','--sequence','-s', required= True, help='sequence file to identify domain baorders')
p.add_argument('-out','--output','-o', required= True, help='output NPX file')
#parser.add_argument('--nargs', nargs='+')
ns = p.parse_args()


bin_step = 0.5
bins = np.array([2.25+bin_step*i for i in range(36)])


rstA = np.load(ns.inputA)
rstB = np.load(ns.inputB)
rstAB = np.load(ns.inputAB)



distA = rstA["dist"].shape[0]
distB = rstB["dist"].shape[0]
distAB = rstAB["dist"].shape[0]


if distA+distB != distAB:
    print ("Not correct sized of distnance matrices",distA,distB,distAB)
    exit(-1)

new_rst = {'dist' : [], 'omega' : [], 'theta' : [], 'phi' : [] } # , 'rep' : []}    
for i in rstAB.files:
    A=rstA[i]
    B=rstB[i]
    AB=rstAB[i]
    AB[0:distA,0:distA]=A    
    AB[distA:,distA:]=B    
    new_rst[i]=AB

np.savez_compressed(ns.output, dist=new_rst['dist'], omega=new_rst['omega'], theta=new_rst['theta'], phi=new_rst['phi']) #, rep=new_rst['rep'])

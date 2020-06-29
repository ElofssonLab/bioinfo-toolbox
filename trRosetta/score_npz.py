#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import argparse
from argparse import RawTextHelpFormatter
        
##args = docopt.docopt(__doc__)
#out_dir = args['--output_folder']
 

p = argparse.ArgumentParser(description = '- Merging three NPZ files (each file + interaction area)',
                            formatter_class=RawTextHelpFormatter)
p.add_argument('-data','--input','-i', required= True, help='Input trRossetta NPZ file')
#p.add_argument('-seq','--sequence','-s', required= True, help='sequence file to identify domain baorders')
p.add_argument('-out','--output','-o', required= True, help='output NPX file')
p.add_argument('-seq','--sequence','-s', required= True, help='sequence file to identify domain baorders')
#parser.add_argument('--nargs', nargs='+')
ns = p.parse_args()


bin_step = 0.5
bins = np.array([2.25+bin_step*i for i in range(36)])

rst = np.load(ns.input)
dist = rst["dist"]
len=




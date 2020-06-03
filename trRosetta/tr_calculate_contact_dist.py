#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import argparse
from argparse import RawTextHelpFormatter
        
##args = docopt.docopt(__doc__)
#out_dir = args['--output_folder']
 

p = argparse.ArgumentParser(description = '- plotting trRosetta maps-',
                            formatter_class=RawTextHelpFormatter)
p.add_argument('-data','--input','-i', required= True, help='Input trRossetta NPZ file')
p.add_argument('-dom','--domain','-d', required= False, help='positions of domain borders', nargs='+')
p.add_argument('-out','--output','-o', required= False, help='output image')
#parser.add_argument('--nargs', nargs='+')
ns = p.parse_args()

input_file = np.load(ns.input)
bin_step = 0.5
bins = np.array([2.25+bin_step*i for i in range(36)])
dist = input_file["dist"]
p_len = dist.shape[0]
res = np.zeros((p_len, p_len))
res.fill(20)
np.fill_diagonal(res, 4)

for i in range(p_len-1):
    for j in range(i+1):
        prob = dist[i, j, 0]
        if prob > 0.5:
            continue
        d_slice = dist[i, j, 1:]
        mean_dist = np.sum(np.multiply(bins, d_slice/np.sum(d_slice)))
        res[i, j] = mean_dist
        res[j, i] = mean_dist

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set(title=ns.input)
cax = ax.matshow(res, cmap="hot")
print (res)
if type(ns.domain) is list:
    for cut in ns.domain:
        x=[0,p_len-1]
        y=[cut,cut]
        ax.plot(x,y,lw=3,c="b")
        ax.plot(y,x,lw=3,c="b")
fig.colorbar(cax)
if ns.output:
    fig.savefig(ns.output)
plt.show()

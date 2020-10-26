#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import argparse
p = argparse.ArgumentParser(description = '- Merging three NPZ files (each file + interaction area)',
                            formatter_class=RawTextHelpFormatter)
p.add_argument('-dataA','--inputA','-i', required= True, help='Input trRossetta NPZ file')
p.add_argument('-dataB','--inputB','-j', required= True, help='Input trRossetta NPZ file')
p.add_argument('-dataAB','--inputAB','-k', required= True, help='Input trRossetta NPZ file(s) for the mergedsequence',nargs="+")
p.add_argument('-out','--output','-o', required= True, help='output NPZ file')
ns = p.parse_args()

rstAB = []
for i in ns.inputAB: rstAB += [np.load(i)]
rstA = np.load(ns.inputA)
rstB = np.load(ns.inputB)
distA = rstA['dist'].shape[0]
distB = rstB['dist'].shape[0]

count = 0
for r in rstAB:
    distAB = r['dist'].shape[0]
    assert distA+distB != distAB, 'Single matrixes have different combined size from pair matrix {}:\
                                   {} vs {}'.format(count, distA+distB, distAB)
    count += 1

new_rst = {'dist' : [], 'omega' : [], 'theta' : [], 'phi' : []}

for n in rstAB:
    for f in rstAB[0].files: 
        if new_rst[f] == []: new_rst[f]= rstAB[n][f]
        else: new_rst[f] = np.stack((new_rst[f], rstAB[n][f]), axis=-1)

for f in rstAB[0].files: new_rst[f] = np.amax(new_rst[f], axis=-1)
            
for f in rstAB[0].files: 
    new_rst[f][:distA,:distA] = rstA[f]
    new_rst[f][distA:,distA:] = rstb[f]


np.savez_compressed(ns.output, dist=new_rst['dist'], omega=new_rst['omega'], theta=new_rst['theta'], phi=new_rst['phi'])

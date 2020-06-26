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
p.add_argument('-seq','--sequence','-s', required= True, help='sequence file to identify domain baorders')
p.add_argument('-out','--output','-o', required= True, help='output NPX file')
#parser.add_argument('--nargs', nargs='+')
ns = p.parse_args()


bin_step = 0.5
bins = np.array([2.25+bin_step*i for i in range(36)])


rstA = np.load(ns.inputA)
rstB = np.load(ns.inputB)
rstAB = np.load(ns.inputAB)

distA = rstA["dist"]
distB = rstB["dist"]
distAB = rstAB["dist"]


if distA+distB != distAB:
    print ("Not correct sized of distnance matrices",distA,distB,distAB)
    exit(-1)

p_lenA = dist.shape[0]
res = np.zeros((p_len, p_len))
res.fill(20)
np.fill_diagonal(res, 4)

#print(rst)



shift=0
new_rst = {'dist' : [], 'omega' : [], 'theta' : [], 'phi' : [], 'rep' : []}
#new_rst=rst
new_rst[x]=np.delete(rst[x],slice(m+shift,m+shift+len(sepseq)),1)


for m in borders:
    print ("Deleting :",m,m+len(sepseq))
    for x in rst.files:
        
        new_rst[x]=np.delete(new_rst[x],slice(m+shift,m+shift+len(sepseq)),0)
    shift+=len(sepseq)






borders=[]
if ns.sequence:
    sepseq=ns.sepseq
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    #import A3MIO
    # We can read it as fasta as we only care about the first sequence withouth gaps
    import re
    with open(ns.sequence, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq=record
            #print (record)
            break
    ns.domain=[]
    for m in re.finditer(sepseq,str(seq.seq)):
        borders+=[m.start()]
        #print(m.start(), m.group())
        for i in range(m.start(),m.start()+len(sepseq)):
            ns.domain+=[i]

        
# o        
#sys.exit()
#print (new_rst)

# Save
np.savez_compressed(ns.output, dist=new_rst['dist'], omega=new_rst['omega'], theta=new_rst['theta'], phi=new_rst['phi'])

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
p.add_argument('-seq','--sequence','-s', required= True, help='sequence file to identify domain baorders')
p.add_argument("--sepseq","-sep","-S",required=False, help='Separation sequence between protein in MSA' ,default="GGGGGGGGGGGGGGGGGGGG")
p.add_argument('-out','--output','-o', required= True, help='output NPX file')
#parser.add_argument('--nargs', nargs='+')
ns = p.parse_args()

rst = np.load(ns.input)
bin_step = 0.5
bins = np.array([2.25+bin_step*i for i in range(36)])
dist = rst["dist"]
p_len = dist.shape[0]
res = np.zeros((p_len, p_len))
res.fill(20)
np.fill_diagonal(res, 4)

#print(rst)

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

shift=0
new_rst = {'dist' : [], 'omega' : [], 'theta' : [], 'phi' : [] } # , 'rep' : []}
#new_rst=rst
for m in borders:
    print ("Deleting :",m,m+len(sepseq))
    for x in rst.files:
        #print (x)
        new_rst[x]=np.delete(rst[x],slice(m+shift,m+shift+len(sepseq)),1)
        new_rst[x]=np.delete(new_rst[x],slice(m+shift,m+shift+len(sepseq)),0)
    shift+=len(sepseq)
        
# o        
#sys.exit()
#print (new_rst)

# Save
np.savez_compressed(ns.output, dist=new_rst['dist'], omega=new_rst['omega'], theta=new_rst['theta'], phi=new_rst['phi'])# , rep=new_rst['rep'])

import pandas as pd
import argparse
import sys

'''
add boundary to already merged MSA. 
'''

parser = argparse.ArgumentParser(description = '''Insert amino acids btw two proteins in an already-merged MSA for obvious visualization. ''')
parser.add_argument('-input', nargs=1, type= str, default=sys.stdin, help ='merged MSA')
parser.add_argument('-seqlen1', nargs=1, type= int, default=sys.stdin, help ='seq length of the first protein.')
parser.add_argument('-aa_pattern', nargs=1, type= str, default=sys.stdin, help = 'Amino acid pattern to insert.')
parser.add_argument('-pattern_length', nargs=1, type= int, default=sys.stdin, help = 'Length of amino acid pattern to insert.')
parser.add_argument('-reverse',nargs=1, type= str, default=sys.stdin, help = 'if reversed alignment is required.')
parser.add_argument('-outdir', nargs=1, type= str, default=sys.stdin, help = 'outdir file to save concatenated alignment.')

args = parser.parse_args()
a3m = pd.read_csv(args.input[0],header=None)

## get all sequences from alignment to a list
sequence=[]
header1=[]
header2=[]
with open(args.input[0],'r') as f:
     for line in f:
         line=line.strip()
         if line.startswith('>'):
            header1.append('>'+line.split('>')[1])
            header2.append('>'+line.split('>')[2])
            continue
         else:
            sequence.append(line)

pattern = args.aa_pattern[0]*args.pattern_length[0]
## two directions: forward merge and reversed merge
sequenceF=[x[:args.seqlen1[0]]+ pattern + x[args.seqlen1[0]:] for x in sequence]
sequenceB=[x[args.seqlen1[0]:]+ pattern + x[:args.seqlen1[0]] for x in sequence]

a3mforward=a3m.copy()
a3mbackward=a3m.copy()

#FORWARD: just make changes to the sequences
a3mforward[0].loc[1::2]=sequenceF

#BACKWARD: make changes to the sequences + header reversed
header_rev=[header2[x]+header1[x] for x in range(len(header1))]
a3mbackward[0].loc[1::2]=sequenceB
a3mbackward[0].loc[0::2]=header_rev

pref=args.input[0].split('/')[-1].split('.')[0]
filename='_'.join(pref.split('_')[:8])

if args.reverse[0]=="True":
   a3mforward.to_csv(args.outdir[0]+'/'+filename+'.a3m',header=False,index=False)
   a3mbackward.to_csv(args.outdir[0]+'/'+filename+'_reversed.a3m',header=False,index=False)

else:
   a3mforward.to_csv(args.outdir[0]+'/'+filename+'.a3m',header=False,index=False)




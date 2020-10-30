import pandas as pd
import sys

##inputs: 
##      i: proteome name
## tmpdir: tmpdir is the root tmpdir of this array job. eg: tmp.ZFkWlVixhY (without '/' in the end)

i=sys.argv[1]
tmpdir=sys.argv[2]

fw=tmpdir+'/forward_trimmed/'+i
bw=tmpdir+'/backward_trimmed/'+i

df_for=pd.read_csv(fw,sep='\t',header=None)
df_for.columns=['qseqid','sseqid','qlen','slen','length','qcovs','pident','qstart','qend','sstart','send','evalue']
df_bac=pd.read_csv(bw,sep='\t',header=None)
df_bac.columns=['qseqid','sseqid','qlen','slen','length','qcovs','pident','qstart','qend','sstart','send','evalue']

#for anna's MI_IPA method 
df_for=df_for.drop_duplicates(subset=['qseqid','sseqid'],keep='last')
df_bac=df_bac.drop_duplicates(subset=['qseqid','sseqid'],keep='last')

pd.merge(df_for, df_bac, left_on=['qseqid','sseqid'], right_on=['sseqid','qseqid']).iloc[:,:12].to_csv(tmpdir+'/rbh_hit/'+i,header=None,index=None,sep='\t')

import pandas as pd
import os,sys
import numpy as np

def main():
    #import the file.
    file = sys.argv[1]
#    direction = sys.argv[2]
    dir = sys.argv[2]

    df = pd.read_csv(file,delimiter='\t',header=None)
    df.columns = ['qseqid', 'sseqid', 'qlen', 'slen', 'length', 'qcovs', 'pident', 'qstart', 'qend', 'sstart', 'send',
                  'evalue']

    #only keep the last iterative result 
    df=df.drop_duplicates(subset=['qseqid','sseqid'],keep='last')
    #apply the first filter:
    #evalue<0.01(already filtered in blast cmd settings) AND hit length>=50 or query coverage>=0.5
    #AND query coverage >=0.5 or hit coverage >=0.5
    df = df[(df.slen >= 50) | (df.qcovs >= 50)]
    df = df[(df.length / df.slen >= 0.5) | (df.length / df.slen >= 0.5)]
    
##    #apply for parallel filtering process
#    if direction == 'forward':
##    #forward:
#	abs_path = '/home/j/juliezhu/pfs/data/ref_proteomes/rbh_blast/blast_forward/'	
#    	forward_hit = filtering(df)
#    	forward_hit.to_csv(abs_path+file+'_fw', header=False, sep='\t', index=False)

#    elif direction == 'backward':
#	abs_path = '/home/j/juliezhu/pfs/data/ref_proteomes/rbh_blast/blast_backward/'
#    	backward_hit = filtering(df)
#	backward_hit.to_csv(abs_path+file+'_bw',header=False,sep='\t',index=False)
    hit = filtering(df)
    hit.to_csv(dir+'/'+file,header=False,sep='\t',index=False)
	
def filtering(df):
    '''
    do the filtering based on the workflow in Song's paper. start from two individual ranking parts
    :param df: data needs trimmed after applying the 3 AND filter(very first step)
    :return: final forward best hits for each ecoli sequence.
    '''

    #group each ecoli sequence, then the query would be a single sequence and target is a complete proteome.
    g1 = df.groupby('qseqid',sort=False)
    g1_name = [name for name,ff in g1]

    #an empty df to save results
    forw_orth=pd.DataFrame()
    for i in g1_name:
        if len(g1.get_group(i)) == 1:                   #if a proteome contains only one hit, take as the result.
            forw_orth = forw_orth.append(g1.get_group(i))
        else:
            g1_df = g1.get_group(i)
            #rank by identity
            ident_result = rank_identity(g1_df)
            eval_result = rank_eval(g1_df)
            forw_orth = forw_orth.append(pd.merge(ident_result, eval_result))


    return forw_orth





def rank_identity(df1):
    '''
    do the trimming for df ranking by identity.
    :param df1: input dataframe after the 1st filter.
    :return:best hits by identity
    '''
    df1_rank = df1.sort_values('pident', ascending=False)
    ident_out = df1_rank.iloc[[0]]
    for i in range(1, len(df1_rank)):
        pre_hit = set(range(int(ident_out.iloc[-1].sstart), int(ident_out.iloc[-1].send + 1)))
        cur_hit = set(range(int(df1_rank.iloc[i].sstart), int(df1_rank.iloc[i].send + 1)))
        overlap = pre_hit & cur_hit
        judgement1 = (len(overlap) >= 0.25 * df1_rank.iloc[i]['length']) & (
                    df1_rank.iloc[i]['pident'] > df1_rank.iloc[0]['pident'] - 10)
        judgement2 = (len(overlap) < 0.25 * df1_rank.iloc[i]['length']) & (
                    df1_rank.iloc[i]['pident'] > df1_rank.iloc[0]['pident'] - 20)
        if judgement1 | judgement2:
            ident_out = ident_out.append(df1_rank.iloc[[i]])
    return ident_out


def rank_eval(df2):
    '''
        do the trimming for df ranking by evalue.
        :param df1: input dataframe after the 1st filter.
        :return:best hits by evalue
    '''

    df2_rank=df2.sort_values('evalue')
    eval_out=df2_rank.iloc[[0]]
    for i in range(1,len(df2_rank)):
        pre_hit=set(range(int(eval_out.iloc[-1].sstart),int(eval_out.iloc[-1].send+1)))
        cur_hit=set(range(int(df2_rank.iloc[i].sstart),int(df2_rank.iloc[i].send+1)))
        overlap=pre_hit & cur_hit
        judgement1 = (len(overlap)>=0.25*df2_rank.iloc[i]['length']) & (df2_rank.iloc[i]['evalue']<df2_rank.iloc[0]['evalue']*1e5)
        judgement2 = len(overlap)<0.25*df2_rank.iloc[i]['length']
        if judgement1 | judgement2:
            eval_out=eval_out.append(df2_rank.iloc[[i]])
    return eval_out


if __name__== "__main__":
    main()

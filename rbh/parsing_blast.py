import pandas as pd
#import numpy as np
#import matplotlib.pyplot as plt
import os
import sys

def main():
    #import the files.	
    file=sys.argv[1]
    #file='test.m8'
    #tmpdir is the root tmpdir of this array job. eg: tmp.ZFkWlVixhY (without '/' in the end)    
    abs_blast=sys.argv[2]
		
    df = pd.read_csv(abs_blast+'/rbh_hit/'+file, delimiter='\t', header=None)
    df.columns = ['qseqid','sseqid','qlen','slen','length','qcovs','pident','qstart','qend','sstart','send','evalue']
    # apply evalue criteria
    df = df[df['evalue'] <= 0.01]

    #set up an empty df for complicated cases: do it manually
    manual_df = pd.DataFrame()

    #if df has multiple records
    if len(df[df.duplicated(subset='sseqid', keep=False)]) != 0:
        df,manual_df = multiTarget_trimming(df, manual_df)

    if len(df[df.duplicated(subset='qseqid', keep=False)]) != 0:
        df, manual_df = multiQuery_trimming(df, manual_df)

    #save new df as the final csv.

    # df.to_csv(proteome+'.m8')
    df.to_csv(abs_blast+'/rbh_hit/'+file,header=False,sep='\t',index=False)

    #check if manual_df is empty or not
    if not manual_df.empty:
       manual_df.to_csv(abs_blast+'/rbh_hit/manual/'+file,header=False,sep='\t',index=False)



def multiTarget_trimming(df,manual_df):
    '''
    remove the multiple records and resulting in a file where each target protein has one ecoli hit.
    '''

    # remove target proteins which has more than 3 hits. since only two records are accepted. Becase if there were still more than oen candidate in a species after the redundancy
    #removal and the merging partial hits, the protein should be excluded.

    df_dropt = df.groupby('sseqid', sort=False).filter(lambda x: len(x) >= 3)
    df = df.drop(index=df_dropt.index)

    #get target proteins with two ecoli hits and groupby 'Tid'.
    g1 = df[df.duplicated(subset='sseqid', keep=False)].sort_values('evalue').groupby('sseqid', sort=False)
    g1_name = [name for name, ff in g1]

    #continue the trimming process only if g1 is not empty.
    if g1_name == []:
        return df,manual_df

    #For each group, apply our criteria: for the target protein, if the alignments of two hits are in the same region,
    # take the ecoli protein with longer alignment. Otherwise, do the trimming manually(manual_df).
    for i in g1_name:
        g1_df = g1.get_group(i)
        sta1, end1 = g1_df.iloc[0]['sstart'], g1_df.iloc[0]['send']
        sta2, end2 = g1_df.iloc[1]['sstart'], g1_df.iloc[1]['send']
        range1, range2 = set(range(sta1, end1)), set(range(sta2, end2))
        idx=[end1-sta1,end2-sta2].index(max([end1-sta1,end2-sta2]))
        # second item includes the first one: drop 1st
        if list(range1 - range2) == []:
            df = df.drop(index=g1_df.index[0])
        # first item includes the second one: drop 2nd
        elif list(range2 - range1) == []:
            df = df.drop(index=g1_df.index[1])
        # two protein not fully overlapped: if overlapped area cover >=0.8 of the shorter protein, drop the shorter protein
        elif len(range1 & range2)/min(len(range1),len(range2))>0.8:
        #drop the shorter one
            df = df.drop(index=g1_df.index[1-idx])

        # intersection or no-overlap at all: drop from the result and do it manually
        else:
            df = df.drop(index=g1_df.index)
            manual_df = manual_df.append(g1_df)



    return df, manual_df


def multiQuery_trimming(df, manual_df):
    '''
    remove the multiple records and resulting in a file where each ecoli query protein has one hit.
    '''
    # remove ecoli proteins which has more than 3 hits. since only two records are accepted. Becase if there were still more than oen candidate in a species after the redundancy
    #removal and the merging partial hits, the protein should be excluded.oved.
    df_dropq = df.groupby('qseqid', sort=False).filter(lambda x: len(x) >= 3)
    df = df.drop(index=df_dropq.index)

    # get ecoli proteins with two target hits and groupby 'Qid'.
    g2 = df[df.duplicated(subset='qseqid', keep=False)].sort_values('evalue').groupby('qseqid', sort=False)
    g2_name = [name for name, ff in g2]

    if g2_name==[]:
        return df,manual_df

    # For each group, apply our criteria: for the ecoli protein, if the alignments of two hits are in the same region,
    # take the target protein with longer alignment. Otherwise, do the trimming manually(manual_df).
    for j in g2_name:
        g2_df = g2.get_group(j)
        sta1, end1 = g2_df.iloc[0]['qstart'], g2_df.iloc[0]['qend']
        sta2, end2 = g2_df.iloc[1]['qstart'], g2_df.iloc[1]['qend']
        range1, range2 = set(range(sta1, end1)), set(range(sta2, end2))
        idx=[end1-sta1,end2-sta2].index(max([end1-sta1,end2-sta2]))
        if list(range1 - range2) == []:
            df = df.drop(index=g2_df.index[0])
        # first item includes the second one: drop 2nd
        elif list(range2 - range1) == []:
            df = df.drop(index=g2_df.index[1])
        # two protein not fully overlapped: if overlapped area cover >=0.8 of the shorter protein, drop the shorter protein
        elif len(range1 & range2)/min(len(range1),len(range2))>0.8:
        #drop the shorter one
            df = df.drop(index=g2_df.index[1-idx])

        # intersection or no-overlap at all: drop from the result and do it manually
        else:
            df = df.drop(index=g2_df.index)
            manual_df = manual_df.append(g2_df)

    return df, manual_df


if __name__== "__main__":
   main()


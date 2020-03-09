#!/usr/bin/env python3 
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from argparse import RawTextHelpFormatter
from keras import backend as K
from keras.models import load_model, Model
import argparse
import csv
import datetime as day
import h5py
import matplotlib.pyplot as plt
import nnlab as nl
import numpy as np
import os
import pandas as pd
import pickle
import random
import re
import seaborn as sb
import tensorflow as tf


codons = [
'ATA', 'ATC', 'ATT', 'ATG', 'ACA', 'ACC', 'ACG', 'ACT', 
'AAC', 'AAT', 'AAA', 'AAG', 'AGC', 'AGT', 'AGA', 'AGG',                  
'CTA', 'CTC', 'CTG', 'CTT', 'CCA', 'CCC', 'CCG', 'CCT', 
'CAC', 'CAT', 'CAA', 'CAG', 'CGA', 'CGC', 'CGG', 'CGT', 
'GTA', 'GTC', 'GTG', 'GTT', 'GCA', 'GCC', 'GCG', 'GCT', 
'GAC', 'GAT', 'GAA', 'GAG', 'GGA', 'GGC', 'GGG', 'GGT', 
'TCA', 'TCC', 'TCG', 'TCT', 'TTC', 'TTT', 'TTA', 'TTG', 
    'TAC', 'TAT', 'TGC', 'TGT', 'TGG',
    "TAG","TAA","TGA"] # These are the stop codons - they will be 0s

nuc_encode = {
'A':[1,0,0,0],
'T':[0,1,0,0],
'C':[0,0,1,0],
'G':[0,0,0,1]}

res_encode = {
'A':[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
'R':[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
'N':[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
'D':[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
'C':[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
'Q':[0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
'E':[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
'G':[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
'H':[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
'I':[0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
'L':[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],
'K':[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],
'M':[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
'F':[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
'P':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],
'S':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
'T':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],
'W':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],
'Y':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
'V':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],
'_':[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
}

codons2aa = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

label = {'S':0, 'X':0.5, 'C':0.5, 'D':1}

##### codon wise one-hot encoding table
pos = 0
cod_encode = {}
for codon in codons: 
    cod_encode[codon] = []
    for n in range(61): 
        if n == pos: cod_encode[codon].append(1)
        else: cod_encode[codon].append(0)
    pos += 1
#####

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description =
                                         '- Running disorder prediction -',
        formatter_class=RawTextHelpFormatter)
    parser.add_argument('-f', required= True, help='Input file')
    #parser.add_argument('-m', required= True, help='Model file')
    parser.add_argument('-gc', required= False, help='GC', action='store_true')
    parser.add_argument('-rna', required= False, help='RNA', action='store_true')
    parser.add_argument('-pro', required= False, help='Protein', action='store_true')
    parser.add_argument('-gcgenomic', required= False, help='GCgenomic', action='store_true')
    parser.add_argument('-gcvalue', required= False,help= "GCgenomic value")
    parser.add_argument('-rnafile',required= False, action='store_true', help=' Specifies that inpout is an rna file (required if -gc or -rna')
    #parser.add_argument('-kingdom', required= False, help='Kingdom', action='store_true')
    ns = parser.parse_args()

    print ("ARGs: ",ns)
    
    if ns.gc or ns.rna: ns.rnafile=True
    
    if (ns.gc): ns.gc='GC'
    else: ns.gc='noGC'

    cutoff =0.4 # Prediction cutoff.
    seed = 42
    tiny=1.e-20
    os.environ['PYTHONHASHSEED'] = '0'
    np.random.seed(seed)
    random.seed(seed)
    config = tf.ConfigProto()
    config.gpu_options.allow_growth = True
    tf.keras.backend.set_session(tf.Session(config=config))
    session_conf = tf.ConfigProto(intra_op_parallelism_threads=1, inter_op_parallelism_threads=1)
    tf.set_random_seed(seed)
    sess = tf.Session(graph=tf.get_default_graph(), config=session_conf)
    K.set_session(sess)

    if ns.pro: feat_len=20
    elif ns.rna: feat_len=61
    else: sys.exit('Unknown ')
    if ns.gc: feat_len+=1
    if ns.gcgenomic:  feat_len+=1

    #   First we need to create the input file from a fasta file (and
    #   potentially other information as well)


    type=''
    extra=0
    if ns.pro:
        extra=20
        type+="pro"
    elif ns.rna:
        extra=21
        type+="rna"

    if ns.gc:
        type+="_GC"
        extra+=1
    if ns.gcgenomic:
        type+="_GCgenomic"
        extra+=1

    #if ns.kingdom:
    #    type+="_kingdom"
    #    extra+=3

    
    # This is GC data for genomes ("we could also just add that in the command line")
    data_dir = "/home/arnee/git/paper_proks_vs_euks_proteins/data/"
    df_reference_file = data_dir  + "df_reference.csv"
    df_reference = pd.read_csv(df_reference_file)
    taxid2gc = df_reference.set_index("TaxID").to_dict()["GC%"]

    fastafile=ns.f
    #GCgenomic=ns.gcvalue

    memo2taxid={}
    with open(data_dir+'speclist.txt') as fp:
        for line in fp:
            if re.match('.*:\sN=',line):
                s=line.split()
                taxid=s[2].replace(':','')
                mnemonic=s[0]
                memo2taxid[mnemonic]=taxid
                
    
    outlist = open('formatted_list','w')

    #dataset=np.array()
    fasta_dir="uniprot/"

    # We skip kingdom for the moment (it did not improve performance)
    #if data[key]['kingdom']=='A':
    #    k=[1,0,0]
    #elif data[key]['kingdom']=='B':
    #    k=[0,1,0]
    #elif data[key]['kingdom']=='E':
    #    k=[0,0,1]
    #else:
    #    k=[0,0,0]
    
    seq=[]
    for record in SeqIO.parse(fastafile, "fasta"):
        #print("%s %i" % (record.id, len(record)))
        bar,id,name=record.id.split("|") 
        #print (record.id,id,name)
        # try:
        if ns.gcvalue:
            GCgenomic=float(ns.gcvalue)
        elif ns.gcgenomic:
            try:
                prot,taxname=name.split("_") 
                GCgenomic=float(taxid2gc[int(memo2taxid[taxname])])
            except:
                continue
        else:
            GCgenomic=0.50
        if ns.rnafile:
            GC=0.
            print("RNARECORD",record)
            for n in record.seq:
                if n in ["G","C"]:
                    GC+=1.
            print ("GCTEST",GC,len(record.seq))
            GC=GC/len(record.seq)
            for i in range(0,len(record.seq), 3):
                triplett=record.seq[i]+record.seq[i+1]+record.seq[i+2]
                print ("TRIPLETT",i,triplett)
                # This assumes we read a DNA file.
                if ns.pro:
                    # This assumed we translat RNA to AA
                    aa=codons2aa[triplett]
                    code=res_encode[aa][:]
                else:
                    code=cod_encode[triplett][:]
                    print ("TEST",triplett,code)
                if ns.gc:
                    code.append(GC)
                if ns.gcgenomic:
                    code.append(GCgenomic)
                    #if ns.kingdom:
                    #    code.append(kingdom)
                seq.append(code)
            print ("RNA: ",seq)
        else:
            print("PRORECORD",record)

            for aa in record.seq:
                code=res_encode[aa][:]
                if ns.gcgenomic:
                    code.append(GCgenomic)
                #if ns.kingdom:
                #    code.append(kingdom)
                seq.append(code)
            print ("PRO: ",seq)    
        print ("FINAL",seq)
        sample = np.array(seq, dtype=np.float64)
        X = sample[:,:-1].reshape(1, len(sample), len(sample[0])-1)
        Y = sample[:,-1]
        #print (protein,X,Y)
        prediction = model.predict_on_batch(X)
        print(prediction)


